# CARE-INDIA (Carbapenem resistance among Gram-negative ESKAPEE pathogens in India, 2017-2024)
# Extended trend tests: weighted regression and Cochran-Armitage across all 47 organism-drug pairs with Holm and Benjamini-Hochberg correction; unweighted sensitivity; pre-/post-2020 interaction model and piecewise grid search (all candidate knots) for the ten organism-carbapenem pairs.
# Authors: Abhishek Akella, Anand Srinivasan  |  Dept of Pharmacology, AIIMS Bhubaneswar, India
# License: MIT (see LICENSE)

suppressPackageStartupMessages(library(tidyverse))

root <- getwd()
if (!file.exists(file.path(root, "R", "amrsn_data_shared.R"))) stop("Run from the repository root; R/amrsn_data_shared.R not found.")
source(file.path(root, "R", "amrsn_data_shared.R"))
OUT <- file.path(root, "output", "trends")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

fl <- full_long %>% arrange(organism, drug, year)
CARB <- c("Imipenem", "Meropenem")

lm_pairs <- fl %>% group_by(organism, drug) %>% group_modify(function(d, k) {
  fw <- stats::lm(res_pct ~ year, data = d, weights = tested)
  fu <- stats::lm(res_pct ~ year, data = d)
  cw <- stats::confint(fw, "year"); cu <- stats::confint(fu, "year")
  tibble(n_years = nrow(d), n_isolates = sum(d$tested),
         slope_weighted = stats::coef(fw)[["year"]], ci_lo_weighted = cw[1], ci_hi_weighted = cw[2],
         r2_weighted = summary(fw)$r.squared, p_weighted = summary(fw)$coefficients["year", "Pr(>|t|)"],
         slope_unweighted = stats::coef(fu)[["year"]], ci_lo_unweighted = cu[1], ci_hi_unweighted = cu[2],
         p_unweighted = summary(fu)$coefficients["year", "Pr(>|t|)"])
}) %>% ungroup() %>%
  mutate(p_weighted_holm = stats::p.adjust(p_weighted, "holm"), p_weighted_bh = stats::p.adjust(p_weighted, "BH")) %>%
  arrange(p_weighted)
write_csv(lm_pairs, file.path(OUT, "table_trend_weighted_lm_all_pairs.csv"))

ca_pairs <- fl %>% group_by(organism, drug) %>% group_modify(function(d, k) {
  tt <- stats::prop.trend.test(x = d$susceptible, n = d$tested, score = d$year)
  tibble(n_years = nrow(d), chi_sq = unname(tt$statistic), p_value = tt$p.value)
}) %>% ungroup() %>%
  mutate(p_holm = stats::p.adjust(p_value, "holm"), p_bh = stats::p.adjust(p_value, "BH")) %>% arrange(p_value)
write_csv(ca_pairs, file.path(OUT, "table_trend_cochran_armitage_all_pairs.csv"))

carb <- fl %>% filter(drug %in% CARB)
inter <- carb %>% group_by(organism, drug) %>% group_modify(function(d, k) {
  d <- d %>% mutate(yc = year - 2020, period = factor(year >= 2020, levels = c(FALSE, TRUE), labels = c("pre", "post")))
  m <- stats::lm(res_pct ~ yc * period, data = d, weights = tested); cf <- stats::coef(m)
  tibble(mean_pre = mean(d$res_pct[d$period == "pre"]), mean_post = mean(d$res_pct[d$period == "post"]),
         slope_pre = cf[["yc"]], slope_post = cf[["yc"]] + cf[["yc:periodpost"]], slope_difference = cf[["yc:periodpost"]],
         p_interaction = summary(m)$coefficients["yc:periodpost", "Pr(>|t|)"])
}) %>% ungroup() %>%
  mutate(p_holm = stats::p.adjust(p_interaction, "holm"), p_bh = stats::p.adjust(p_interaction, "BH"))
write_csv(inter, file.path(OUT, "table_covid_interaction_carbapenems.csv"))

grid <- carb %>% group_by(organism, drug) %>% group_modify(function(d, k) {
  fit0 <- stats::lm(res_pct ~ year, data = d, weights = tested)
  map_dfr(2019:2022, function(bp) {
    d$pw <- pmax(0, d$year - bp)
    fit1 <- stats::lm(res_pct ~ year + pw, data = d, weights = tested); a <- stats::anova(fit0, fit1)
    tibble(candidate_knot = bp, aic_null = stats::AIC(fit0), aic_piecewise = stats::AIC(fit1), delta_aic = stats::AIC(fit1) - stats::AIC(fit0),
           f_statistic = a$F[2], f_p_nominal = a$`Pr(>F)`[2],
           slope_before = stats::coef(fit1)[["year"]], slope_after = stats::coef(fit1)[["year"]] + stats::coef(fit1)[["pw"]],
           r2_null = summary(fit0)$r.squared, r2_piecewise = summary(fit1)$r.squared)
  })
}) %>% ungroup()
write_csv(grid, file.path(OUT, "table_piecewise_grid_all_knots.csv"))
write_csv(grid %>% group_by(organism, drug) %>% slice_min(aic_piecewise, n = 1, with_ties = FALSE) %>% ungroup(),
          file.path(OUT, "table_piecewise_grid_best_knot.csv"))

cat(sprintf("\nPairs tested: %d (weighted regression); Holm-significant: %d; BH-significant: %d\n",
            nrow(lm_pairs), sum(lm_pairs$p_weighted_holm < 0.05), sum(lm_pairs$p_weighted_bh < 0.05)))
print(as.data.frame(lm_pairs %>% filter(drug %in% c(CARB, "Ertapenem"), organism %in% c("K. pneumoniae", "E. coli")) %>%
  transmute(organism, drug, slope_weighted = round(slope_weighted, 2), ci = sprintf("%.2f to %.2f", ci_lo_weighted, ci_hi_weighted),
            p_weighted = signif(p_weighted, 2), p_holm = signif(p_weighted_holm, 2), p_bh = signif(p_weighted_bh, 2),
            slope_unweighted = round(slope_unweighted, 2), ci_unw = sprintf("%.2f to %.2f", ci_lo_unweighted, ci_hi_unweighted))), row.names = FALSE)
cat("\nInteraction model (period split at 2020):\n")
print(as.data.frame(inter %>% mutate(across(where(is.numeric), ~ round(., 3)))), row.names = FALSE)
cat("\nGrid search, K. pneumoniae meropenem, all candidate knots:\n")
print(as.data.frame(grid %>% filter(organism == "K. pneumoniae", drug == "Meropenem") %>% mutate(across(where(is.numeric), ~ round(., 3)))), row.names = FALSE)
cat("\nWrote output/trends/: table_trend_weighted_lm_all_pairs.csv, table_trend_cochran_armitage_all_pairs.csv, table_covid_interaction_carbapenems.csv, table_piecewise_grid_all_knots.csv, table_piecewise_grid_best_knot.csv\n")
