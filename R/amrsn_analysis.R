# CARE-INDIA (Carbapenem resistance among Gram-negative ESKAPEE pathogens in India, 2017-2024)
# Primary trend analysis: Cochran-Armitage trend tests and weighted linear regression of carbapenem resistance (2017-2024), with a segmented-regression/Davies breakpoint sensitivity, main-text resistance figures and the isolate-burden figure.
# Authors: Abhishek Akella, Anand Srinivasan  |  Dept of Pharmacology, AIIMS Bhubaneswar, India
# License: MIT (see LICENSE)

set.seed(42)

# Setup: portable root detection, output directory, and dataset loading
suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

detect_root <- function() {
  candidates <- c(getwd(),
                  dirname(sys.frame(1)$ofile %||% ""),
                  ".")
  for (p in candidates) {
    if (nzchar(p) && file.exists(file.path(p, "R", "amrsn_data_shared.R"))) {
      return(normalizePath(p, winslash = "/", mustWork = FALSE))
    }
  }
  stop("Could not locate R/amrsn_data_shared.R. Run from the project root.")
}
`%||%` <- function(a, b) if (is.null(a) || !nzchar(a)) b else a

root <- detect_root()
setwd(root)

dir.create(file.path(root, "output"), showWarnings = FALSE, recursive = TRUE)

source(file.path(root, "R", "amrsn_data_shared.R"))


# 1. Cochran-Armitage trend tests (resistance proportion vs year)
run_trend_test <- function(data, organism_name, abx, res_col, tested_col) {
  d <- data %>% filter(organism == organism_name)
  tt <- prop.trend.test(
    x     = d[[res_col]],
    n     = d[[tested_col]],
    score = d$year
  )
  tibble(
    organism   = organism_name,
    antibiotic = abx,
    chi_sq     = unname(tt$statistic),
    df         = unname(tt$parameter),
    p_value    = tt$p.value
  )
}

# Full ESKAPEE panel (five organisms) for the trend tests, matching the manuscript
organisms <- c("K. pneumoniae", "A. baumannii", "P. aeruginosa", "Enterobacter spp.", "E. coli")

trend_results <- bind_rows(
  map_dfr(organisms, ~run_trend_test(amrsn_ee, .x, "Imipenem",  "imi_res", "imi_tested")),
  map_dfr(organisms, ~run_trend_test(amrsn_ee, .x, "Meropenem", "mer_res", "mer_tested"))
) %>%
  mutate(significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE            ~ "ns"
  ))

cat("\n═══ COCHRAN–ARMITAGE TREND TEST RESULTS ═══\n\n")
print(trend_results, n = Inf)


# 2. Absolute and annualised resistance change, 2017 to 2024
summary_stats <- amrsn_ee %>%
  group_by(organism) %>%
  summarise(
    imi_res_2017   = imi_res_pct[year == 2017],
    imi_res_2024   = imi_res_pct[year == 2024],
    imi_abs_change = imi_res_2024 - imi_res_2017,
    imi_annualised = imi_abs_change / 7,
    mer_res_2017   = mer_res_pct[year == 2017],
    mer_res_2024   = mer_res_pct[year == 2024],
    mer_abs_change = mer_res_2024 - mer_res_2017,
    mer_annualised = mer_abs_change / 7,
    total_isolates = sum(total_n),
    .groups = "drop"
  )

cat("\n═══ SUMMARY: RESISTANCE CHANGE 2017→2024 ═══\n\n")
print(summary_stats, width = Inf)


# 3. Weighted linear regression: slope and 95% CI
run_weighted_lm <- function(data, organism_name, abx, res_pct_col, tested_col) {
  d   <- data %>% filter(organism == organism_name)
  fit <- lm(as.formula(paste(res_pct_col, "~ year")), data = d, weights = d[[tested_col]])
  ci  <- confint(fit, "year", level = 0.95)
  tibble(
    organism         = organism_name,
    antibiotic       = abx,
    slope_ppt_per_yr = unname(coef(fit)["year"]),
    ci_lower         = ci[1],
    ci_upper         = ci[2],
    r_squared        = summary(fit)$r.squared,
    p_value          = summary(fit)$coefficients["year", "Pr(>|t|)"]
  )
}

lm_results <- bind_rows(
  map_dfr(organisms, ~run_weighted_lm(amrsn_ee, .x, "Imipenem",  "imi_res_pct", "imi_tested")),
  map_dfr(organisms, ~run_weighted_lm(amrsn_ee, .x, "Meropenem", "mer_res_pct", "mer_tested"))
)

cat("\n═══ WEIGHTED LINEAR REGRESSION: SLOPE (PPT/YEAR) ═══\n\n")
print(lm_results, width = Inf)


# 3b. Breakpoint sensitivity: segmented regression + Davies test
#     Weighted lm per organism-carbapenem pair (weights = isolates tested), Davies
#     test for the existence of a breakpoint, and a segmented fit anchored at psi=2020.5.
#     With only eight annual points, segmented() often does not converge (reported as
#     seg_converged = FALSE); the Davies p-value is the primary breakpoint statistic.
if (!requireNamespace("segmented", quietly = TRUE)) {
  tryCatch(install.packages("segmented", repos = "https://cloud.r-project.org", quiet = TRUE),
           error = function(e) invisible(NULL))
}
library(segmented)

run_davies_segmented <- function(data, organism_name, abx, res_pct_col, tested_col) {
  d    <- data %>% filter(organism == organism_name)
  fit0 <- lm(as.formula(paste(res_pct_col, "~ year")), data = d, weights = d[[tested_col]])
  dav  <- tryCatch(
    suppressWarnings(davies.test(fit0, seg.Z = ~year, k = 10)$p.value),
    error = function(e) NA_real_
  )
  seg  <- tryCatch({
    fs <- segmented(fit0, seg.Z = ~year, psi = list(year = 2020.5),
                    control = seg.control(it.max = 100, tol = 1e-3, n.boot = 20))
    list(bp = unname(summary(fs)$psi[, "Est."]), converged = TRUE)
  }, error = function(e) list(bp = NA_real_, converged = FALSE))
  tibble(
    organism      = organism_name,
    antibiotic    = abx,
    davies_p      = round(dav, 4),
    breakpoint    = round(seg$bp, 1),
    seg_converged = seg$converged
  )
}

davies_results <- bind_rows(
  map_dfr(organisms, ~run_davies_segmented(amrsn_ee, .x, "Imipenem",  "imi_res_pct", "imi_tested")),
  map_dfr(organisms, ~run_davies_segmented(amrsn_ee, .x, "Meropenem", "mer_res_pct", "mer_tested"))
)

cat("\n═══ DAVIES BREAKPOINT TEST (segmented regression) ═══\n\n")
print(davies_results, n = Inf)
write_csv(davies_results, file.path(root, "output", "table_davies_segmented_breakpoint.csv"))


# 4. Figures
theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(fill = NA, colour = "grey70"),
    legend.position  = "bottom",
    strip.text       = element_text(face = "bold", size = 11),
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle    = element_text(hjust = 0.5, size = 10, colour = "grey40")
  )

organism_colours <- c(
  "K. pneumoniae"     = "#E41A1C",
  "A. baumannii"      = "#377EB8",
  "P. aeruginosa"     = "#4DAF4A",
  "Enterobacter spp." = "#FF7F00"
)

fig1 <- amrsn_long %>%
  filter(antibiotic == "Meropenem") %>%
  ggplot(aes(x = year, y = resistance_pct, colour = organism, shape = organism)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_colour_manual(values = organism_colours) +
  scale_x_continuous(breaks = 2017:2024) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    title    = "Meropenem Resistance Trends — Gram-Negative ESKAPE",
    subtitle = "ICMR-AMRSN, 2017–2024",
    x = "Year", y = "Resistance (%)",
    colour = "Organism", shape = "Organism"
  ) +
  theme_pub

ggsave(file.path(root, "output", "fig1_meropenem_trends.png"), fig1,
       width = 10, height = 6, dpi = 300)
ggsave(file.path(root, "output", "fig1_meropenem_trends.pdf"), fig1,
       width = 10, height = 6)

fig2 <- amrsn_long %>%
  filter(antibiotic == "Imipenem") %>%
  ggplot(aes(x = year, y = resistance_pct, colour = organism, shape = organism)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_colour_manual(values = organism_colours) +
  scale_x_continuous(breaks = 2017:2024) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    title    = "Imipenem Resistance Trends — Gram-Negative ESKAPE",
    subtitle = "ICMR-AMRSN, 2017–2024",
    x = "Year", y = "Resistance (%)",
    colour = "Organism", shape = "Organism"
  ) +
  theme_pub

ggsave(file.path(root, "output", "fig2_imipenem_trends.png"), fig2,
       width = 10, height = 6, dpi = 300)
ggsave(file.path(root, "output", "fig2_imipenem_trends.pdf"), fig2,
       width = 10, height = 6)

fig3 <- amrsn_long %>%
  ggplot(aes(x = year, y = resistance_pct, colour = antibiotic, linetype = antibiotic)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~organism, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c("Imipenem" = "#D95F02", "Meropenem" = "#1B9E77")) +
  scale_x_continuous(breaks = seq(2017, 2024, 2)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    title    = "Carbapenem Resistance by Organism and Antibiotic",
    subtitle = "ICMR-AMRSN, 2017–2024",
    x = "Year", y = "Resistance (%)",
    colour = "Carbapenem", linetype = "Carbapenem"
  ) +
  theme_pub

ggsave(file.path(root, "output", "fig3_faceted_organism.png"), fig3,
       width = 10, height = 8, dpi = 300)

fig4 <- amrsn %>%
  ggplot(aes(x = year, y = total_n / 1000, fill = organism)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = organism_colours) +
  scale_x_continuous(breaks = 2017:2024) +
  labs(
    title    = "Number of Isolates per Year by Organism",
    subtitle = "ICMR-AMRSN, 2017–2024",
    x = "Year", y = "Isolates (×1000)", fill = "Organism"
  ) +
  theme_pub

ggsave(file.path(root, "output", "fig4_isolate_counts.png"), fig4,
       width = 10, height = 6, dpi = 300)


# 5. Exports
annualised_summary <- amrsn_ee %>%
  group_by(organism) %>%
  summarise(
    imi_total_change_pp = imi_res_pct[year == 2024] - imi_res_pct[year == 2017],
    imi_pp_per_year     = imi_total_change_pp / 7,
    mer_total_change_pp = mer_res_pct[year == 2024] - mer_res_pct[year == 2017],
    mer_pp_per_year     = mer_total_change_pp / 7,
    .groups = "drop"
  )

write_csv(annualised_summary, file.path(root, "output", "table_annualised_pp_change.csv"))
write_csv(trend_results,       file.path(root, "output", "table_cochran_armitage_results.csv"))
write_csv(lm_results,          file.path(root, "output", "table_weighted_lm_slopes.csv"))


# 6. Session info (reproducibility)
writeLines(capture.output(sessionInfo()),
           file.path(root, "output", "sessionInfo_amrsn_analysis.txt"))

cat("\n═══ ALL OUTPUTS SAVED TO ", file.path(root, "output"), " ═══\n", sep = "")
