# CARE-INDIA (Carbapenem resistance among Gram-negative ESKAPEE pathogens in India, 2017-2024)
# Forecast sensitivity: 2028 projections for the ten organism-carbapenem pairs under alternative fitting windows and forms
# (primary 2017-2024 logit trend; 2020-2024 refit; leave-one-year-out refits; piecewise trend with a 2019 knot), with 95% prediction intervals.
# Authors: Abhishek Akella, Anand Srinivasan  |  Dept of Pharmacology, AIIMS Bhubaneswar, India
# License: MIT (see LICENSE)

suppressPackageStartupMessages(library(tidyverse))
root <- getwd()
if (!file.exists(file.path(root, "R", "amrsn_data_shared.R"))) stop("Run from the repository root; R/amrsn_data_shared.R not found.")
source(file.path(root, "R", "amrsn_data_shared.R"))
OUT <- file.path(root, "output", "forecast"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

EPS <- 5e-4
lgt <- function(p) { pp <- pmax(EPS, pmin(1 - EPS, p / 100)); log(pp / (1 - pp)) }
inv <- function(x) 100 / (1 + exp(-x))
CARB <- c("Imipenem", "Meropenem"); TARGET <- 2028

proj <- function(d, knot = NULL) {
  if (is.null(knot)) { m <- stats::lm(lgt(res_pct) ~ year, data = d, weights = tested); nd <- data.frame(year = TARGET) }
  else { d$pw <- pmax(0, d$year - knot); m <- stats::lm(lgt(res_pct) ~ year + pw, data = d, weights = tested)
         nd <- data.frame(year = TARGET, pw = TARGET - knot) }
  p <- stats::predict(m, nd, interval = "prediction", level = 0.95, weights = mean(d$tested))
  tibble(point = inv(p[1, "fit"]), lo95 = inv(p[1, "lwr"]), hi95 = inv(p[1, "upr"]), n_years = nrow(d),
         slope_logit = unname(stats::coef(m)[["year"]] + if (!is.null(knot)) stats::coef(m)[["pw"]] else 0))
}

carb <- full_long %>% filter(drug %in% CARB) %>% arrange(organism, drug, year)
rows <- carb %>% group_by(organism, drug) %>% group_modify(function(d, k) {
  bind_rows(
    proj(d) %>% mutate(scenario = "primary: 2017-2024 logit trend"),
    proj(d %>% filter(year >= 2020)) %>% mutate(scenario = "refit on 2020-2024 only"),
    proj(d, knot = 2019) %>% mutate(scenario = "piecewise trend, knot 2019 (post-2019 slope extrapolated)"),
    map_dfr(2017:2024, function(y) proj(d %>% filter(year != y)) %>% mutate(scenario = sprintf("leave-one-year-out: drop %d", y))))
}) %>% ungroup() %>% select(organism, drug, scenario, n_years, point, lo95, hi95, slope_logit)
write_csv(rows, file.path(OUT, "table_forecast_scenarios_2028.csv"))

summ <- rows %>% group_by(organism, drug) %>% summarise(
  primary_2028 = point[scenario == "primary: 2017-2024 logit trend"],
  primary_pi95 = sprintf("%.1f-%.1f", lo95[scenario == "primary: 2017-2024 logit trend"], hi95[scenario == "primary: 2017-2024 logit trend"]),
  refit_2020_2024 = point[scenario == "refit on 2020-2024 only"],
  piecewise_2019 = point[scenario == "piecewise trend, knot 2019 (post-2019 slope extrapolated)"],
  loyo_min = min(point[grepl("leave-one", scenario)]), loyo_max = max(point[grepl("leave-one", scenario)]),
  all_scenarios_min = min(point), all_scenarios_max = max(point), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ round(., 1)))
write_csv(summ, file.path(OUT, "table_forecast_scenarios_2028_summary.csv"))
cat("2028 projections across scenarios (ten carbapenem pairs):\n"); print(as.data.frame(summ), row.names = FALSE)
cat("\nWrote output/forecast/: table_forecast_scenarios_2028.csv, table_forecast_scenarios_2028_summary.csv\n")
