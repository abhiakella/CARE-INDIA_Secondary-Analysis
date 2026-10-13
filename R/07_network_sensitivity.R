# CARE-INDIA (Carbapenem resistance among Gram-negative ESKAPEE pathogens in India, 2017-2024)
# Network-composition sensitivity: weighted carbapenem trends restricted to the years with a stable, RC-coded network (2020-2024) and to
# 2019-2024 (21 participating centres), and with annual network throughput (log total isolates reported) or listed regional centres as a covariate.
# Authors: Abhishek Akella, Anand Srinivasan  |  Dept of Pharmacology, AIIMS Bhubaneswar, India
# License: MIT (see LICENSE)

suppressPackageStartupMessages(library(tidyverse))
root <- getwd()
if (!file.exists(file.path(root, "R", "amrsn_data_shared.R"))) stop("Run from the repository root; R/amrsn_data_shared.R not found.")
source(file.path(root, "R", "amrsn_data_shared.R"))
OUT <- file.path(root, "output", "network"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
net <- read_csv(file.path(root, "data", "network_by_year.csv"), show_col_types = FALSE) %>%
  mutate(log_total = log(total_isolates_reported))
carb <- full_long %>% filter(drug %in% c("Imipenem", "Meropenem")) %>% left_join(net, by = "year") %>% arrange(organism, drug, year)
slope <- function(d, formula, label) {
  m <- stats::lm(formula, data = d, weights = tested); ci <- stats::confint(m, "year")
  tibble(analysis = label, n_years = nrow(d), slope = unname(stats::coef(m)[["year"]]), ci_lo = ci[1], ci_hi = ci[2],
         p = summary(m)$coefficients["year", "Pr(>|t|)"])
}
res <- carb %>% group_by(organism, drug) %>% group_modify(function(d, k) bind_rows(
  slope(d, res_pct ~ year, "primary: 2017-2024"),
  slope(d %>% filter(year >= 2019), res_pct ~ year, "restricted: 2019-2024 (21 participating centres)"),
  slope(d %>% filter(year >= 2020), res_pct ~ year, "restricted: 2020-2024 (RC-coded tables, 20-21 RCs)"),
  slope(d, res_pct ~ year + log_total, "adjusted for log(total isolates reported)"),
  slope(d, res_pct ~ year + regional_centres_listed, "adjusted for listed regional centres"))) %>% ungroup()
write_csv(res, file.path(OUT, "table_network_sensitivity.csv"))
kp <- res %>% filter(organism %in% c("K. pneumoniae", "A. baumannii")) %>% mutate(across(c(slope, ci_lo, ci_hi), ~ round(., 2)), p = signif(p, 2))
print(as.data.frame(kp), row.names = FALSE)
cat("\nWrote output/network/table_network_sensitivity.csv\n")
