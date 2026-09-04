# CARE-INDIA (Carbapenem resistance among Gram-negative ESKAPEE pathogens in India, 2017-2024)
# Loads the ICMR-AMRSN antibiogram from CSV and builds shared long/wide surveillance objects.
# Authors: Abhishek Akella, Anand Srinivasan  |  Dept of Pharmacology, AIIMS Bhubaneswar, India
# License: MIT (see LICENSE)

suppressPackageStartupMessages(library(tidyverse))

years <- 2017:2024

dd <- Sys.getenv("AMRSN_DATA_DIR", "data")

full_long <- read.csv(file.path(dd, "data_full_long.csv"), stringsAsFactors = FALSE) %>%
  as_tibble()
if (!"resistant" %in% names(full_long)) {
  full_long <- full_long %>% mutate(resistant = tested - susceptible)
}
if (!"res_pct" %in% names(full_long)) {
  full_long <- full_long %>% mutate(res_pct = 100 - susc_pct)
}
full_long <- full_long %>% arrange(organism, drug, year)

totals_path <- file.path(dd, "organism_totals.csv")
if (file.exists(totals_path)) {
  totals <- read.csv(totals_path, stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    arrange(organism, year)
} else {
  warning("organism_totals.csv not found; total_n set to NA")
  totals <- full_long %>%
    distinct(organism, year) %>%
    mutate(total_n = NA_real_) %>%
    arrange(organism, year)
}

build_carb_wide <- function(drug_data, total_data) {
  org <- unique(drug_data$organism)
  imi <- drug_data %>% filter(drug == "Imipenem")
  mer <- drug_data %>% filter(drug == "Meropenem")
  tibble(
    organism     = org,
    year         = years,
    total_n      = total_data$total_n,
    imi_tested   = imi$tested,
    imi_susc     = imi$susceptible,
    imi_susc_pct = imi$susc_pct,
    mer_tested   = mer$tested,
    mer_susc     = mer$susceptible,
    mer_susc_pct = mer$susc_pct
  ) %>%
    mutate(
      imi_res     = imi_tested - imi_susc,
      imi_res_pct = 100 - imi_susc_pct,
      mer_res     = mer_tested - mer_susc,
      mer_res_pct = 100 - mer_susc_pct
    )
}

amrsn <- bind_rows(
  build_carb_wide(full_long %>% filter(organism == "K. pneumoniae"),
                  totals %>% filter(organism == "K. pneumoniae")),
  build_carb_wide(full_long %>% filter(organism == "A. baumannii"),
                  totals %>% filter(organism == "A. baumannii")),
  build_carb_wide(full_long %>% filter(organism == "P. aeruginosa"),
                  totals %>% filter(organism == "P. aeruginosa")),
  build_carb_wide(full_long %>% filter(organism == "Enterobacter spp."),
                  totals %>% filter(organism == "Enterobacter spp."))
)

amrsn_ee <- bind_rows(
  amrsn,
  build_carb_wide(full_long %>% filter(organism == "E. coli"),
                  totals %>% filter(organism == "E. coli"))
)

amrsn_long <- amrsn %>%
  pivot_longer(
    cols = c(imi_res_pct, mer_res_pct),
    names_to = "antibiotic", values_to = "resistance_pct"
  ) %>%
  mutate(antibiotic = recode(antibiotic,
    "imi_res_pct" = "Imipenem", "mer_res_pct" = "Meropenem"
  ))

amrsn_ee_long <- amrsn_ee %>%
  pivot_longer(
    cols = c(imi_res_pct, mer_res_pct),
    names_to = "antibiotic", values_to = "resistance_pct"
  ) %>%
  mutate(antibiotic = recode(antibiotic,
    "imi_res_pct" = "Imipenem", "mer_res_pct" = "Meropenem"
  ))

cat(sprintf("Full dataset: %s organism-drug-years | organisms: %d\n",
            format(nrow(full_long), big.mark = ","),
            length(unique(full_long$organism))))
cat(sprintf("Total isolates (all organisms): %s\n",
            format(sum(totals$total_n, na.rm = TRUE), big.mark = ",")))
cat(sprintf("ESKAPE carbapenem wide: %d rows | ESKAPEE: %d rows\n",
            nrow(amrsn), nrow(amrsn_ee)))
