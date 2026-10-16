# CARE-INDIA (Carbapenem resistance among Gram-negative ESKAPEE pathogens in India, 2017-2024)
# PELT mean-variance changepoint sensitivity for all ten organism-carbapenem pairs (ESKAPEE, including E. coli),
# matching Module 10 of amrsn_ml_enhancements.R but run on amrsn_ee so that Supplementary Table S4 covers all five organisms.
# Authors: Abhishek Akella, Anand Srinivasan  |  Dept of Pharmacology, AIIMS Bhubaneswar, India
# License: MIT (see LICENSE)

suppressPackageStartupMessages(library(tidyverse))
root <- getwd()
if (!file.exists(file.path(root, "R", "amrsn_data_shared.R"))) stop("Run from the repository root; R/amrsn_data_shared.R not found.")
source(file.path(root, "R", "amrsn_data_shared.R"))
if (!requireNamespace("changepoint", quietly = TRUE)) install.packages("changepoint", repos = "https://cloud.r-project.org", quiet = TRUE)
library(changepoint)
OUT <- file.path(root, "output", "ml"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

organisms <- c("K. pneumoniae", "A. baumannii", "P. aeruginosa", "Enterobacter spp.", "E. coli")
rows <- list()
for (org in organisms) {
  d <- amrsn_ee %>% filter(organism == org) %>% arrange(year)
  for (abx in c("Meropenem", "Imipenem")) {
    x <- if (abx == "Meropenem") d$mer_res_pct / 100 else d$imi_res_pct / 100
    cp <- tryCatch(changepoint::cpt.meanvar(x, method = "PELT"), error = function(e) NULL)
    if (is.null(cp)) { rows[[length(rows) + 1]] <- tibble(organism = org, antibiotic = abx, n_cpt = NA_integer_, cpt_years = NA_character_); next }
    ix <- changepoint::cpts(cp)
    rows[[length(rows) + 1]] <- tibble(organism = org, antibiotic = abx, n_cpt = length(ix),
                                       cpt_years = if (length(ix)) paste(d$year[ix], collapse = ";") else "")
  }
}
tbl <- bind_rows(rows)
print(as.data.frame(tbl))
write_csv(tbl, file.path(OUT, "table_changepoint_pelt_carbapenem_eskapee.csv"))
cat("Wrote output/ml/table_changepoint_pelt_carbapenem_eskapee.csv\n")
