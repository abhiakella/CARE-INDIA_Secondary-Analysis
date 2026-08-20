# CARE-INDIA (Carbapenem resistance among Gram-negative ESKAPEE pathogens in India, 2017-2024)
# Builds India-vs-global genetic feature tables and a genotype-informed therapy decision matrix from the curated ESKAPEE feature workbook.
# Authors: Abhishek Akella, Anand Srinivasan  |  Dept of Pharmacology, AIIMS Bhubaneswar, India
# License: MIT (see LICENSE)

if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl", repos = "https://cloud.r-project.org")
}
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

xl_path <- "data/ESKAPEE_Genetic_Feature_Table.xlsx"

if (!file.exists(xl_path)) {
  warning(sprintf("Genetic features file %s not found. Skipping genetic modules.", xl_path))
} else {

  # Read genetic feature sheets and normalize column names
  global_features <- read_excel(xl_path, sheet = "Global Genetic Features")
  india_features  <- read_excel(xl_path, sheet = "India Genetic Features")
  global_vs_india <- read_excel(xl_path, sheet = "Global vs India Comparison")

  clean_names <- function(df) {
    names(df) <- names(df) %>%
      str_to_lower() %>%
      str_replace_all("[^a-z0-9]+", "_") %>%
      str_remove_all("^_|_$")
    df
  }

  global_features <- clean_names(global_features)
  india_features  <- clean_names(india_features)
  global_vs_india <- clean_names(global_vs_india)

  # Genotype-informed therapeutic decision matrix (India-focused)
  therapy_matrix <- tibble::tribble(
    ~organism, ~dominant_mechanism, ~phenotype_driver, ~first_line_therapy, ~salvage_therapy, ~india_context_alert,
    "K. pneumoniae", "blaNDM-1 / blaNDM-5", "MBL (Zinc-dependent)", "Aztreonam + Ceftazidime-Avibactam", "Colistin + Fosfomycin / Cefiderocol", "Ceftazidime-Avibactam alone is INEFFECTIVE in India due to NDM dominance.",
    "K. pneumoniae", "blaOXA-48 / blaOXA-181", "Oxacillinase", "Ceftazidime-Avibactam", "Colistin + Tigecycline", " blaOXA-181 is more common than blaOXA-48 in parts of India.",
    "A. baumannii", "blaOXA-23", "Class D CHDL", "Colistin/Polymyxin B + Minocycline", "High-dose Ampicillin-Sulbactam + Colistin", "Extremely high prevalence (>90%). Nearly pan-resistant.",
    "A. baumannii", "blaOXA-23 + blaNDM-1", "Dual Carbapenemase", "Colistin + Minocycline", "Cefiderocol (if available)", "Co-carriage rates up to 29% in India; negates novel beta-lactamase inhibitors.",
    "P. aeruginosa", "MexAB-OprM + OprD loss", "Efflux + Porin loss", "High-dose extended infusion Meropenem", "Colistin + Aminoglycoside", "Non-enzymatic resistance is a major driver of CRPA in India.",
    "P. aeruginosa", "blaNDM-1", "MBL", "Aztreonam + Ceftazidime-Avibactam", "Colistin", "blaNDM-1 is more prevalent than blaVIM in India, unlike global trends.",
    "Enterobacter spp.", "derepressed AmpC + blaCTX-M", "AmpC + ESBL", "Cefepime (if MIC low) or Carbapenem", "Tigecycline / Fosfomycin", "High empiric carbapenem use due to ESBLs drives CRE emergence.",
    "Enterobacter spp.", "blaNDM-1", "MBL", "Aztreonam + Ceftazidime-Avibactam", "Colistin", "blaNDM dominant in India vs blaKPC globally.",
    "E. coli", "blaCTX-M-15", "ESBL", "Ertapenem / Meropenem / Nitrofurantoin (UTI)", "Fosfomycin", "Globally dominant ESBL. Resistance rate >60% in India.",
    "E. coli", "blaNDM-5 (ST167/ST410)", "MBL", "Aztreonam + Ceftazidime-Avibactam", "Colistin based combos", "Rapidly rising. Always test for susceptibility to novel combinations."
  )

  # Organism-specific molecular profiles (India)
  molecular_profile_india <- india_features %>%
    group_by(pathogen) %>%
    summarise(
      key_genes = paste(unique(key_amr_genes_india), collapse = " | "),
      dominant_alleles = paste(unique(dominant_alleles_in_india), collapse = " | "),
      key_sts = paste(unique(na.omit(dominant_indian_sts)), collapse = " | ")
    )

  # Save output tables
  dir.create("output/genetic", showWarnings = FALSE)
  write.csv(global_features, "output/genetic/table_genetic_features_global.csv", row.names = FALSE)
  write.csv(india_features, "output/genetic/table_genetic_features_india.csv", row.names = FALSE)
  write.csv(global_vs_india, "output/genetic/table_global_vs_india_comparison.csv", row.names = FALSE)
  write.csv(therapy_matrix, "output/genetic/table_genotype_therapy_decision_matrix.csv", row.names = FALSE)

  cat("Genetic features parsed and decision matrices generated.\n")
}
