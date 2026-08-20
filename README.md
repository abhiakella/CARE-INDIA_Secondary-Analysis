# CARE-INDIA analysis code

Analysis code for **CARE-INDIA** — *Carbapenem Resistance among Gram-negative ESKAPEE Pathogens in India: An Eight-Year Secondary Analysis of ICMR-AMRSN Surveillance Data (2017-2024) with Resistance Forecasting.*

This repository contains **code only**. The surveillance data are a secondary analysis of the publicly available ICMR-AMRSN annual reports and are not redistributed here (see **Data availability**).

## Authors
- Abhishek Akella — Department of Pharmacology, AIIMS Bhubaneswar, Odisha, India
- Anand Srinivasan (corresponding; anandsrinivasan@aiimsbhubaneswar.edu.in) — Department of Pharmacology, AIIMS Bhubaneswar, Odisha, India

## Repository layout
```
R/         Trend analysis, forecast/changepoint ML modules, clustering/archetypes, genetic-feature tables (and the shared data loader)
python/    Supervised empiric-failure classifiers, calibration/baselines
scripts/   Supplementary figure/table modules and bespoke manuscript figures
data/      NOT included - place the input files here (see Data availability)
output/    Created by the scripts (figures, tables, metrics)
LICENSE    MIT License
requirements.txt   Python dependencies
```

## Scripts

### R/
| Script | What it does | Reads | Writes |
|---|---|---|---|
| `amrsn_data_shared.R` | Shared data loader (Module 0): reads the antibiogram and builds the long/wide objects used by the other R scripts. `source()`d by the two scripts below. | `data/data_full_long.csv`, `data/organism_totals.csv` | (in-memory objects) |
| `amrsn_analysis.R` | Primary trend analysis: Cochran-Armitage trend tests, weighted linear-regression slopes/CIs, segmented-regression/Davies breakpoint sensitivity, main-text resistance figures, isolate-burden figure. | (sources `amrsn_data_shared.R`) | `output/fig1-4*`, `output/table_cochran_armitage_results.csv`, `output/table_weighted_lm_slopes.csv`, `output/table_davies_segmented_breakpoint.csv`, `output/table_annualised_pp_change.csv` |
| `amrsn_ml_enhancements.R` | R ML modules: auto.arima + linear-mixed-model forecast ensemble, hierarchical/k-medoids clustering with NMF archetypes, and PELT changepoint sensitivity. | (sources `amrsn_data_shared.R`) | `output/ml/*` (forecast table+figs, clustering, changepoint/regime map) |
| `forecast_calibration.R` | Forecast calibration: rolling-origin backtest + bias-correction + logit split-conformal + Jeffreys floor; produces the calibrated 2025-2028 projection (headline K. pneumoniae meropenem ~74% by 2028). Run **after** `amrsn_ml_enhancements.R`. | (sources `amrsn_data_shared.R`; reads `output/ml/table_forecast_resistance_2025_2028.csv`) | `output/ml/table_forecast_resistance_2025_2028_tier1logit.csv`, backtest/metrics CSVs, `RUN_SUMMARY.md` |
| `05_trajectory_clustering.R` | k-medoids (PAM) trajectory clustering of organisms on per-drug resistance-slope features; slope heatmap and silhouette plot (appendix SR.1 resistance archetypes). | `data/AMRSN_GN_ESKAPEE_Susceptibility_Trends.csv` | `output/table_organism_drug_slopes.csv`, `output/table_pam_trajectory_clusters.csv`, `output/fig_trajectory_*`, `output/fig_silhouette_*` |
| `genetic_features.R` | India-vs-global genetic-feature tables and a genotype-to-therapy decision matrix from the curated feature workbook. | `data/ESKAPEE_Genetic_Feature_Table.xlsx` | `output/genetic/*.csv` |

### python/
| Script | What it does | Reads | Writes |
|---|---|---|---|
| `consolidated_ml.py` | Supervised high-risk (susceptibility < 60%) classifier (logistic regression + random forest) with 75/25 holdout, external, and rolling-origin validation; headline metrics. | `data/data_full_long.csv`, `data/atlas_external_india.csv` (external) | `output/headline_*` (JSON, CSVs, figures) |
| `02_classifier_baselines_calibration.py` | Naive / LR / RandomForest baselines with bootstrap-CI AUC, Brier score, and calibration slope/intercept. | `data/data_full_long.csv`, `data/atlas_external_india.csv` | `output/table_baselines_*`, `output/classifier_baselines_metrics.json`, calibration/AUC figures |
| `empiric_classifier_shap.py` | Random-forest empiric-failure classifier with permutation and SHAP interpretability and external validation. | `data/data_full_long.csv` (`--csv`), `data/atlas_external_india.csv` (`--external_csv`) | `output/ml/empiric_*`, SHAP summary figures |

### scripts/
| Script | What it does | Reads | Writes |
|---|---|---|---|
| `new_modules.R` | Supplementary modules: specimen-stratified meropenem trends + 2024 antibiogram, organism isolation share, carbapenemase gene trends, K. pneumoniae genotype-phenotype. | `data/specimen_stratified.csv`, `data/isolation_trends.csv`, `data/gene_trends.csv`, `data/data_full_long.csv` | `output/specimen/*`, `output/burden/*`, `output/genes/*` |
| `bespoke_figs.R` | Bespoke manuscript figures: 2024 antibiogram heatmap, therapeutic-tier map, carbapenemase gene-prevalence panels. | `data/data_full_long.csv`, `data/gene_trends.csv` | `output/fig2_antibiogram_heatmap.png`, `output/fig3_tier_map.png`, `output/fig4_india_carbapenemase.png` |

> **Model parameters, not data.** The three Python classifiers contain a small fixed `GENETIC_MAP` — per-organism literature-derived feature constants (an India CRE-rate plus binary NDM / OXA-23 indicators for each of the five organisms). These are model input **parameters**, not the surveillance dataset (which is read from CSV), and are kept inline so the classifiers reproduce without an extra file.

## Requirements
- **R** >= 4.6.0 with: `tidyverse`, `scales`, `cluster`, `forecast`, `changepoint`, `lme4`, `segmented`, `readxl` (and, optionally, `NMF` + `BiocManager` for the archetype module).
- **Python** 3.12 with the packages in `requirements.txt` (`numpy`, `pandas`, `scikit-learn`, `shap`, `matplotlib`, `scipy`, `openpyxl`). Install with `pip install -r requirements.txt`.

## How to run
1. Place the input files (see **Data availability**) in `data/`, and create an `output/` directory.
2. **R:** run `R/amrsn_analysis.R` and `R/amrsn_ml_enhancements.R` (both `source()` `R/amrsn_data_shared.R`), then `R/forecast_calibration.R` (uses the forecast CSV written by `amrsn_ml_enhancements.R`), then `R/05_trajectory_clustering.R`, `R/genetic_features.R`, and the `scripts/` modules, from the repository root.
3. **Python:** `pip install -r requirements.txt`, then run each classifier, e.g. `python python/empiric_classifier_shap.py --csv data/data_full_long.csv --external_csv data/atlas_external_india.csv`.

Scripts read from `data/` and write to `output/` using paths relative to the repository root; several accept explicit path arguments. Randomised steps are seeded (`set.seed(42)` / `random_state=42`) for reproducibility.

## Data availability
The analysis is a secondary analysis of the publicly available **ICMR-AMRSN annual surveillance reports** (https://iamrsn.icmr.org.in). The derived analysis dataset used by these scripts is deposited separately (DOI assigned on acceptance); external validation uses **Pfizer ATLAS** India data obtained through the **Vivli AMR Register** (https://amr.vivli.org) and is not redistributed here. Expected files in `data/`:

- `data_full_long.csv` — long-format antibiogram (organism, year, drug, tested, susceptible, susc_pct); primary input.
- `organism_totals.csv` — organism-year total isolate counts (organism, year, total_n), used by the loader.
- `AMRSN_GN_ESKAPEE_Susceptibility_Trends.csv` — full-panel susceptibility trends (trajectory clustering).
- `gene_trends.csv` — carbapenemase gene prevalence by organism and year.
- `specimen_stratified.csv`, `isolation_trends.csv` — specimen-stratified antibiogram and isolation share.
- `ESKAPEE_Genetic_Feature_Table.xlsx` — curated genetic-feature workbook.
- `atlas_external_india.csv` — external validation set (the Pfizer ATLAS India export).

## License
Code released under the **MIT License** (see `LICENSE`). Copyright (c) 2026 Abhishek Akella and Anand Srinivasan. If you also use the deposited dataset, note it is licensed separately (CC-BY-4.0).

## Citation
If you use this code, please cite the CARE-INDIA article (citation to be added on publication).
