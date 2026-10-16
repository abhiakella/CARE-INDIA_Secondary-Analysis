# CARE-INDIA analysis code

Analysis code for **CARE-INDIA** -  *Carbapenem Resistance among Gram-negative ESKAPEE Pathogens in India: An Eight-Year Secondary Analysis of ICMR-AMRSN Surveillance Data (2017-2024) with Resistance Forecasting.*

This repository contains **code only**. The surveillance data are a secondary analysis of the publicly available ICMR-AMRSN annual reports and are not redistributed here (see **Data availability**).

## Authors
- Abhishek Akella -  Department of Pharmacology, AIIMS Bhubaneswar, Odisha, India
- Anand Srinivasan (corresponding; anandsrinivasan@aiimsbhubaneswar.edu.in) -  Department of Pharmacology, AIIMS Bhubaneswar, Odisha, India

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
| `02_trend_tests_extended.R` | Extended trend tests: weighted regression (primary) and Cochran-Armitage across all 47 organism-drug pairs with Holm and Benjamini-Hochberg correction; unweighted sensitivity regression; pre-/post-2020 interaction model and piecewise grid search over candidate knots 2019-2022 for the ten organism-carbapenem pairs (Supplementary Tables S1-S3, S15, S17). | (sources `amrsn_data_shared.R`) | `output/trends/*.csv` |
| `amrsn_ml_enhancements.R` | R ML modules: hierarchical/k-medoids clustering with NMF archetypes, and PELT changepoint sensitivity. | (sources `amrsn_data_shared.R`) | `output/ml/*` (clustering, changepoint/regime map) |
| `03_forecast_trend.R` | Resistance forecast 2025-2028: per organism-drug pair weighted logit regression (weights = isolates tested) with 80%/95% prediction intervals; rolling-origin backtest (origins 2021-2023) against naive persistence and an auto.arima + mixed-model ensemble comparator; Figure 5 and Supplementary Table S8. | (sources `amrsn_data_shared.R`) | `output/forecast/table_forecast_trend_2025_2028.csv`, `table_S8_forecast_carbapenems.csv`, backtest CSVs, `fig5_forecast_trend.png/.pdf`, `RUN_SUMMARY.md` |
| `05_trajectory_clustering.R` | k-medoids (PAM) trajectory clustering of organisms on per-drug resistance-slope features; slope heatmap and silhouette plot (appendix SR.1 resistance archetypes). | `data/AMRSN_GN_ESKAPEE_Susceptibility_Trends.csv` | `output/table_organism_drug_slopes.csv`, `output/table_pam_trajectory_clusters.csv`, `output/fig_trajectory_*`, `output/fig_silhouette_*` |
| `genetic_features.R` | India-vs-global genetic-feature tables and a genotype-to-therapy decision matrix from the curated feature workbook. | `data/ESKAPEE_Genetic_Feature_Table.xlsx` | `output/genetic/*.csv` |
| `04_forecast_scenarios.R` | Forecast sensitivity: 2028 projections for the ten organism-carbapenem pairs under alternative fitting windows and forms (2017-2024 logit trend; 2020-2024 refit; leave-one-year-out refits; piecewise trend with a 2019 knot) with 95% prediction intervals (Supplementary Table S23). | (sources `amrsn_data_shared.R`) | `output/forecast/*scenario*` |
| `06_revision_figures.R` | Journal-style main Figures 1-6 (TIFF 300 dpi, PNG, JPEG and vector PDF; lowercase panel tags; E. coli in Figure 1; colour-blind-safe red-yellow-blue antibiogram scale with plain grey blocks for untested pairs) and the year-by-year antibiogram supplementary figure (S13). | (sources `amrsn_data_shared.R`), `data/gene_trends.csv`, `data/fig4_panel_b_curated_literature.csv`, `data/syndrome_key_antibiogram_2024.csv`, `output/forecast/table_forecast_trend_2025_2028.csv` | `output/figures_revision/*` |
| `07_network_sensitivity.R` | Network-composition sensitivity: weighted carbapenem trends restricted to the stable-network years (2020-2024 and 2019-2024) and with annual network throughput or the number of listed regional centres as a covariate (Supplementary Table S22). | (sources `amrsn_data_shared.R`), `data/network_by_year.csv` | `output/network/table_network_sensitivity.csv` |
| `08_pelt_eskapee.R` | PELT mean-variance changepoint sensitivity for all ten organism-carbapenem pairs including E. coli (Supplementary Table S4). | (sources `amrsn_data_shared.R`) | `output/ml/table_changepoint_pelt_carbapenem_eskapee.csv` |

### python/
| Script | What it does | Reads | Writes |
|---|---|---|---|
| `consolidated_ml.py` | Supervised high-risk (susceptibility < 60%) classifier (logistic regression + random forest) with 75/25 holdout, external, and rolling-origin validation; headline metrics. | `data/data_full_long.csv`, `data/atlas_external_india.csv` (external) | `output/headline_*` (JSON, CSVs, figures) |
| `02_classifier_baselines_calibration.py` | Naive / LR / RandomForest baselines with bootstrap-CI AUC, Brier score, and calibration slope/intercept. | `data/data_full_long.csv`, `data/atlas_external_india.csv` | `output/table_baselines_*`, `output/classifier_baselines_metrics.json`, calibration/AUC figures |
| `empiric_classifier_shap.py` | Random-forest empiric-failure classifier with permutation and SHAP interpretability and external validation. | `data/data_full_long.csv` (`--csv`), `data/atlas_external_india.csv` (`--external_csv`) | `output/ml/empiric_*`, SHAP summary figures |
| `03_classifier_cv.py` | 5-fold stratified cross-validation of the same classifiers (random forest, full and no-lag logistic regression); per-fold ROC-AUC and Brier score with fold means and SDs. | `data/data_full_long.csv` | `output/ml_cv/classifier_cv_metrics.json`, `output/ml_cv/table_classifier_cv_folds.csv` |

> **Lag feature.** In all four classifier scripts the lagged susceptibility is the previous calendar year's value, computed on the complete organism-drug series *before* any train/test partitioning; first-year lags are imputed with the training-partition (or training-fold) mean only. External (ATLAS) rows take their lag from the preceding ATLAS year, or from the AMRSN training data where no ATLAS history exists.

### scripts/
| Script | What it does | Reads | Writes |
|---|---|---|---|
| `new_modules.R` | Supplementary modules: specimen-stratified meropenem trends + 2024 antibiogram, organism isolation share, carbapenemase gene trends, K. pneumoniae genotype-phenotype. | `data/specimen_stratified.csv`, `data/isolation_trends.csv`, `data/gene_trends.csv`, `data/data_full_long.csv` | `output/specimen/*`, `output/burden/*`, `output/genes/*` |
| `bespoke_figs.R` | Bespoke manuscript figures: 2024 antibiogram heatmap, therapeutic-tier map, carbapenemase gene-prevalence panels. | `data/data_full_long.csv`, `data/gene_trends.csv` | `output/fig2_antibiogram_heatmap.png`, `output/fig3_tier_map.png`, `output/fig4_india_carbapenemase.png` |

> **Genetic constants are annotation only.** The Python classifiers carry a small fixed `GENETIC_MAP` (per-organism literature-derived constants: an India CRE-rate plus binary NDM / OXA-23 indicators). Since revision 1 these columns are merged for annotation only and are **not** model features, because they are constant within organism and therefore collinear with organism identity; the classifiers use organism, drug, year, log(isolates tested) and lagged susceptibility. The surveillance data themselves are read from CSV.

## Requirements
- **R** >= 4.6.0 with: `tidyverse`, `scales`, `cluster`, `forecast`, `changepoint`, `lme4`, `segmented`, `readxl` (and, optionally, `NMF` + `BiocManager` for the archetype module).
- **Python** 3.12 with the packages in `requirements.txt` (`numpy`, `pandas`, `scikit-learn`, `shap`, `matplotlib`, `scipy`, `openpyxl`). Install with `pip install -r requirements.txt`.

## How to run
1. Place the input files (see **Data availability**) in `data/`, and create an `output/` directory.
2. **R:** run `R/amrsn_analysis.R`, `R/02_trend_tests_extended.R`, `R/amrsn_ml_enhancements.R`, `R/03_forecast_trend.R`, `R/04_forecast_scenarios.R`, `R/05_trajectory_clustering.R`, `R/genetic_features.R`, `R/07_network_sensitivity.R`, `R/08_pelt_eskapee.R`, `R/06_revision_figures.R` (all `source()` `R/amrsn_data_shared.R`; run `03_forecast_trend.R` before `06_revision_figures.R`), and the `scripts/` modules, from the repository root.
3. **Python:** `pip install -r requirements.txt`, then run each classifier, e.g. `python python/empiric_classifier_shap.py --csv data/data_full_long.csv --external_csv data/atlas_external_india.csv`.

Scripts read from `data/` and write to `output/` using paths relative to the repository root; several accept explicit path arguments. Randomised steps are seeded (`set.seed(42)` / `random_state=42`) for reproducibility.

## Data availability
The analysis is a secondary analysis of the publicly available **ICMR-AMRSN annual surveillance reports** (https://iamrsn.icmr.org.in). The derived analysis dataset used by these scripts is deposited separately (DOI assigned on acceptance); external validation uses **Pfizer ATLAS** India data obtained through the **Vivli AMR Register** (https://amr.vivli.org) and is not redistributed here. Expected files in `data/`:

- `data_full_long.csv` -  long-format antibiogram (organism, year, drug, tested, susceptible, susc_pct); primary input.
- `organism_totals.csv` -  organism-year total isolate counts (organism, year, total_n), used by the loader.
- `AMRSN_GN_ESKAPEE_Susceptibility_Trends.csv` -  full-panel susceptibility trends (trajectory clustering).
- `gene_trends.csv` -  carbapenemase gene prevalence by organism and year.
- `specimen_stratified.csv`, `isolation_trends.csv` -  specimen-stratified antibiogram and isolation share.
- `ESKAPEE_Genetic_Feature_Table.xlsx` -  curated genetic-feature workbook.
- `atlas_external_india.csv` -  external validation set (the Pfizer ATLAS India export).
- `network_by_year.csv` -  annual network composition (listed regional centres, total isolates reported), used by the network sensitivity.
- `fig4_panel_b_curated_literature.csv`, `syndrome_key_antibiogram_2024.csv` -  curated literature values for Figure 4b and the 2024 syndrome antibiogram for Figure 6.

## License
Code released under the **MIT License** (see `LICENSE`). Copyright (c) 2026 Abhishek Akella and Anand Srinivasan. If you also use the deposited dataset, note it is licensed separately (CC-BY-4.0).

## Citation
If you use this code, please cite the CARE-INDIA article (citation to be added on publication).
