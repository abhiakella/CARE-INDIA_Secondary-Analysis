# CARE-INDIA analysis code

Code for *Carbapenem Resistance among Gram-negative ESKAPEE Pathogens in India (CARE-INDIA): An Eight-Year Secondary Analysis of ICMR-AMRSN Surveillance Data (2017-2024) with Resistance Forecasting*.

Abhishek Akella, Anand Srinivasan (corresponding: anandsrinivasan@aiimsbhubaneswar.edu.in). Department of Pharmacology, AIIMS Bhubaneswar, Odisha, India.

Code only. Data are not redistributed: the analysis is a secondary analysis of the publicly available ICMR-AMRSN annual reports (https://iamrsn.icmr.org.in); the derived analysis dataset will be deposited separately with a DOI on acceptance; external validation data (Pfizer ATLAS, India) are available to qualified researchers through the Vivli AMR Register (https://amr.vivli.org).

## Layout

| Path | Produces |
|---|---|
| `R/amrsn_data_shared.R` | shared data loader (sourced by the R scripts) |
| `R/amrsn_analysis.R` | trend tests, weighted-regression slopes, Davies breakpoint test, main trend figures |
| `R/amrsn_ml_enhancements.R` | ARIMA + mixed-model forecast ensemble, clustering/NMF archetypes, PELT changepoints |
| `R/forecast_calibration.R` | rolling-origin backtest and calibrated 2025-2028 projection (run after `amrsn_ml_enhancements.R`) |
| `R/05_trajectory_clustering.R` | k-medoids trajectory clustering |
| `R/genetic_features.R` | genetic-feature tables and genotype-therapy matrix |
| `python/consolidated_ml.py` | empiric-failure classifier: holdout, rolling-origin and external validation |
| `python/02_classifier_baselines_calibration.py` | baseline comparators, calibration, paired bootstrap |
| `python/03_classifier_cv.py` | 5-fold cross-validation |
| `python/empiric_classifier_shap.py` | random-forest classifier with permutation importance and SHAP |
| `scripts/new_modules.R` | specimen-stratified, isolation-share and gene-trend modules |
| `scripts/bespoke_figs.R` | antibiogram heatmap, tier map, carbapenemase panels |

## Running

R >= 4.6.0 (`tidyverse`, `scales`, `cluster`, `forecast`, `changepoint`, `lme4`, `segmented`, `readxl`; optionally `NMF`) and Python 3.12 (`pip install -r requirements.txt`). Place the input files in `data/` and run each script from the repository root; outputs are written to `output/`. Python scripts accept `--csv`, `--external_csv` and `--out`. Random seeds are fixed (42).

Expected inputs in `data/`: `data_full_long.csv`, `organism_totals.csv`, `AMRSN_GN_ESKAPEE_Susceptibility_Trends.csv`, `gene_trends.csv`, `specimen_stratified.csv`, `isolation_trends.csv`, `ESKAPEE_Genetic_Feature_Table.xlsx`, `atlas_external_india.csv`.

## License and citation

MIT License (see `LICENSE`). Copyright (c) 2026 Abhishek Akella and Anand Srinivasan. Please cite the CARE-INDIA article (citation to be added on publication). Questions about the methods or reuse: contact the corresponding author.
