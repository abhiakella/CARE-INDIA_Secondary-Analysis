# CARE-INDIA (Carbapenem resistance among Gram-negative ESKAPEE pathogens in India, 2017-2024)
# 5-fold stratified cross-validation of the empiric-failure classifiers (random forest; full and no-lag logistic regression).
# Authors: Abhishek Akella, Anand Srinivasan  |  Dept of Pharmacology, AIIMS Bhubaneswar, India
# License: MIT (see LICENSE)

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from sklearn.compose import ColumnTransformer
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import brier_score_loss, roc_auc_score
    from sklearn.model_selection import StratifiedKFold
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import OneHotEncoder, StandardScaler
except ImportError as e:
    print("Missing dependency:", e, file=sys.stderr)
    sys.exit(1)

RANDOM_STATE     = 42
HIGH_RISK_CUTOFF = 60
N_SPLITS         = 5

GENETIC_MAP = {
    "K. pneumoniae":     {"india_cre_rate": 58.0, "has_ndm": 1, "is_oxa23_driven": 0},
    "A. baumannii":      {"india_cre_rate": 92.0, "has_ndm": 1, "is_oxa23_driven": 1},
    "P. aeruginosa":     {"india_cre_rate": 30.0, "has_ndm": 1, "is_oxa23_driven": 0},
    "Enterobacter spp.": {"india_cre_rate": 20.0, "has_ndm": 1, "is_oxa23_driven": 0},
    "E. coli":           {"india_cre_rate": 34.0, "has_ndm": 1, "is_oxa23_driven": 0},
}
NO_LAG_FEATURES = ["organism", "drug", "year", "log_tested"]
FULL_FEATURES   = ["organism", "drug", "year", "log_tested", "lag_susc_pct",
                   "india_cre_rate", "has_ndm", "is_oxa23_driven"]
MODELS = [("no_lag_LR", "lr", NO_LAG_FEATURES), ("full_LR", "lr", FULL_FEATURES), ("full_RF", "rf", FULL_FEATURES)]

def load_required(csv_path: Path) -> pd.DataFrame:
    if not csv_path.exists():
        raise SystemExit(f"CSV not found: {csv_path}")
    df = pd.read_csv(csv_path)
    miss = {"organism", "drug", "year", "tested", "susc_pct", "res_pct"} - set(df.columns)
    if miss:
        raise SystemExit(f"CSV missing columns: {miss}")
    return df

def add_series_lag(df: pd.DataFrame) -> pd.DataFrame:
    df = df.reset_index(drop=True).copy()
    prev = df[["organism", "drug", "year", "susc_pct"]].copy()
    prev["year"] = prev["year"] + 1
    prev = prev.rename(columns={"susc_pct": "lag_susc_pct"})
    df["lag_susc_pct"] = df[["organism", "drug", "year"]].merge(prev, on=["organism", "drug", "year"], how="left")["lag_susc_pct"].to_numpy()
    df["log_tested"] = np.log1p(df["tested"].astype(float))
    gen = pd.DataFrame.from_dict(GENETIC_MAP, orient="index").reset_index().rename(columns={"index": "organism"})
    return df.merge(gen, on="organism", how="left")

def make_pipe(kind: str, features: list[str]) -> Pipeline:
    cat = [c for c in ["organism", "drug"] if c in features]
    num = [c for c in ["year", "log_tested", "lag_susc_pct", "india_cre_rate"] if c in features]
    bin_ = [c for c in ["has_ndm", "is_oxa23_driven"] if c in features]
    prep = ColumnTransformer([("cat", OneHotEncoder(handle_unknown="ignore"), cat),
                              ("num", StandardScaler(), num), ("bin", "passthrough", bin_)])
    clf = (LogisticRegression(max_iter=1000, class_weight="balanced", random_state=RANDOM_STATE) if kind == "lr"
           else RandomForestClassifier(n_estimators=300, max_depth=8, min_samples_leaf=2, class_weight="balanced",
                                       random_state=RANDOM_STATE, n_jobs=-1))
    return Pipeline([("prep", prep), ("clf", clf)])

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=Path, default=Path("data/data_full_long.csv"))
    ap.add_argument("--out", type=Path, default=Path("output/ml_cv"))
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    df = add_series_lag(load_required(args.csv))
    y = (df["susc_pct"] < HIGH_RISK_CUTOFF).astype(int).to_numpy()
    skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE)
    rows: list[dict] = []
    for fold, (tr, te) in enumerate(skf.split(df, y), start=1):
        train, test = df.iloc[tr].copy(), df.iloc[te].copy()
        fill = float(train["susc_pct"].mean())
        train["lag_susc_pct"] = train["lag_susc_pct"].fillna(fill)
        test["lag_susc_pct"]  = test["lag_susc_pct"].fillna(fill)
        for name, kind, feats in MODELS:
            pipe = make_pipe(kind, feats).fit(train[feats], y[tr])
            p = pipe.predict_proba(test[feats])[:, 1]
            rows.append({"model": name, "fold": fold, "n_train": int(len(tr)), "n_test": int(len(te)),
                         "auc": float(roc_auc_score(y[te], p)), "brier": float(brier_score_loss(y[te], p))})
    folds = pd.DataFrame(rows)
    folds.to_csv(args.out / "table_classifier_cv_folds.csv", index=False)
    summ = folds.groupby("model").agg(auc_mean=("auc", "mean"), auc_sd=("auc", "std"),
                                      brier_mean=("brier", "mean"), brier_sd=("brier", "std")).reset_index()
    print(f"\n{N_SPLITS}-fold stratified CV (n={len(df)}, lag computed on the full series before folding):")
    for r in summ.itertuples():
        print(f"  {r.model:>10}: AUC {r.auc_mean:.3f} +/- {r.auc_sd:.3f}   Brier {r.brier_mean:.3f} +/- {r.brier_sd:.3f}")
    (args.out / "classifier_cv_metrics.json").write_text(json.dumps({
        "method": f"{N_SPLITS}-fold stratified CV, shuffle=True, RANDOM_STATE={RANDOM_STATE}",
        "high_risk_cutoff": HIGH_RISK_CUTOFF, "n": int(len(df)),
        "lag_feature": "previous calendar year's susceptibility, computed on the full series before folding; "
                       "first-year lags imputed with the training-fold mean only",
        "summary": summ.to_dict(orient="records"), "folds": rows}, indent=2), encoding="utf-8")
    print("Wrote:", args.out / "classifier_cv_metrics.json", "and", args.out / "table_classifier_cv_folds.csv")

if __name__ == "__main__":
    main()
