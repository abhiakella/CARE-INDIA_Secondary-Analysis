# CARE-INDIA (Carbapenem resistance among Gram-negative ESKAPEE pathogens in India, 2017-2024)
# Empiric-failure random-forest classifier with SHAP interpretability and external validation.
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
    import matplotlib.pyplot as plt
    import shap
    from sklearn.compose import ColumnTransformer
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.inspection import permutation_importance
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import classification_report, roc_auc_score
    from sklearn.model_selection import train_test_split
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import OneHotEncoder, StandardScaler
except ImportError as e:
    print("Missing dependency:", e, file=sys.stderr)
    print("Install: pip install -r python/requirements.txt", file=sys.stderr)
    sys.exit(1)

RANDOM_STATE = 42

GENETIC_MAP = {
    "K. pneumoniae":     {"india_cre_rate": 58.0, "has_ndm": 1, "is_oxa23_driven": 0},
    "A. baumannii":      {"india_cre_rate": 92.0, "has_ndm": 1, "is_oxa23_driven": 1},
    "P. aeruginosa":     {"india_cre_rate": 30.0, "has_ndm": 1, "is_oxa23_driven": 0},
    "Enterobacter spp.": {"india_cre_rate": 20.0, "has_ndm": 1, "is_oxa23_driven": 0},
    "E. coli":           {"india_cre_rate": 34.0, "has_ndm": 1, "is_oxa23_driven": 0},
}

FEATURE_COLS = [
    "organism", "drug", "year", "log_tested", "lag_susc_pct",
    "india_cre_rate", "has_ndm", "is_oxa23_driven",
]
CAT_COLS = ["organism", "drug"]
NUM_COLS = ["year", "log_tested", "lag_susc_pct", "india_cre_rate"]
BIN_COLS = ["has_ndm", "is_oxa23_driven"]

def series_lag(df: pd.DataFrame) -> np.ndarray:
    prev = df[["organism", "drug", "year", "susc_pct"]].copy()
    prev["year"] = prev["year"] + 1
    prev = prev.rename(columns={"susc_pct": "lag_susc_pct"})
    return df[["organism", "drug", "year"]].merge(prev, on=["organism", "drug", "year"], how="left")["lag_susc_pct"].to_numpy()

def add_series_lag(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["lag_susc_pct"] = series_lag(df)
    return df

def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]

def load_amrsn(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    need = {"organism", "drug", "year", "tested", "susc_pct", "res_pct"}
    miss = need - set(df.columns)
    if miss:
        raise SystemExit(f"AMRSN CSV missing columns: {miss}")
    return df

def load_external(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    need = {"organism", "drug", "year", "tested", "susc_pct", "res_pct"}
    miss = need - set(df.columns)
    if miss:
        raise SystemExit(f"External CSV missing columns: {miss}")
    return df

def add_lag_and_genetic(df: pd.DataFrame, lag_fill: float, lag_lookup: pd.DataFrame | None = None) -> pd.DataFrame:
    df = df.copy()
    df["log_tested"] = np.log1p(df["tested"].astype(float))
    df.sort_values(["organism", "drug", "year"], inplace=True)
    if "lag_susc_pct" not in df.columns:
        df["lag_susc_pct"] = series_lag(df)

    if lag_lookup is not None:
        need_lag = df["lag_susc_pct"].isna()
        if need_lag.any():
            lag_lookup = lag_lookup.copy()
            lag_lookup["join_year"] = lag_lookup["year"] + 1
            merged = df.loc[need_lag].merge(
                lag_lookup[["organism", "drug", "join_year", "susc_pct"]]
                    .rename(columns={"susc_pct": "lag_from_train", "join_year": "year"}),
                on=["organism", "drug", "year"], how="left",
            )
            df.loc[need_lag, "lag_susc_pct"] = merged["lag_from_train"].values

    df["lag_susc_pct"] = df["lag_susc_pct"].fillna(lag_fill)

    gen_df = pd.DataFrame.from_dict(GENETIC_MAP, orient="index").reset_index()
    gen_df.rename(columns={"index": "organism"}, inplace=True)
    df = df.merge(gen_df, on="organism", how="left")
    return df

def make_preprocessor() -> ColumnTransformer:
    return ColumnTransformer([
        ("cat", OneHotEncoder(handle_unknown="ignore"), CAT_COLS),
        ("num", StandardScaler(),                       NUM_COLS),
        ("bin", "passthrough",                          BIN_COLS),
    ])

def make_rf() -> RandomForestClassifier:
    return RandomForestClassifier(
        n_estimators=300,
        max_depth=8,
        min_samples_leaf=2,
        class_weight="balanced",
        random_state=RANDOM_STATE,
        n_jobs=-1,
    )

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv",          type=Path,
                        default=repo_root() / "data" / "data_full_long.csv",
                        help="Path to AMRSN data_full_long.csv")
    parser.add_argument("--external_csv", type=Path,
                        default=repo_root() / "data" / "atlas_external_india.csv",
                        help="Path to external (non-AMRSN) test CSV")
    parser.add_argument("--out",          type=Path, default=repo_root() / "output" / "ml")
    args = parser.parse_args()
    out: Path = args.out
    out.mkdir(parents=True, exist_ok=True)

    if not args.csv.exists():
        print(
            f"Input not found: {args.csv}\n"
            "Place data_full_long.csv in the data/ directory (see README, Data availability),\n"
            "or pass --csv to point at it.",
            file=sys.stderr,
        )
        sys.exit(1)

    df_raw = add_series_lag(load_amrsn(args.csv))

    y_all = (df_raw["susc_pct"] < 60).astype(int)
    train_df, test_df = train_test_split(
        df_raw, test_size=0.25, stratify=y_all, random_state=RANDOM_STATE
    )

    train_lag_fill = float(train_df["susc_pct"].mean())
    train_feat = add_lag_and_genetic(train_df, lag_fill=train_lag_fill)
    test_feat  = add_lag_and_genetic(test_df,  lag_fill=train_lag_fill)

    X_train, y_train = train_feat[FEATURE_COLS], (train_feat["susc_pct"] < 60).astype(int)
    X_test,  y_test  = test_feat[FEATURE_COLS],  (test_feat["susc_pct"]  < 60).astype(int)

    pipe = Pipeline([("prep", make_preprocessor()), ("rf", make_rf())])
    pipe.fit(X_train, y_train)
    proba = pipe.predict_proba(X_test)[:, 1]
    pred  = (proba >= 0.5).astype(int)
    holdout_auc = float(roc_auc_score(y_test, proba))

    lr_pipe = Pipeline([
        ("prep", make_preprocessor()),
        ("lr",   LogisticRegression(max_iter=1000, class_weight="balanced",
                                    random_state=RANDOM_STATE)),
    ])
    lr_pipe.fit(X_train, y_train)
    lr_auc = float(roc_auc_score(y_test, lr_pipe.predict_proba(X_test)[:, 1]))

    external_metrics = {}
    if args.external_csv.exists():
        ext_raw = load_external(args.external_csv)
        ext_feat = add_lag_and_genetic(
            ext_raw, lag_fill=train_lag_fill, lag_lookup=train_df,
        )
        X_ext = ext_feat[FEATURE_COLS]
        y_ext = (ext_feat["susc_pct"] < 60).astype(int)
        proba_ext = pipe.predict_proba(X_ext)[:, 1]

        external_metrics = {
            "n_external":           int(len(y_ext)),
            "external_prevalence":  float(y_ext.mean()),
        }
        if y_ext.nunique() > 1:
            external_metrics["external_roc_auc"] = float(roc_auc_score(y_ext, proba_ext))

        ext_scored = ext_feat.assign(
            empiric_failure_prob=proba_ext,
            high_risk_label=y_ext.values,
            predicted_label=(proba_ext >= 0.5).astype(int),
        )
        ext_scored.to_csv(out / "table_external_validation_predictions.csv", index=False)
        print(f"\nExternal validation: n={len(y_ext)} from "
              f"{ext_raw['source_first_author'].nunique() if 'source_first_author' in ext_raw else '?'} papers")
        if "external_roc_auc" in external_metrics:
            print(f"  External ROC-AUC: {external_metrics['external_roc_auc']:.3f}")
        else:
            print("  Single class in external set — AUC not computed; "
                  "see table_external_validation_predictions.csv.")
    else:
        print(f"External CSV not found at {args.external_csv} — skipping external validation.")

    metrics = {
        "holdout_roc_auc":              holdout_auc,
        "holdout_logreg_baseline":      lr_auc,
        "n_train":                      int(len(y_train)),
        "n_test":                       int(len(y_test)),
        "prevalence_high_risk_train":   float(y_train.mean()),
        "prevalence_high_risk_test":    float(y_test.mean()),
        **external_metrics,
    }
    (out / "empiric_classifier_metrics.json").write_text(
        json.dumps(metrics, indent=2), encoding="utf-8"
    )
    (out / "empiric_classifier_report.txt").write_text(
        classification_report(y_test, pred, digits=3), encoding="utf-8"
    )

    perm = permutation_importance(
        pipe, X_test, y_test, n_repeats=20, random_state=RANDOM_STATE, n_jobs=-1
    )
    pd.DataFrame({
        "feature":         FEATURE_COLS,
        "importance_mean": perm.importances_mean,
        "importance_sd":   perm.importances_std,
    }).sort_values("importance_mean", ascending=False).to_csv(
        out / "table_permutation_importance.csv", index=False
    )

    prep = pipe.named_steps["prep"]
    rf   = pipe.named_steps["rf"]

    def _dense(X):
        if hasattr(X, "toarray"):
            return X.toarray().astype(np.float64)
        return np.asarray(X, dtype=np.float64)

    Xtr_enc = _dense(prep.transform(X_train))
    Xte_enc = _dense(prep.transform(X_test))

    explainer = shap.TreeExplainer(rf)
    bg = shap.sample(Xtr_enc, min(100, Xtr_enc.shape[0]), random_state=RANDOM_STATE)
    ex = shap.sample(Xte_enc, min(150, Xte_enc.shape[0]), random_state=RANDOM_STATE + 1)
    sv = explainer.shap_values(ex, check_additivity=False)
    if isinstance(sv, list):
        sv = sv[1]
    elif getattr(sv, "ndim", 2) == 3:
        sv = sv[:, :, 1]

    feature_names = prep.get_feature_names_out()
    shap.summary_plot(sv, ex, feature_names=feature_names, show=False, max_display=20)
    plt.tight_layout()
    plt.savefig(out / "fig_shap_summary_empiric_risk.png", dpi=200)
    plt.savefig(out / "fig_shap_summary_empiric_risk_light.png", dpi=200,
                facecolor="white", edgecolor="white")
    plt.close()

    test_feat.assign(
        empiric_failure_prob=proba,
        high_risk_label=y_test.values,
    ).to_csv(out / "table_empiric_risk_scores_holdout.csv", index=False)

    print("\nWrote:", out / "empiric_classifier_metrics.json")
    print("      ", out / "fig_shap_summary_empiric_risk.png")
    print("      ", out / "table_empiric_risk_scores_holdout.csv")
    print("      ", out / "table_permutation_importance.csv")
    if external_metrics:
        print("      ", out / "table_external_validation_predictions.csv")

if __name__ == "__main__":
    main()
