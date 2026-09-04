# CARE-INDIA (Carbapenem resistance among Gram-negative ESKAPEE pathogens in India, 2017-2024)
# Supervised high-risk classifier (logistic regression + random forest): 75/25 holdout, external, and rolling-origin validation
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
    from sklearn.compose import ColumnTransformer
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import brier_score_loss, roc_auc_score
    from sklearn.model_selection import train_test_split
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import OneHotEncoder, StandardScaler
except ImportError as e:
    print("Missing dependency:", e, file=sys.stderr)
    sys.exit(1)

RANDOM_STATE        = 42
HIGH_RISK_CUTOFF    = 60
N_BOOTSTRAP         = 1000
TEST_YEARS          = (2022, 2023, 2024)
THIN_FOLD_THRESHOLD = 200

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

def series_lag(df: pd.DataFrame) -> np.ndarray:
    prev = df[["organism", "drug", "year", "susc_pct"]].copy()
    prev["year"] = prev["year"] + 1
    prev = prev.rename(columns={"susc_pct": "lag_susc_pct"})
    return df[["organism", "drug", "year"]].merge(prev, on=["organism", "drug", "year"], how="left")["lag_susc_pct"].to_numpy()

def add_series_lag(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["lag_susc_pct"] = series_lag(df)
    return df

def load_required(csv_path: Path, label: str) -> pd.DataFrame:
    if not csv_path.exists():
        raise SystemExit(f"{label} CSV not found: {csv_path}")
    df = pd.read_csv(csv_path)
    need = {"organism", "drug", "year", "tested", "susc_pct", "res_pct"}
    miss = need - set(df.columns)
    if miss:
        raise SystemExit(f"{label} CSV missing columns: {miss}")
    return df

def add_features(
    df: pd.DataFrame,
    lag_fill: float,
    lag_lookup: pd.DataFrame | None = None,
) -> pd.DataFrame:
    df = df.copy()
    df["log_tested"] = np.log1p(df["tested"].astype(float))
    df.sort_values(["organism", "drug", "year"], inplace=True)
    if "lag_susc_pct" not in df.columns:
        df["lag_susc_pct"] = series_lag(df)

    if lag_lookup is not None:
        need_lag = df["lag_susc_pct"].isna()
        if need_lag.any():
            lookup = lag_lookup.copy()
            lookup["join_year"] = lookup["year"] + 1
            merged = df.loc[need_lag].merge(
                lookup[["organism", "drug", "join_year", "susc_pct"]]
                    .rename(columns={"susc_pct": "lag_from_train", "join_year": "year"}),
                on=["organism", "drug", "year"], how="left",
            )
            df.loc[need_lag, "lag_susc_pct"] = merged["lag_from_train"].values

    df["lag_susc_pct"] = df["lag_susc_pct"].fillna(lag_fill)

    gen_df = pd.DataFrame.from_dict(GENETIC_MAP, orient="index").reset_index()
    gen_df.rename(columns={"index": "organism"}, inplace=True)
    df = df.merge(gen_df, on="organism", how="left")
    for col in ["india_cre_rate", "has_ndm", "is_oxa23_driven"]:
        if col not in df.columns or df[col].isna().any():
            df[col] = df[col].fillna(0)
    return df

def make_pipe(model: str, features: list[str]) -> Pipeline:
    cat = [c for c in ["organism", "drug"] if c in features]
    num = [c for c in ["year", "log_tested", "lag_susc_pct", "india_cre_rate"] if c in features]
    bin_ = [c for c in ["has_ndm", "is_oxa23_driven"] if c in features]
    prep = ColumnTransformer([
        ("cat", OneHotEncoder(handle_unknown="ignore"), cat),
        ("num", StandardScaler(),                       num),
        ("bin", "passthrough",                          bin_),
    ])
    if model == "lr":
        clf = LogisticRegression(max_iter=1000, class_weight="balanced",
                                 random_state=RANDOM_STATE)
    elif model == "rf":
        clf = RandomForestClassifier(
            n_estimators=300, max_depth=8, min_samples_leaf=2,
            class_weight="balanced", random_state=RANDOM_STATE, n_jobs=-1,
        )
    else:
        raise ValueError(f"unknown model {model}")
    return Pipeline([("prep", prep), ("clf", clf)])

MODELS = [
    ("no_lag_LR",  "lr", NO_LAG_FEATURES, "PRIMARY"),
    ("full_LR",    "lr", FULL_FEATURES,   "Sensitivity (best-calibrated)"),
    ("full_RF",    "rf", FULL_FEATURES,   "Sensitivity (legacy production model)"),
]

def bootstrap_auc_ci(y_true: np.ndarray, y_proba: np.ndarray,
                     n_bootstrap: int = N_BOOTSTRAP,
                     random_state: int = RANDOM_STATE) -> tuple[float, float]:
    if len(np.unique(y_true)) < 2:
        return float("nan"), float("nan")
    rng = np.random.default_rng(random_state)
    n = len(y_true)
    aucs: list[float] = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, n, n)
        if len(np.unique(y_true[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y_true[idx], y_proba[idx]))
    if not aucs:
        return float("nan"), float("nan")
    return float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))

def calibration_slope_intercept(y_true: np.ndarray, y_proba: np.ndarray) -> tuple[float, float]:
    if len(np.unique(y_true)) < 2 or len(np.unique(y_proba)) < 3:
        return float("nan"), float("nan")
    eps = 1e-6
    p = np.clip(y_proba, eps, 1 - eps)
    logit = np.log(p / (1 - p)).reshape(-1, 1)
    fit = LogisticRegression(max_iter=1000).fit(logit, y_true)
    return float(fit.coef_[0, 0]), float(fit.intercept_[0])

def evaluate(name: str, role: str, y_true: np.ndarray, y_proba: np.ndarray) -> dict:
    y_true = np.asarray(y_true, dtype=int)
    y_proba = np.asarray(y_proba, dtype=float)
    has_both = len(np.unique(y_true)) > 1
    auc = float(roc_auc_score(y_true, y_proba)) if has_both else float("nan")
    auc_lo, auc_hi = bootstrap_auc_ci(y_true, y_proba)
    cal_s, cal_i = calibration_slope_intercept(y_true, y_proba)
    pred_class = (y_proba >= 0.5).astype(int)
    return {
        "model":         name,
        "role":          role,
        "n":             int(len(y_true)),
        "prevalence":    float(y_true.mean()),
        "auc":           auc,
        "auc_ci_lo":     auc_lo,
        "auc_ci_hi":     auc_hi,
        "brier":         float(brier_score_loss(y_true, y_proba)),
        "cal_slope":     cal_s,
        "cal_intercept": cal_i,
        "accuracy":      float((pred_class == y_true).mean()),
    }

def fit_predict_all(train_feat: pd.DataFrame, eval_feat: pd.DataFrame,
                    y_train: np.ndarray) -> dict:
    out: dict = {}
    for name, kind, feats, _role in MODELS:
        pipe = make_pipe(kind, feats)
        pipe.fit(train_feat[feats], y_train)
        proba = pipe.predict_proba(eval_feat[feats])[:, 1]
        out[name] = proba
    return out

def plot_calibration(y_true: np.ndarray, preds: dict, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot([0, 1], [0, 1], "k--", linewidth=1, label="perfect")
    for name in preds:
        df = pd.DataFrame({"p": preds[name], "y": y_true})
        try:
            df["bin"] = pd.qcut(df["p"], q=5, duplicates="drop")
        except ValueError:
            continue
        binned = df.groupby("bin", observed=True).agg(
            mean_pred=("p", "mean"),
            mean_obs =("y", "mean"),
            n        =("y", "size"),
        ).reset_index()
        ax.plot(binned["mean_pred"], binned["mean_obs"], marker="o", label=name)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel("Mean predicted probability")
    ax.set_ylabel("Observed fraction high-risk")
    ax.set_title("Calibration on 25% holdout")
    ax.legend(loc="lower right", fontsize=9)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, facecolor="white", edgecolor="white")
    plt.close()

def plot_auc_forest(rows: list[dict], title: str, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 0.6 * len(rows) + 1.5))
    ys = np.arange(len(rows))[::-1]
    for r, y in zip(rows, ys):
        if not np.isnan(r["auc"]):
            ax.errorbar(r["auc"], y,
                        xerr=[[r["auc"] - r["auc_ci_lo"]], [r["auc_ci_hi"] - r["auc"]]],
                        fmt="o", capsize=4, color="black")
    ax.set_yticks(ys); ax.set_yticklabels([f"{r['model']}\n({r['role']})" for r in rows])
    ax.set_xlim(0.3, 1.0)
    ax.axvline(0.5, color="grey", linestyle=":", linewidth=1)
    ax.set_xlabel("ROC-AUC (95% bootstrap CI)")
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, facecolor="white", edgecolor="white")
    plt.close()

def plot_rolling_origin_auc(rolling_rows: list[dict], out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    by_model: dict[str, list[tuple[int, float]]] = {}
    for r in rolling_rows:
        by_model.setdefault(r["model"], []).append((r["test_year"], r["auc"]))
    for name, points in by_model.items():
        points = sorted(points)
        xs = [p[0] for p in points]; ys = [p[1] for p in points]
        ax.plot(xs, ys, marker="o", linewidth=2, label=name)
    ax.set_xlabel("Test year (rolling-origin)")
    ax.set_ylabel("ROC-AUC")
    ax.set_ylim(0.5, 1.0)
    ax.set_title("Rolling-origin AUC across test years")
    ax.legend(loc="lower left", fontsize=9)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, facecolor="white", edgecolor="white")
    plt.close()

def run_rolling_origin(df_raw: pd.DataFrame) -> list[dict]:
    rows: list[dict] = []
    for test_year in TEST_YEARS:
        train_df = df_raw[df_raw["year"] <  test_year].copy()
        test_df  = df_raw[df_raw["year"] == test_year].copy()
        if train_df.empty or test_df.empty:
            continue
        lag_fill = float(train_df["susc_pct"].mean())
        train_feat = add_features(train_df, lag_fill=lag_fill)
        test_feat  = add_features(test_df,  lag_fill=lag_fill, lag_lookup=train_df)

        y_tr = (train_feat["susc_pct"] < HIGH_RISK_CUTOFF).astype(int).values
        y_te = (test_feat["susc_pct"]  < HIGH_RISK_CUTOFF).astype(int).values

        preds = fit_predict_all(train_feat, test_feat, y_tr)
        for name, _kind, _feats, role in MODELS:
            metrics = evaluate(name, role, y_te, preds[name])
            metrics.update({
                "test_year":       int(test_year),
                "n_train":         int(len(y_tr)),
                "n_test":          int(len(y_te)),
                "thin_fold":       bool(len(y_tr) < THIN_FOLD_THRESHOLD),
                "prevalence_test": float(y_te.mean()),
            })
            rows.append(metrics)
    return rows

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=Path,
                        default=Path("data/data_full_long.csv"))
    parser.add_argument("--external_csv", type=Path,
                        default=Path("data/atlas_external_india.csv"))
    parser.add_argument("--out", type=Path,
                        default=Path("output"))
    parser.add_argument("--label", type=str, default="CARE-INDIA carbapenem-resistance classifier")
    parser.add_argument("--prefix", type=str, default="headline",
                        help="Filename prefix for outputs (e.g. 'headline' or 'pathb').")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    px = args.prefix

    print(f"\n=== Consolidated headline ML - {args.label} ===")
    print(f"Data:     {args.csv}")
    print(f"External: {args.external_csv}\n")

    df_raw = add_series_lag(load_required(args.csv, "Source"))
    y_all  = (df_raw["susc_pct"] < HIGH_RISK_CUTOFF).astype(int)

    train_df, test_df = train_test_split(
        df_raw, test_size=0.25, stratify=y_all, random_state=RANDOM_STATE
    )
    train_lag_fill = float(train_df["susc_pct"].mean())
    train_feat = add_features(train_df, lag_fill=train_lag_fill)
    test_feat  = add_features(test_df,  lag_fill=train_lag_fill)

    y_train = (train_feat["susc_pct"] < HIGH_RISK_CUTOFF).astype(int).values
    y_test  = (test_feat["susc_pct"]  < HIGH_RISK_CUTOFF).astype(int).values

    preds_holdout = fit_predict_all(train_feat, test_feat, y_train)
    holdout_rows = [
        evaluate(name, role, y_test, preds_holdout[name])
        for name, _kind, _feats, role in MODELS
    ]

    print("Holdout (n={}):".format(len(y_test)))
    print(f"  {'Model':>12}  {'Role':<40} {'AUC':>5}  {'95% CI':>14}  {'Brier':>6}  {'CalSlope':>9}")
    for r in holdout_rows:
        ci = f"[{r['auc_ci_lo']:.2f}, {r['auc_ci_hi']:.2f}]"
        cs = f"{r['cal_slope']:.2f}" if not np.isnan(r['cal_slope']) else "N/A"
        print(f"  {r['model']:>12}  {r['role']:<40} {r['auc']:>5.3f}  {ci:>14}  {r['brier']:>6.3f}  {cs:>9}")

    pd.DataFrame(holdout_rows).to_csv(args.out / f"{px}_holdout_metrics.csv", index=False)
    plot_calibration(y_test, preds_holdout, args.out / f"fig_{px}_calibration.png")
    plot_auc_forest(holdout_rows, "Holdout AUC (95% bootstrap CI)",
                    args.out / f"fig_{px}_auc_forest.png")

    external_rows: list[dict] = []
    if args.external_csv.exists():
        ext_raw = load_required(args.external_csv, "External")
        ext_feat = add_features(ext_raw, lag_fill=train_lag_fill, lag_lookup=train_df)
        y_ext = (ext_feat["susc_pct"] < HIGH_RISK_CUTOFF).astype(int).values
        preds_ext = fit_predict_all(train_feat, ext_feat, y_train)
        external_rows = [
            evaluate(name, role, y_ext, preds_ext[name])
            for name, _kind, _feats, role in MODELS
        ]
        print(f"\nExternal (n={len(y_ext)}):")
        print(f"  {'Model':>12}  {'AUC':>5}  {'95% CI':>14}  {'Acc':>5}")
        for r in external_rows:
            ci = f"[{r['auc_ci_lo']:.2f}, {r['auc_ci_hi']:.2f}]"
            print(f"  {r['model']:>12}  {r['auc']:>5.3f}  {ci:>14}  {r['accuracy']:>5.3f}")
        pd.DataFrame(external_rows).to_csv(args.out / f"{px}_external_metrics.csv", index=False)
    else:
        print(f"\nNo external CSV at {args.external_csv} - skipped.")

    print("\nRolling-origin walk-forward (test years 2022, 2023, 2024):")
    rolling_rows = run_rolling_origin(df_raw)
    if rolling_rows:
        df_roll = pd.DataFrame(rolling_rows)
        df_roll.to_csv(args.out / f"{px}_rolling_origin_metrics.csv", index=False)
        print(f"  {'Year':>5}  {'Model':>12}  {'AUC':>5}  {'95% CI':>14}  {'Brier':>6}")
        for r in rolling_rows:
            ci = f"[{r['auc_ci_lo']:.2f}, {r['auc_ci_hi']:.2f}]"
            print(f"  {r['test_year']:>5}  {r['model']:>12}  {r['auc']:>5.3f}  {ci:>14}  {r['brier']:>6.3f}")

        print("\n  Mean rolling-origin AUC by model:")
        for name, _kind, _feats, _role in MODELS:
            sub = df_roll[df_roll["model"] == name]
            if not sub.empty:
                print(f"    {name:>12}: {sub['auc'].mean():.3f} +/- {sub['auc'].std():.3f}")

        plot_rolling_origin_auc(rolling_rows,
                                args.out / f"fig_{px}_rolling_origin_auc.png")

    summary = {
        "label":              args.label,
        "data_source":        str(args.csv),
        "external_source":    str(args.external_csv) if args.external_csv.exists() else None,
        "primary_model":      "no_lag_LR",
        "high_risk_cutoff":   HIGH_RISK_CUTOFF,
        "split":              "75/25 stratified, RANDOM_STATE=42",
        "lag_feature":        "previous calendar year's susceptibility, computed on the full series before partitioning; "
                              "first-year lags imputed with the training-partition mean only",
        "n_bootstrap":        N_BOOTSTRAP,
        "n_train":            int(len(y_train)),
        "n_holdout":          int(len(y_test)),
        "test_years_rolling": list(TEST_YEARS),
        "holdout_metrics":    holdout_rows,
        "external_metrics":   external_rows,
        "rolling_origin":     rolling_rows,
    }
    (args.out / f"{px}_classifier_metrics.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    if holdout_rows:
        nl = next(r for r in holdout_rows if r["model"] == "no_lag_LR")
        rf = next(r for r in holdout_rows if r["model"] == "full_RF")
        print(f"\nHeadline:")
        print(f"  Primary No-lag LR  : AUC {nl['auc']:.3f} [{nl['auc_ci_lo']:.2f}, {nl['auc_ci_hi']:.2f}]"
              f"  Brier {nl['brier']:.3f}  CalSlope {nl['cal_slope']:.2f}")
        print(f"  Sensitivity Full RF: AUC {rf['auc']:.3f} [{rf['auc_ci_lo']:.2f}, {rf['auc_ci_hi']:.2f}]"
              f"  Brier {rf['brier']:.3f}  CalSlope {rf['cal_slope']:.2f}")
        print(f"  Margin RF over No-lag LR: {rf['auc'] - nl['auc']:+.3f}")

    print("\nWrote:")
    for f in [
        f"{px}_classifier_metrics.json",
        f"{px}_holdout_metrics.csv",
        f"{px}_external_metrics.csv",
        f"{px}_rolling_origin_metrics.csv",
        f"fig_{px}_calibration.png",
        f"fig_{px}_auc_forest.png",
        f"fig_{px}_rolling_origin_auc.png",
    ]:
        print(" ", args.out / f)

if __name__ == "__main__":
    main()
