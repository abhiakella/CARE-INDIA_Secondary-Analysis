#!/usr/bin/env python3
# CARE-INDIA (Carbapenem resistance among Gram-negative ESKAPEE pathogens in India, 2017-2024)
# Classifier baseline comparators with calibration and bootstrap-CI diagnostics.
# Authors: Abhishek Akella, Anand Srinivasan  |  Dept of Pharmacology, AIIMS Bhubaneswar, India

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
    print("Install: pip install -r requirements.txt", file=sys.stderr)
    sys.exit(1)

RANDOM_STATE     = 42
HIGH_RISK_CUTOFF = 60
N_BOOTSTRAP      = 1000

GENETIC_MAP = {
    "K. pneumoniae":     {"india_cre_rate": 58.0, "has_ndm": 1, "is_oxa23_driven": 0},
    "A. baumannii":      {"india_cre_rate": 92.0, "has_ndm": 1, "is_oxa23_driven": 1},
    "P. aeruginosa":     {"india_cre_rate": 30.0, "has_ndm": 1, "is_oxa23_driven": 0},
    "Enterobacter spp.": {"india_cre_rate": 20.0, "has_ndm": 1, "is_oxa23_driven": 0},
    "E. coli":           {"india_cre_rate": 34.0, "has_ndm": 1, "is_oxa23_driven": 0},
}

FULL_FEATURES = [
    "organism", "drug", "year", "log_tested", "lag_susc_pct",
    "india_cre_rate", "has_ndm", "is_oxa23_driven",
]
NO_LAG_FEATURES = ["organism", "drug", "year", "log_tested"]
NO_VOLUME_FEATURES = ["organism", "drug", "year"]

def series_lag(df: pd.DataFrame) -> np.ndarray:
    prev = df[["organism", "drug", "year", "susc_pct"]].copy()
    prev["year"] = prev["year"] + 1
    prev = prev.rename(columns={"susc_pct": "lag_susc_pct"})
    return df[["organism", "drug", "year"]].merge(prev, on=["organism", "drug", "year"], how="left")["lag_susc_pct"].to_numpy()

def add_series_lag(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["lag_susc_pct"] = series_lag(df)
    return df

def aas_root() -> Path:
    return Path(__file__).resolve().parents[1]

def load_required(csv_path: Path, label: str) -> pd.DataFrame:
    if not csv_path.exists():
        raise SystemExit(f"{label} CSV not found: {csv_path}")
    df = pd.read_csv(csv_path)
    need = {"organism", "drug", "year", "tested", "susc_pct", "res_pct"}
    miss = need - set(df.columns)
    if miss:
        raise SystemExit(f"{label} CSV missing columns: {miss}")
    return df

def add_lag_and_genetic(
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
    return df

def make_lr_pipe(features: list[str]) -> Pipeline:
    cat = [c for c in ["organism", "drug"] if c in features]
    num = [c for c in ["year", "log_tested", "lag_susc_pct", "india_cre_rate"] if c in features]
    bin_ = [c for c in ["has_ndm", "is_oxa23_driven"] if c in features]
    return Pipeline([
        ("prep", ColumnTransformer([
            ("cat", OneHotEncoder(handle_unknown="ignore"), cat),
            ("num", StandardScaler(),                       num),
            ("bin", "passthrough",                          bin_),
        ])),
        ("lr", LogisticRegression(max_iter=1000, class_weight="balanced",
                                  random_state=RANDOM_STATE)),
    ])

def make_rf_pipe(features: list[str]) -> Pipeline:
    cat = [c for c in ["organism", "drug"] if c in features]
    num = [c for c in ["year", "log_tested", "lag_susc_pct", "india_cre_rate"] if c in features]
    bin_ = [c for c in ["has_ndm", "is_oxa23_driven"] if c in features]
    return Pipeline([
        ("prep", ColumnTransformer([
            ("cat", OneHotEncoder(handle_unknown="ignore"), cat),
            ("num", StandardScaler(),                       num),
            ("bin", "passthrough",                          bin_),
        ])),
        ("rf", RandomForestClassifier(
            n_estimators=300, max_depth=8, min_samples_leaf=2,
            class_weight="balanced", random_state=RANDOM_STATE, n_jobs=-1,
        )),
    ])

def bootstrap_auc_ci(y_true: np.ndarray, y_proba: np.ndarray,
                     n_bootstrap: int = N_BOOTSTRAP,
                     random_state: int = RANDOM_STATE) -> tuple[float, float]:
    if len(np.unique(y_true)) < 2:
        return float("nan"), float("nan")
    rng = np.random.default_rng(random_state)
    n = len(y_true)
    aucs = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, n, n)
        if len(np.unique(y_true[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y_true[idx], y_proba[idx]))
    if not aucs:
        return float("nan"), float("nan")
    return float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))

def calibration_slope_intercept(y_true: np.ndarray, y_proba: np.ndarray) -> tuple[float, float]:
    unique_p = np.unique(y_proba)
    if len(unique_p) < 3 or len(np.unique(y_true)) < 2:
        return float("nan"), float("nan")
    eps   = 1e-6
    p     = np.clip(y_proba, eps, 1 - eps)
    logit = np.log(p / (1 - p)).reshape(-1, 1)
    fit   = LogisticRegression(max_iter=1000).fit(logit, y_true)
    return float(fit.coef_[0, 0]), float(fit.intercept_[0])

def evaluate(name: str, y_true: np.ndarray, y_proba: np.ndarray,
             y_class: np.ndarray) -> dict:
    y_true = np.asarray(y_true, dtype=int)
    y_proba = np.asarray(y_proba, dtype=float)
    y_class = np.asarray(y_class, dtype=int)

    has_both = len(np.unique(y_true)) > 1
    auc      = float(roc_auc_score(y_true, y_proba)) if has_both else float("nan")
    auc_lo, auc_hi = bootstrap_auc_ci(y_true, y_proba)
    brier    = float(brier_score_loss(y_true, y_proba))
    cal_s, cal_i = calibration_slope_intercept(y_true, y_proba)
    accuracy = float((y_class == y_true).mean())

    return {
        "model":            name,
        "n":                int(len(y_true)),
        "prevalence":       float(y_true.mean()),
        "auc":              auc,
        "auc_ci_lo":        auc_lo,
        "auc_ci_hi":        auc_hi,
        "brier":            brier,
        "cal_slope":        cal_s,
        "cal_intercept":    cal_i,
        "accuracy":         accuracy,
    }

def predictions_for_all_models(
    train_feat: pd.DataFrame, eval_feat: pd.DataFrame,
    y_train: pd.Series, y_eval: pd.Series,
) -> dict:
    out = {}

    naive_class = (eval_feat["lag_susc_pct"] < HIGH_RISK_CUTOFF).astype(int).values
    out["naive_lag_rule"] = {
        "y_proba": naive_class.astype(float),
        "y_class": naive_class,
    }

    lag_score = (1.0 - eval_feat["lag_susc_pct"].clip(0, 100) / 100.0).values
    out["lag_as_score"] = {
        "y_proba": lag_score,
        "y_class": (lag_score >= 0.5).astype(int),
    }

    pair_mean   = train_feat.groupby(["organism", "drug"])["susc_pct"].mean()
    global_mean = float(train_feat["susc_pct"].mean())
    cm_susc = eval_feat.apply(
        lambda r: pair_mean.get((r["organism"], r["drug"]), global_mean), axis=1
    ).astype(float).values
    out["cell_mean_lookup"] = {
        "y_proba": (1.0 - np.clip(cm_susc, 0, 100) / 100.0),
        "y_class": (cm_susc < HIGH_RISK_CUTOFF).astype(int),
    }

    lr_no_lag = make_lr_pipe(NO_LAG_FEATURES)
    lr_no_lag.fit(train_feat[NO_LAG_FEATURES], y_train)
    no_lag_proba = lr_no_lag.predict_proba(eval_feat[NO_LAG_FEATURES])[:, 1]
    out["no_lag_LR"] = {
        "y_proba": no_lag_proba,
        "y_class": (no_lag_proba >= 0.5).astype(int),
    }

    lr_no_vol = make_lr_pipe(NO_VOLUME_FEATURES)
    lr_no_vol.fit(train_feat[NO_VOLUME_FEATURES], y_train)
    no_vol_proba = lr_no_vol.predict_proba(eval_feat[NO_VOLUME_FEATURES])[:, 1]
    out["no_lag_LR_no_volume"] = {
        "y_proba": no_vol_proba,
        "y_class": (no_vol_proba >= 0.5).astype(int),
    }

    lr_full = make_lr_pipe(FULL_FEATURES)
    lr_full.fit(train_feat[FULL_FEATURES], y_train)
    full_lr_proba = lr_full.predict_proba(eval_feat[FULL_FEATURES])[:, 1]
    out["full_LR"] = {
        "y_proba": full_lr_proba,
        "y_class": (full_lr_proba >= 0.5).astype(int),
    }

    rf_full = make_rf_pipe(FULL_FEATURES)
    rf_full.fit(train_feat[FULL_FEATURES], y_train)
    full_rf_proba = rf_full.predict_proba(eval_feat[FULL_FEATURES])[:, 1]
    out["full_RF"] = {
        "y_proba": full_rf_proba,
        "y_class": (full_rf_proba >= 0.5).astype(int),
    }

    return out

def plot_calibration(preds: dict, y_true: np.ndarray, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot([0, 1], [0, 1], "k--", linewidth=1, label="perfect")
    for name in ["lag_as_score", "no_lag_LR", "full_LR", "full_RF"]:
        proba = preds[name]["y_proba"]
        df = pd.DataFrame({"p": proba, "y": y_true})
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
    plt.savefig(out_path, dpi=200)
    plt.close()

def plot_auc_forest(rows: list[dict], title: str, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 0.6 * len(rows) + 1.5))
    ys = np.arange(len(rows))[::-1]
    aucs    = [r["auc"]        for r in rows]
    los     = [r["auc_ci_lo"]  for r in rows]
    his     = [r["auc_ci_hi"]  for r in rows]
    labels  = [r["model"]      for r in rows]
    for y, a, lo, hi in zip(ys, aucs, los, his):
        if not np.isnan(a):
            ax.errorbar(a, y, xerr=[[a - lo], [hi - a]], fmt="o",
                        capsize=4, color="black")
    ax.set_yticks(ys); ax.set_yticklabels(labels)
    ax.set_xlim(0.4, 1.0)
    ax.axvline(0.5, color="grey", linestyle=":", linewidth=1)
    ax.set_xlabel("ROC-AUC (95% bootstrap CI)")
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=Path,
                        default=aas_root() / "data" / "data_full_long.csv")
    parser.add_argument("--external_csv", type=Path,
                        default=aas_root() / "data" / "atlas_external_india.csv")
    parser.add_argument("--out", type=Path,
                        default=aas_root() / "output")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    df_raw = add_series_lag(load_required(args.csv, "AMRSN"))
    y_all  = (df_raw["susc_pct"] < HIGH_RISK_CUTOFF).astype(int)

    train_df, test_df = train_test_split(
        df_raw, test_size=0.25, stratify=y_all, random_state=RANDOM_STATE
    )
    train_lag_fill = float(train_df["susc_pct"].mean())
    train_feat = add_lag_and_genetic(train_df, lag_fill=train_lag_fill)
    test_feat  = add_lag_and_genetic(test_df,  lag_fill=train_lag_fill)

    y_train = (train_feat["susc_pct"] < HIGH_RISK_CUTOFF).astype(int).values
    y_test  = (test_feat["susc_pct"]  < HIGH_RISK_CUTOFF).astype(int).values

    preds_holdout = predictions_for_all_models(train_feat, test_feat,
                                               y_train, y_test)
    holdout_rows = [
        evaluate(name, y_test, p["y_proba"], p["y_class"])
        for name, p in preds_holdout.items()
    ]

    print("\n25% Holdout -  model comparison")
    print(f"{'Model':>20}  {'AUC':>5}  {'95% CI':>15}  {'Brier':>6}  "
          f"{'CalSlope':>9}  {'CalInt':>7}  {'Acc':>5}")
    for r in holdout_rows:
        ci  = f"[{r['auc_ci_lo']:.2f}, {r['auc_ci_hi']:.2f}]" if not np.isnan(r["auc_ci_lo"]) else " - "
        cal_s = f"{r['cal_slope']:.2f}" if not np.isnan(r['cal_slope']) else "N/A"
        cal_i = f"{r['cal_intercept']:.2f}" if not np.isnan(r['cal_intercept']) else "N/A"
        print(f"{r['model']:>20}  {r['auc']:>5.3f}  {ci:>15}  "
              f"{r['brier']:>6.3f}  {cal_s:>9}  {cal_i:>7}  {r['accuracy']:>5.3f}")

    pd.DataFrame(holdout_rows).to_csv(args.out / "table_baselines_holdout.csv", index=False)
    plot_calibration(preds_holdout, y_test,
                     args.out / "fig_calibration_holdout.png")
    plot_auc_forest(holdout_rows, "AUC on 25% holdout (95% bootstrap CI)",
                    args.out / "fig_auc_forest_holdout.png")

    external_rows: list[dict] = []
    if args.external_csv.exists():
        ext_raw = load_required(args.external_csv, "External")
        ext_feat = add_lag_and_genetic(
            ext_raw, lag_fill=train_lag_fill, lag_lookup=train_df,
        )
        y_ext = (ext_feat["susc_pct"] < HIGH_RISK_CUTOFF).astype(int).values

        preds_external = predictions_for_all_models(train_feat, ext_feat,
                                                    y_train, y_ext)
        external_rows = [
            evaluate(name, y_ext, p["y_proba"], p["y_class"])
            for name, p in preds_external.items()
        ]

        print("\nExternal (n={}) -  model comparison".format(len(y_ext)))
        for r in external_rows:
            ci = f"[{r['auc_ci_lo']:.2f}, {r['auc_ci_hi']:.2f}]" if not np.isnan(r["auc_ci_lo"]) else " - "
            print(f"  {r['model']:>20}  AUC={r['auc']:.3f}  CI={ci}  "
                  f"Brier={r['brier']:.3f}  Acc={r['accuracy']:.3f}")

        pd.DataFrame(external_rows).to_csv(
            args.out / "table_baselines_external.csv", index=False
        )
        plot_auc_forest(external_rows,
                        f"AUC on external set (n={len(y_ext)}, 95% bootstrap CI)",
                        args.out / "fig_auc_forest_external.png")

        def paired_auc_diff(a: np.ndarray, b: np.ndarray, n_boot: int = N_BOOTSTRAP,
                            seed: int = RANDOM_STATE) -> dict:
            rng = np.random.default_rng(seed); n = len(y_ext); diffs = []
            for _ in range(n_boot):
                idx = rng.integers(0, n, n)
                if len(np.unique(y_ext[idx])) < 2:
                    continue
                diffs.append(roc_auc_score(y_ext[idx], a[idx]) - roc_auc_score(y_ext[idx], b[idx]))
            d = np.asarray(diffs)
            return {"diff": float(roc_auc_score(y_ext, a) - roc_auc_score(y_ext, b)),
                    "ci_lo": float(np.percentile(d, 2.5)), "ci_hi": float(np.percentile(d, 97.5)),
                    "p_two_sided": float(2 * min((d <= 0).mean(), (d >= 0).mean()))}
        pairs = [("full_RF", "lag_as_score"), ("full_RF", "naive_lag_rule"), ("full_LR", "lag_as_score"),
                 ("no_lag_LR", "cell_mean_lookup"), ("no_lag_LR", "no_lag_LR_no_volume")]
        pb_rows = [{"model": a, "baseline": b,
                    **paired_auc_diff(preds_external[a]["y_proba"], preds_external[b]["y_proba"])} for a, b in pairs]
        pd.DataFrame(pb_rows).to_csv(args.out / "table_paired_bootstrap_external.csv", index=False)
        print("\nPaired bootstrap AUC differences (external):")
        for r in pb_rows:
            print(f"  {r['model']:>10} - {r['baseline']:<20} {r['diff']:+.3f} "
                  f"[{r['ci_lo']:+.3f}, {r['ci_hi']:+.3f}]  p={r['p_two_sided']:.2f}")
    else:
        print(f"\nExternal CSV not found at {args.external_csv} -  skipped.")

    summary = {
        "method":               "Seven-model comparator: naive_lag_rule, lag_as_score, cell_mean_lookup, no_lag_LR, "
                                "no_lag_LR_no_volume, full_LR, full_RF",
        "lag_feature":          "previous calendar year's susceptibility, computed on the full series before partitioning; "
                                "first-year lags imputed with the training-partition mean only",
        "split":                "75/25 stratified, RANDOM_STATE=42",
        "high_risk_cutoff":     HIGH_RISK_CUTOFF,
        "n_bootstrap":          N_BOOTSTRAP,
        "n_train":              int(len(y_train)),
        "n_holdout":            int(len(y_test)),
        "prevalence_train":     float(y_train.mean()),
        "prevalence_holdout":   float(y_test.mean()),
        "holdout_metrics":      holdout_rows,
        "external_metrics":     external_rows,
    }
    (args.out / "classifier_baselines_metrics.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    if holdout_rows:
        rf  = next(r for r in holdout_rows if r["model"] == "full_RF")
        nl  = next(r for r in holdout_rows if r["model"] == "naive_lag_rule")
        nol = next(r for r in holdout_rows if r["model"] == "no_lag_LR")
        print("\nHeadline:")
        print(f"  Full RF AUC                : {rf['auc']:.3f} "
              f"[{rf['auc_ci_lo']:.2f}, {rf['auc_ci_hi']:.2f}]")
        print(f"  Naive lag rule AUC         : {nl['auc']:.3f}  "
              f"(margin RF − naive: {rf['auc'] - nl['auc']:+.3f})")
        print(f"  No-lag LR AUC              : {nol['auc']:.3f}  "
              f"(margin RF − no-lag LR: {rf['auc'] - nol['auc']:+.3f})")
        if rf["auc"] - nl["auc"] < 0.02:
            print("  RF margin over naive-lag rule < 0.02 -  "
                  "the model is essentially memorising last year's susceptibility.")
        if rf["auc"] - nol["auc"] > 0.10:
            print("  RF beats no-lag LR by > 0.10 -  lag feature is doing real work.")

    print("\nWrote:")
    for f in [
        "classifier_baselines_metrics.json",
        "table_baselines_holdout.csv",
        "table_baselines_external.csv",
        "fig_calibration_holdout.png",
        "fig_auc_forest_holdout.png",
        "fig_auc_forest_external.png",
        "table_paired_bootstrap_external.csv",
    ]:
        print(" ", args.out / f)

if __name__ == "__main__":
    main()
