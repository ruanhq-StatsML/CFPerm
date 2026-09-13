#!/usr/bin/env python3
"""Flatten TencentGR seq → user behavior counts, StandardScaler, plot, detect.

Not unique CS/CD. Not CATE. Time-late vs early is a batch-shift probe.
Feature ranking is localization of sequence summaries, not a unique decomposition.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import IsolationForest
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from tencentgr.behavior_features import (  # noqa: E402
    STRUCTURAL_COLS,
    explode_seq,
    rank_against_binary,
    standard_scale_behavior,
    synthesize_user_behavior,
    try_featuretools_dfs,
    z_matrix,
)
from tencentgr.config import DEFAULT_CFG  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--cache-dir", default=DEFAULT_CFG["cache_dir"])
    p.add_argument("--out-dir", default="experiments/tencentgr/behavior")
    p.add_argument("--events-out", default="", help="optional parquet path for exploded events")
    p.add_argument("--max-users", type=int, default=0, help="0 = all users in cache")
    return p.parse_args()


def _savefig(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def _boxplot(ax, data, names, ylabel, title):
    ax.boxplot(data)
    ax.set_xticklabels(names)
    ax.set_ylabel(ylabel)
    ax.set_title(title)


def _cv_auc(Z: np.ndarray, y: np.ndarray, seed: int = 0) -> tuple[list[float], LogisticRegression]:
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    aucs = []
    for tr, te in skf.split(Z, y):
        clf = LogisticRegression(max_iter=400, solver="lbfgs")
        clf.fit(Z[tr], y[tr])
        p = clf.predict_proba(Z[te])[:, 1]
        aucs.append(float(roc_auc_score(y[te], p)))
    clf_all = LogisticRegression(max_iter=400, solver="lbfgs")
    clf_all.fit(Z, y)
    return aucs, clf_all


def _rank_bar(path: Path, table: pd.DataFrame, title: str, value_col: str, n: int = 12) -> None:
    top = table.head(n)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    names = list(reversed(top["feature"].tolist()))
    vals = list(reversed(top[value_col].tolist()))
    ax.barh(names, vals, color="#4C78A8")
    ax.set_xlabel(value_col)
    ax.set_title(title)
    _savefig(fig, path)


def main() -> None:
    args = parse_args()
    cache = Path(args.cache_dir)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    seq = pd.read_parquet(cache / "seq_df.parquet")
    if args.max_users and args.max_users > 0:
        seq = seq.head(int(args.max_users))
    user_feat = pd.read_parquet(cache / "user_feat.parquet") if (cache / "user_feat.parquet").exists() else None
    item_feat = pd.read_parquet(cache / "item_feat.parquet") if (cache / "item_feat.parquet").exists() else None

    events = explode_seq(seq)
    users = synthesize_user_behavior(events, user_feat=user_feat, item_feat=item_feat)
    scaled, scaler, cols = standard_scale_behavior(users)
    Z = z_matrix(scaled, cols)

    events_path = Path(args.events_out) if args.events_out else cache / "events.parquet"
    events_path.parent.mkdir(parents=True, exist_ok=True)
    events.to_parquet(events_path, index=False)
    users.to_parquet(out / "user_behavior.parquet", index=False)
    users.to_csv(out / "user_behavior.csv", index=False)

    ft_fm = try_featuretools_dfs(events.head(min(len(events), 20000)) if len(events) > 20000 else events)
    featuretools_ok = ft_fm is not None
    if featuretools_ok:
        ft_fm.to_csv(out / "featuretools_dfs_sample.csv", index=False)

    # --- plots ---
    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    hist_cols = ["n_click", "n_conversion", "engage_rate", "span_days", "decay_engage", "n_exp_to_click"]
    for ax, col in zip(axes.ravel(), hist_cols):
        ax.hist(users[col].to_numpy(), bins=40, color="#4C78A8", edgecolor="white")
        ax.set_title(col)
        ax.set_xlabel(col)
        ax.set_ylabel("users")
    _savefig(fig, out / "hist_counts.png")

    late = users["time_late"].to_numpy()
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    _boxplot(
        ax,
        [users.loc[late == 0, "engage_rate"], users.loc[late == 1, "engage_rate"]],
        ["early tmax", "late tmax"],
        "engage_rate (click+conv)/n",
        "window engage rate by last-timestamp batch",
    )
    _savefig(fig, out / "box_engage_by_time.png")

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    _boxplot(
        ax,
        [users.loc[late == 0, "n_click"], users.loc[late == 1, "n_click"]],
        ["early tmax", "late tmax"],
        "n_click",
        "click counts by last-timestamp batch",
    )
    _savefig(fig, out / "box_nclick_by_time.png")

    pca = PCA(n_components=2, random_state=0)
    xy = pca.fit_transform(Z)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, color_by, title in (
        (axes[0], users["any_click"].to_numpy(), "any_click in window"),
        (axes[1], users["time_late"].to_numpy(), "time_late (tmax ≥ median)"),
    ):
        sc = ax.scatter(xy[:, 0], xy[:, 1], c=color_by, cmap="coolwarm", s=8, alpha=0.7)
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.2f})")
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.2f})")
        ax.set_title(title)
        fig.colorbar(sc, ax=ax, fraction=0.046)
    _savefig(fig, out / "pca_scaled_behavior.png")

    z_names = [f"z_{c}" for c in cols]
    corr = scaled[z_names].corr().to_numpy()
    fig, ax = plt.subplots(figsize=(9.5, 8))
    im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(cols)))
    ax.set_yticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=90, fontsize=7)
    ax.set_yticklabels(cols, fontsize=7)
    ax.set_title("corr of StandardScaled behavior features")
    fig.colorbar(im, ax=ax, fraction=0.046)
    _savefig(fig, out / "corr_scaled.png")

    iso = IsolationForest(n_estimators=200, contamination=0.05, random_state=0)
    iso_pred = iso.fit_predict(Z)
    iso_score = iso.decision_function(Z)
    scaled["iso_outlier"] = (iso_pred == -1).astype(np.int8)
    scaled["iso_score"] = iso_score
    users["iso_outlier"] = scaled["iso_outlier"]
    users["iso_score"] = iso_score

    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.scatter(xy[:, 0], xy[:, 1], c=scaled["iso_outlier"], cmap="coolwarm", s=8, alpha=0.75)
    ax.set_title("IsolationForest outliers on scaled counts (5%)")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    _savefig(fig, out / "pca_isoforest.png")

    y_time = users["time_late"].to_numpy()
    aucs_time, clf_time = _cv_auc(Z, y_time)
    coef = {c: float(v) for c, v in zip(cols, clf_time.coef_.ravel())}
    top = sorted(coef.items(), key=lambda kv: abs(kv[1]), reverse=True)[:12]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    names = [k for k, _ in reversed(top)]
    vals = [coef[k] for k in names]
    ax.barh(names, vals, color="#4C78A8")
    ax.set_xlabel("logistic coef (z-scored features → time_late)")
    ax.set_title("which counts move with the time batch")
    _savefig(fig, out / "logit_time_late_coefs.png")

    perm = permutation_importance(clf_time, Z, y_time, n_repeats=10, random_state=0, scoring="roc_auc")
    perm_rows = [
        {"feature": c, "import_mean": float(m), "import_std": float(s)}
        for c, m, s in zip(cols, perm.importances_mean, perm.importances_std)
    ]
    perm_tbl = pd.DataFrame(perm_rows).sort_values("import_mean", ascending=False)
    perm_tbl.to_csv(out / "perm_importance_time_late.csv", index=False)
    _rank_bar(out / "perm_importance_time_late.png", perm_tbl, "permutation importance → time_late (AUC drop)", "import_mean")

    rank_time = rank_against_binary(scaled, cols, y_time)
    rank_click = rank_against_binary(scaled, cols, users["any_click"].to_numpy())
    struct_cols = [c for c in STRUCTURAL_COLS if c in cols]
    Z_struct = z_matrix(scaled, struct_cols)
    aucs_click_struct, clf_struct = _cv_auc(Z_struct, users["any_click"].to_numpy())
    rank_click_struct = rank_against_binary(scaled, struct_cols, users["any_click"].to_numpy())
    rank_time.to_csv(out / "rank_vs_time_late.csv", index=False)
    rank_click.to_csv(out / "rank_vs_any_click.csv", index=False)
    rank_click_struct.to_csv(out / "rank_structural_vs_any_click.csv", index=False)
    _rank_bar(out / "rank_auc_time_late.png", rank_time, "univariate |AUC| vs time_late (ranking, not CS/CD)", "auc_abs")
    _rank_bar(
        out / "rank_auc_click_structural.png",
        rank_click_struct,
        "structural features |AUC| vs any_click (no n_click leak)",
        "auc_abs",
    )

    mean_z_early = scaled.loc[late == 0, z_names].mean()
    mean_z_late = scaled.loc[late == 1, z_names].mean()
    delta = (mean_z_late - mean_z_early).sort_values(key=np.abs, ascending=False)
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    names = list(reversed([c[2:] if c.startswith("z_") else c for c in delta.head(14).index]))
    vals = list(reversed(delta.head(14).to_numpy().tolist()))
    ax.barh(names, vals, color="#4C78A8")
    ax.axvline(0.0, color="black", linewidth=0.8)
    ax.set_xlabel("mean z (late) − mean z (early)")
    ax.set_title("batch shift in scaled sequence summaries")
    _savefig(fig, out / "mean_z_shift_time.png")

    scaled.to_parquet(out / "user_behavior_scaled.parquet", index=False)

    summary = {
        "n_users": int(len(users)),
        "n_events": int(len(events)),
        "events_parquet": str(events_path),
        "scale_cols": list(cols),
        "scaler_mean": {c: float(m) for c, m in zip(cols, scaler.mean_)},
        "scaler_scale": {c: float(s) for c, s in zip(cols, scaler.scale_)},
        "pca_var": [float(x) for x in pca.explained_variance_ratio_],
        "mean_n_click": float(users["n_click"].mean()),
        "mean_n_conversion": float(users["n_conversion"].mean()),
        "mean_engage_rate": float(users["engage_rate"].mean()),
        "frac_any_click": float(users["any_click"].mean()),
        "frac_any_conversion": float(users["any_conversion"].mean()),
        "frac_has_user_feat": float(users["has_user_feat"].mean()),
        "frac_iso_outlier": float(users["iso_outlier"].mean()),
        "engage_rate_early": float(users.loc[users["time_late"] == 0, "engage_rate"].mean()),
        "engage_rate_late": float(users.loc[users["time_late"] == 1, "engage_rate"].mean()),
        "n_click_early": float(users.loc[users["time_late"] == 0, "n_click"].mean()),
        "n_click_late": float(users.loc[users["time_late"] == 1, "n_click"].mean()),
        "time_late_cv_auc": aucs_time,
        "time_late_cv_auc_mean": float(np.mean(aucs_time)),
        "time_late_logit_top": [{"feature": k, "coef": v} for k, v in top],
        "perm_importance_top": perm_tbl.head(8).to_dict(orient="records"),
        "rank_time_late_top": rank_time.head(8).to_dict(orient="records"),
        "any_click_structural_cv_auc": aucs_click_struct,
        "any_click_structural_cv_auc_mean": float(np.mean(aucs_click_struct)),
        "featuretools_available": bool(featuretools_ok),
        "other_sequence_summaries": [
            "counts and rates (n_click, click_rate, engage_rate)",
            "recency (recency_to_end_days, pos_last_click, decay_engage, last-20)",
            "pacing (span_days, mean/median inter-event gap)",
            "diversity (n_unique_items, repeat_rate, n_item_format)",
            "transitions (exposure→click)",
            "item-side mix (mean_item_log_freq)",
            "missingness (has_user_feat; mm missing stays on the tower path)",
        ],
        "note": (
            "Counts/rates StandardScaled for tabular/DR. IsolationForest flags unusual "
            "sequences. Logistic AUC and permutation importance are a batch-shift probe "
            "(time_late from counts), not CATE. Rank tables localize sequence summaries; "
            "they are not a unique CS/CD decomposition. featuretools is optional: the "
            "default path is pandas DFS (COUNT/NUM_UNIQUE/MEAN/LAST/window/decay)."
        ),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps({k: summary[k] for k in (
        "n_users", "n_events", "mean_engage_rate", "frac_any_click",
        "engage_rate_early", "engage_rate_late", "time_late_cv_auc_mean",
        "any_click_structural_cv_auc_mean", "frac_iso_outlier", "featuretools_available",
    )}, indent=2))
    print("wrote", out)


if __name__ == "__main__":
    main()
