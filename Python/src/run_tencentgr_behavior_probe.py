#!/usr/bin/env python3
"""Flatten TencentGR seq → user behavior counts, StandardScaler, plot, detect.

Not unique CS/CD. Not CATE. Time-late vs early is a batch shift probe.
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
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from tencentgr.behavior_features import (  # noqa: E402
    SCALE_COLS,
    explode_seq,
    standard_scale_behavior,
    synthesize_user_behavior,
    z_matrix,
)
from tencentgr.config import DEFAULT_CFG  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--cache-dir", default=DEFAULT_CFG["cache_dir"])
    p.add_argument("--out-dir", default="experiments/tencentgr/behavior")
    p.add_argument("--events-out", default="", help="optional parquet path for exploded events")
    return p.parse_args()


def _savefig(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    cache = Path(args.cache_dir)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    seq = pd.read_parquet(cache / "seq_df.parquet")
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
    scaled.to_parquet(out / "user_behavior_scaled.parquet", index=False)
    users.to_csv(out / "user_behavior.csv", index=False)

    # --- plots ---
    fig, axes = plt.subplots(2, 2, figsize=(9, 7))
    for ax, col in zip(axes.ravel(), ["n_click", "n_conversion", "engage_rate", "span_days"]):
        ax.hist(users[col].to_numpy(), bins=40, color="#4C78A8", edgecolor="white")
        ax.set_title(col)
        ax.set_xlabel(col)
        ax.set_ylabel("users")
    _savefig(fig, out / "hist_counts.png")

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    late = users["time_late"].to_numpy()
    ax.boxplot(
        [users.loc[late == 0, "engage_rate"], users.loc[late == 1, "engage_rate"]],
        labels=["early tmax", "late tmax"],
    )
    ax.set_ylabel("engage_rate (click+conv)/n")
    ax.set_title("window engage rate by last-timestamp batch")
    _savefig(fig, out / "box_engage_by_time.png")

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.boxplot(
        [users.loc[late == 0, "n_click"], users.loc[late == 1, "n_click"]],
        labels=["early tmax", "late tmax"],
    )
    ax.set_ylabel("n_click")
    ax.set_title("click counts by last-timestamp batch")
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
    fig, ax = plt.subplots(figsize=(8.5, 7))
    im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(cols)))
    ax.set_yticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=90, fontsize=7)
    ax.set_yticklabels(cols, fontsize=7)
    ax.set_title("corr of StandardScaled behavior features")
    fig.colorbar(im, ax=ax, fraction=0.046)
    _savefig(fig, out / "corr_scaled.png")

    # IsolationForest: unusual sequences in z-space (detection, not treatment)
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

    # Can behavior counts recover the time-batch label? Diagnostic AUC, not CATE.
    y = users["time_late"].to_numpy()
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
    aucs = []
    for tr, te in skf.split(Z, y):
        clf = LogisticRegression(max_iter=400, solver="lbfgs")
        clf.fit(Z[tr], y[tr])
        p = clf.predict_proba(Z[te])[:, 1]
        aucs.append(float(roc_auc_score(y[te], p)))
    clf_all = LogisticRegression(max_iter=400, solver="lbfgs")
    clf_all.fit(Z, y)
    coef = {c: float(v) for c, v in zip(cols, clf_all.coef_.ravel())}
    top = sorted(coef.items(), key=lambda kv: abs(kv[1]), reverse=True)[:12]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    names = [k for k, _ in reversed(top)]
    vals = [coef[k] for k in names]
    ax.barh(names, vals, color="#4C78A8")
    ax.set_xlabel("logistic coef (z-scored features → time_late)")
    ax.set_title("which counts move with the time batch")
    _savefig(fig, out / "logit_time_late_coefs.png")

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
        "time_late_cv_auc": aucs,
        "time_late_cv_auc_mean": float(np.mean(aucs)),
        "time_late_logit_top": [{"feature": k, "coef": v} for k, v in top],
        "note": (
            "Counts/rates StandardScaled for tabular/DR. IsolationForest flags unusual "
            "sequences. Logistic AUC is a batch-shift probe (time_late from counts), not CATE. "
            "featuretools not required: this is pandas DFS (COUNT/NUM_UNIQUE/MEAN/LAST/window)."
        ),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    scaled.to_parquet(out / "user_behavior_scaled.parquet", index=False)

    print(json.dumps({k: summary[k] for k in (
        "n_users", "n_events", "mean_engage_rate", "frac_any_click",
        "engage_rate_early", "engage_rate_late", "time_late_cv_auc_mean", "frac_iso_outlier",
    )}, indent=2))
    print("wrote", out)


if __name__ == "__main__":
    main()
