#!/usr/bin/env python3
"""Parsimonious CausalForest empirical tool on TencentGR W1/W2 graph feats.

Empirical stacking aid — NOT a causal ATE claim:
  W = period indicator (0=W1, 1=W2); Y = continuous proxy (PCA rank of X).
  Use CF for (i) feature_importances_  (ii) |τ̂| / τ̂² aggregated by entity
  as an alternate drift score next to MMD — then compare rankings.

  PYTHONPATH=. python3 scripts/tencent_gr/run_cf_parsimonious_smoke.py \\
    --root data/tencent_subset --max-users 15000 --gap-days 30
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LogisticRegression

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
_PY = Path(__file__).resolve().parents[2] / "Python" / "src"
if str(_PY) not in sys.path:
    sys.path.insert(0, str(_PY))

from time_window_feats import (  # noqa: E402
    DAY,
    feature_engineer,
    fsds_feature_columns,
    propose_two_windows,
    scan_time_range,
)
from run_standardize_mmd_fsds import (  # noqa: E402
    _matrix,
    fit_standardizer,
    rbf_mmd2,
    standardize,
    w1w2_candidate_columns,
)
from run_three_step_subset_localize import (  # noqa: E402
    attach_merchant,
    load_item_merchant_map,
)

ROOT = Path(__file__).resolve().parents[2]


def _pca_rank_y(X: np.ndarray, seed: int) -> np.ndarray:
    z = PCA(n_components=1, random_state=seed).fit_transform(X).ravel()
    order = np.argsort(z, kind="mergesort")
    y = np.empty(len(z), float)
    y[order] = np.linspace(0.0, 1.0, len(z))
    return y


def fit_causal_forest(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    *,
    seed: int,
    n_estimators: int,
    min_samples_leaf: int,
):
    from econml.dml import CausalForestDML

    model_y = RandomForestRegressor(
        n_estimators=80,
        max_depth=6,
        min_samples_leaf=8,
        random_state=seed,
        n_jobs=1,
    )
    model_t = LogisticRegression(max_iter=500, random_state=seed)
    cf = CausalForestDML(
        model_y=model_y,
        model_t=model_t,
        n_estimators=n_estimators,
        min_samples_leaf=min_samples_leaf,
        max_depth=6,
        discrete_treatment=True,
        random_state=seed,
    )
    cf.fit(Y=Y, T=W.astype(int), X=X)
    return cf


def entity_score_from_tau(
    df: pd.DataFrame,
    tau2: np.ndarray,
    *,
    entity_col: str,
    min_n: int,
) -> pd.DataFrame:
    d = df[[entity_col]].copy()
    d["tau2"] = tau2
    g = d.groupby(entity_col, sort=False)["tau2"].agg(["mean", "count"])
    g = g[g["count"] >= min_n].reset_index()
    g = g.rename(columns={"mean": "cf_tau2_mean", "count": "n"})
    g = g.sort_values("cf_tau2_mean", ascending=False).reset_index(drop=True)
    g["rank"] = np.arange(1, len(g) + 1)
    return g


def entity_mmd_scores(
    g1: pd.DataFrame,
    g2: pd.DataFrame,
    cols: Sequence[str],
    sc,
    *,
    entity_col: str,
    seed: int,
    max_cand: int,
    min_edges: int,
    mmd_max_n: int,
) -> pd.DataFrame:
    c1 = g1[entity_col].value_counts()
    c2 = g2[entity_col].value_counts()
    common = sorted(set(c1[c1 >= min_edges].index) & set(c2[c2 >= min_edges].index))
    vol = {e: int(c1.get(e, 0) + c2.get(e, 0)) for e in common}
    cand = sorted(common, key=lambda e: -vol[e])[:max_cand]
    rng = np.random.default_rng(seed)
    rows = []
    for e in cand:
        a = standardize(sc, _matrix(g1[g1[entity_col] == e], cols))
        b = standardize(sc, _matrix(g2[g2[entity_col] == e], cols))
        rows.append(
            {
                entity_col: int(e),
                "mmd2": rbf_mmd2(a, b, max_n=mmd_max_n, rng=rng),
                "n_W1": int(len(a)),
                "n_W2": int(len(b)),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values("mmd2", ascending=False).reset_index(drop=True)


def rank_overlap(a: Sequence, b: Sequence, k: int) -> float:
    sa, sb = set(list(a)[:k]), set(list(b)[:k])
    if not sa or not sb:
        return float("nan")
    return float(len(sa & sb) / k)


def spearman_on_common(df_a: pd.DataFrame, df_b: pd.DataFrame, key: str, col_a: str, col_b: str) -> float:
    m = df_a[[key, col_a]].merge(df_b[[key, col_b]], on=key, how="inner")
    if len(m) < 5:
        return float("nan")
    ra = m[col_a].rank()
    rb = m[col_b].rank()
    return float(np.corrcoef(ra, rb)[0, 1])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT / "data" / "tencent_subset")
    ap.add_argument("--max-users", type=int, default=15000)
    ap.add_argument("--window-days", type=int, default=45)
    ap.add_argument("--gap-days", type=int, default=30)
    ap.add_argument("--merchant-col", type=str, default="122")
    ap.add_argument("--max-rows", type=int, default=12000)
    ap.add_argument("--n-estimators", type=int, default=80)
    ap.add_argument("--min-leaf", type=int, default=20)
    ap.add_argument("--min-edges", type=int, default=3)
    ap.add_argument("--max-cand", type=int, default=200)
    ap.add_argument("--topk", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "results" / "tencent_gr_cf_parsimonious",
    )
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print("FE + merchant map ...", flush=True)
    imap = load_item_merchant_map(args.root, merchant_col=args.merchant_col)
    t_min, t_max, _ = scan_time_range(args.root, max_users=args.max_users)
    w1, w2 = propose_two_windows(
        t_min, t_max, window_days=args.window_days, gap_days=args.gap_days
    )
    gap_d = (w2.t_start - w1.t_end) / DAY
    p1 = feature_engineer(
        args.root, w1.t_start, w1.t_end,
        max_users=args.max_users, window_name=w1.name, terminal_action=1,
    )
    p2 = feature_engineer(
        args.root, w2.t_start, w2.t_end,
        max_users=args.max_users, window_name=w2.name, terminal_action=1,
    )
    g1 = attach_merchant(p1["grid"], imap)
    g2 = attach_merchant(p2["grid"], imap)
    cols = w1w2_candidate_columns([c for c in fsds_feature_columns(g1) if c in g2.columns])
    cols = [c for c in cols if c not in ("merchant_id", "merchant_mapped")]
    print(f"gap={gap_d:.0f}d edges W1={len(g1)} W2={len(g2)} feats={len(cols)}", flush=True)

    sc = fit_standardizer(_matrix(g1, cols))
    rng = np.random.default_rng(args.seed)
    n1 = min(len(g1), args.max_rows // 2)
    n2 = min(len(g2), args.max_rows // 2)
    i1 = rng.choice(len(g1), n1, replace=False)
    i2 = rng.choice(len(g2), n2, replace=False)
    g1s, g2s = g1.iloc[i1].copy(), g2.iloc[i2].copy()
    X1 = standardize(sc, _matrix(g1s, cols))
    X2 = standardize(sc, _matrix(g2s, cols))
    X = np.vstack([X1, X2])
    W = np.concatenate([np.zeros(len(X1), int), np.ones(len(X2), int)])
    Y = _pca_rank_y(X, args.seed)
    pooled = pd.concat([g1s, g2s], ignore_index=True)

    print("fit CausalForestDML (parsimonious) ...", flush=True)
    t0 = time.time()
    cf = fit_causal_forest(
        X, Y, W, seed=args.seed,
        n_estimators=args.n_estimators, min_samples_leaf=args.min_leaf,
    )
    sec = time.time() - t0
    tau = np.asarray(cf.effect(X)).reshape(-1)
    tau2 = tau ** 2
    if hasattr(cf, "feature_importances_"):
        imp = np.asarray(cf.feature_importances_, float).ravel()[: len(cols)]
    else:
        imp = np.zeros(len(cols))
    feat_rank = (
        pd.DataFrame({"feature": cols, "cf_importance": imp})
        .sort_values("cf_importance", ascending=False)
        .reset_index(drop=True)
    )
    feat_rank["rank"] = np.arange(1, len(feat_rank) + 1)
    print(f"  CF fit sec={sec:.1f}  mean(tau^2)={float(tau2.mean()):.6f}", flush=True)

    # entity scores: CF tau2 vs MMD
    print("entity scores CF vs MMD ...", flush=True)
    cf_mer = entity_score_from_tau(pooled, tau2, entity_col="merchant_id", min_n=args.min_edges)
    cf_usr = entity_score_from_tau(pooled, tau2, entity_col="user_id", min_n=max(2, args.min_edges - 1))
    mmd_mer = entity_mmd_scores(
        g1, g2, cols, sc, entity_col="merchant_id", seed=args.seed,
        max_cand=args.max_cand, min_edges=args.min_edges, mmd_max_n=64,
    )
    mmd_usr = entity_mmd_scores(
        g1, g2, cols, sc, entity_col="user_id", seed=args.seed + 1,
        max_cand=args.max_cand, min_edges=2, mmd_max_n=64,
    )

    k = args.topk
    cmp = {
        "merchant_spearman_cf_vs_mmd": spearman_on_common(
            cf_mer, mmd_mer, "merchant_id", "cf_tau2_mean", "mmd2"
        ),
        "merchant_top_overlap": rank_overlap(
            cf_mer["merchant_id"].tolist(), mmd_mer["merchant_id"].tolist(), k
        ),
        "user_spearman_cf_vs_mmd": spearman_on_common(
            cf_usr, mmd_usr, "user_id", "cf_tau2_mean", "mmd2"
        ),
        "user_top_overlap": rank_overlap(
            cf_usr["user_id"].tolist(), mmd_usr["user_id"].tolist(), k
        ),
        "global_po_risk_proxy": float(tau2.mean()),
        "cf_fit_sec": sec,
        "disclaimer": (
            "W=period indicator; CF is a parsimonious empirical heterogeneity tool, "
            "NOT a causal ATE claim. Compare rankings to MMD localization."
        ),
    }

    feat_rank.to_csv(args.out_dir / "cf_feature_importance.csv", index=False)
    cf_mer.to_csv(args.out_dir / "cf_merchant_tau2.csv", index=False)
    cf_usr.to_csv(args.out_dir / "cf_user_tau2.csv", index=False)
    mmd_mer.to_csv(args.out_dir / "mmd_merchant.csv", index=False)
    mmd_usr.to_csv(args.out_dir / "mmd_user.csv", index=False)
    (args.out_dir / "summary.json").write_text(json.dumps({**cmp, "gap_days": gap_d, "n_rows": int(len(W)), "n_features": len(cols)}, indent=2))

    # plot
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    ax = axes[0]
    top = feat_rank.head(12)
    ax.barh(range(len(top)), top["cf_importance"][::-1], color="#4C78A8")
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["feature"][::-1], fontsize=8)
    ax.set_xlabel("CF feature_importances_")
    ax.set_title("CausalForest (empirical) feature rank")

    ax = axes[1]
    m = cf_mer.merge(mmd_mer, on="merchant_id", how="inner")
    if len(m):
        ax.scatter(m["mmd2"], m["cf_tau2_mean"], s=14, alpha=0.7, color="#F58518")
    ax.set_xlabel("MMD² (merchant)")
    ax.set_ylabel("CF mean(τ̂²)")
    ax.set_title(
        f"Merchant rank agree  Spearman={cmp['merchant_spearman_cf_vs_mmd']:.2f}  "
        f"overlap@{k}={cmp['merchant_top_overlap']:.2f}"
    )
    fig.suptitle(f"CF parsimonious vs MMD  gap={gap_d:.0f}d  (W=period, non-causal claim)", fontsize=10)
    fig.tight_layout()
    fig.savefig(args.out_dir / "cf_parsimonious_smoke.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    md = [
        "# CausalForest parsimonious empirical smoke",
        "",
        "> **Not a causal ATE claim.** `W` = period (W1/W2). CF is a lean empirical tool",
        "> for heterogeneity / feature rank next to **MMD localization**.",
        "",
        f"- gap={gap_d:.0f}d | rows={len(W)} | feats={len(cols)} | CF trees={args.n_estimators}",
        f"- global mean(τ̂²) ≈ **{cmp['global_po_risk_proxy']:.6f}** (PO-risk-like scalar)",
        f"- merchant Spearman(CF τ̂², MMD) = **{cmp['merchant_spearman_cf_vs_mmd']:.3f}**",
        f"- merchant top-{k} overlap = **{cmp['merchant_top_overlap']:.3f}**",
        f"- user Spearman = **{cmp['user_spearman_cf_vs_mmd']:.3f}** | overlap@{k} = **{cmp['user_top_overlap']:.3f}**",
        "",
        "## CF feature importance (head)",
        "| rank | feature | importance |",
        "|---:|---|---:|",
    ]
    for _, r in feat_rank.head(12).iterrows():
        md.append(f"| {int(r['rank'])} | `{r['feature']}` | {r['cf_importance']:.4f} |")
    md += [
        "",
        "## Role in the pipeline",
        "- **MMD path**: subset localization (who drifted)",
        "- **CF path**: parsimonious tau-hat / feature_importances_ (alternate drift score + ranking)",
        "- **Stop-drill / business logic**: unchanged — CF does not replace drill triggers",
        "",
    ]
    (args.out_dir / "CF_PARSIMONIOUS_REPORT.md").write_text("\n".join(md))
    print("wrote", args.out_dir)
    print(json.dumps(cmp, indent=2))


if __name__ == "__main__":
    main()
