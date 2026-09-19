#!/usr/bin/env python3
"""Hierarchical *dissolve* of leaf attribution scores (讲武德).

Note on naming
--------------
PyPI package ``dissolve`` (jelmer/dissolve) is an **API-deprecation migrator** —
not for attribution. Here "dissolve" means the **geopandas-style aggregate**:
group leaf rows by a parent key and roll scores up the hierarchy.

Hierarchy (TencentGR drill):
  order/edge (leaf)  →  user  →  merchant

Use: take a leaf-level shift score (τ̂² proxy, ‖x−μ‖₂, |δ|-energy, …),
dissolve upward, compare to top-down MMD drill for landability.

  from attribution_dissolve import dissolve_scores, drill_hierarchy
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


HIERARCHY_DEFAULT: Tuple[str, ...] = ("order_id", "user_id", "merchant_id")


def dissolve_scores(
    leaf: pd.DataFrame,
    *,
    score_col: str,
    by: str,
    how: str = "sum",
    weight_col: Optional[str] = None,
) -> pd.DataFrame:
    """Aggregate leaf scores to a parent key (dissolve up one level).

    Parameters
    ----------
    leaf : DataFrame with ``by`` and ``score_col``
    how : 'sum' | 'mean' | 'max' | 'l2' (sqrt of sum of squares)
    weight_col : optional mass for weighted mean
    """
    if by not in leaf.columns:
        raise KeyError(f"missing parent key {by}")
    if score_col not in leaf.columns:
        raise KeyError(f"missing score {score_col}")
    g = leaf.dropna(subset=[by]).copy()
    s = g[score_col].astype(float)
    if how == "sum":
        agg = g.groupby(by, sort=False)[score_col].sum()
    elif how == "mean":
        if weight_col and weight_col in g.columns:
            w = g[weight_col].astype(float).clip(lower=0)
            num = g.assign(_ws=s * w).groupby(by, sort=False)["_ws"].sum()
            den = g.groupby(by, sort=False)[weight_col].sum().clip(lower=1e-12)
            agg = num / den
        else:
            agg = g.groupby(by, sort=False)[score_col].mean()
    elif how == "max":
        agg = g.groupby(by, sort=False)[score_col].max()
    elif how == "l2":
        agg = np.sqrt(g.assign(_sq=s**2).groupby(by, sort=False)["_sq"].sum())
    else:
        raise ValueError(f"unknown how={how}")
    n = g.groupby(by, sort=False).size()
    out = (
        pd.DataFrame({by: agg.index, "score": agg.to_numpy(), "n_leaf": n.to_numpy()})
        .sort_values("score", ascending=False)
        .reset_index(drop=True)
    )
    out["rank"] = np.arange(1, len(out) + 1)
    out["dissolve_how"] = how
    out["score_col"] = score_col
    return out


def drill_hierarchy(
    leaf: pd.DataFrame,
    *,
    score_col: str,
    levels: Sequence[str] = ("user_id", "merchant_id"),
    how: str = "sum",
    top_k: Optional[Dict[str, int]] = None,
) -> Dict[str, pd.DataFrame]:
    """Dissolve leaf → each parent level; optionally keep TopK per level.

    Returns dict level_name → ranked table.
    """
    top_k = top_k or {}
    out: Dict[str, pd.DataFrame] = {"leaf": leaf.copy()}
    for level in levels:
        tab = dissolve_scores(leaf, score_col=score_col, by=level, how=how)
        k = top_k.get(level)
        if k is not None:
            tab["selected"] = (tab.index < int(k)).astype(int)
        else:
            tab["selected"] = 0
        out[level] = tab
    return out


def jaccard_topk(a: Sequence, b: Sequence, k: int) -> float:
    sa, sb = set(list(a)[:k]), set(list(b)[:k])
    if not sa and not sb:
        return 1.0
    return float(len(sa & sb) / len(sa | sb))


def rank_overlap_table(
    dissolve_rank: pd.DataFrame,
    reference_rank: pd.DataFrame,
    *,
    key: str,
    ks: Sequence[int] = (5, 10, 20, 50),
) -> pd.DataFrame:
    """Compare dissolve TopK vs a reference entity ranking (e.g. MMD drill)."""
    a = dissolve_rank[key].tolist()
    b = reference_rank[key].tolist()
    rows = []
    for k in ks:
        if k > max(len(a), len(b)):
            continue
        rows.append(
            {
                "k": k,
                "jaccard": jaccard_topk(a, b, k),
                "n_dissolve": min(k, len(a)),
                "n_ref": min(k, len(b)),
            }
        )
    return pd.DataFrame(rows)


def leaf_shift_scores(
    g_w1: pd.DataFrame,
    g_w2: pd.DataFrame,
    cols: Sequence[str],
    *,
    entity_keys: Sequence[str] = ("user_id", "item_id", "merchant_id"),
    seed: int = 0,
) -> pd.DataFrame:
    """Build order-like leaves with a simple period-shift score.

    Score = ‖x_W2 − μ_W1‖₂ on W1-standardized features (same spirit as
    three-step order score). Always keep entity keys for dissolve.
    """
    from sklearn.preprocessing import StandardScaler

    cols = list(cols)
    keys = [k for k in entity_keys if k in g_w1.columns or k in g_w2.columns]
    sc = StandardScaler().fit(g_w1.loc[:, cols].to_numpy(float))
    X1 = sc.transform(g_w1.loc[:, cols].to_numpy(float))
    mu = X1.mean(axis=0)
    # prefer W2 rows as leaves (late-window mass)
    g = g_w2.copy()
    X2 = sc.transform(g.loc[:, cols].to_numpy(float))
    score = np.linalg.norm(X2 - mu, axis=1)
    out = pd.DataFrame({"shift_l2": score})
    for k in keys:
        if k in g.columns:
            out[k] = g[k].to_numpy()
    # synthetic order id
    uid = out["user_id"].astype(str) if "user_id" in out.columns else "0"
    iid = out["item_id"].astype(str) if "item_id" in out.columns else np.arange(len(out)).astype(str)
    out["order_id"] = uid + "_" + iid + "_" + pd.Series(np.arange(len(out))).astype(str)
    out["y_convert"] = g["y_convert"].to_numpy(int) if "y_convert" in g.columns else 0
    return out
