#!/usr/bin/env python3
"""Pseudo-outcome (PO) risk helpers for graph-feature FSDS.

Locked framing (讲武德):
  - W = period (W1/W2), not treatment — no ATE claim
  - Fit PO **once** on a support; use τ̂² / VIMP to help selection
  - Does not replace official FSDS; feeds it

Usage from iterate / localize scripts:

  from po_risk_fsds import fit_period_po, po_feature_table, po_help_select

DS one-liner (period shift → ranking prior → FSDS pool):

  report = po_help_select(g_w1, g_w2, cols, k=15)
  # report["pool"] → hand to official FSDS; report["feature_table"] for notebooks
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler


def _pca_rank_y(X: np.ndarray, seed: int) -> np.ndarray:
    z = PCA(n_components=1, random_state=seed).fit_transform(X).ravel()
    order = np.argsort(z, kind="mergesort")
    y = np.empty(len(z), float)
    y[order] = np.linspace(0.0, 1.0, len(z))
    return y


def fit_period_po(
    X: np.ndarray,
    W: np.ndarray,
    *,
    seed: int = 0,
    n_trees: int = 40,
    y: Optional[np.ndarray] = None,
) -> Dict:
    """DR-style PO → τ̂(X); risk = mean(τ̂²). Same pattern as repo localize scripts.

    Parameters
    ----------
    X : (n, d) features (prefer W1-standardized graph feats)
    W : (n,) period in {0,1}
    y : optional continuous outcome; default = PCA-rank of X (rare-label safe)
    """
    X = np.asarray(X, float)
    W = np.asarray(W, int).ravel()
    n = len(W)
    if y is None:
        Yf = _pca_rank_y(X, seed)
    else:
        Yf = np.asarray(y, float).ravel()

    m_hat = np.zeros(n)
    e_hat = np.zeros(n)
    n_splits = 3 if n >= 60 else 2
    # Stratified on W (period), not on rare convert labels
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for fold, (tr, te) in enumerate(cv.split(X, W)):
        m = RandomForestRegressor(
            n_estimators=n_trees,
            max_depth=6,
            min_samples_leaf=4,
            random_state=seed + fold,
            n_jobs=1,
        )
        e = RandomForestClassifier(
            n_estimators=n_trees,
            max_depth=6,
            min_samples_leaf=4,
            random_state=seed + 40 + fold,
            n_jobs=1,
        )
        m.fit(X[tr], Yf[tr])
        e.fit(X[tr], W[tr])
        m_hat[te] = m.predict(X[te])
        e_hat[te] = e.predict_proba(X[te])[:, 1]
    e_hat = np.clip(e_hat, 0.05, 0.95)
    po = (Yf - m_hat) * (W - e_hat)
    tau = RandomForestRegressor(
        n_estimators=max(60, n_trees),
        max_depth=8,
        min_samples_leaf=4,
        random_state=seed + 7,
        n_jobs=1,
    )
    tau.fit(X, po)
    tau_hat = tau.predict(X)
    tau2 = tau_hat ** 2
    return {
        "risk": float(np.mean(tau2)),
        "tau_hat": tau_hat.astype(float),
        "tau2": tau2.astype(float),
        "po": po.astype(float),
        "vimp": tau.feature_importances_.astype(float),
        "tau_model": tau,
        "n": int(n),
        "note": "W=period; PO-risk is a shift proxy, not an ATE",
    }


def po_feature_table(cols: Sequence[str], vimp: np.ndarray) -> pd.DataFrame:
    tab = (
        pd.DataFrame({"feature": list(cols), "po_vimp": np.asarray(vimp, float)})
        .sort_values("po_vimp", ascending=False)
        .reset_index(drop=True)
    )
    tab["rank"] = np.arange(1, len(tab) + 1)
    s = tab["po_vimp"].to_numpy()
    tab["po_share"] = s / (s.sum() + 1e-12)
    return tab


def fit_po_on_windows(
    g_w1: pd.DataFrame,
    g_w2: pd.DataFrame,
    cols: Sequence[str],
    *,
    seed: int = 0,
    max_n: int = 8000,
) -> Dict:
    """Pool W1/W2 rows, W1-fit scaler, fit period-PO once."""
    rng = np.random.default_rng(seed)
    n1 = min(len(g_w1), max_n // 2)
    n2 = min(len(g_w2), max_n // 2)
    i1 = rng.choice(len(g_w1), size=n1, replace=False)
    i2 = rng.choice(len(g_w2), size=n2, replace=False)
    X1 = g_w1.iloc[i1].loc[:, list(cols)].to_numpy(float)
    X2 = g_w2.iloc[i2].loc[:, list(cols)].to_numpy(float)
    sc = StandardScaler().fit(X1)
    Xp = np.vstack([sc.transform(X1), sc.transform(X2)])
    Wp = np.concatenate([np.zeros(n1, dtype=int), np.ones(n2, dtype=int)])
    fit = fit_period_po(Xp, Wp, seed=seed)
    fit["scaler"] = sc
    fit["cols"] = list(cols)
    fit["feature_table"] = po_feature_table(cols, fit["vimp"])
    fit["n_w1"] = int(n1)
    fit["n_w2"] = int(n2)
    return fit


def tau2_on_frame(
    g: pd.DataFrame,
    cols: Sequence[str],
    po_fit: Dict,
) -> np.ndarray:
    """Score rows with τ̂² using a fitted period-PO (scaler + τ forest)."""
    X = g.loc[:, list(cols)].to_numpy(float)
    Xs = po_fit["scaler"].transform(X)
    model = po_fit.get("tau_model")
    if model is not None:
        return (model.predict(Xs) ** 2).astype(float)
    # fallback: VIMP-weighted squared standardized features
    v = np.asarray(po_fit["vimp"], float)
    v = v / (v.sum() + 1e-12)
    return (Xs ** 2 @ v).astype(float)


def filter_by_tau2_quantile(
    g: pd.DataFrame,
    cols: Sequence[str],
    po_fit: Dict,
    *,
    q: float = 0.5,
    y_col: str = "y_convert",
) -> pd.DataFrame:
    """Keep rows with τ̂²-proxy ≥ quantile q; always keep all positives.

    Period-shift mass filter for FSDS selection (讲武德: not an ATE weight).
    """
    scores = tau2_on_frame(g, cols, po_fit)
    thr = float(np.quantile(scores, q))
    y = g[y_col].to_numpy(int) if y_col in g.columns else np.zeros(len(g), int)
    keep = (scores >= thr) | (y > 0)
    out = g.loc[keep].copy()
    out.attrs["tau2_q"] = q
    out.attrs["tau2_thr"] = thr
    out.attrs["tau2_kept"] = int(keep.sum())
    return out


def blend_cmean_po_scores(
    cols: Sequence[str],
    cmean_abs: np.ndarray,
    po_vimp: np.ndarray,
    *,
    alpha: float = 0.5,
) -> pd.DataFrame:
    """Blend normalized |δ| and PO-VIMP for feature ranking (α in [0,1])."""
    c = np.asarray(cmean_abs, float)
    p = np.asarray(po_vimp, float)
    c = c / (c.max() + 1e-12)
    p = p / (p.max() + 1e-12)
    score = alpha * p + (1.0 - alpha) * c
    tab = (
        pd.DataFrame(
            {
                "feature": list(cols),
                "cmean_abs_norm": c,
                "po_vimp_norm": p,
                "blend": score,
            }
        )
        .sort_values("blend", ascending=False)
        .reset_index(drop=True)
    )
    tab["rank"] = np.arange(1, len(tab) + 1)
    return tab


def rare_pos_n_splits(y: np.ndarray, *, prefer: int = 5) -> int:
    """Cap CV folds by positive count (rare-convert regime)."""
    n_pos = int(np.sum(np.asarray(y).ravel() > 0))
    if n_pos <= 1:
        return 2
    return max(2, min(prefer, n_pos))


def bootstrap_pi_select(
    X: np.ndarray,
    y: np.ndarray,
    cols: Sequence[str],
    *,
    k: int,
    n_boot: int = 40,
    seed: int = 0,
    score_fn=f_classif,
) -> Tuple[List[str], pd.DataFrame]:
    """Stratified bootstrap π-stability (rare-pos friendly vs K-fold).

    Each draw keeps class balance; π_j = fraction of boots where j ∈ TopK.
    With tiny n_pos this still runs (unlike 5-fold when a fold has 0 pos).
    """
    X = np.asarray(X, float)
    y = np.asarray(y, int).ravel()
    cols = list(cols)
    n, d = X.shape
    assert d == len(cols)
    rng = np.random.default_rng(seed)
    pos = np.where(y == 1)[0]
    neg = np.where(y == 0)[0]
    hits = np.zeros(d, dtype=float)
    score_sum = np.zeros(d, dtype=float)
    n_ok = 0
    kk = min(k, d, max(1, n - 1))
    for b in range(n_boot):
        if len(pos) == 0 or len(neg) == 0:
            break
        # draw all positives with replacement; same count of negatives
        i_pos = rng.choice(pos, size=len(pos), replace=True)
        i_neg = rng.choice(neg, size=len(neg), replace=True)
        idx = np.concatenate([i_pos, i_neg])
        yb = y[idx]
        if len(np.unique(yb)) < 2:
            continue
        sel = SelectKBest(score_fn, k=min(kk, len(idx) - 1))
        try:
            sel.fit(X[idx], yb)
        except Exception:
            continue
        hits += sel.get_support().astype(float)
        sc = np.nan_to_num(sel.scores_, nan=0.0)
        score_sum += sc
        n_ok += 1
    if n_ok == 0:
        # fallback: single F on full data
        sel = SelectKBest(score_fn, k=kk)
        sel.fit(X, y)
        selected = [cols[j] for j, m in enumerate(sel.get_support()) if m]
        tab = (
            pd.DataFrame(
                {
                    "feature": cols,
                    "pi": sel.get_support().astype(float),
                    "mean_score": np.nan_to_num(sel.scores_, nan=0.0),
                }
            )
            .sort_values(["pi", "mean_score"], ascending=False)
            .reset_index(drop=True)
        )
        tab["rank"] = np.arange(1, len(tab) + 1)
        tab["selected"] = tab["feature"].isin(selected).astype(int)
        tab["n_boot_ok"] = 0
        return selected, tab
    pi = hits / n_ok
    mean_sc = score_sum / n_ok
    order = sorted(range(d), key=lambda j: (-pi[j], -mean_sc[j]))
    selected = [cols[j] for j in order[:k]]
    tab = (
        pd.DataFrame({"feature": cols, "pi": pi, "mean_score": mean_sc})
        .sort_values(["pi", "mean_score"], ascending=False)
        .reset_index(drop=True)
    )
    tab["rank"] = np.arange(1, len(tab) + 1)
    tab["selected"] = tab["feature"].isin(selected).astype(int)
    tab["n_boot_ok"] = int(n_ok)
    return selected, tab


def po_help_select(
    g_w1: pd.DataFrame,
    g_w2: pd.DataFrame,
    cols: Sequence[str],
    *,
    k: int = 15,
    seed: int = 0,
    alpha: float = 0.5,
    max_n: int = 6000,
    pool_extra: int = 3,
) -> Dict:
    """DS entrypoint: fit PO once, blend with |δ|, return FSDS pool + tables.

    Returns
    -------
    dict with keys:
      risk, feature_table, blend_table, pool, note
    """
    cols = list(cols)
    po = fit_po_on_windows(g_w1, g_w2, cols, seed=seed, max_n=max_n)
    X1 = g_w1.loc[:, cols].to_numpy(float)
    X2 = g_w2.loc[:, cols].to_numpy(float)
    dlt = np.abs(X2.mean(axis=0) - X1.mean(axis=0))
    blend = blend_cmean_po_scores(cols, dlt, po["vimp"], alpha=alpha)
    pre_n = min(len(cols), max(k, k + pool_extra))
    pool = list(blend["feature"].head(pre_n))
    return {
        "risk": po["risk"],
        "feature_table": po["feature_table"],
        "blend_table": blend,
        "pool": pool,
        "k": int(k),
        "n_w1": po["n_w1"],
        "n_w2": po["n_w2"],
        "note": (
            "W=period; PO-risk/VIMP is a shift ranking prior for FSDS — "
            "not an ATE. Hand `pool` to official Scaler→Var→SelectKBest."
        ),
    }
