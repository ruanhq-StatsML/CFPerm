#!/usr/bin/env python3
"""Pseudo-outcome (PO) risk helpers for graph-feature FSDS.

Locked framing (讲武德):
  - W = period (W1/W2), not treatment — no ATE claim
  - Fit PO **once** on a support; use τ̂² / VIMP to help selection
  - Does not replace official FSDS; feeds it

Usage from iterate / localize scripts:

  from po_risk_fsds import fit_period_po, po_feature_table
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
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
