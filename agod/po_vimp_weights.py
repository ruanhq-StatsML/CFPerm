"""Tunable PO weighting: quantile IPTW + LOCO-PO VIMP feature sampling.

Default streaming strategy remains **uniform**.
After OnlineRFPerm significance only, optional post-hoc knobs:

1. Instance weights from PO-risk **quantiles** (empirical CDF rank).
2. Feature-specific LOCO PO-risk → probability a feature is sampled
   in each RF tree (VIMP-reweighted RF).
"""
from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.tree import DecisionTreeRegressor

from agod.po_iptw import po_iptw_weights
from agod.po_refit import build_t01_windows, refit_po_risk_t01


def quantile_po_weights(
    po: np.ndarray,
    *,
    scheme: str = "cdf",
    floor: float = 0.25,
    ceil: float = 4.0,
    power: float = 1.0,
    eps: float = 1e-6,
) -> np.ndarray:
    """Map PO-risk → sample weights via within-batch quantiles.

    Schemes
    -------
    cdf     : w_i = floor + (ceil-floor) * F̂(PO_i)^power
              (empirical rank quantile; soft, tunable)
    topq    : w_i = ceil if PO_i >= q_{1-α} else floor  (α=power if in (0,1) else 0.2)
    softcap : w_i = 1 + (ceil-1) * σ( z_i ) with z = rank-quantile centered
    """
    po = np.asarray(po, float).ravel()
    po = np.maximum(po, eps)
    n = len(po)
    # mid-rank empirical CDF in (0,1]
    order = np.argsort(np.argsort(po))
    q = (order + 1.0) / (n + 1.0)

    if scheme == "cdf":
        w = floor + (ceil - floor) * np.power(q, power)
    elif scheme == "topq":
        alpha = power if 0.0 < power < 1.0 else 0.2
        thr = np.quantile(po, 1.0 - alpha)
        w = np.where(po >= thr, ceil, floor).astype(float)
    elif scheme == "softcap":
        z = (q - 0.5) * 6.0
        sig = 1.0 / (1.0 + np.exp(-z))
        w = 1.0 + (ceil - 1.0) * sig
    else:
        raise ValueError(f"unknown quantile scheme {scheme!r}")

    w = w / (w.mean() + eps)
    return np.clip(w, floor, ceil)


def loco_po_vimp(
    X0: np.ndarray,
    y0: np.ndarray,
    X1: np.ndarray,
    y1: np.ndarray,
    *,
    seed: int = 0,
    n_estimators: int = 40,
) -> np.ndarray:
    """Feature-specific LOCO PO-risk on T=1 under μ0 fit on T=0.

    LOCO_j = mean|Y1 − μ0^{−j}(X1)| − mean|Y1 − μ0(X1)|
    (positive ⇒ feature j helps reduce PO-risk on the shifted batch).

    Returns non-negative importance vector (len = d), floored at eps then
    normalized to a probability simplex for RF feature sampling.
    """
    X0 = np.asarray(X0, float)
    X1 = np.asarray(X1, float)
    y0 = np.asarray(y0, float).ravel()
    y1 = np.asarray(y1, float).ravel()
    n0, d = X0.shape
    rng = np.random.default_rng(seed)

    def _fit(X, y, s):
        rf = RandomForestRegressor(
            n_estimators=n_estimators,
            max_depth=6,
            min_samples_leaf=2,
            random_state=s,
            n_jobs=1,
        )
        rf.fit(X, y)
        return rf

    mu_full = _fit(X0, y0, seed)
    base = float(np.mean(np.abs(y1 - mu_full.predict(X1))))
    imp = np.zeros(d, float)
    for j in range(d):
        # LOCO: drop column j (set to train mean → remove signal)
        X0m = X0.copy()
        X1m = X1.copy()
        fill = float(X0[:, j].mean())
        X0m[:, j] = fill
        X1m[:, j] = fill
        mu_j = _fit(X0m, y0, seed + 17 * (j + 1))
        risk_j = float(np.mean(np.abs(y1 - mu_j.predict(X1m))))
        imp[j] = max(risk_j - base, 0.0)

    # if all zero (no measurable LOCO), fall back to flat
    if imp.sum() <= 1e-12:
        imp = np.ones(d, float)
    # mild smoothing so no feature is hard-zeroed
    imp = imp + 0.05 * imp.mean()
    return imp / imp.sum()


class VimpFeatureRF:
    """Bagged trees that sample features with P(j) ∝ LOCO-PO VIMP."""

    def __init__(
        self,
        *,
        n_estimators: int = 40,
        max_depth: int = 8,
        min_samples_leaf: int = 2,
        max_features: float = 0.5,
        seed: int = 0,
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.min_samples_leaf = min_samples_leaf
        self.max_features = max_features
        self.seed = seed
        self.trees_: List[DecisionTreeRegressor] = []
        self.feat_idx_: List[np.ndarray] = []
        self.n_features_in_: int = 0

    def fit(
        self,
        X: np.ndarray,
        y: np.ndarray,
        feature_probs: np.ndarray,
        sample_weight: Optional[np.ndarray] = None,
    ) -> "VimpFeatureRF":
        X = np.asarray(X, float)
        y = np.asarray(y, float).ravel()
        n, d = X.shape
        self.n_features_in_ = d
        p = np.asarray(feature_probs, float).ravel()
        p = np.maximum(p, 1e-12)
        p = p / p.sum()
        mtry = max(1, int(round(self.max_features * d))) if self.max_features <= 1 else int(self.max_features)
        mtry = min(mtry, d)
        rng = np.random.default_rng(self.seed)
        sw = None if sample_weight is None else np.asarray(sample_weight, float).ravel()

        self.trees_, self.feat_idx_ = [], []
        for b in range(self.n_estimators):
            # bootstrap rows
            idx = rng.integers(0, n, size=n)
            # sample features w/o replacement ∝ VIMP
            feats = rng.choice(d, size=mtry, replace=False, p=p)
            feats = np.sort(feats)
            tree = DecisionTreeRegressor(
                max_depth=self.max_depth,
                min_samples_leaf=self.min_samples_leaf,
                random_state=self.seed + b,
            )
            xb = X[idx][:, feats]
            yb = y[idx]
            if sw is None:
                tree.fit(xb, yb)
            else:
                tree.fit(xb, yb, sample_weight=sw[idx])
            self.trees_.append(tree)
            self.feat_idx_.append(feats)
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        X = np.asarray(X, float)
        preds = np.zeros(len(X), float)
        for tree, feats in zip(self.trees_, self.feat_idx_):
            preds += tree.predict(X[:, feats])
        return preds / max(len(self.trees_), 1)


def gated_quantile_weights_from_stream(
    stream: list,
    t: int,
    *,
    seed: int = 0,
    n_control: int = 1,
    window_mode: str = "recent_ood",
    scheme: str = "cdf",
    floor: float = 0.25,
    ceil: float = 4.0,
    power: float = 1.0,
) -> np.ndarray:
    """Recent/OOD re-fit PO on current batch → quantile weights."""
    X0, y0, X1, y1 = build_t01_windows(
        stream, t, n_control=n_control, window_mode=window_mode  # type: ignore[arg-type]
    )
    po1 = refit_po_risk_t01(X0, y0, X1, y1, seed=seed)
    n_cur = len(stream[t][1])
    po_cur = po1[-n_cur:]
    return quantile_po_weights(po_cur, scheme=scheme, floor=floor, ceil=ceil, power=power)


def gated_loco_vimp_from_stream(
    stream: list,
    t: int,
    *,
    seed: int = 0,
    n_control: int = 1,
    window_mode: str = "recent_ood",
) -> np.ndarray:
    """LOCO-PO feature probabilities from recent/OOD cut at time t."""
    X0, y0, X1, y1 = build_t01_windows(
        stream, t, n_control=n_control, window_mode=window_mode  # type: ignore[arg-type]
    )
    return loco_po_vimp(X0, y0, X1, y1, seed=seed)
