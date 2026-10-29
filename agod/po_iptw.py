"""Green IPTW-style sample weights from continuous-batch PO-risk.

PO-risk as OOD score (preferred empirically vs DRE):

  w = 1              # uniform
  w ∝ PO             # prop  — hard upweight high-risk
  w ∝ sqrt(PO)       # soft  — w_i = √PO(X_i,Y_i,T_i=1)
  w ∝ 1 / PO         # inv
  w ∝ p/(1-p)        # dre  — logistic density-ratio baseline

Normalize to mean 1 so RF / LR scales stay comparable.
"""
from __future__ import annotations

from typing import Literal

import numpy as np

WeightMode = Literal["uniform", "prop", "sqrt", "inv", "dre"]


def po_iptw_weights(
    po: np.ndarray,
    mode: WeightMode = "sqrt",
    *,
    eps: float = 1e-6,
    clip: tuple[float, float] = (0.05, 20.0),
    treated: np.ndarray | None = None,
) -> np.ndarray:
    """Map PO-risk → sample weights.

    If ``treated`` is given (T_i∈{0,1}), only T_i=1 rows get PO-shaped
    weights; untreated stay at 1 (then re-normalize).
    """
    po = np.asarray(po, dtype=float).ravel()
    po = np.maximum(po, eps)
    if mode == "uniform":
        w = np.ones_like(po)
    elif mode == "prop":
        w = po.copy()
    elif mode == "sqrt":
        # w_i = sqrt(PO-risk(X_i, Y_i, T_i=1))
        w = np.sqrt(po)
    elif mode == "inv":
        w = 1.0 / po
    elif mode == "dre":
        raise ValueError("mode='dre' needs dre_weights(X_ref, X_cur), not PO")
    else:
        raise ValueError(f"unknown mode {mode!r}")
    if treated is not None:
        t = np.asarray(treated, dtype=float).ravel()
        w = np.where(t > 0.5, w, 1.0)
    w = w / (w.mean() + eps)
    lo, hi = clip
    return np.clip(w, lo, hi)


def dre_weights(
    X_ref: np.ndarray,
    X_cur: np.ndarray,
    *,
    eps: float = 1e-6,
    clip: tuple[float, float] = (0.05, 20.0),
    seed: int = 0,
) -> np.ndarray:
    """Logistic density-ratio baseline: w ∝ p(cur|x) / p(ref|x)."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    Xr = np.asarray(X_ref, float)
    Xc = np.asarray(X_cur, float)
    X = np.vstack([Xr, Xc])
    Xs = StandardScaler().fit_transform(X)
    y = np.concatenate([np.zeros(len(Xr)), np.ones(len(Xc))])
    clf = LogisticRegression(max_iter=400, random_state=seed)
    clf.fit(Xs, y)
    p = clf.predict_proba(Xs[len(Xr) :])[:, 1]
    w = p / np.maximum(1.0 - p, eps)
    w = w / (w.mean() + eps)
    lo, hi = clip
    return np.clip(w, lo, hi)


def instance_po_risk(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    *,
    batch_po: float | None = None,
    mix: float = 0.5,
) -> np.ndarray:
    """Observation PO proxy = |residual| (optionally blended with batch PO)."""
    resid = np.abs(np.asarray(y_true, float) - np.asarray(y_pred, float)).ravel()
    if batch_po is None:
        return resid
    bp = max(float(batch_po), 0.0)
    return (1.0 - mix) * resid + mix * (resid * (1.0 + bp))
