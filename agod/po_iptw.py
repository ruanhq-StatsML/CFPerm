"""Green IPTW-style sample weights from continuous-batch PO-risk.

No density-ratio / DGA — just:

  w ∝ PO        # upweight high concept-risk rows
  w ∝ 1 / PO    # inverse
  w = 1         # uniform

Normalize to mean 1 so RF / LR scales stay comparable.
"""
from __future__ import annotations

from typing import Literal

import numpy as np

WeightMode = Literal["uniform", "prop", "inv"]


def po_iptw_weights(
    po: np.ndarray,
    mode: WeightMode = "prop",
    *,
    eps: float = 1e-6,
    clip: tuple[float, float] = (0.05, 20.0),
) -> np.ndarray:
    po = np.asarray(po, dtype=float).ravel()
    po = np.maximum(po, eps)
    if mode == "uniform":
        w = np.ones_like(po)
    elif mode == "prop":
        w = po.copy()
    elif mode == "inv":
        w = 1.0 / po
    else:
        raise ValueError(f"unknown mode {mode!r}")
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
