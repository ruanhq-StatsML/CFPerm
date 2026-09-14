"""Observation-level PO-risk → sample weights (hard-reweight).

After OnlineRFPerm reject on batch O:

  1. Fit μ0 on recent control R
  2. PO_i = |Y_i − μ0(X_i)|  (+ optional μ-gap blend) on O
  3. w_i = transform(PO_i)   — prop / √PO / ∛PO / quantile / hybrid
  4. Re-fit downstream learner on O with sample_weight=w

Non-reject → w=1. Hard-reweight, not an OOD detector.
"""
from __future__ import annotations

from typing import Literal

import numpy as np

from agod.online_rfperm import gate_blend
from agod.po_iptw import po_iptw_weights
from agod.po_vimp_weights import quantile_po_weights

ObsWeightMode = Literal[
    "uniform",
    "prop",
    "sqrt",
    "cbrt",
    "quantile",
    "hybrid",
]


def obs_po_to_weights(
    po: np.ndarray,
    mode: ObsWeightMode = "sqrt",
    *,
    power: float | None = None,
    q_floor: float = 0.25,
    q_ceil: float = 4.0,
    q_power: float = 1.0,
    hybrid_power: float = 0.5,
    clip: tuple[float, float] = (0.05, 20.0),
    eps: float = 1e-6,
) -> np.ndarray:
    """Map observation-level PO-risk → mean-1 sample weights.

    Modes
    -----
    uniform  : w = 1
    prop     : w ∝ PO
    sqrt     : w ∝ √PO
    cbrt     : w ∝ PO^{1/3}
    quantile : w from within-batch CDF rank (soft, bounded)
    hybrid   : w ∝ F̂(PO)^{q_power} · PO^{hybrid_power}
    """
    po = np.maximum(np.asarray(po, float).ravel(), eps)
    if mode == "uniform":
        return np.ones_like(po)
    if mode in ("prop", "sqrt", "cbrt"):
        return po_iptw_weights(po, mode=mode, power=power, clip=clip, eps=eps)  # type: ignore[arg-type]
    if mode == "quantile":
        return quantile_po_weights(
            po, scheme="cdf", floor=q_floor, ceil=q_ceil, power=q_power, eps=eps
        )
    if mode == "hybrid":
        q = quantile_po_weights(
            po, scheme="cdf", floor=q_floor, ceil=q_ceil, power=q_power, eps=eps
        )
        mag = po_iptw_weights(po, mode="sqrt", power=hybrid_power, clip=clip, eps=eps)
        w = q * mag
        w = w / (w.mean() + eps)
        lo, hi = clip
        return np.clip(w, lo, hi)
    raise ValueError(f"unknown obs PO weight mode {mode!r}")


def gated_obs_po_weights(
    po: np.ndarray,
    *,
    reject: bool,
    mode: ObsWeightMode = "sqrt",
    soft: bool = False,
    p: float = 1.0,
    alpha: float = 0.05,
    **kwargs,
) -> np.ndarray:
    """Uniform unless reject (or soft-blend by p-value)."""
    ones = np.ones(len(np.asarray(po).ravel()), float)
    if not reject and not soft:
        return ones
    w = obs_po_to_weights(po, mode=mode, **kwargs)
    return gate_blend(reject, w, soft=soft, p=p, alpha=alpha)
