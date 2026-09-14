"""Observation-level PO-risk → sample weights (hard-reweight).

After OnlineRFPerm reject on batch O:

  1. Fit μ0 on recent control R
  2. PO_i = |Y_i − μ0(X_i)|  (+ optional μ-gap blend) on O
  3. w_i = transform(PO_i)
  4. Re-fit downstream learner on O with sample_weight=w

Non-reject → w=1. Hard-reweight, not an OOD detector.

Iteration note
--------------
Raw IPTW (esp. prop / √PO) often hurts next-MSE vs uniform while still
ranking hard rows well. Prefer soft maps + tempering / top-k support.
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
    "qrt",  # PO^{1/4} — softer than cbrt
    "quantile",
    "hybrid",
    "topk",  # boost only top-q hard rows
]


def _mean1_clip(w: np.ndarray, clip: tuple[float, float], eps: float) -> np.ndarray:
    w = np.asarray(w, float).ravel()
    w = w / (w.mean() + eps)
    lo, hi = clip
    return np.clip(w, lo, hi)


def temper_weights(
    w: np.ndarray,
    *,
    lam: float = 0.5,
    clip: tuple[float, float] = (0.05, 20.0),
    eps: float = 1e-6,
) -> np.ndarray:
    """Mix PO weights with uniform: w' = (1−λ)·1 + λ·w  (then mean-1).

    λ=0 → uniform; λ=1 → full PO map. Softens MSE blow-ups while keeping
    relative hard-row emphasis.
    """
    lam = float(np.clip(lam, 0.0, 1.0))
    w = np.asarray(w, float).ravel()
    mixed = (1.0 - lam) * np.ones_like(w) + lam * w
    return _mean1_clip(mixed, clip, eps)


def topk_boost_weights(
    po: np.ndarray,
    *,
    frac: float = 0.2,
    boost: float = 2.0,
    base: float = 1.0,
    clip: tuple[float, float] = (0.05, 20.0),
    eps: float = 1e-6,
) -> np.ndarray:
    """Only the hardest top-``frac`` rows get ``boost``; others stay ``base``."""
    po = np.asarray(po, float).ravel()
    n = len(po)
    k = max(1, int(round(n * frac)))
    w = np.full(n, base, float)
    top = np.argpartition(po, -k)[-k:]
    w[top] = boost
    return _mean1_clip(w, clip, eps)


def obs_po_to_weights(
    po: np.ndarray,
    mode: ObsWeightMode = "sqrt",
    *,
    power: float | None = None,
    q_floor: float = 0.25,
    q_ceil: float = 4.0,
    q_power: float = 1.0,
    hybrid_power: float = 0.5,
    temper: float = 1.0,
    topk_frac: float = 0.2,
    topk_boost: float = 2.0,
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
    qrt      : w ∝ PO^{1/4}
    quantile : within-batch CDF rank (soft, bounded)
    hybrid   : F̂(PO)^{q_power} · PO^{hybrid_power}
    topk     : only top-``topk_frac`` hard rows get ``topk_boost``

    ``temper`` ∈ [0,1] mixes the chosen map with uniform after shaping
    (1 = full map, 0 = uniform). Default 1 keeps backward compatibility.
    """
    po = np.maximum(np.asarray(po, float).ravel(), eps)
    if mode == "uniform":
        return np.ones_like(po)
    if mode == "qrt":
        w = po_iptw_weights(po, mode="sqrt", power=0.25 if power is None else power, clip=clip, eps=eps)
    elif mode in ("prop", "sqrt", "cbrt"):
        w = po_iptw_weights(po, mode=mode, power=power, clip=clip, eps=eps)  # type: ignore[arg-type]
    elif mode == "quantile":
        w = quantile_po_weights(
            po, scheme="cdf", floor=q_floor, ceil=q_ceil, power=q_power, eps=eps
        )
    elif mode == "hybrid":
        q = quantile_po_weights(
            po, scheme="cdf", floor=q_floor, ceil=q_ceil, power=q_power, eps=eps
        )
        mag = po_iptw_weights(po, mode="sqrt", power=hybrid_power, clip=clip, eps=eps)
        w = _mean1_clip(q * mag, clip, eps)
    elif mode == "topk":
        w = topk_boost_weights(
            po, frac=topk_frac, boost=topk_boost, clip=clip, eps=eps
        )
    else:
        raise ValueError(f"unknown obs PO weight mode {mode!r}")

    if temper < 1.0 - 1e-12:
        w = temper_weights(w, lam=temper, clip=clip, eps=eps)
    return w


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
