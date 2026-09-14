"""Observation-level PO-risk → sample weights (hard-reweight).

After OnlineRFPerm reject on batch O:

  1. Fit μ0 on recent control R
  2. PO_i = |Y_i − μ0(X_i)|  (+ optional μ-gap blend) on O
  3. w_i = transform(PO_i)
  4. Re-fit downstream learner on O with sample_weight=w

Non-reject → w=1. Hard-reweight, not an OOD detector.

Logic we buy (v4)
-----------------
- **Hard-rank always**: obs PO is a hardness score (Spearman / P@k).
  Mechanism that matches this claim: put mass on the hard support
  (``hard_support`` / top-k), not diffuse soft IPTW.
- **Pack MSE only under beijing-class drift**: when hard-tail ≈ shift
  signal, soft temper can lift all-row next-MSE; on calm packs hard ≈
  noise and uniform wins. Encode with a *higher* drift gate for the
  pack-MSE path; do not chase pack MSE on mild rejects.
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
    "hard_support",  # top-k boost scaled by λ (hard-claim mechanism)
]

# Beijing-class drift gate for pack-MSE path (stricter than mild reject).
BEIJING_DRIFT_GATE = 0.45


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


def hard_support_weights(
    po: np.ndarray,
    *,
    frac: float = 0.2,
    boost_max: float = 3.0,
    lam: float = 1.0,
    clip: tuple[float, float] = (0.05, 20.0),
    eps: float = 1e-6,
) -> np.ndarray:
    """Hard-claim mechanism: mass only on top-``frac`` PO rows.

    ``boost = 1 + λ·(boost_max − 1)``. λ=0 → uniform; λ=1 → full hard boost.
    Easy rows stay at base weight 1 (before mean-1 renormalization).
    """
    lam = float(np.clip(lam, 0.0, 1.0))
    boost = 1.0 + lam * (float(boost_max) - 1.0)
    if boost <= 1.0 + 1e-12:
        return np.ones(len(np.asarray(po).ravel()), float)
    return topk_boost_weights(po, frac=frac, boost=boost, clip=clip, eps=eps)


def is_beijing_class_drift(drift: float, *, gate: float = BEIJING_DRIFT_GATE) -> bool:
    """True when drift looks like shift-signal (beijing-class), not mild noise."""
    return float(drift) > float(gate)


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
    boost_max: float = 3.0,
    clip: tuple[float, float] = (0.05, 20.0),
    eps: float = 1e-6,
) -> np.ndarray:
    """Map observation-level PO-risk → mean-1 sample weights.

    Modes
    -----
    uniform      : w = 1
    prop         : w ∝ PO
    sqrt         : w ∝ √PO
    cbrt         : w ∝ PO^{1/3}
    qrt          : w ∝ PO^{1/4}
    quantile     : within-batch CDF rank (soft, bounded)
    hybrid       : F̂(PO)^{q_power} · PO^{hybrid_power}
    topk         : only top-``topk_frac`` hard rows get ``topk_boost``
    hard_support : top-k boost scaled by ``temper`` as λ (hard-claim path)

    ``temper`` ∈ [0,1] mixes soft maps with uniform after shaping
    (1 = full map, 0 = uniform). For ``hard_support``, ``temper`` is λ
    on the hard boost (no second temper mix).
    """
    po = np.maximum(np.asarray(po, float).ravel(), eps)
    if mode == "uniform":
        return np.ones_like(po)
    if mode == "hard_support":
        return hard_support_weights(
            po, frac=topk_frac, boost_max=boost_max, lam=temper, clip=clip, eps=eps
        )
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


def drift_intensity(
    po: np.ndarray,
    *,
    control_resid: np.ndarray | None = None,
    p: float | None = None,
    T: float | None = None,
    alpha: float = 0.05,
    eps: float = 1e-6,
) -> float:
    """Scalar drift / shift intensity in [0, 1].

    Primary signal (beijing-like shift)
    -----------------------------------
    ``mean(PO_ood) / mean(|Y−μ0| on recent control) − 1``
    — hard mass relative to the control regime. Calm packs sit near 0;
    drifted packs go large.

    Secondary: within-OOD hard-tail q90/q50, RFPerm p-strength, T>0.
    High ⇒ hard rows look like **shift signal**; low ⇒ hard ≈ noise.
    """
    po = np.asarray(po, float).ravel()
    if len(po) == 0:
        return 0.0

    if control_resid is not None and len(control_resid) > 0:
        cr = np.asarray(control_resid, float).ravel()
        gap = float(np.mean(po) / (np.mean(np.abs(cr)) + eps) - 1.0)
        # gap≈0 → matched; gap≳1 → OOD twice as hard as control
        gap_n = float(np.clip(gap / 1.0, 0.0, 1.0))
    else:
        q50 = float(np.quantile(po, 0.5))
        q90 = float(np.quantile(po, 0.9))
        gap_n = float(np.clip((q90 / (q50 + eps) - 1.0) / 2.0, 0.0, 1.0))

    p_strength = 0.0
    if p is not None:
        p_strength = float(np.clip((alpha - float(p)) / max(alpha, eps), 0.0, 1.0))

    t_strength = 0.0
    if T is not None:
        t_strength = float(np.tanh(max(float(T), 0.0) / 10.0))

    # gap dominates; p/T are tie-breakers
    return float(np.clip(0.70 * gap_n + 0.20 * p_strength + 0.10 * t_strength, 0.0, 1.0))


def adaptive_temper(
    drift: float,
    *,
    lam_max: float = 0.75,
    drift_gate: float = 0.20,
) -> float:
    """Map drift intensity → temper λ.

    Below ``drift_gate`` → λ=0 (keep uniform even on reject — mild/noise).
    Above gate → ramp to ``lam_max``.

    Use ``drift_gate≈0.2`` for the hard-support path (mild fire OK).
    Use ``BEIJING_DRIFT_GATE≈0.45`` for the pack-MSE path (only shift packs).
    """
    d = float(np.clip(drift, 0.0, 1.0))
    if d <= drift_gate:
        return 0.0
    return float(lam_max * (d - drift_gate) / max(1.0 - drift_gate, 1e-6))


def hard_subset_mask(score: np.ndarray, *, frac: float = 0.2) -> np.ndarray:
    """Boolean mask of the hardest top-``frac`` rows by ``score``."""
    s = np.asarray(score, float).ravel()
    n = len(s)
    k = max(1, int(round(n * frac)))
    idx = np.argpartition(s, -k)[-k:]
    m = np.zeros(n, dtype=bool)
    m[idx] = True
    return m
