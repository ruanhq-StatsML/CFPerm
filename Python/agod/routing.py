"""Clever-covariate routing: MSG state → dynamic modality weights."""

from __future__ import annotations

from typing import Dict, Optional, Sequence, Tuple

import numpy as np

from .chronoberg import MODALITIES
from .msg import MSGState


def softmax_weights(gaps: np.ndarray, tau: float = 0.35) -> np.ndarray:
    """alpha = softmax(g / tau). Smaller tau → sharper routing."""
    g = np.asarray(gaps, dtype=np.float64)
    tau = max(float(tau), 1e-6)
    z = g / tau
    z = z - np.max(z)
    e = np.exp(np.clip(z, -40.0, 40.0))
    s = float(e.sum())
    if not np.isfinite(s) or s <= 0:
        return np.ones_like(g) / len(g)
    return e / s


def smoothness_penalty(alpha: np.ndarray, alpha_prev: Optional[np.ndarray]) -> float:
    if alpha_prev is None:
        return 0.0
    return float(np.sum((np.asarray(alpha) - np.asarray(alpha_prev)) ** 2))


def gate_weights(
    alpha: np.ndarray,
    aucs: np.ndarray,
    *,
    theta_low: float = 0.55,
    theta_high: float = 0.80,
    decay: float = 0.45,
) -> Tuple[np.ndarray, Tuple[int, ...]]:
    """Decay stable modalities; flag high-AUC ones for post-hoc localization."""
    a = np.asarray(alpha, dtype=np.float64).copy()
    aucs = np.asarray(aucs, dtype=np.float64)
    localize: list[int] = []
    for i, auc in enumerate(aucs):
        if auc < theta_low:
            a[i] *= decay
        if auc > theta_high:
            localize.append(i)
    s = float(a.sum())
    if s <= 0:
        a = np.ones_like(a) / len(a)
    else:
        a = a / s
    return a, tuple(localize)


def ema_smooth(alpha: np.ndarray, alpha_prev: Optional[np.ndarray], momentum: float) -> np.ndarray:
    if alpha_prev is None or momentum <= 0:
        return np.asarray(alpha, dtype=np.float64)
    m = float(np.clip(momentum, 0.0, 0.95))
    blended = (1.0 - m) * np.asarray(alpha, dtype=np.float64) + m * np.asarray(alpha_prev, dtype=np.float64)
    s = float(blended.sum())
    return blended / s if s > 0 else blended


def route_from_state(
    state: MSGState,
    *,
    tau: float = 0.35,
    theta_low: float = 0.55,
    theta_high: float = 0.80,
    decay: float = 0.45,
    momentum: float = 0.35,
    alpha_prev: Optional[np.ndarray] = None,
    modalities: Sequence[str] = MODALITIES,
    baseline: str = "B3",
) -> Tuple[Dict[str, float], Tuple[str, ...], float]:
    """Map an MSG state to gated, temporally smoothed modality weights.

    Baselines
    ---------
    B1 : static uniform weights.
    B2 : covariate-only, alpha ∝ AUC.
    B3 : AGOD, alpha ∝ MSG gap.
    """
    mods = tuple(modalities)
    if baseline == "B1":
        alpha = np.ones(len(mods), dtype=np.float64) / len(mods)
        localize: Tuple[str, ...] = ()
        omega = smoothness_penalty(alpha, alpha_prev)
        return {m: float(alpha[i]) for i, m in enumerate(mods)}, localize, omega

    if baseline == "B2":
        aucs = np.array([state.details[m].auc for m in mods], dtype=np.float64)
        scores = np.clip(2.0 * aucs - 1.0, 0.0, None) + 1e-3
        alpha = softmax_weights(scores, tau=tau)
    elif baseline == "B3":
        alpha = softmax_weights(state.vector(mods), tau=tau)
    else:
        raise ValueError(f"Unknown baseline {baseline!r}; expected B1, B2 or B3.")

    aucs = np.array([state.details[m].auc for m in mods], dtype=np.float64)
    alpha, loc_idx = gate_weights(
        alpha, aucs, theta_low=theta_low, theta_high=theta_high, decay=decay
    )
    alpha = ema_smooth(alpha, alpha_prev, momentum)
    omega = smoothness_penalty(alpha, alpha_prev)
    localize_names = tuple(mods[i] for i in loc_idx)
    return {m: float(alpha[i]) for i, m in enumerate(mods)}, localize_names, omega
