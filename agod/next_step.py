"""Proportion prediction + next-step adapt planning.

Given α history and current sensors, forecast next-window attribution
proportions and fold landscape gains into the LR actuator.

  â_{t+1} = EMA(α_t) + residual from sensor delta
  LR_m    = α_to_lr(â) ⊙ landscape_gain_m
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from .lr_controller import alpha_to_lr, softmax_scores, z_norm


def predict_next_alpha(
    alpha_hist: Sequence[Mapping[str, float]],
    mods: Sequence[str],
    *,
    sensors_now: Mapping[str, Mapping[str, float]] | None = None,
    sensors_prev: Mapping[str, Mapping[str, float]] | None = None,
    ema: float = 0.55,
    tau: float = 0.35,
    sensor_key: str = "g",
) -> dict[str, float]:
    """Predict next α proportions from history + optional sensor delta.

    If only one snapshot: return that α. Else EMA of last two, with a
    Softmax residual from sensor deltas (PO/disc/proto score ``g``).
    """
    mods = list(mods)
    if not alpha_hist:
        return {m: 1.0 / len(mods) for m in mods}
    cur = {m: float(alpha_hist[-1][m]) for m in mods}
    if len(alpha_hist) == 1:
        base = cur
    else:
        prev = {m: float(alpha_hist[-2][m]) for m in mods}
        base = {m: ema * prev[m] + (1.0 - ema) * cur[m] for m in mods}

    if sensors_now is None or sensor_key not in sensors_now:
        s = sum(base.values())
        return {m: base[m] / s for m in mods}

    now = sensors_now[sensor_key]
    if sensors_prev is not None and sensor_key in sensors_prev:
        prv = sensors_prev[sensor_key]
        delta = {m: float(now.get(m, 0.0)) - float(prv.get(m, 0.0)) for m in mods}
    else:
        delta = {m: float(now.get(m, 0.0)) for m in mods}

    # residual routing mass from sensor delta
    res = softmax_scores(delta, mods, tau)
    mix = {m: 0.70 * base[m] + 0.30 * res[m] for m in mods}
    s = sum(mix.values())
    return {m: float(mix[m] / s) for m in mods}


def next_step_lr(
    alpha_hat: Mapping[str, float],
    mods: Sequence[str],
    *,
    landscape_gain: Mapping[str, float] | None = None,
    beta: float = 0.10,
) -> dict[str, float]:
    """Actuator: predicted α → LR, optionally × landscape gain."""
    gain = None
    if landscape_gain is not None:
        gain = {m: float(landscape_gain.get(m, 1.0)) for m in mods}
    return alpha_to_lr(alpha_hat, mods, beta=beta, gain=gain)


def proportion_report(
    alpha: Mapping[str, float],
    alpha_hat: Mapping[str, float],
    mods: Sequence[str],
) -> dict:
    """Telemetry for 3-mod attribution proportions + forecast error."""
    mods = list(mods)
    a = z_norm(alpha, mods)
    h = z_norm(alpha_hat, mods)
    mae = float(np.mean([abs(a[m] - h[m]) for m in mods]))
    return {
        "alpha": a,
        "alpha_hat_next": h,
        "prop_mae": mae,
        "top_mod": max(mods, key=lambda m: a[m]),
        "top_hat": max(mods, key=lambda m: h[m]),
    }
