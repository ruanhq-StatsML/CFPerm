"""LR actuator: Softmax routing + optional concept-intensity gain."""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np


def z_norm(d: Mapping[str, float], mods: Sequence[str]) -> dict[str, float]:
    """Non-negative then simplex-normalize across modalities."""
    v = np.array([max(float(d[m]), 0.0) for m in mods], float)
    if v.sum() <= 1e-12:
        return {m: 1.0 / len(mods) for m in mods}
    v = v / v.sum()
    return {m: float(v[i]) for i, m in enumerate(mods)}


def softmax_scores(
    d: Mapping[str, float], mods: Sequence[str], tau: float
) -> dict[str, float]:
    z = np.array([d[m] for m in mods], float) / max(tau, 1e-6)
    z = z - z.max()
    e = np.exp(z)
    e = e / e.sum()
    return {m: float(e[i]) for i, m in enumerate(mods)}


def alpha_to_lr(
    alpha: Mapping[str, float],
    mods: Sequence[str],
    *,
    beta: float = 0.10,
    gain: Mapping[str, float] | None = None,
) -> dict[str, float]:
    """Map routing mass → per-modality LR multipliers (shared = mean)."""
    inv = float(len(mods))
    gain = gain or {m: 1.0 for m in mods}
    out = {
        m: float((beta + (1.0 - beta) * float(alpha[m]) * inv) * float(gain[m]))
        for m in mods
    }
    out["shared"] = float(np.mean([out[m] for m in mods]))
    return out


def intensity_gain(
    decomp: Mapping[str, Mapping[str, float]],
    mods: Sequence[str],
    *,
    kappa: float = 1.25,
) -> dict[str, float]:
    """Optional B4: larger steps when concept dominates covariate."""
    return {
        m: float(1.0 + kappa * max(decomp["con"][m] - decomp["cov"][m], 0.0))
        for m in mods
    }


class EMARouter:
    """Online EMA over Softmax α (control-plane state)."""

    def __init__(self, mods: Sequence[str], *, ema: float = 0.40):
        self.mods = list(mods)
        self.ema = float(ema)
        self.state = {m: 1.0 / len(self.mods) for m in self.mods}

    def update(self, raw: Mapping[str, float]) -> dict[str, float]:
        for m in self.mods:
            self.state[m] = self.ema * self.state[m] + (1.0 - self.ema) * float(raw[m])
        s = sum(self.state.values())
        return {m: self.state[m] / s for m in self.mods}
