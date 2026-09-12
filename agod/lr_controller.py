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


# ---------------------------------------------------------------------------
# Per-modality LR *schedulers* (L2 maps over the online-window stream)
# All keep FWD on; they only reshape next-stage step sizes.
# ---------------------------------------------------------------------------

SCHEDULER_NAMES = (
    "equal",        # dense equal LR (B1)
    "soft",         # α → LR_m (current default)
    "soft_cosine",  # soft × shared cosine decay over windows
    "soft_budget",  # soft with fixed excess-LR budget Σ(LR-β)
    "soft_warmup",  # soft × linear warmup over early windows
    "damp",         # high-α → *lower* LR (covariate damp / inverse)
    "soft_entropy", # blend soft↔equal by α-entropy (uncertain → flatter)
    "soft_gradcos", # soft × gradient-alignment gain (cosine similarity)
)


def _with_shared(lr: dict[str, float], mods: Sequence[str]) -> dict[str, float]:
    out = {m: float(lr[m]) for m in mods}
    out["shared"] = float(np.mean([out[m] for m in mods]))
    return out


def equal_lr(mods: Sequence[str]) -> dict[str, float]:
    return _with_shared({m: 1.0 for m in mods}, mods)


def soft_cosine_lr(
    alpha: Mapping[str, float],
    mods: Sequence[str],
    *,
    t: int,
    t_max: int,
    beta: float = 0.10,
    floor: float = 0.35,
    gain: Mapping[str, float] | None = None,
) -> dict[str, float]:
    """Soft LR × shared cosine envelope (same attribution mass, less late overshoot)."""
    base = alpha_to_lr(alpha, mods, beta=beta, gain=gain)
    # t in [0, t_max-1] → cos from 1 → floor
    u = float(t) / max(int(t_max) - 1, 1)
    env = floor + (1.0 - floor) * 0.5 * (1.0 + np.cos(np.pi * u))
    return _with_shared({m: float(base[m] * env) for m in mods}, mods)


def soft_budget_lr(
    alpha: Mapping[str, float],
    mods: Sequence[str],
    *,
    beta: float = 0.10,
    budget: float | None = None,
    gain: Mapping[str, float] | None = None,
) -> dict[str, float]:
    """Keep Σ_m (LR_m − β) = budget (default = |M|·(1−β) = equal total excess).

    Redistributes a *fixed* adapt spend according to α — fairer FLOPs-matched
    comparison than unconstrained soft LR.
    """
    inv = float(len(mods))
    budget = float(inv * (1.0 - beta) if budget is None else budget)
    gain = gain or {m: 1.0 for m in mods}
    # unnormalized excess ∝ α · gain
    raw = np.array([max(float(alpha[m]), 0.0) * float(gain[m]) for m in mods], float)
    if raw.sum() <= 1e-12:
        excess = np.full(len(mods), budget / inv)
    else:
        excess = budget * (raw / raw.sum())
    out = {m: float(beta + excess[i]) for i, m in enumerate(mods)}
    return _with_shared(out, mods)


def soft_warmup_lr(
    alpha: Mapping[str, float],
    mods: Sequence[str],
    *,
    t: int,
    warmup: int = 2,
    beta: float = 0.10,
    gain: Mapping[str, float] | None = None,
) -> dict[str, float]:
    """Soft LR × min(1, (t+1)/warmup) — blunt early noisy α."""
    base = alpha_to_lr(alpha, mods, beta=beta, gain=gain)
    w = min(1.0, float(t + 1) / max(int(warmup), 1))
    # warmup blends toward equal rather than scaling all to 0
    return _with_shared(
        {m: float((1.0 - w) * 1.0 + w * base[m]) for m in mods},
        mods,
    )


def damp_lr(
    alpha: Mapping[str, float],
    mods: Sequence[str],
    *,
    beta: float = 0.10,
    power: float = 2.0,
    gain: Mapping[str, float] | None = None,
) -> dict[str, float]:
    """Inverse actuator: high routing mass → *smaller* steps (covariate damp).

    LR_m = β + (1−β)·(1−α_m)^p · |M| · gain_m
    Use when α mostly tracks covariate / nuisance rather than concept.
    """
    inv = float(len(mods))
    gain = gain or {m: 1.0 for m in mods}
    out = {
        m: float(
            (beta + (1.0 - beta) * ((1.0 - float(alpha[m])) ** power) * inv)
            * float(gain[m])
        )
        for m in mods
    }
    return _with_shared(out, mods)


def soft_entropy_lr(
    alpha: Mapping[str, float],
    mods: Sequence[str],
    *,
    beta: float = 0.10,
    gain: Mapping[str, float] | None = None,
) -> dict[str, float]:
    """Blend soft↔equal by normalized α-entropy (uncertain routing → flatter LR)."""
    soft = alpha_to_lr(alpha, mods, beta=beta, gain=gain)
    p = np.array([max(float(alpha[m]), 1e-12) for m in mods], float)
    p = p / p.sum()
    h = float(-(p * np.log(p)).sum())
    h_max = float(np.log(len(mods)))
    # u=0 peaked → full soft; u=1 uniform → equal
    u = float(np.clip(h / max(h_max, 1e-8), 0.0, 1.0))
    return _with_shared(
        {m: float((1.0 - u) * soft[m] + u * 1.0) for m in mods},
        mods,
    )


def cos_sim(a: np.ndarray, b: np.ndarray) -> float:
    """Cosine similarity of two flat vectors; 0 if either is ~0."""
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    if na < 1e-12 or nb < 1e-12:
        return 0.0
    return float(np.dot(a, b) / (na * nb))


def flatten_grads(params) -> np.ndarray:
    """Concatenate parameter gradients (zeros if missing)."""
    chunks = []
    for p in params:
        if p is None:
            continue
        g = getattr(p, "grad", None)
        if g is None:
            chunks.append(np.zeros(p.numel(), dtype=np.float64))
        else:
            chunks.append(g.detach().float().cpu().numpy().ravel())
    if not chunks:
        return np.zeros(1, dtype=np.float64)
    return np.concatenate(chunks)


def common_dim_grad_signature(params) -> np.ndarray:
    """Build a common-length grad signature across differently shaped modality towers.

    For each parameter tensor, reduce to a 1-D signature by averaging over all
    but the leading output dimension (bias kept as-is). Signatures are then
    L2-normalized and concatenated. This lets cos(g_m, g_m') be well-defined
    when input dims differ (e.g. video 768 vs audio 512 → shared FUSE).
    """
    parts = []
    for p in params:
        if p is None:
            continue
        g = getattr(p, "grad", None)
        if g is None:
            arr = np.zeros(p.shape, dtype=np.float64)
        else:
            arr = g.detach().float().cpu().numpy()
        if arr.ndim == 0:
            parts.append(np.array([float(arr)], dtype=np.float64))
        elif arr.ndim == 1:
            parts.append(arr.astype(np.float64))
        else:
            # weight (out, in, ...): mean over non-leading dims → (out,)
            axes = tuple(range(1, arr.ndim))
            parts.append(arr.mean(axis=axes).astype(np.float64))
    if not parts:
        return np.zeros(1, dtype=np.float64)
    return np.concatenate(parts)


def aligned_cos_sim(a: np.ndarray, b: np.ndarray) -> float:
    """Cosine sim with zero-pad to equal length (after common-dim signatures)."""
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    n = max(len(a), len(b))
    if len(a) < n:
        a = np.pad(a, (0, n - len(a)))
    if len(b) < n:
        b = np.pad(b, (0, n - len(b)))
    return cos_sim(a, b)


def modality_grad_cosine(
    grads: Mapping[str, np.ndarray],
    mods: Sequence[str],
    *,
    shared: np.ndarray | None = None,
) -> dict:
    """Characterize per-modality gradient geometry via cosine similarity.

    ``grads[m]`` / ``shared`` should preferably be *common-dim signatures*
    (see ``common_dim_grad_signature``) so towers with different input sizes
    remain comparable. Falls back to zero-pad if lengths still differ.
    """
    vecs = {m: np.asarray(grads[m], dtype=float).ravel() for m in mods}
    pair = {}
    vals = []
    for i, mi in enumerate(mods):
        for mj in mods[i + 1 :]:
            c = aligned_cos_sim(vecs[mi], vecs[mj])
            pair[f"{mi}|{mj}"] = c
            vals.append(c)
    mean_pair = float(np.mean(vals)) if vals else 0.0
    frac_neg = float(np.mean([v < 0.0 for v in vals])) if vals else 0.0

    cos_shared = {}
    if shared is not None:
        sh = np.asarray(shared, dtype=float).ravel()
        for m in mods:
            cos_shared[m] = aligned_cos_sim(vecs[m], sh)
        align = {m: float(0.5 * (1.0 + cos_shared[m])) for m in mods}
    else:
        for m in mods:
            others = [aligned_cos_sim(vecs[m], vecs[o]) for o in mods if o != m]
            cos_shared[m] = float(np.mean(others)) if others else 0.0
        align = {m: float(0.5 * (1.0 + cos_shared[m])) for m in mods}

    return {
        "pair_cos": pair,
        "mean_pair_cos": mean_pair,
        "frac_conflict": frac_neg,
        "cos_to_shared": cos_shared,
        "align_gain": align,
    }


def soft_gradcos_lr(
    alpha: Mapping[str, float],
    mods: Sequence[str],
    *,
    align_gain: Mapping[str, float],
    beta: float = 0.10,
    lambda_align: float = 0.75,
    gain: Mapping[str, float] | None = None,
) -> dict[str, float]:
    """Soft α→LR modulated by gradient cosine alignment.

    gain_m = (1-λ) + λ · align_gain_m,  align_gain_m = ½(1+cos(g_m, g_shared))

    Aligned modality grads get larger next-stage steps; conflicting ones shrink
    toward the β floor — without hard-dropping FWD.
    """
    base = alpha_to_lr(alpha, mods, beta=beta, gain=gain)
    out = {}
    for m in mods:
        a = float(np.clip(align_gain.get(m, 0.5), 0.0, 1.0))
        scale = (1.0 - lambda_align) + lambda_align * a
        out[m] = float(base[m] * scale)
    return _with_shared(out, mods)


def schedule_modality_lr(
    name: str,
    alpha: Mapping[str, float],
    mods: Sequence[str],
    *,
    t: int = 0,
    t_max: int = 6,
    beta: float = 0.10,
    gain: Mapping[str, float] | None = None,
    warmup: int = 2,
    cosine_floor: float = 0.35,
    damp_power: float = 2.0,
    align_gain: Mapping[str, float] | None = None,
    lambda_align: float = 0.75,
) -> dict[str, float]:
    """Dispatch per-modality LR scheduler by name."""
    if name not in SCHEDULER_NAMES:
        raise ValueError(f"unknown scheduler {name!r}; expected {SCHEDULER_NAMES}")
    if name == "equal":
        return equal_lr(mods)
    if name == "soft":
        return alpha_to_lr(alpha, mods, beta=beta, gain=gain)
    if name == "soft_cosine":
        return soft_cosine_lr(
            alpha, mods, t=t, t_max=t_max, beta=beta, floor=cosine_floor, gain=gain
        )
    if name == "soft_budget":
        return soft_budget_lr(alpha, mods, beta=beta, gain=gain)
    if name == "soft_warmup":
        return soft_warmup_lr(alpha, mods, t=t, warmup=warmup, beta=beta, gain=gain)
    if name == "damp":
        return damp_lr(alpha, mods, beta=beta, power=damp_power, gain=gain)
    if name == "soft_gradcos":
        if align_gain is None:
            align_gain = {m: 0.5 for m in mods}
        return soft_gradcos_lr(
            alpha,
            mods,
            align_gain=align_gain,
            beta=beta,
            lambda_align=lambda_align,
            gain=gain,
        )
    # soft_entropy
    return soft_entropy_lr(alpha, mods, beta=beta, gain=gain)


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
