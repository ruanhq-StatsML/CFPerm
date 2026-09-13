"""Gradient-projection memory + landscape eigenstructure.

Store last-K common-dim grad signatures per modality. Current step:
  proj residual energy = ||g − P_span(g)|| / ||g||
  landscape = Gram eigenstructure over the memory bank (erank, top_frac)

Fusion with FSDS:
  - high residual ⇒ new direction worth learning (boost adapt)
  - low residual + low erank ⇒ revisiting old basin (damp / adapter-only)

Reuses ``common_dim_grad_signature`` / ``gram_effective_rank`` geometry.
"""
from __future__ import annotations

from collections import deque
from typing import Mapping, Sequence

import numpy as np

from .ensemble_decorr import gram_effective_rank
from .lr_controller import cos_sim


def _as_vec(v) -> np.ndarray:
    return np.asarray(v, float).ravel()


class GradProjMemory:
    """Per-modality FIFO of grad signatures + projection residual."""

    def __init__(self, mods: Sequence[str], *, capacity: int = 8):
        self.mods = list(mods)
        self.capacity = int(capacity)
        self.bank: dict[str, deque[np.ndarray]] = {
            m: deque(maxlen=self.capacity) for m in self.mods
        }

    def push(self, sigs: Mapping[str, np.ndarray]) -> None:
        for m in self.mods:
            if m not in sigs:
                continue
            v = _as_vec(sigs[m])
            if v.size == 0 or not np.isfinite(v).all():
                continue
            n = float(np.linalg.norm(v))
            if n < 1e-12:
                continue
            self.bank[m].append(v / n)

    def residual_energy(self, sigs: Mapping[str, np.ndarray]) -> dict[str, float]:
        """||g − Proj_span(bank) g|| / ||g|| ∈ [0,1]; 1 if empty bank."""
        out = {}
        for m in self.mods:
            g = _as_vec(sigs.get(m, np.zeros(1)))
            ng = float(np.linalg.norm(g))
            if ng < 1e-12:
                out[m] = 0.0
                continue
            g = g / ng
            hist = list(self.bank[m])
            if len(hist) < 1:
                out[m] = 1.0
                continue
            # columns = past directions (same dim; truncate/pad)
            d = g.size
            cols = []
            for h in hist:
                if h.size != d:
                    # align by truncate/pad
                    vv = np.zeros(d, float)
                    n = min(d, h.size)
                    vv[:n] = h[:n]
                    h = vv
                    nh = float(np.linalg.norm(h))
                    if nh < 1e-12:
                        continue
                    h = h / nh
                cols.append(h)
            if not cols:
                out[m] = 1.0
                continue
            A = np.stack(cols, axis=1)  # d × k
            # least-squares proj onto col(A)
            try:
                coef, *_ = np.linalg.lstsq(A, g, rcond=None)
                proj = A @ coef
            except Exception:
                out[m] = 1.0
                continue
            r = g - proj
            out[m] = float(np.clip(np.linalg.norm(r), 0.0, 1.0))
        return out

    def landscape(self) -> dict:
        """Eigenstructure of cross-mod cos Gram from *latest* bank tips."""
        tips = {}
        for m in self.mods:
            if self.bank[m]:
                tips[m] = self.bank[m][-1]
        if len(tips) < 2:
            return {
                "effective_rank": float(len(self.mods)),
                "participation_ratio": float(len(self.mods)),
                "top_frac": 0.0,
                "pair_cos": {},
            }
        pair = {}
        mods = [m for m in self.mods if m in tips]
        for i, a in enumerate(mods):
            for b in mods[i + 1 :]:
                pair[f"{a}||{b}"] = float(cos_sim(tips[a], tips[b]))
        spec = gram_effective_rank(pair, mods)
        spec["pair_cos"] = pair
        return spec


def landscape_adapt_gain(
    residual: Mapping[str, float],
    landscape: Mapping[str, float],
    mods: Sequence[str],
    *,
    erank_floor: float = 1.4,
    residual_boost: float = 0.75,
    collinear_damp: float = 0.55,
) -> dict[str, float]:
    """Per-mod LR gain from residual energy × global eigenstructure.

    New residual directions → boost; collinear low-erank landscape → damp.
    """
    mods = list(mods)
    erank = float(landscape.get("effective_rank", len(mods)))
    # scale: erank near 1 ⇒ damp; erank near |M| ⇒ allow boosts
    span = max(float(len(mods)) - 1.0, 1e-6)
    erank_n = float(np.clip((erank - 1.0) / span, 0.0, 1.0))
    collinear = erank <= float(erank_floor)
    out = {}
    for m in mods:
        r = float(residual.get(m, 0.5))
        g = 1.0 + residual_boost * (r - 0.5)
        if collinear:
            g *= collinear_damp + (1.0 - collinear_damp) * erank_n
        else:
            g *= 0.85 + 0.30 * erank_n
        out[m] = float(np.clip(g, 0.35, 1.75))
    out["shared"] = float(np.mean([out[m] for m in mods]))
    out["erank"] = erank
    out["erank_n"] = erank_n
    return out
