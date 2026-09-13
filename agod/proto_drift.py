"""Modality prototype bank + prototype-conditioned drift.

Prototypes = EMA centroids of modality features (optionally class-conditional).
Drift logic:
  d_m = 1 − cos(proto_m, mean(batch_m))     # how far current cloud moved
  Coupled with FSDS: high PO + high proto-drift ⇒ true concept move;
  high MMD + low proto-drift ⇒ covariate mush around a stable centroid.

L0 sensor; does not skip FWD — only reshapes adapt spend via α.
"""
from __future__ import annotations

from typing import Mapping, MutableMapping, Sequence

import numpy as np


def _mean_row(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, float)
    if X.ndim == 1:
        return X
    return X.mean(axis=0)


def _cos(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float).ravel()
    b = np.asarray(b, float).ravel()
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    if na < 1e-12 or nb < 1e-12:
        return 0.0
    return float(np.dot(a, b) / (na * nb))


class ModalityPrototypeBank:
    """EMA centroids per modality (and optional class labels)."""

    def __init__(self, mods: Sequence[str], *, ema: float = 0.85):
        self.mods = list(mods)
        self.ema = float(ema)
        self.proto: dict[str, np.ndarray | None] = {m: None for m in self.mods}
        self.class_proto: dict[str, dict[int, np.ndarray]] = {m: {} for m in self.mods}
        self.n_updates = 0

    def update(
        self,
        blocks: Mapping[str, np.ndarray],
        y: np.ndarray | None = None,
    ) -> None:
        for m in self.mods:
            mu = _mean_row(blocks[m])
            if self.proto[m] is None:
                self.proto[m] = mu.copy()
            else:
                self.proto[m] = self.ema * self.proto[m] + (1.0 - self.ema) * mu
            if y is not None:
                y = np.asarray(y).astype(int).ravel()
                X = np.asarray(blocks[m], float)
                for c in np.unique(y):
                    mask = y == c
                    if mask.sum() < 2:
                        continue
                    cmu = X[mask].mean(axis=0)
                    bank = self.class_proto[m]
                    if c not in bank:
                        bank[c] = cmu.copy()
                    else:
                        bank[c] = self.ema * bank[c] + (1.0 - self.ema) * cmu
        self.n_updates += 1

    def drift(
        self,
        blocks: Mapping[str, np.ndarray],
        y: np.ndarray | None = None,
    ) -> dict[str, float]:
        """Per-mod prototype drift in [0, 2] typically; 0 = aligned."""
        out = {}
        for m in self.mods:
            if self.proto[m] is None:
                out[m] = 0.0
                continue
            mu = _mean_row(blocks[m])
            d = 1.0 - _cos(self.proto[m], mu)
            # class-conditional add-on: mean drift of present classes
            if y is not None and self.class_proto[m]:
                yv = np.asarray(y).astype(int).ravel()
                X = np.asarray(blocks[m], float)
                cds = []
                for c in np.unique(yv):
                    if c not in self.class_proto[m]:
                        continue
                    mask = yv == c
                    if mask.sum() < 2:
                        continue
                    cds.append(1.0 - _cos(self.class_proto[m][c], X[mask].mean(0)))
                if cds:
                    d = 0.5 * d + 0.5 * float(np.mean(cds))
            out[m] = float(max(d, 0.0))
        return out


def compose_proto_fsds(
    po: Mapping[str, float],
    mmd: Mapping[str, float],
    proto_d: Mapping[str, float],
    disc: Mapping[str, float],
    mods: Sequence[str],
    *,
    w_po: float = 1.0,
    w_mmd: float = 0.75,
    w_proto: float = 0.75,
    w_disc: float = 0.50,
) -> dict[str, float]:
    """FSDS-style score: concept-ish − covariate-ish, with proto/disc.

    g_m ∝ w_po·PO + w_disc·disc + w_proto·proto_drift − w_mmd·MMD
    (then caller Softmax / EMA → α).
    """
    mods = list(mods)
    raw = {
        m: float(w_po) * max(float(po.get(m, 0.0)), 0.0)
        + float(w_disc) * max(float(disc.get(m, 0.0)), 0.0)
        + float(w_proto) * max(float(proto_d.get(m, 0.0)), 0.0)
        - float(w_mmd) * max(float(mmd.get(m, 0.0)), 0.0)
        for m in mods
    }
    # shift to non-neg then normalize
    mn = min(raw.values())
    raw = {m: raw[m] - mn for m in mods}
    s = sum(raw.values())
    if s <= 1e-12:
        return {m: 1.0 / len(mods) for m in mods}
    return {m: float(raw[m] / s) for m in mods}
