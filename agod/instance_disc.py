"""Nonparametric instance discrimination (feature-learning sensor).

Cheap, model-free discrimination between reference and current batches
without a deep NCE head:

  score_m = NN purity / LOO domain accuracy on modality block m

High score ⇒ instances are separable across the window boundary
(feature-level shift worth adapting). Low score ⇒ overlapping clouds
(noise / covariate mush — do not boost).

Fits L0 telemetry: returns a per-modality [0,1] sensor usable like PO.
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler


def _l2_rows(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, float)
    n = np.linalg.norm(X, axis=1, keepdims=True)
    return X / np.clip(n, 1e-8, None)


def nn_domain_purity(
    X_ref: np.ndarray,
    X_cur: np.ndarray,
    *,
    k: int = 5,
    seed: int = 0,
) -> float:
    """kNN domain purity: fraction of neighbors sharing the query's domain.

    Ref=0, Cur=1. Returns mean purity in [0.5, 1] roughly (chance≈0.5 when
    clouds overlap; →1 when domains are locally pure).
    """
    del seed  # sklearn NN is deterministic given data
    X0 = np.asarray(X_ref, float)
    X1 = np.asarray(X_cur, float)
    if len(X0) < k + 1 or len(X1) < k + 1:
        return 0.5
    X = np.vstack([X0, X1])
    y = np.array([0] * len(X0) + [1] * len(X1), int)
    Xs = StandardScaler().fit_transform(X)
    Xs = _l2_rows(Xs)
    nn = NearestNeighbors(n_neighbors=k + 1, algorithm="auto")
    nn.fit(Xs)
    idx = nn.kneighbors(Xs, return_distance=False)[:, 1:]  # drop self
    same = (y[idx] == y[:, None]).mean(axis=1)
    return float(np.clip(same.mean(), 0.0, 1.0))


def nn_discrimination_scores(
    blocks_ref: Mapping[str, np.ndarray],
    blocks_cur: Mapping[str, np.ndarray],
    mods: Sequence[str],
    *,
    k: int = 5,
    seed: int = 0,
) -> dict[str, float]:
    """Per-modality nonparametric discrimination sensor."""
    return {
        m: nn_domain_purity(blocks_ref[m], blocks_cur[m], k=k, seed=seed + i)
        for i, m in enumerate(mods)
    }


def disc_as_concept(
    disc: Mapping[str, float],
    mods: Sequence[str],
    *,
    floor: float = 0.5,
) -> dict[str, float]:
    """Map purity → non-negative concept-like mass (0 if ≤ chance)."""
    return {m: float(max(float(disc[m]) - floor, 0.0)) for m in mods}
