"""Recall@K, forgetting, and attribution-consistency metrics."""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
from scipy.stats import spearmanr


def _l2_normalize(x: np.ndarray, eps: float = 1e-8) -> np.ndarray:
    nrm = np.linalg.norm(x, axis=1, keepdims=True)
    return x / np.clip(nrm, eps, None)


def recall_at_k(
    query: np.ndarray,
    gallery: np.ndarray,
    *,
    k: int = 5,
    labels: np.ndarray | None = None,
) -> float:
    """Self-retrieval against a gallery, or same-label retrieval when labels are given."""
    q = _l2_normalize(np.asarray(query, dtype=np.float64))
    g = _l2_normalize(np.asarray(gallery, dtype=np.float64))
    sim = q @ g.T
    n = q.shape[0]
    k = max(1, min(int(k), n if labels is None else g.shape[0]))
    if labels is None:
        # Recover the diagonal identity: each query's true neighbor is its index.
        ranked = np.argsort(-sim, axis=1)[:, :k]
        hits = np.any(ranked == np.arange(n)[:, None], axis=1)
        return float(np.mean(hits))
    labels = np.asarray(labels)
    ranked = np.argsort(-sim, axis=1)[:, : k + 1]
    hits = []
    for i in range(n):
        neighbors = [j for j in ranked[i] if j != i][:k]
        if not neighbors:
            hits.append(0.0)
            continue
        hits.append(float(np.mean(labels[neighbors] == labels[i])))
    return float(np.mean(hits))


def forgetting_rate(scores_ref: Sequence[float], scores_later: Sequence[float]) -> float:
    """Mean drop on the stable reference task after later updates."""
    a = np.asarray(scores_ref, dtype=np.float64)
    b = np.asarray(scores_later, dtype=np.float64)
    return float(np.mean(a - b))


def attribution_consistency(
    estimated_rank: Mapping[str, float],
    true_rank: Mapping[str, float],
) -> float:
    """Spearman correlation between estimated MSG ranks and ground-truth shift mass."""
    keys = list(estimated_rank.keys())
    est = np.array([estimated_rank[k] for k in keys], dtype=np.float64)
    tru = np.array([true_rank[k] for k in keys], dtype=np.float64)
    if np.allclose(est, est[0]) or np.allclose(tru, tru[0]):
        return 0.0
    rho, _ = spearmanr(est, tru)
    return float(rho) if np.isfinite(rho) else 0.0


def cosine_alignment(a: np.ndarray, b: np.ndarray) -> float:
    a = _l2_normalize(np.asarray(a, dtype=np.float64))
    b = _l2_normalize(np.asarray(b, dtype=np.float64))
    return float(np.mean(np.sum(a * b, axis=1)))


def discretize_y(y: np.ndarray, q: float = 0.5) -> np.ndarray:
    thr = np.quantile(y, q)
    return (np.asarray(y) >= thr).astype(int)
