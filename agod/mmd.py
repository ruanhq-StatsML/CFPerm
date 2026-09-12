"""Unbiased RBF MMD² sensor (FSDS / AGOD covariate channel)."""
from __future__ import annotations

import numpy as np

DEFAULT_MAX_N = 96


def whiten_pair(X0: np.ndarray, X1: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Column-standardize on pooled ref+cur (stable RBF bandwidth)."""
    X0 = np.asarray(X0, float)
    X1 = np.asarray(X1, float)
    pooled = np.concatenate([X0, X1], 0)
    mu = pooled.mean(0, keepdims=True)
    sd = pooled.std(0, keepdims=True) + 1e-6
    return (X0 - mu) / sd, (X1 - mu) / sd


def rbf_mmd2(
    X0: np.ndarray,
    X1: np.ndarray,
    *,
    max_n: int = DEFAULT_MAX_N,
    seed: int = 0,
) -> float:
    """Unbiased squared MMD with RBF kernel (median heuristic).

    Same estimator family as FSDS MMD-LOCO: measures RKHS mean-embedding
    distance between two batches — the natural online ``P(X)`` sensor.
    """
    rng = np.random.default_rng(seed)
    X0 = np.asarray(X0, float)
    X1 = np.asarray(X1, float)
    if len(X0) > max_n:
        X0 = X0[rng.choice(len(X0), max_n, replace=False)]
    if len(X1) > max_n:
        X1 = X1[rng.choice(len(X1), max_n, replace=False)]
    if len(X0) < 4 or len(X1) < 4:
        return 0.0

    Z = np.vstack([X0, X1])
    idx = rng.choice(len(Z), size=min(256, len(Z)), replace=False)
    S = Z[idx]
    d2 = np.sum((S[:, None, :] - S[None, :, :]) ** 2, axis=-1)
    med = float(np.median(d2[d2 > 0])) if np.any(d2 > 0) else 1.0
    gamma = 1.0 / max(med, 1e-8)

    def k(A, B):
        dd = np.sum((A[:, None, :] - B[None, :, :]) ** 2, axis=-1)
        return np.exp(-gamma * dd)

    K00, K11, K01 = k(X0, X0), k(X1, X1), k(X0, X1)
    n0, n1 = len(X0), len(X1)
    mmd2 = (
        (K00.sum() - np.trace(K00)) / max(n0 * (n0 - 1), 1)
        + (K11.sum() - np.trace(K11)) / max(n1 * (n1 - 1), 1)
        - 2.0 * K01.mean()
    )
    return float(max(mmd2, 0.0))
