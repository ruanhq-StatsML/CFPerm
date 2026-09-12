"""LOCO-MMD utilities (restored for real-data overlap benchmarks)."""
from __future__ import annotations

import numpy as np


def _rank_positions(vimp):
    order = np.argsort(-np.asarray(vimp, dtype=float))
    rank_pos = np.empty_like(order)
    rank_pos[order] = np.arange(len(vimp))
    return rank_pos


def _subsample_batch(X, max_n, rng):
    X = np.asarray(X, dtype=float)
    if X.shape[0] <= max_n:
        return X
    idx = rng.choice(X.shape[0], max_n, replace=False)
    return X[idx]


def _median_bandwidth(Z, max_n=80, seed=0):
    Z = np.asarray(Z, dtype=float)
    if Z.shape[0] > max_n:
        rng = np.random.default_rng(seed)
        Z = Z[rng.choice(Z.shape[0], max_n, replace=False)]
    if Z.shape[0] < 2:
        return 1.0
    # subsample pairs rather than materializing (n, n, d)
    A2 = np.sum(Z * Z, axis=1, keepdims=True)
    dist2 = np.maximum(A2 + A2.T - 2.0 * (Z @ Z.T), 0.0)
    dist = np.sqrt(dist2)
    med = np.median(dist[dist > 0])
    return float(med) if med > 0 else 1.0


def _rbf_mmd2_unbiased(X, Y, gamma):
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    n = X.shape[0]
    m = Y.shape[0]
    if n < 2 or m < 2:
        return 0.0

    def _sq_dists(A, B):
        A2 = np.sum(A * A, axis=1, keepdims=True)
        B2 = np.sum(B * B, axis=1, keepdims=True)
        return np.maximum(A2 + B2.T - 2.0 * (A @ B.T), 0.0)

    Kxx = np.exp(-gamma * _sq_dists(X, X))
    Kyy = np.exp(-gamma * _sq_dists(Y, Y))
    Kxy = np.exp(-gamma * _sq_dists(X, Y))
    np.fill_diagonal(Kxx, 0.0)
    np.fill_diagonal(Kyy, 0.0)
    term_xx = Kxx.sum() / (n * (n - 1))
    term_yy = Kyy.sum() / (m * (m - 1))
    term_xy = Kxy.mean()
    return float(max(term_xx + term_yy - 2.0 * term_xy, 0.0))


class MMD:
    def __init__(self, compute_kernel="rbf", gamma=None):
        if compute_kernel != "rbf":
            raise ValueError("unsupported kernel: %s" % compute_kernel)
        self.compute_kernel = compute_kernel
        self.gamma = gamma

    def __call__(self, X, Y):
        X = np.asarray(X, dtype=float)
        Y = np.asarray(Y, dtype=float)
        gamma = self.gamma
        if gamma is None:
            pooled = np.vstack([X, Y])
            med = _median_bandwidth(pooled)
            gamma = 1.0 / (2.0 * med * med + 1e-12)
        return _rbf_mmd2_unbiased(X, Y, gamma), float("nan")


def compute_loco_mmd_batches(df1_X, df2_X, max_n=80, seed=0):
    rng = np.random.default_rng(seed)
    X_exist = _subsample_batch(np.asarray(df1_X, dtype=float), max_n, rng)
    X_new = _subsample_batch(np.asarray(df2_X, dtype=float), max_n, rng)
    mmd_test = MMD(compute_kernel="rbf")
    mmd_orig, _ = mmd_test(X_exist, X_new)
    p = X_exist.shape[1]
    mmd_dic = np.zeros(p, dtype=float)
    gamma = mmd_test.gamma
    if gamma is None:
        pooled = np.vstack([X_exist, X_new])
        med = _median_bandwidth(pooled, seed=seed)
        gamma = 1.0 / (2.0 * med * med + 1e-12)
    for j in range(p):
        mask = np.ones(p, dtype=bool)
        mask[j] = False
        mmd_dic[j] = _rbf_mmd2_unbiased(X_exist[:, mask], X_new[:, mask], gamma)
    mmd_vimp = mmd_orig - mmd_dic
    return mmd_vimp, _rank_positions(mmd_vimp)
