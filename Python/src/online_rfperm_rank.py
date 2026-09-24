"""OnlineRFPerm ranking probe. Information retrieval, not quantile regression.

Frozen scorer on D_ref. Each batch is a set of slates. The rent check is
NDCG@k vs the hold-out slice of D_ref, not MSE, not pinball.

T is the batch label. Y is graded relevance, last column, never a feature.
`groups` is the slate id. It is not an X column.

    rec = onlinePermOOB_rank(df, groups, k=5, ref_batch_size=..., batch_size=...)
"""
from __future__ import annotations

import numpy as np

from online_rfperm import fit_frozen_rf, probe_mse
from streaming_po_risk import large_deviation, mmd_vs_reference, rbf_bandwidth


def ndcg_at_k(y, score, k: int) -> float:
    y = np.asarray(y, dtype=float).ravel()
    score = np.asarray(score, dtype=float).ravel()
    n = len(y)
    if n == 0:
        return 0.0
    kk = min(int(k), n)
    order = np.argsort(-score, kind="mergesort")
    gain = (2.0 ** y[order][:kk] - 1.0) / np.log2(np.arange(kk) + 2.0)
    dcg = float(gain.sum())
    ideal = np.sort(y)[::-1][:kk]
    idcg = float(((2.0 ** ideal - 1.0) / np.log2(np.arange(kk) + 2.0)).sum())
    if idcg <= 1e-12:
        return 0.0
    return dcg / idcg


def mean_ndcg(y, score, groups, k: int) -> float:
    y = np.asarray(y, dtype=float).ravel()
    score = np.asarray(score, dtype=float).ravel()
    groups = np.asarray(groups)
    vals = []
    for g in np.unique(groups):
        m = groups == g
        if int(m.sum()) < 2:
            continue
        vals.append(ndcg_at_k(y[m], score[m], k))
    if not vals:
        return float("nan")
    return float(np.mean(vals))


def _first(flags) -> int:
    for t, v in enumerate(flags):
        if v:
            return int(t)
    return -1


def onlinePermOOB_rank(
    df,
    groups,
    k=5,
    ref_batch_size=160,
    batch_size=40,
    seed=2026,
    gate=2.0,
):
    """Frozen RF scorer. Loss = 1 − NDCG@k, gated vs hold-out D_ref.

    df last column is Y. groups aligns with rows and is not in X.
    """
    df = np.asarray(df, dtype=float)
    groups = np.asarray(groups)
    if df.ndim != 2 or df.shape[1] < 2:
        raise ValueError("df needs X columns and Y last")
    if len(groups) != len(df):
        raise ValueError("groups must align with df rows")
    X, Y = df[:, :-1], df[:, -1]
    n_ref = int(ref_batch_size)
    bs = int(batch_size)
    kk = int(k)
    cut = max(int(0.7 * n_ref), 16)
    scorer = fit_frozen_rf(X[:cut], Y[:cut], seed=seed)
    ndcg_ref = mean_ndcg(Y[cut:n_ref], scorer.predict(X[cut:n_ref]), groups[cut:n_ref], kk)
    loss_ref = 1.0 - float(ndcg_ref)
    mse_ref = probe_mse(scorer, X[cut:n_ref], Y[cut:n_ref])
    sigma = rbf_bandwidth(X[:n_ref], seed=seed)
    rng = np.random.default_rng(int(seed) + 5)
    i1 = rng.choice(n_ref, size=min(40, n_ref), replace=True)
    i2 = rng.choice(n_ref, size=min(40, n_ref), replace=True)
    mmd_base = max(float(mmd_vs_reference(X[i1], X[i2], sigma=sigma, seed=seed)), 0.015)

    ndcg, mse, mmd = [], [], []
    hop_n, hop_m, hop_x = [], [], []
    i = n_ref
    while i < len(X):
        Xb, Yb, gb = X[i : i + bs], Y[i : i + bs], groups[i : i + bs]
        if len(Xb) < 8:
            break
        pred = np.asarray(scorer.predict(Xb), dtype=float)
        n_t = mean_ndcg(Yb, pred, gb, kk)
        loss = 1.0 - float(n_t)
        m_t = probe_mse(scorer, Xb, Yb)
        x_t = float(mmd_vs_reference(X[:n_ref], Xb, sigma=sigma, seed=seed))
        hop_n.append(bool(large_deviation(loss, max(loss_ref, 1e-3), ratio=float(gate))))
        hop_m.append(bool(large_deviation(m_t, mse_ref, ratio=float(gate))))
        hop_x.append(bool(large_deviation(x_t, mmd_base, ratio=float(gate))))
        ndcg.append(float(n_t))
        mse.append(float(m_t))
        mmd.append(x_t)
        i += bs

    return {
        "k": kk,
        "ndcg_ref": float(ndcg_ref),
        "mse_ref": float(mse_ref),
        "NDCG_list": np.asarray(ndcg, dtype=float),
        "MSE_list": np.asarray(mse, dtype=float),
        "MMD_list": np.asarray(mmd, dtype=float),
        "hop_ndcg": np.asarray(hop_n, dtype=bool),
        "hop_mse": np.asarray(hop_m, dtype=bool),
        "hop_x": np.asarray(hop_x, dtype=bool),
        "ndcg_1": _first(hop_n),
        "mse_1": _first(hop_m),
        "x_1": _first(hop_x),
        "y_in_x": False,
    }


def make_rank_df(n_queries=60, slate=8, p=6, kind="quiet", onset_q=30, seed=0):
    """Slates of graded relevance. Last column Y. groups returned separately.

    quiet      : frozen ranking
    rank_flip  : after onset, X0 sign flips in the relevance rule
    covariate  : unused column walks; ranking rule frozen
    """
    rng = np.random.default_rng(int(seed))
    n_q = int(n_queries)
    L = int(slate)
    p = int(p)
    X = rng.normal(size=(n_q * L, p))
    groups = np.repeat(np.arange(n_q), L)
    w = np.zeros(p)
    w[0], w[1] = 1.6, 0.7
    qid = groups
    latent = X @ w
    if kind == "rank_flip":
        w2 = w.copy()
        w2[0] = -1.6
        latent = np.where(qid >= int(onset_q), X @ w2, X @ w)
    if kind == "covariate":
        X[qid >= int(onset_q), p - 1] += 2.8
        latent = X @ w
    Y = np.zeros(n_q * L, dtype=float)
    for q in range(n_q):
        m = np.where(qid == q)[0]
        order = np.argsort(-latent[m], kind="mergesort")
        rel = np.zeros(L, dtype=float)
        rel[order[:2]] = 2.0
        rel[order[2:4]] = 1.0
        Y[m] = rel
    return np.column_stack([X, Y]), groups
