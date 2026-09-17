#!/usr/bin/env python3
"""Post-hoc subset localization: not just conditional means.

Pull subset indices, then call the two ready statistics on every pair:

    MMD(X[I_s], X[I_t])
    PO-risk on the union, with W = pair membership

That is the localization efficacy / significance, not mean(Y | subset).
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
VENDOR = ROOT / "vendor"
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))

from fsds.VIMP_mmd_benchmark import MMD  # noqa: E402

MIN_N = 12
MAX_MMD_N = 250
N_PERM = 30
CLIP = 1e-3


def zscore(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    sd = X.std(axis=0)
    sd = np.where(sd < 1e-8, 1.0, sd)
    return (X - X.mean(axis=0)) / sd


def subset_indices_from_labels(labels: np.ndarray, *, prefix: str = "S") -> dict[str, np.ndarray]:
    """Pull index lists for each label. That is the only grouping step."""
    labels = np.asarray(labels)
    out = {}
    for g in np.unique(labels):
        try:
            key = f"{prefix}{int(g)}"
        except (TypeError, ValueError):
            key = f"{prefix}{g}"
        out[key] = np.flatnonzero(labels == g)
    return out


def subset_indices_from_feature(x: np.ndarray, *, n_bins: int = 4, prefix: str = "Q") -> dict[str, np.ndarray]:
    """Quartile (or n-bin) split of one coordinate — same discretization as the board."""
    x = np.asarray(x, dtype=float).ravel()
    qs = np.quantile(x, np.linspace(0.0, 1.0, n_bins + 1)[1:-1])
    qs = np.unique(qs)
    if qs.size == 0:
        return {f"{prefix}0": np.arange(len(x))}
    bins = np.digitize(x, qs, right=True)
    return subset_indices_from_labels(bins, prefix=prefix)


def _subsample(X: np.ndarray, max_n: int, rng: np.random.Generator) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    if len(X) <= max_n:
        return X
    return X[rng.choice(len(X), max_n, replace=False)]


def mmd_pair(Xa: np.ndarray, Xb: np.ndarray, *, max_n: int = MAX_MMD_N, seed: int = 0) -> float:
    """Ready MMD() on two index-pulled blocks."""
    rng = np.random.default_rng(seed)
    mmd = MMD(compute_kernel="rbf")
    stat, _ = mmd(_subsample(Xa, max_n, rng), _subsample(Xb, max_n, rng))
    return float(stat)


def po_risk(X: np.ndarray, Y: np.ndarray, W: np.ndarray) -> float:
    """Ready PO-risk on a pair: φ = (Y-μ)(W-e), τ̂(X) ≈ φ, risk = mean(τ̂²)."""
    Xd = np.column_stack([np.ones(len(X)), zscore(X)])
    Y = np.asarray(Y, dtype=float).ravel()
    W = np.asarray(W, dtype=float).ravel()
    mu, *_ = np.linalg.lstsq(Xd, Y, rcond=None)
    e, *_ = np.linalg.lstsq(Xd, W, rcond=None)
    mu_hat = Xd @ mu
    e_hat = np.clip(Xd @ e, CLIP, 1.0 - CLIP)
    phi = (Y - mu_hat) * (W - e_hat)
    tau, *_ = np.linalg.lstsq(Xd, phi, rcond=None)
    tau_hat = Xd @ tau
    return float(np.mean(tau_hat ** 2))


def _perm_p(obs: float, null: np.ndarray) -> float:
    return float((1.0 + np.sum(null >= obs - 1e-15)) / (len(null) + 1.0))


def pairwise_subset_mmd(
    X: np.ndarray,
    subsets: dict[str, np.ndarray],
    *,
    n_perm: int = N_PERM,
    seed: int = 0,
    max_n: int = MAX_MMD_N,
) -> list[dict]:
    """MMD(X[I_s], X[I_t]) for every subset pair, with a label-permutation p-value."""
    names = [k for k, idx in subsets.items() if len(idx) >= MIN_N]
    rows = []
    rng = np.random.default_rng(seed)
    for a, b in combinations(names, 2):
        ia, ib = np.asarray(subsets[a]), np.asarray(subsets[b])
        Xa, Xb = X[ia], X[ib]
        obs = mmd_pair(Xa, Xb, max_n=max_n, seed=seed)
        pooled_idx = np.concatenate([ia, ib])
        n0 = len(ia)
        null = np.zeros(n_perm, dtype=float)
        for t in range(n_perm):
            perm = rng.permutation(len(pooled_idx))
            take = pooled_idx[perm]
            null[t] = mmd_pair(X[take[:n0]], X[take[n0:]], max_n=max_n, seed=seed + 1 + t)
        rows.append(
            {
                "a": a,
                "b": b,
                "n_a": int(len(ia)),
                "n_b": int(len(ib)),
                "mmd": obs,
                "mmd_p": _perm_p(obs, null),
            }
        )
    return rows


def pairwise_subset_po_risk(
    X: np.ndarray,
    Y: np.ndarray,
    subsets: dict[str, np.ndarray],
    *,
    n_perm: int = N_PERM,
    seed: int = 0,
) -> list[dict]:
    """PO-risk on I_s ∪ I_t with W = pair membership. Permute W for the p-value."""
    names = [k for k, idx in subsets.items() if len(idx) >= MIN_N]
    Y = np.asarray(Y, dtype=float).ravel()
    rows = []
    rng = np.random.default_rng(seed)
    for a, b in combinations(names, 2):
        ia, ib = np.asarray(subsets[a]), np.asarray(subsets[b])
        idx = np.concatenate([ia, ib])
        W = np.concatenate([np.zeros(len(ia), dtype=int), np.ones(len(ib), dtype=int)])
        Xi, Yi = X[idx], Y[idx]
        obs = po_risk(Xi, Yi, W)
        null = np.zeros(n_perm, dtype=float)
        for t in range(n_perm):
            null[t] = po_risk(Xi, Yi, rng.permutation(W))
        rows.append(
            {
                "a": a,
                "b": b,
                "n_a": int(len(ia)),
                "n_b": int(len(ib)),
                "po_risk": obs,
                "po_p": _perm_p(obs, null),
                "mean_Y_a": float(Y[ia].mean()),
                "mean_Y_b": float(Y[ib].mean()),
            }
        )
    return rows


def merge_pairwise(mmd_rows: list[dict], po_rows: list[dict]) -> list[dict]:
    by_po = {(r["a"], r["b"]): r for r in po_rows}
    out = []
    for m in mmd_rows:
        rec = dict(m)
        rec.update({k: v for k, v in by_po.get((m["a"], m["b"]), {}).items() if k not in rec})
        out.append(rec)
    return out


def localize(
    X: np.ndarray,
    Y: np.ndarray,
    *,
    group_labels: np.ndarray | None = None,
    feature: np.ndarray | None = None,
    n_bins: int = 4,
    n_perm: int = N_PERM,
    seed: int = 0,
) -> dict:
    """Means are kept. The decision numbers are pairwise MMD and PO-risk."""
    blocks = {}
    if group_labels is not None:
        blocks["groups"] = subset_indices_from_labels(group_labels, prefix="T")
    if feature is not None:
        blocks["bins"] = subset_indices_from_feature(feature, n_bins=n_bins, prefix="Q")
    out = {}
    for name, subsets in blocks.items():
        means = {
            k: {
                "n": int(len(idx)),
                "mean_Y": float(np.mean(Y[idx])) if len(idx) else None,
                "mean_X": float(np.mean(X[idx])) if len(idx) and X.ndim == 1 else None,
            }
            for k, idx in subsets.items()
        }
        mmd_rows = pairwise_subset_mmd(X if X.ndim == 2 else X.reshape(-1, 1), subsets, n_perm=n_perm, seed=seed)
        po_rows = pairwise_subset_po_risk(X if X.ndim == 2 else X.reshape(-1, 1), Y, subsets, n_perm=n_perm, seed=seed + 7)
        pairs = merge_pairwise(mmd_rows, po_rows)
        n_sig = sum(int(r.get("mmd_p", 1) <= 0.05 or r.get("po_p", 1) <= 0.05) for r in pairs)
        out[name] = {
            "n_subsets": int(len(subsets)),
            "subset_sizes": {k: int(len(v)) for k, v in subsets.items()},
            "conditional_mean": means,
            "pairwise": pairs,
            "n_sig_pairs": int(n_sig),
        }
    return out
