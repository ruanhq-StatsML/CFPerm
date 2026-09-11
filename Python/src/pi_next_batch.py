"""Modality π as next-batch group learning rates.

Not online learning. Not a stream. Not regret.

Batch A (two-sample on W) → freeze π.
Batch B (the *next* batch) → finite-step GD with group step sizes
    η_m = η_0 · M · π̃_m     (mean-preserving reallocation)

High-π blocks drifted, so their coordinates are more stale: give them a
larger step on B. Low-π blocks stay slow. The overall mean step is η_0,
so this is not “train harder”, it is *where* the next batch spends its
budget of steps.

The same simplex already enters packing, π-RF, and adaptive group-ridge.
Group LR is the GD/linear twin of those procedures (trees ignore monotone
column scale; GD does not).

Inverse allocation η_m ∝ 1/π_m is the ablation that *damps* the drifted
block — useful only if that block is treated as nuisance.
"""
from __future__ import annotations

from typing import Dict, List, Sequence

import numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from clever_covariate_gap import (
    CLIP_E,
    EPS,
    ModalitySpec,
    _as_1d,
    decompose_modality_gap,
    fit_vimp,
    relu_normalize,
    smoothed_pi,
    vimp_mass_share,
)


def block_learning_rates(
    spec: ModalitySpec,
    pi: np.ndarray,
    *,
    eta0: float = 0.4,
    floor: float = 0.05,
    mode: str = "boost",
) -> np.ndarray:
    """Mean-preserving group steps. boost: η_m = η0 M π_m; damp: ∝ (1/M − π)⁺."""
    pi = smoothed_pi(np.asarray(pi, dtype=float), floor=floor)
    m = max(len(pi), 1)
    if mode == "boost":
        eta = eta0 * m * pi
    elif mode == "damp":
        inv = relu_normalize(1.0 / np.maximum(pi, EPS))
        eta = eta0 * m * inv
    elif mode == "uniform":
        eta = np.full(m, float(eta0))
    else:
        raise ValueError(f"unknown mode {mode}")
    return eta.astype(float)


def coordinate_learning_rates(spec: ModalitySpec, eta_block: np.ndarray, n_features: int) -> np.ndarray:
    eta_block = np.asarray(eta_block, dtype=float).reshape(-1)
    out = np.zeros(int(n_features), dtype=float)
    for m, sl in enumerate(spec.slices):
        out[sl] = eta_block[m]
    return out


def _sigmoid(z: np.ndarray) -> np.ndarray:
    z = np.clip(np.asarray(z, dtype=float), -30.0, 30.0)
    return 1.0 / (1.0 + np.exp(-z))


def gd_logistic_block_lr(
    X: np.ndarray,
    y: np.ndarray,
    spec: ModalitySpec,
    eta_block: np.ndarray,
    *,
    n_steps: int = 25,
    ridge: float = 0.05,
    record: bool = False,
) -> dict:
    """Finite-step logistic GD. Path matters; do not run to convergence."""
    X = np.asarray(X, dtype=float)
    y = _as_1d(y)
    n, d = X.shape
    eta_c = coordinate_learning_rates(spec, eta_block, d)
    eta0 = float(np.mean(eta_block))
    beta = np.zeros(d)
    b0 = 0.0
    path: List[tuple] = []
    for _ in range(int(n_steps)):
        p = _sigmoid(X @ beta + b0)
        resid = (p - y) / n
        g = X.T @ resid + (ridge / n) * beta
        beta = beta - eta_c * g
        b0 = b0 - eta0 * float(resid.sum())
        if record:
            path.append((beta.copy(), float(b0)))
    return {"beta": beta, "intercept": b0, "path": path}


def gd_square_block_lr(
    X: np.ndarray,
    y: np.ndarray,
    spec: ModalitySpec,
    eta_block: np.ndarray,
    *,
    n_steps: int = 25,
    ridge: float = 0.05,
) -> dict:
    X = np.asarray(X, dtype=float)
    y = _as_1d(y)
    n, d = X.shape
    eta_c = coordinate_learning_rates(spec, eta_block, d)
    eta0 = float(np.mean(eta_block))
    beta = np.zeros(d)
    b0 = 0.0
    for _ in range(int(n_steps)):
        pred = X @ beta + b0
        resid = (pred - y) / n
        g = X.T @ resid + (ridge / n) * beta
        beta = beta - eta_c * g
        b0 = b0 - eta0 * float(resid.sum())
    return {"beta": beta, "intercept": b0}


def _predict_logit(model: dict, X: np.ndarray) -> np.ndarray:
    return _sigmoid(np.asarray(X, dtype=float) @ model["beta"] + model["intercept"])


def _auc(y, s) -> float:
    y = _as_1d(y).astype(int)
    if y.min() == y.max():
        return float("nan")
    return float(roc_auc_score(y, s))


def next_batch_lr_eval(
    X: np.ndarray,
    W: np.ndarray,
    Y: np.ndarray | None,
    spec: ModalitySpec,
    *,
    gt: str | None = None,
    seed: int = 0,
    n_steps: int = 20,
    eta0: float = 0.45,
    n_estimators: int = 40,
    light: bool = True,
) -> dict:
    """Batch A: freeze π from W. Batch B: finite-step GD. Batch C: holdout.

    Three disjoint slices — not a stream, not an online update of π.
    """
    X = np.asarray(X, dtype=float)
    W = _as_1d(W).astype(int)
    idx = np.arange(len(W))
    A, rest = train_test_split(idx, test_size=0.6, random_state=seed, stratify=W)
    B, C = train_test_split(rest, test_size=0.45, random_state=seed + 1, stratify=W[rest])
    gap = decompose_modality_gap(
        X[A], W[A], spec, seed=seed, n_splits=3, n_estimators=n_estimators, light=light,
    )
    pi = np.asarray(gap.pi_consensus, dtype=float)
    vimp = fit_vimp(X[A], W[A], seed=seed, n_estimators=n_estimators)
    pi_vimp = np.array([vimp_mass_share(vimp, spec)[n] for n in spec.names], dtype=float)
    scaler = StandardScaler()
    XB = scaler.fit_transform(X[B])
    XC = scaler.transform(X[C])
    names = list(spec.names)

    def eta_for(kind: str) -> np.ndarray:
        if kind == "uniform":
            return block_learning_rates(spec, pi, eta0=eta0, mode="uniform")
        if kind == "pi_boost":
            return block_learning_rates(spec, pi, eta0=eta0, mode="boost")
        if kind == "vimp_boost":
            return block_learning_rates(spec, pi_vimp, eta0=eta0, mode="boost")
        if kind == "pi_damp":
            return block_learning_rates(spec, pi, eta0=eta0, mode="damp")
        if kind == "oracle" and gt is not None and gt in names:
            one = np.zeros(len(names))
            one[names.index(gt)] = 1.0
            return block_learning_rates(spec, one, eta0=eta0, mode="boost", floor=0.0)
        raise KeyError(kind)

    kinds = ["uniform", "pi_boost", "vimp_boost", "pi_damp"]
    if gt is not None and gt in names:
        kinds.append("oracle")
    domain_auc: Dict[str, float] = {}
    mass_on_gt: Dict[str, float] = {}
    for kind in kinds:
        eta = eta_for(kind)
        fit = gd_logistic_block_lr(XB, W[B], spec, eta, n_steps=n_steps)
        domain_auc[kind] = round(_auc(W[C], _predict_logit(fit, XC)), 4)
        abs_b = np.abs(fit["beta"])
        share = vimp_mass_share(abs_b / (abs_b.sum() + EPS), spec)
        if gt is not None and gt in names:
            mass_on_gt[kind] = round(float(share[gt]), 4)

    y_mse: Dict[str, float] = {}
    if Y is not None:
        Y = _as_1d(Y)
        yb = (Y[B] - Y[B].mean()) / (Y[B].std() + EPS)
        yc = (Y[C] - Y[B].mean()) / (Y[B].std() + EPS)
        for kind in kinds:
            eta = eta_for(kind) * 0.12
            fit = gd_square_block_lr(XB, yb, spec, eta, n_steps=n_steps, ridge=0.2)
            pred = XC @ fit["beta"] + fit["intercept"]
            y_mse[kind] = round(float(np.mean((pred - yc) ** 2)), 4)

    return {
        "n_A": int(len(A)),
        "n_B": int(len(B)),
        "n_C": int(len(C)),
        "n_steps": int(n_steps),
        "pi": {n: round(float(v), 4) for n, v in zip(names, pi)},
        "pi_vimp": {n: round(float(v), 4) for n, v in zip(names, pi_vimp)},
        "eta_pi_boost": {
            n: round(float(v), 4)
            for n, v in zip(names, block_learning_rates(spec, pi, eta0=eta0, mode="boost"))
        },
        "domain_auc": domain_auc,
        "mass_on_gt": mass_on_gt,
        "y_mse": y_mse,
        "gt": gt,
        "note": (
            "Finite-step GD on the next batch. π is frozen from batch A. "
            "Not an online update."
        ),
    }


def gd_domain_path(
    X: np.ndarray,
    W: np.ndarray,
    spec: ModalitySpec,
    eta_block: np.ndarray,
    *,
    n_steps: int = 30,
    ridge: float = 0.05,
    test: tuple | None = None,
) -> np.ndarray:
    """Holdout AUC after each GD step (learning-curve for the step allocation)."""
    fit = gd_logistic_block_lr(
        X, W, spec, eta_block, n_steps=n_steps, ridge=ridge, record=True,
    )
    if test is None:
        return np.asarray([np.nan] * n_steps)
    Xt, Wt = test
    aucs = []
    for beta, b0 in fit["path"]:
        aucs.append(_auc(Wt, _sigmoid(Xt @ beta + b0)))
    return np.asarray(aucs, dtype=float)
