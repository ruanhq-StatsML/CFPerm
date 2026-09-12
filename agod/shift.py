"""Concept / covariate shift decomposition for AGOD LR control."""
from __future__ import annotations

from typing import Any, Mapping, Sequence

import numpy as np
from sklearn.ensemble import RandomForestRegressor

from .lr_controller import z_norm
from .mmd import rbf_mmd2, whiten_pair

RESID_TREES = 40
MMD_MAX_N = 96


def residual_concept(
    X0: np.ndarray,
    Y0: np.ndarray,
    X1: np.ndarray,
    Y1: np.ndarray,
    *,
    seed: int,
) -> float:
    """P(Y|X) proxy: f_ref fit on ref; error gap + residual MMD on cur."""
    if len(X0) < 12 or len(X1) < 12:
        return 0.0
    rf = RandomForestRegressor(
        n_estimators=RESID_TREES,
        max_depth=5,
        min_samples_leaf=3,
        random_state=seed,
        n_jobs=1,
    )
    rf.fit(X0, Y0)
    e0 = float(np.mean(np.abs(Y0 - rf.predict(X0))))
    e1 = float(np.mean(np.abs(Y1 - rf.predict(X1))))
    r0 = (Y0 - rf.predict(X0)).reshape(-1, 1)
    r1 = (Y1 - rf.predict(X1)).reshape(-1, 1)
    mmd_r = rbf_mmd2(r0, r1, max_n=min(MMD_MAX_N, 64), seed=seed + 3)
    return max(e1 - e0, 0.0) + 0.5 * mmd_r


def _score(
    con: Mapping[str, float],
    cov: Mapping[str, float],
    mods: Sequence[str],
    *,
    lam_concept: float = 1.0,
    lam_cov: float = 1.0,
) -> dict[str, float]:
    return {m: lam_concept * con[m] - lam_cov * cov[m] for m in mods}


def decompose_mmd(
    b0: Mapping[str, np.ndarray],
    b1: Mapping[str, np.ndarray],
    y0: np.ndarray,
    y1: np.ndarray,
    mods: Sequence[str],
    *,
    seed: int,
    lam_concept: float = 1.0,
    lam_cov: float = 1.0,
) -> dict[str, Any]:
    """Pure MMD: cov=MMD²(X); concept=joint-excess + residual gap."""
    y0 = np.asarray(y0, float)
    y1 = np.asarray(y1, float)
    cov_raw, con_raw, joint_raw, resid_raw = {}, {}, {}, {}
    y_all = np.concatenate([y0, y1])
    ys = y_all.std() + 1e-6
    y0n = ((y0 - y_all.mean()) / ys).reshape(-1, 1)
    y1n = ((y1 - y_all.mean()) / ys).reshape(-1, 1)

    for i, m in enumerate(mods):
        X0, X1 = whiten_pair(b0[m], b1[m])
        cov = rbf_mmd2(X0, X1, seed=seed + i)
        joint = rbf_mmd2(
            np.hstack([X0, y0n]), np.hstack([X1, y1n]), seed=seed + 17 + i
        )
        excess = max(joint - cov, 0.0)
        resid = residual_concept(X0, y0, X1, y1, seed=seed + 31 + i)
        cov_raw[m] = cov
        joint_raw[m] = joint
        resid_raw[m] = resid
        con_raw[m] = excess + resid

    cov, con = z_norm(cov_raw, mods), z_norm(con_raw, mods)
    return {
        "cov_raw": cov_raw,
        "con_raw": con_raw,
        "joint_raw": joint_raw,
        "resid_raw": resid_raw,
        "cov": cov,
        "con": con,
        "score": _score(con, cov, mods, lam_concept=lam_concept, lam_cov=lam_cov),
    }


def decompose_rf(
    msg: Any,
    mods: Sequence[str],
    *,
    lam_concept: float = 1.0,
    lam_cov: float = 1.0,
) -> dict[str, Any]:
    """RF Domain AUC/VIMP covariate + PO concept (legacy Amazon rule)."""
    cov_raw = {
        m: max(float(msg.auc[m]) - 0.5, 0.0) * (1.0 + float(msg.vimp[m]))
        for m in mods
    }
    con_raw = {m: max(float(msg.po[m]), 0.0) for m in mods}
    cov, con = z_norm(cov_raw, mods), z_norm(con_raw, mods)
    return {
        "cov_raw": cov_raw,
        "con_raw": con_raw,
        "cov": cov,
        "con": con,
        "score": _score(con, cov, mods, lam_concept=lam_concept, lam_cov=lam_cov),
    }


def decompose_hybrid(
    msg: Any,
    mmd_d: Mapping[str, Any],
    mods: Sequence[str],
    *,
    lam_concept: float = 1.0,
    lam_cov: float = 1.0,
) -> dict[str, Any]:
    """FSDS coupling: MMD² for P(X), PO for P(Y|X) — Acc-winning Amazon rule."""
    cov_raw = dict(mmd_d["cov_raw"])
    con_raw = {m: max(float(msg.po[m]), 0.0) for m in mods}
    cov, con = z_norm(cov_raw, mods), z_norm(con_raw, mods)
    return {
        "cov_raw": cov_raw,
        "con_raw": con_raw,
        "cov": cov,
        "con": con,
        "score": _score(con, cov, mods, lam_concept=lam_concept, lam_cov=lam_cov),
    }
