"""Ensemble decorrelation under high modality correlation.

When pairwise grad cosines are high, modalities are *not* independent voters —
they form a near-rank-1 ensemble. Skipping alone is blunt; the control plane
should **decompose roles** and **decorrelate** step sizes:

  Gram G[i,j] = cos(g_i, g_j)           # alignment geometry
  effective_rank(G) ↓  ⇒  collinear ensemble
  leader      = argmax_m α_m · (1−ρ_m)_+  # owns shared update direction
  diversifier = high residual uniqueness after removing leader direction
  redundant   = collinear with leader / low residual → damp adapt LR

LR mapping (FWD always on; only BWD step sizes reshape):
  LR_m ∝ α_m · role_gain_m
  role_gain: leader↑, diversifier mid (residual), redundant→β floor

Justification: paying proj-BWD on a redundant tower retrains the *same*
direction; diversifiers buy orthogonal residual mass — that is the only
ensemble member worth keeping when correlation is high.
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from .compute_ratio import modality_redundancy, pair_cos_matrix, unique_adapt_mass


ROLE_NAMES = ("leader", "diversifier", "redundant")


def gram_effective_rank(pair_cos: Mapping[str, float], mods: Sequence[str]) -> dict:
    """Spectral summary of the pairwise cosine Gram matrix.

    effective_rank = exp(H(λ/Σλ)) ∈ [1, |M|]  (Roy & Vetterli).
    Near 1 ⇒ one shared direction dominates ⇒ decorrelate / role-split.
    """
    mods = list(mods)
    g = pair_cos_matrix(pair_cos, mods)
    # symmetrize + clip for numeric PSD-ish
    g = 0.5 * (g + g.T)
    np.fill_diagonal(g, 1.0)
    evals = np.linalg.eigvalsh(g)
    evals = np.clip(evals, 0.0, None)
    s = float(evals.sum())
    if s <= 1e-12:
        return {
            "evals": evals.tolist(),
            "effective_rank": float(len(mods)),
            "participation_ratio": float(len(mods)),
            "top_frac": 0.0,
        }
    p = evals / s
    p = p[p > 1e-12]
    h = float(-np.sum(p * np.log(p)))
    erank = float(np.exp(h))
    # participation ratio (alternative rank proxy)
    pr = float((s * s) / max(float(np.sum(evals * evals)), 1e-12))
    top_frac = float(evals.max() / s)
    return {
        "evals": [float(x) for x in evals],
        "effective_rank": erank,
        "participation_ratio": pr,
        "top_frac": top_frac,
        "n_mods": len(mods),
    }


def residual_uniqueness(
    pair_cos: Mapping[str, float],
    mods: Sequence[str],
    leader: str,
) -> dict[str, float]:
    """How much of each modality is *not* explained by the leader direction.

    r_m = (1 − max(cos(g_m, g_leader), 0))_+
    leader itself gets r=0 (it *is* the shared direction).
    """
    mods = list(mods)
    mat = pair_cos_matrix(pair_cos, mods)
    idx = {m: i for i, m in enumerate(mods)}
    li = idx[leader]
    out = {}
    for m in mods:
        if m == leader:
            out[m] = 0.0
            continue
        c = float(mat[idx[m], li])
        out[m] = float(max(1.0 - max(c, 0.0), 0.0))
    return out


def assign_ensemble_roles(
    alpha: Mapping[str, float],
    pair_cos: Mapping[str, float],
    mods: Sequence[str],
    *,
    eta: float = 0.75,
    rho_high: float = 0.55,
    residual_div: float = 0.25,
    erank_collinear: float = 1.55,
) -> dict:
    """Decompose modalities into ensemble roles under correlation geometry.

    Trigger: mean redundancy ≥ rho_high **or** effective_rank ≤ erank_collinear.
    When *not* triggered (low correlation / full-rank ensemble), every modality
    is treated as an independent voter → role=diversifier, no damp.

    Returns roles, unique/residual masses, spectral stats, and whether
    decorrelation mode is active.
    """
    mods = list(mods)
    rho = modality_redundancy(pair_cos, mods)
    uniq = unique_adapt_mass(alpha, rho, mods, eta=eta)
    spec = gram_effective_rank(pair_cos, mods)
    mean_rho = float(np.mean([rho[m] for m in mods]))
    decorr = bool(mean_rho >= float(rho_high) or spec["effective_rank"] <= float(erank_collinear))

    # leader = highest unique mass (attribution after redundancy discount)
    leader = max(mods, key=lambda m: uniq[m])
    resid = residual_uniqueness(pair_cos, mods, leader)
    # residual *mass* for budgeting: α · residual uniqueness
    resid_mass_raw = {m: float(alpha[m]) * float(resid[m]) for m in mods}
    s = sum(resid_mass_raw.values())
    resid_mass = (
        {m: resid_mass_raw[m] / s for m in mods}
        if s > 1e-12
        else {m: 0.0 for m in mods}
    )

    roles: dict[str, str] = {}
    if not decorr:
        roles = {m: "diversifier" for m in mods}
    else:
        roles[leader] = "leader"
        for m in mods:
            if m == leader:
                continue
            # diversifier if residual uniqueness clears threshold
            if resid[m] >= float(residual_div):
                roles[m] = "diversifier"
            else:
                roles[m] = "redundant"

    return {
        "decorr_active": decorr,
        "leader": leader,
        "roles": roles,
        "redundancy": rho,
        "unique_mass": uniq,
        "residual": resid,
        "residual_mass": resid_mass,
        "mean_redundancy": mean_rho,
        "spectral": spec,
        "rho_high": float(rho_high),
        "residual_div": float(residual_div),
        "erank_collinear": float(erank_collinear),
        "eta": float(eta),
    }


def role_lr_gains(
    roles: Mapping[str, str],
    mods: Sequence[str],
    *,
    residual: Mapping[str, float] | None = None,
    leader_gain: float = 1.35,
    diversifier_gain: float = 1.00,
    redundant_gain: float = 0.25,
    residual_boost: float = 0.50,
) -> dict[str, float]:
    """Map ensemble roles → multiplicative LR gains (clipped ≥0).

    Diversifiers get an extra residual_boost · r_m so more-orthogonal members
    step harder than near-collinear ones still labeled diversifier.
    """
    mods = list(mods)
    residual = residual or {m: 0.0 for m in mods}
    out = {}
    for m in mods:
        role = roles.get(m, "diversifier")
        if role == "leader":
            out[m] = float(leader_gain)
        elif role == "redundant":
            out[m] = float(redundant_gain)
        else:
            out[m] = float(diversifier_gain + residual_boost * float(residual.get(m, 0.0)))
    return out


def soft_decorr_lr(
    alpha: Mapping[str, float],
    pair_cos: Mapping[str, float],
    mods: Sequence[str],
    *,
    beta: float = 0.10,
    eta: float = 0.75,
    rho_high: float = 0.55,
    residual_div: float = 0.25,
    leader_gain: float = 1.35,
    diversifier_gain: float = 1.00,
    redundant_gain: float = 0.25,
    residual_boost: float = 0.50,
    gain: Mapping[str, float] | None = None,
) -> dict:
    """Soft α→LR with ensemble decorrelation gains when correlation is high.

    When decorr is inactive (low ρ / high effective rank), reduces to soft LR
    (role gains ≈ 1 with mild residual boost only if roles are all diversifier
    and residual_boost is applied — here we force gain=1 when inactive).
    """
    from .lr_controller import alpha_to_lr  # local import avoids cycle

    mods = list(mods)
    decomp = assign_ensemble_roles(
        alpha,
        pair_cos,
        mods,
        eta=eta,
        rho_high=rho_high,
        residual_div=residual_div,
    )
    if decomp["decorr_active"]:
        rg = role_lr_gains(
            decomp["roles"],
            mods,
            residual=decomp["residual"],
            leader_gain=leader_gain,
            diversifier_gain=diversifier_gain,
            redundant_gain=redundant_gain,
            residual_boost=residual_boost,
        )
    else:
        rg = {m: 1.0 for m in mods}

    # fold optional external gain
    gain = gain or {m: 1.0 for m in mods}
    combined = {m: float(rg[m]) * float(gain[m]) for m in mods}
    lr = alpha_to_lr(alpha, mods, beta=beta, gain=combined)
    return {
        "lr": lr,
        "role_gain": rg,
        "decomp": decomp,
    }


def characterize_ensemble_traj(
    traj: Sequence[Mapping],
    mods: Sequence[str],
    *,
    eta: float = 0.75,
    rho_high: float = 0.55,
    residual_div: float = 0.25,
) -> dict:
    """Summarize role assignments / spectral rank over a window trajectory."""
    mods = list(mods)
    if not traj:
        return {}
    eranks, mean_rhos, decorr_flags = [], [], []
    role_counts = {r: 0 for r in ROLE_NAMES}
    leaders: dict[str, int] = {}
    for r in traj:
        alpha = r.get("alpha") or {}
        pair = r.get("pair_cos") or {}
        if not pair and "mean_pair_cos" in r:
            mpc = float(r["mean_pair_cos"])
            pair = {
                f"{a}|{b}": mpc
                for i, a in enumerate(mods)
                for b in mods[i + 1 :]
            }
        decomp = assign_ensemble_roles(
            alpha, pair, mods, eta=eta, rho_high=rho_high, residual_div=residual_div
        )
        eranks.append(decomp["spectral"]["effective_rank"])
        mean_rhos.append(decomp["mean_redundancy"])
        decorr_flags.append(1.0 if decomp["decorr_active"] else 0.0)
        leaders[decomp["leader"]] = leaders.get(decomp["leader"], 0) + 1
        for m in mods:
            role_counts[decomp["roles"][m]] = role_counts.get(decomp["roles"][m], 0) + 1
    n = max(len(traj), 1)
    return {
        "n_windows": len(traj),
        "mean_redundancy": float(np.nanmean(mean_rhos)),
        "mean_effective_rank": float(np.nanmean(eranks)),
        "frac_decorr_active": float(np.nanmean(decorr_flags)),
        "leader_counts": leaders,
        "role_counts": role_counts,
        "role_frac": {k: float(v) / (n * len(mods)) for k, v in role_counts.items()},
        "rho_high": float(rho_high),
        "residual_div": float(residual_div),
        "eta": float(eta),
    }
