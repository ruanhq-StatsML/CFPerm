"""Computational-ratio characterization & prediction under modality correlation.

Green control target: FWD inference is always paid; only *adapt BWD* is
optional. When modality gradients are highly correlated (collinear updates),
paying a second tower's proj-BWD buys little new direction — so the skip is
justifiable and the FLOPs ratio is predictable from cosine geometry.

Accounting (matches ``flops_rel_proj``):
  C_fwd = |M|·c_pf + c_sf          # always on
  C_bwd(A) = |A|·c_pb + c_sb       # A = active adapt set
  flops_rel(A) = (C_fwd + C_bwd(A)) / (C_fwd + C_bwd(M))

Correlation → redundancy:
  ρ_m = mean_{m'≠m} max(cos(g_m, g_{m'}), 0)     # aligned mass only
  unique_m = α_m · (1 − η·ρ_m)_+                 # attribution after discount
  skip low-unique mods → predict flops_rel without running the gate
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from .adapter import flops_rel_proj


def pair_cos_matrix(
    pair_cos: Mapping[str, float],
    mods: Sequence[str],
) -> np.ndarray:
    """Build a |M|×|M| cosine matrix from ``modality_grad_cosine`` pair keys."""
    mods = list(mods)
    n = len(mods)
    mat = np.eye(n, dtype=float)
    idx = {m: i for i, m in enumerate(mods)}
    for key, val in pair_cos.items():
        # keys look like "video|audio" or "text|image"
        if "|" not in key:
            continue
        a, b = key.split("|", 1)
        if a in idx and b in idx:
            i, j = idx[a], idx[b]
            c = float(val)
            mat[i, j] = mat[j, i] = c
    return mat


def modality_redundancy(
    pair_cos: Mapping[str, float],
    mods: Sequence[str],
    *,
    only_positive: bool = True,
) -> dict[str, float]:
    """Per-modality redundancy ρ_m ∈ [0,1] from pairwise grad cosines.

    Justification: if cos(g_m, g_{m'})≈1, the two updates are nearly collinear,
    so m's adapt-BWD is largely redundant given m' (and vice versa).
    Negative cosines are *conflict*, not redundancy — they do not justify skip
    (you may still want both, or damp differently). Default: only average the
    positive part.
    """
    mods = list(mods)
    mat = pair_cos_matrix(pair_cos, mods)
    out = {}
    for i, m in enumerate(mods):
        others = []
        for j, _ in enumerate(mods):
            if i == j:
                continue
            c = float(mat[i, j])
            others.append(max(c, 0.0) if only_positive else c)
        out[m] = float(np.mean(others)) if others else 0.0
    return out


def unique_adapt_mass(
    alpha: Mapping[str, float],
    redundancy: Mapping[str, float],
    mods: Sequence[str],
    *,
    eta: float = 0.75,
) -> dict[str, float]:
    """Attribution mass after correlation discount: u_m = α_m · (1 − η·ρ_m)_+.

    η=0 → ignore correlation (pure α). η=1 → fully discount collinear mass.
    Renormalized to a simplex for budgeting / keep-set selection.
    """
    mods = list(mods)
    raw = {
        m: max(float(alpha[m]) * (1.0 - float(eta) * float(redundancy.get(m, 0.0))), 0.0)
        for m in mods
    }
    s = sum(raw.values())
    if s <= 1e-12:
        return {m: 1.0 / len(mods) for m in mods}
    return {m: float(raw[m] / s) for m in mods}


def predict_keep_set(
    unique_mass: Mapping[str, float],
    mods: Sequence[str],
    *,
    keep_frac: float = 0.67,
    ensure_one: bool = True,
    rho: Mapping[str, float] | None = None,
    rho_skip: float = 0.70,
    pair_cos: Mapping[str, float] | None = None,
) -> dict[str, bool]:
    """Predict which modalities to keep for adapt-BWD (green skip set).

    Justifiable rules:
      1) keep modalities in descending unique mass until cumulative mass ≥ keep_frac
         (or at least one)
      2) if a pair has cos ≥ rho_skip, drop the *lower*-unique member of that pair
         (collinear update ⇒ one adapt-BWD suffices)

    FWD stays on for all mods; this only predicts adapt-BWD activity.
    """
    mods = list(mods)
    u = {m: float(unique_mass[m]) for m in mods}
    order = sorted(mods, key=lambda m: u[m], reverse=True)
    keep = {m: False for m in mods}
    cum = 0.0
    for m in order:
        keep[m] = True
        cum += u[m]
        if cum >= float(keep_frac) and sum(keep.values()) >= 1:
            break
    # collinear-pair skip: drop lower-unique when highly correlated
    if pair_cos:
        for key, c in pair_cos.items():
            if "|" not in key or float(c) < float(rho_skip):
                continue
            a, b = key.split("|", 1)
            if a not in keep or b not in keep:
                continue
            if not (keep[a] and keep[b]):
                continue
            # drop the less unique one
            loser = a if u[a] < u[b] else b
            # never drop the global top-unique
            if loser != order[0]:
                keep[loser] = False
    elif rho is not None and order:
        top = order[0]
        for m in mods:
            if m == top:
                continue
            if float(rho.get(m, 0.0)) >= float(rho_skip) and u[m] < u[top]:
                keep[m] = False
    if ensure_one and not any(keep.values()):
        keep[order[0]] = True
    return keep


def predict_flops_rel(
    active: Mapping[str, bool],
    mods: Sequence[str],
    *,
    c_pf: float = 1.0,
    c_pb: float = 2.0,
    c_sf: float = 0.5,
    c_sb: float = 1.0,
) -> float:
    """Predicted relative adapt FLOPs given a keep/active set (FWD always on)."""
    return flops_rel_proj(
        active, mods, c_pf=c_pf, c_pb=c_pb, c_sf=c_sf, c_sb=c_sb
    )


def predict_flops_from_correlation(
    alpha: Mapping[str, float],
    pair_cos: Mapping[str, float],
    mods: Sequence[str],
    *,
    eta: float = 0.75,
    keep_frac: float = 0.67,
    rho_skip: float = 0.70,
    c_pf: float = 1.0,
    c_pb: float = 2.0,
    c_sf: float = 0.5,
    c_sb: float = 1.0,
) -> dict:
    """End-to-end green predictor: correlation geometry → flops_rel.

    Returns redundancy, unique mass, predicted keep set, flops_rel, and a
    closed-form *savings proxy* that does not need the keep set:

      savings_proxy = (c_pb / C_full) · Σ_m ρ_m / |M|

    Interpretation: if every modality were on average ρ-fraction redundant,
    expected skippable proj-BWD mass is ρ̄·|M|·c_pb, hence relative savings
    ≈ ρ̄ · c_pb·|M| / C_full. Keep-set prediction refines this into a concrete
    flops_rel.
    """
    mods = list(mods)
    rho = modality_redundancy(pair_cos, mods)
    uniq = unique_adapt_mass(alpha, rho, mods, eta=eta)
    keep = predict_keep_set(
        uniq,
        mods,
        keep_frac=keep_frac,
        rho=rho,
        rho_skip=rho_skip,
        pair_cos=pair_cos,
    )
    flops = predict_flops_rel(
        keep, mods, c_pf=c_pf, c_pb=c_pb, c_sf=c_sf, c_sb=c_sb
    )
    c_fwd = len(mods) * c_pf + c_sf
    c_full = c_fwd + len(mods) * c_pb + c_sb
    mean_rho = float(np.mean([rho[m] for m in mods]))
    savings_proxy = float((c_pb * len(mods) * mean_rho) / max(c_full, 1e-12))
    flops_proxy = float(max(1.0 - savings_proxy, (c_fwd + c_sb) / c_full))
    n_keep = int(sum(1 for m in mods if keep[m]))
    return {
        "redundancy": rho,
        "unique_mass": uniq,
        "keep": keep,
        "n_keep": n_keep,
        "flops_rel_pred": flops,
        "mean_redundancy": mean_rho,
        "savings_proxy": savings_proxy,
        "flops_rel_proxy": flops_proxy,
        "flops_rel_full": 1.0,
        "eta": float(eta),
        "keep_frac": float(keep_frac),
        "rho_skip": float(rho_skip),
    }


def characterize_compute_ratio(
    traj: Sequence[Mapping],
    mods: Sequence[str],
    *,
    eta: float = 0.75,
    keep_frac: float = 0.67,
    rho_skip: float = 0.70,
) -> dict:
    """Summarize predicted vs realized compute ratios over a window trajectory.

    Each traj row should carry:
      - alpha, pair_cos (or mean_pair_cos + pair_cos), optional flops_rel / active
    """
    mods = list(mods)
    if not traj:
        return {}
    preds, reals, proxies, mean_rhos = [], [], [], []
    keep_fracs = []
    for r in traj:
        alpha = r.get("alpha") or {}
        pair = r.get("pair_cos") or r.get("pair_cos") or {}
        # tolerate soft_gradcos traj that only logged mean_pair_cos
        if not pair and "mean_pair_cos" in r:
            # synthesize uniform pair matrix from mean (conservative)
            mpc = float(r["mean_pair_cos"])
            pair = {
                f"{a}|{b}": mpc
                for i, a in enumerate(mods)
                for b in mods[i + 1 :]
            }
        pred = predict_flops_from_correlation(
            alpha, pair, mods, eta=eta, keep_frac=keep_frac, rho_skip=rho_skip
        )
        preds.append(pred["flops_rel_pred"])
        proxies.append(pred["flops_rel_proxy"])
        mean_rhos.append(pred["mean_redundancy"])
        keep_fracs.append(pred["n_keep"] / max(len(mods), 1))
        if "flops_rel" in r:
            reals.append(float(r["flops_rel"]))
        elif "active" in r:
            reals.append(float(flops_rel_proj(r["active"], mods)))
    out = {
        "n_windows": len(traj),
        "mean_redundancy": float(np.nanmean(mean_rhos)),
        "mean_flops_rel_pred": float(np.nanmean(preds)),
        "mean_flops_rel_proxy": float(np.nanmean(proxies)),
        "mean_keep_frac": float(np.nanmean(keep_fracs)),
        "eta": float(eta),
        "keep_frac": float(keep_frac),
        "rho_skip": float(rho_skip),
    }
    if reals:
        reals_a = np.asarray(reals, float)
        preds_a = np.asarray(preds[: len(reals)], float)
        out["mean_flops_rel_real"] = float(np.nanmean(reals_a))
        out["mae_pred_vs_real"] = float(np.nanmean(np.abs(preds_a - reals_a)))
        out["mean_signed_err"] = float(np.nanmean(preds_a - reals_a))
    return out
