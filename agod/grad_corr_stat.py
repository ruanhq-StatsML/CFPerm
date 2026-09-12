"""Modality-gradient correlation → next-stage LR.

Two complementary views of the same Gram geometry ``R_ij = cos(g_i, g_j)``:

Ensemble (voters)
    Modalities are members of a gradient ensemble. High pairwise correlation
    collapses diversity (Krogh–Vedelsby ambiguity → 0); the fused step should
    *not* treat them as independent. Precision-weighted stacking
    ``π* ∝ R^{-1} α`` is the continuous analogue of leader / diversifier /
    redundant roles.

Statistics (GLS / matched filter)
    If votes are correlated, Gauss–Markov / SNR maximisation gives the same
    ``π* ∝ R^{-1} α``. Effective ensemble size ``N_eff = 1ᵀ R^{-1} 1`` and
    variance-stabilising scale ``η ∝ 1/√(πᵀ R π)`` set the *next* global
    step. Shrinkage + PSD projection keep R invertible when |M| is tiny.

FWD stays on; only next-window adapt LR is reshaped.
"""
from __future__ import annotations

import math
from typing import Mapping, Sequence

import numpy as np

from .compute_ratio import pair_cos_matrix


EPS = 1e-8


def correlation_from_pairs(
    pair_cos: Mapping[str, float],
    mods: Sequence[str],
) -> np.ndarray:
    """|M|×|M| cosine Gram, diag=1, entries clipped to [-1, 1]."""
    r = pair_cos_matrix(pair_cos, mods)
    r = 0.5 * (r + r.T)
    np.fill_diagonal(r, 1.0)
    return np.clip(r, -1.0, 1.0)


def psd_project(r: np.ndarray, *, eps: float = EPS) -> np.ndarray:
    """Higham-style PSD projection (clip eigenvalues)."""
    r = 0.5 * (np.asarray(r, float) + np.asarray(r, float).T)
    evals, evecs = np.linalg.eigh(r)
    evals = np.clip(evals, float(eps), None)
    out = (evecs * evals) @ evecs.T
    # restore correlation-scale diagonal
    d = np.sqrt(np.clip(np.diag(out), eps, None))
    out = out / np.outer(d, d)
    np.fill_diagonal(out, 1.0)
    return 0.5 * (out + out.T)


def _inv_psd(r: np.ndarray, *, eps: float = EPS) -> np.ndarray:
    evals, evecs = np.linalg.eigh(0.5 * (r + r.T))
    evals = np.clip(evals, float(eps), None)
    return (evecs * (1.0 / evals)) @ evecs.T


def mean_offdiag(r: np.ndarray) -> float:
    n = r.shape[0]
    if n <= 1:
        return 0.0
    mask = ~np.eye(n, dtype=bool)
    return float(np.clip(r[mask].mean(), -1.0, 1.0))


def equicorr_target(r: np.ndarray) -> np.ndarray:
    """Equicorrelation target T = (1-ρ̄)I + ρ̄ 11ᵀ (diag = 1)."""
    n = r.shape[0]
    rho = mean_offdiag(r)
    t = (1.0 - rho) * np.eye(n) + rho * np.ones((n, n))
    np.fill_diagonal(t, 1.0)
    return t


def shrink_lambda(r: np.ndarray, *, cond_cap: float = 12.0, lam_max: float = 0.45) -> float:
    """Larger shrinkage when R is ill-conditioned (tiny |M| sample Gram)."""
    evals = np.clip(np.linalg.eigvalsh(0.5 * (r + r.T)), EPS, None)
    cond = float(evals.max() / evals.min())
    u = float(np.clip(np.log(max(cond, 1.0)) / np.log(max(cond_cap, 1.01)), 0.0, 1.0))
    return float(u * lam_max)


def shrink_correlation(
    r: np.ndarray,
    *,
    lam: float | None = None,
    target: str = "equicorr",
) -> dict:
    """Ledoit–Wolf-style shrink of the cosine Gram.

    ``R̃ = (1-λ) R_psd + λ T``, T ∈ {equicorr, identity}.
    Auto-λ grows with log(cond(R)).
    """
    r_psd = psd_project(r)
    if target == "identity":
        t = np.eye(r.shape[0], dtype=float)
    else:
        t = equicorr_target(r_psd)
    lam_auto = shrink_lambda(r_psd)
    lam_used = float(lam_auto if lam is None else np.clip(lam, 0.0, 1.0))
    r_tilde = (1.0 - lam_used) * r_psd + lam_used * t
    r_tilde = psd_project(r_tilde)
    evals = np.clip(np.linalg.eigvalsh(r_tilde), EPS, None)
    return {
        "R": r_psd,
        "R_shrink": r_tilde,
        "target": t,
        "lam": lam_used,
        "lam_auto": lam_auto,
        "cond": float(evals.max() / evals.min()),
        "mean_rho": mean_offdiag(r_tilde),
        "evals": [float(x) for x in evals],
    }


def effective_ensemble_size(r: np.ndarray) -> dict:
    """Kish N_eff of the equal-weight average vs GLS N_eff of the BLUE.

    Equicorrelated case: both equal ``M / (1+(M-1)ρ)``.
    GLS is weakly larger when correlations are heterogeneous.
    """
    n = r.shape[0]
    ones = np.ones(n)
    kish = float((n * n) / max(float(ones @ r @ ones), EPS))
    prec = _inv_psd(r)
    gls = float(ones @ prec @ ones)
    # participation / spectral rank of R (Roy–Vetterli)
    evals = np.clip(np.linalg.eigvalsh(r), 0.0, None)
    s = float(evals.sum())
    p = evals / max(s, EPS)
    p = p[p > 1e-12]
    erank = float(np.exp(-np.sum(p * np.log(p)))) if len(p) else float(n)
    return {
        "n_mods": n,
        "n_eff_kish": kish,
        "n_eff_gls": gls,
        "n_eff": gls,
        "effective_rank": erank,
        "equicorr_neff": float(n / (1.0 + (n - 1) * mean_offdiag(r))),
    }


def gls_weights(
    r: np.ndarray,
    alpha: Mapping[str, float],
    mods: Sequence[str],
    *,
    clip_negative: bool = True,
) -> dict:
    """Matched-filter / GLS stacking weights ``π* ∝ R^{-1} α``.

    SNR = (πᵀα)² / (πᵀ R π) is maximised by π ∝ R^{-1} α.
    Negative weights = short-sale / opposing vote → clip then renormalise
    (no negative learning rate).
    """
    mods = list(mods)
    a = np.array([max(float(alpha[m]), 0.0) for m in mods], float)
    if a.sum() <= EPS:
        a = np.full(len(mods), 1.0 / len(mods))
    else:
        a = a / a.sum()
    prec = _inv_psd(r)
    raw = prec @ a
    unconstrained = {m: float(raw[i]) for i, m in enumerate(mods)}
    w = np.array(raw, float)
    if clip_negative:
        w = np.clip(w, 0.0, None)
    if w.sum() <= EPS:
        w = a.copy()
    pi = w / w.sum()
    var = float(pi @ r @ pi)
    snr = float((pi @ a) ** 2 / max(var, EPS))
    return {
        "pi": {m: float(pi[i]) for i, m in enumerate(mods)},
        "alpha": {m: float(a[i]) for i, m in enumerate(mods)},
        "unconstrained": unconstrained,
        "var_combo": var,
        "snr": snr,
        "clipped": bool(clip_negative and np.any(raw < -1e-9)),
    }


def partial_uniqueness(r: np.ndarray, mods: Sequence[str]) -> dict:
    """Precision-matrix uniqueness and partial correlations.

    Residual variance of modality m given the others is ``1 / P_mm``.
    Independent R=I → uniqueness=1; collinear ρ→1 → uniqueness→0.
    """
    mods = list(mods)
    n = len(mods)
    p = _inv_psd(r)
    uniq = {}
    for i, m in enumerate(mods):
        uniq[m] = float(1.0 / max(float(p[i, i]), EPS))
    partial = {}
    for i, mi in enumerate(mods):
        for j in range(i + 1, n):
            mj = mods[j]
            den = float(np.sqrt(max(p[i, i] * p[j, j], EPS)))
            partial[f"{mi}|{mj}"] = float(-p[i, j] / den)
    return {"uniqueness": uniq, "partial_corr": partial, "precision_diag": {m: float(p[i, i]) for i, m in enumerate(mods)}}


def fisher_z_test(rho: float, n_obs: float, *, alternative: str = "greater") -> dict:
    """Fisher z-interval / test for a correlation (one window or a trajectory).

    ``n_obs`` is the effective observation count: gradient-signature length
    for a single Gram, or number of windows for a pooled ρ̄.
    """
    rho = float(np.clip(rho, -0.999, 0.999))
    n = max(float(n_obs), 4.0)
    z = float(np.arctanh(rho))
    se = float(1.0 / np.sqrt(n - 3.0))
    z_stat = z / max(se, EPS)
    # one-sided Φ̄(|z|) via erfc for a normal tail (no scipy)
    def _sf(x: float) -> float:
        return 0.5 * float(math.erfc(x / math.sqrt(2.0)))

    if alternative == "greater":
        p = _sf(z_stat)
    elif alternative == "less":
        p = _sf(-z_stat)
    else:
        p = 2.0 * _sf(abs(z_stat))
    z_lo, z_hi = z - 1.96 * se, z + 1.96 * se
    return {
        "rho": rho,
        "n_obs": n,
        "z": z,
        "se": se,
        "z_stat": float(z_stat),
        "p_value": float(np.clip(p, 0.0, 1.0)),
        "ci95": [float(np.tanh(z_lo)), float(np.tanh(z_hi))],
        "significant_pos": bool(alternative == "greater" and p < 0.05 and rho > 0.0),
    }


def temporal_gain(
    phi: float,
    *,
    floor: float = 0.55,
    cap: float = 1.25,
) -> float:
    """AR(1)-style gain from ``φ = cos(g_t, g_{t-1})``.

    φ→+1 (stable direction) → expand the trust region; φ→−1 (chatter) → shrink.
    """
    phi = float(np.clip(phi, -1.0, 1.0))
    mid = 0.5 * (cap + floor)
    half = 0.5 * (cap - floor)
    return float(np.clip(mid + half * phi, floor, cap))


def variance_stabilizing_scale(
    pi: np.ndarray,
    r: np.ndarray,
    *,
    eta_min: float = 0.45,
    eta_max: float = 1.35,
) -> float:
    """Keep Var(πᵀg) ≈ Var(equal, R=I) = 1/M.

    Independent equal-weight: η=1. Fully collinear equal-weight: η=1/√M.
    """
    n = r.shape[0]
    v0 = 1.0 / max(n, 1)
    v = float(pi @ r @ pi)
    eta = float(np.sqrt(v0 / max(v, EPS)))
    return float(np.clip(eta, eta_min, eta_max))


def conflict_damp(
    mods: Sequence[str],
    *,
    cos_to_shared: Mapping[str, float] | None = None,
    pair_cos: Mapping[str, float] | None = None,
    floor: float = 0.35,
) -> dict[str, float]:
    """Multiplicative damp when a modality fights the shared head (or peers)."""
    mods = list(mods)
    gain = {m: 1.0 for m in mods}
    if cos_to_shared:
        for m in mods:
            c = float(cos_to_shared.get(m, 0.0))
            if c < 0.0:
                # linear map c∈[-1,0] → [floor, 1]
                gain[m] *= float(floor + (1.0 - floor) * (1.0 + c))
    if pair_cos:
        for m in mods:
            peers = [
                float(v)
                for k, v in pair_cos.items()
                if m in k.split("|") and k.count("|") == 1
            ]
            if peers and float(np.mean(peers)) < 0.0:
                gain[m] *= 0.85
    return gain


def characterize_corr(
    pair_cos: Mapping[str, float],
    mods: Sequence[str],
    *,
    alpha: Mapping[str, float] | None = None,
    n_obs: float = 32.0,
    lam: float | None = None,
) -> dict:
    """Full statistical readout of a modality-grad Gram (no LR yet)."""
    mods = list(mods)
    alpha = alpha or {m: 1.0 / len(mods) for m in mods}
    r_raw = correlation_from_pairs(pair_cos, mods)
    shrunk = shrink_correlation(r_raw, lam=lam)
    r = shrunk["R_shrink"]
    neff = effective_ensemble_size(r)
    gls = gls_weights(r, alpha, mods)
    part = partial_uniqueness(r, mods)
    rho = shrunk["mean_rho"]
    ztest = fisher_z_test(rho, n_obs)
    return {
        "mods": mods,
        "mean_rho": rho,
        "R_shrink": r,
        "shrink": {k: shrunk[k] for k in ("lam", "lam_auto", "cond", "mean_rho")},
        "n_eff": neff,
        "gls": gls,
        "partial": part,
        "fisher_z": ztest,
        "ambiguity": float(max(1.0 - max(rho, 0.0), 0.0)),
    }


def stat_corr_lr(
    alpha: Mapping[str, float],
    pair_cos: Mapping[str, float],
    mods: Sequence[str],
    *,
    beta: float = 0.10,
    lam: float | None = None,
    n_obs: float = 32.0,
    temporal_cos: Mapping[str, float] | None = None,
    cos_to_shared: Mapping[str, float] | None = None,
    use_temporal: bool = True,
    use_var_scale: bool = True,
    use_conflict: bool = True,
    gain: Mapping[str, float] | None = None,
    gls_only: bool = False,
) -> dict:
    """Next-stage LR from gradient correlation (ensemble GLS + statistical scale).

    ``gls_only=True`` → ``soft_gls``: precision weights only, η=1.
    Otherwise ``soft_stat``:

        π_gls = relu(R̃^{{-1}} α) / 1ᵀrelu(·)       # stacking / matched filter
        γ     = clip(1 − N_eff/|M|, 0, 1)          # trust GLS only if collinear
        π*    = (1−γ) α + γ π_gls
        η_g   = √( (1/M) / (π*ᵀ R̃ π*) ) · φ_t     # variance + AR(1)
        LR_m  = η_g · (β + (1−β) |M| π*_m) · c_m

    Reduces to soft α→LR when R̃ ≈ I (independent voters).
    """
    from .lr_controller import alpha_to_lr  # local import avoids cycle

    mods = list(mods)
    info = characterize_corr(pair_cos, mods, alpha=alpha, n_obs=n_obs, lam=lam)
    r = info["R_shrink"]
    gls = info["gls"]
    pi_gls = np.array([gls["pi"][m] for m in mods], float)
    a = np.array([gls["alpha"][m] for m in mods], float)
    n_mods = float(len(mods))
    n_eff = float(info["n_eff"]["n_eff"])
    # Empirical-Bayes blend: trust GLS concentration only as N_eff drops.
    # Independent voters (N_eff≈|M|) → π=α; collinear (N_eff≈1) → full GLS.
    gamma = float(np.clip(1.0 - n_eff / max(n_mods, 1.0), 0.0, 1.0))
    if gls_only:
        pi = pi_gls
        gamma = 1.0
    else:
        pi = (1.0 - gamma) * a + gamma * pi_gls
        s = float(pi.sum())
        pi = pi / s if s > EPS else a.copy()

    eta = 1.0
    if use_var_scale and not gls_only:
        eta *= variance_stabilizing_scale(pi, r)
    phi = float("nan")
    if use_temporal and temporal_cos and not gls_only:
        vals = [float(temporal_cos[m]) for m in mods if m in temporal_cos]
        if vals:
            phi = float(np.mean(vals))
            eta *= temporal_gain(phi)
    eta = float(eta)

    c_gain = {m: 1.0 for m in mods}
    if use_conflict and not gls_only:
        c_gain = conflict_damp(mods, cos_to_shared=cos_to_shared, pair_cos=pair_cos)
    ext = gain or {m: 1.0 for m in mods}
    combined = {m: float(c_gain[m]) * float(ext[m]) * eta for m in mods}

    # map blended stacking masses through the same β-floor actuator as soft LR
    pi_map = {m: float(pi[i]) for i, m in enumerate(mods)}
    lr = alpha_to_lr(pi_map, mods, beta=beta, gain=combined)

    # discrete ensemble roles (companion view)
    from .ensemble_decorr import assign_ensemble_roles, role_lr_gains

    roles = assign_ensemble_roles(alpha, pair_cos, mods)
    rg = (
        role_lr_gains(roles["roles"], mods, residual=roles["residual"])
        if roles["decorr_active"]
        else {m: 1.0 for m in mods}
    )

    return {
        "lr": lr,
        "pi": pi_map,
        "pi_gls": gls["pi"],
        "gamma_gls": gamma,
        "eta_global": eta,
        "phi_temporal": phi,
        "conflict_gain": c_gain,
        "stats": info,
        "roles": roles["roles"],
        "leader": roles["leader"],
        "decorr_active": roles["decorr_active"],
        "role_gain": rg,
        "gls_only": bool(gls_only),
    }


def characterize_stat_traj(
    traj: Sequence[Mapping],
    mods: Sequence[str],
    *,
    n_obs: float = 32.0,
) -> dict:
    """Pool window-level Grams into a trajectory-level statistical summary."""
    mods = list(mods)
    if not traj:
        return {}
    rhos, neffs, etas, conds, lams, ambs, gammas = [], [], [], [], [], [], []
    pi_stack = {m: [] for m in mods}
    for r in traj:
        pair = r.get("pair_cos") or {}
        if not pair:
            mpc = float(r.get("mean_pair_cos", 0.0))
            pair = {
                f"{a}|{b}": mpc
                for i, a in enumerate(mods)
                for b in mods[i + 1 :]
            }
        alpha = r.get("alpha") or {m: 1.0 / len(mods) for m in mods}
        packed = stat_corr_lr(
            alpha,
            pair,
            mods,
            n_obs=n_obs,
            temporal_cos=r.get("temporal_cos"),
            cos_to_shared=r.get("cos_to_shared"),
        )
        st = packed["stats"]
        rhos.append(st["mean_rho"])
        neffs.append(st["n_eff"]["n_eff"])
        etas.append(packed["eta_global"])
        conds.append(st["shrink"]["cond"])
        lams.append(st["shrink"]["lam"])
        ambs.append(st["ambiguity"])
        gammas.append(packed["gamma_gls"])
        for m in mods:
            pi_stack[m].append(packed["pi"][m])
    rho_bar = float(np.nanmean(rhos))
    ztest = fisher_z_test(rho_bar, max(len(traj), 4.0))
    return {
        "n_windows": len(traj),
        "mean_rho": rho_bar,
        "mean_n_eff": float(np.nanmean(neffs)),
        "mean_eta_global": float(np.nanmean(etas)),
        "mean_cond": float(np.nanmean(conds)),
        "mean_lam": float(np.nanmean(lams)),
        "mean_ambiguity": float(np.nanmean(ambs)),
        "mean_gamma": float(np.nanmean(gammas)),
        "mean_pi": {m: float(np.nanmean(v)) for m, v in pi_stack.items()},
        "fisher_z_pooled": ztest,
    }
