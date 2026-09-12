"""Other stacking logics: not MoE, finite-step, entropy, MGDA contrast.

MoE (Jacobs et al. 1991) is an *input-conditional gate* π(x)=softmax(g(x)).
Nothing here is that. π is a global (or slowly time-varying) combination of
expert *votes*, trained on out-of-fold / holdout loss — Super Learner, not a
mixture-of-experts layer.

Finite-step Taylor of holdout loss along the stacked direction:

    L(θ − η Gπ) = L − η πᵀs + (η²/2) πᵀ H π + O(η³)

H ≈ GᵀG = R (Gauss–Newton). Dropping the remainder, simplex stacking is

    min_{π∈Δ}  −πᵀs + (η/2) πᵀ R π

η → 0  → vertex (linear gain / discrete OSL)
η → ∞  → min-var (πᵀ R π), i.e. Bates–Granger / MGDA-like
Training LR (~1e-3) lives in the linear regime, so a *real* SGD step is
not enough curvature to combine; you must use the quadratic meta-loss
(or direction-match) explicitly.

Other logics in this file are *also not MoE*:
    entropy     — π = softmax(s/τ); τ-path vertex→equal (no x-gate)
    mgda        — min ||Gπ||², *no holdout target* (Pareto / Sener–Koltun)
    mixloss     — Vovk aggregating / log-sum-exp of expert losses
    pseudo_bma  — π ∝ exp(−L_m); Yao et al. 2018: not stacking
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from .online_stacking import EPS, _simplex, project_simplex


def _vec(d: Mapping[str, float], experts: Sequence[str]) -> np.ndarray:
    return np.array([float(d[e]) for e in experts], float)


def _pgd_simplex(
    grad_fn,
    w0: np.ndarray,
    lip: float,
    *,
    n_iter: int = 250,
) -> np.ndarray:
    """Projected gradient on the probability simplex; step = 1/Lipschitz."""
    w = project_simplex(np.asarray(w0, float))
    step = 1.0 / max(float(lip), 1e-8)
    for _ in range(n_iter):
        w = project_simplex(w - step * np.asarray(grad_fn(w), float))
    return w


def quadratic_finite_step(
    scores: Mapping[str, float],
    r: np.ndarray,
    experts: Sequence[str],
    *,
    eta: float,
    n_iter: int = 250,
) -> dict[str, float]:
    """min_π −sᵀπ + (η/2) πᵀ R π on the simplex (2nd-order virtual step)."""
    experts = list(experts)
    s = _vec(scores, experts)
    rr = np.asarray(r, float)
    eta = float(max(eta, 0.0))
    if eta <= 1e-12:
        w = np.zeros(len(experts))
        w[int(np.argmax(s))] = 1.0
        return {e: float(w[i]) for i, e in enumerate(experts)}
    z = s - s.max()
    w0 = _simplex(np.exp(z))
    lip = eta * max(float(np.linalg.norm(rr, 2)), 1e-8)
    w = _pgd_simplex(lambda p: -s + eta * (rr @ p), w0, lip, n_iter=n_iter)
    return {e: float(w[i]) for i, e in enumerate(experts)}


def entropy_regularized(
    scores: Mapping[str, float],
    experts: Sequence[str],
    *,
    tau: float,
) -> dict[str, float]:
    """π = softmax(s/τ). τ→0 vertex, τ→∞ equal. Not an x-conditional MoE gate."""
    experts = list(experts)
    s = _vec(scores, experts)
    tau = float(max(tau, 1e-8))
    z = s / tau
    z = z - z.max()
    return {e: float(p) for e, p in zip(experts, _simplex(np.exp(z)))}


def mgda_min_norm(
    r: np.ndarray,
    experts: Sequence[str],
    *,
    n_iter: int = 250,
) -> dict[str, float]:
    """min_π ||Gπ||² = πᵀ R π on Δ. No holdout target — not stacking.

    Sener & Koltun NeurIPS 2018 (MGDA): min-norm point of conv{g_m}.
    Contrast: stacking uses g_hold / y; MGDA only uses conflict geometry.
    """
    experts = list(experts)
    rr = np.asarray(r, float)
    w0 = np.full(len(experts), 1.0 / len(experts))
    lip = max(float(np.linalg.norm(rr, 2)), 1e-8)
    w = _pgd_simplex(lambda p: rr @ p, w0, lip, n_iter=n_iter)
    return {e: float(w[i]) for i, e in enumerate(experts)}


def mixloss_weights(
    losses: Mapping[str, float],
    experts: Sequence[str],
    *,
    eta: float = 1.0,
) -> dict[str, float]:
    """Vovk aggregating / mixloss: π_m ∝ exp(−η L_m).

    Proper scoring of the mixture, not of a linear combo of votes.
    """
    experts = list(experts)
    L = _vec(losses, experts)
    z = -float(eta) * (L - L.min())
    return {e: float(p) for e, p in zip(experts, _simplex(np.exp(z)))}


def pseudo_bma_weights(
    losses: Mapping[str, float],
    experts: Sequence[str],
    *,
    n_eff: float = 1.0,
) -> dict[str, float]:
    """Pseudo-BMA: π ∝ exp(−n_eff L). Yao, Vehtari, Simpson, Gelman (BA 2018):

    BMA / pseudo-BMA put mass on posterior model probability. With
    ``n_eff`` = window size (ELPD scale) this concentrates on one
    misspecified model (M-open). Stacking targets predictive risk and
    does *not* multiply the softmax by n. η=1 on per-sample loss is mild;
    η=n is the collapse.
    """
    return mixloss_weights(losses, experts, eta=float(max(n_eff, 0.0)))


def frank_wolfe_linear(
    scores: Mapping[str, float],
    experts: Sequence[str],
) -> dict[str, float]:
    """One Frank–Wolfe step on max πᵀs is a vertex. Interior needs a quadratic."""
    experts = list(experts)
    s = _vec(scores, experts)
    w = np.zeros(len(experts))
    w[int(np.argmax(s))] = 1.0
    return {e: float(w[i]) for i, e in enumerate(experts)}


def cauchy_step(
    pi: Mapping[str, float],
    scores: Mapping[str, float],
    r: np.ndarray,
    experts: Sequence[str],
) -> float:
    """Exact line-search η along a *fixed* stacked direction (1-D quadratic).

    L(η) ≈ −η πᵀs + (η²/2) πᵀRπ  ⇒  η* = (πᵀs) / (πᵀRπ)  if πᵀs > 0.
    This is the Cauchy step. Training LR ≪ η* stays in the linear regime.
    """
    experts = list(experts)
    p = _vec(pi, experts)
    s = _vec(scores, experts)
    num = float(p @ s)
    den = float(p @ np.asarray(r, float) @ p)
    if den < 1e-12 or num <= 0.0:
        return 0.0
    return float(num / den)


def rayleigh_gls(
    scores: Mapping[str, float],
    r: np.ndarray,
    experts: Sequence[str],
    *,
    clip_negative: bool = True,
) -> dict:
    """Joint min over (η, π) of the 2nd-order Taylor = GLS / matched filter.

    Sequential: freeze tiny training η, then min_π −η πᵀs + (η²/2) πᵀRπ
    → vertex (linear gain). Joint: for fixed π, η* = (πᵀs)/(πᵀRπ); plug in

        min_π  −(πᵀs)² / (2 πᵀRπ)

    which is max Rayleigh / SNR. Unconstrained π ∝ R⁻¹s — the same
    object as ``direction_match_weights`` / static GLS. Finite step size
    is *how* you leave the vertex, not a different estimator.
    """
    from .grad_corr_stat import gls_weights

    experts = list(experts)
    s_pos = {e: max(float(scores[e]), 0.0) for e in experts}
    gls = gls_weights(r, s_pos, experts, clip_negative=clip_negative)
    pi = gls["pi"]
    eta_star = cauchy_step(pi, scores, r, experts)
    p = _vec(pi, experts)
    s = _vec(scores, experts)
    return {
        "pi": pi,
        "eta_star": eta_star,
        "snr": float((p @ s) ** 2 / max(float(p @ np.asarray(r, float) @ p), 1e-18)),
        "var_combo": gls["var_combo"],
    }


def taylor_remainder_scale(eta: float, g_norm: float = 1.0) -> dict:
    """How much the O(η³) remainder can matter vs linear / quadratic terms.

    linear ~ η, quadratic ~ η², remainder ~ η³ ||∇³L|| ||Gπ||³.
    With typical adapt LR η~3e-3 the quadratic is already ~η times smaller
    than linear, so a *real* SGD step ≈ linear gain ≈ vertex.
    """
    eta = float(eta)
    lin = abs(eta) * g_norm
    quad = 0.5 * (eta**2) * (g_norm**2)
    cub = (abs(eta) ** 3) * (g_norm**3) / 6.0
    return {
        "eta": eta,
        "linear": lin,
        "quadratic": quad,
        "cubic": cub,
        "quad_over_lin": float(quad / max(lin, 1e-18)),
        "cubic_over_lin": float(cub / max(lin, 1e-18)),
        "linear_regime": bool(quad / max(lin, 1e-18) < 0.05),
    }


def eta_path(
    scores: Mapping[str, float],
    r: np.ndarray,
    experts: Sequence[str],
    etas: Sequence[float],
) -> list[dict]:
    rows = []
    for eta in etas:
        pi = quadratic_finite_step(scores, r, experts, eta=float(eta))
        p = np.array([pi[e] for e in experts])
        h = float(-(p[p > 1e-12] * np.log(p[p > 1e-12])).sum())
        rows.append(
            {
                "eta": float(eta),
                "pi": pi,
                "n_eff": float(np.exp(h)),
                "top": max(experts, key=lambda e: pi[e]),
            }
        )
    return rows


def n_path(
    losses: Mapping[str, float],
    experts: Sequence[str],
    ns: Sequence[float],
) -> list[dict]:
    rows = []
    for n in ns:
        pi = pseudo_bma_weights(losses, experts, n_eff=float(n))
        p = np.array([pi[e] for e in experts])
        h = float(-(p[p > 1e-12] * np.log(p[p > 1e-12])).sum())
        rows.append(
            {
                "n_eff": float(n),
                "pi": pi,
                "n_eff_pi": float(np.exp(h)),
                "top": max(experts, key=lambda e: pi[e]),
            }
        )
    return rows


def tau_path(
    scores: Mapping[str, float],
    experts: Sequence[str],
    taus: Sequence[float],
) -> list[dict]:
    rows = []
    for tau in taus:
        pi = entropy_regularized(scores, experts, tau=float(tau))
        p = np.array([pi[e] for e in experts])
        h = float(-(p[p > 1e-12] * np.log(p[p > 1e-12])).sum())
        rows.append(
            {
                "tau": float(tau),
                "pi": pi,
                "n_eff": float(np.exp(h)),
                "top": max(experts, key=lambda e: pi[e]),
            }
        )
    return rows


LOGIC_NAMES = (
    "linear_vertex",
    "finite_step_quadratic",
    "joint_rayleigh",
    "entropy",
    "direction_match",
    "mgda",
    "mixloss",
    "pseudo_bma",
    "frank_wolfe",
)
