"""Streaming CFPerm subset for AGOD L0 gate + L1 intensity.

Subset used here (from ``Python/src`` / package README):

* **DRPerm** — batch-level distribution-shift test via PO-risk + permute-``W``.
  This is the L0 trigger: recent control = ``W=0``, current batch = ``W=1``.
* **RRPerm** — R-risk analogue (optional; same permute-``W`` skeleton).
* **CFPerm-VIMP** (permuCATE / LOCO / GRF) — feature-level attribution after
  reject; *not* the primary stream gate (too heavy), but usable as an L1
  secondary signal (``vimp_focus``) when enabled.

OnlineRFPerm is **not** used here. Estimation and evaluation hooks live in
this module so dual/blend can sit strictly post-CFPerm.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Literal, Optional, Sequence, Tuple

import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import KFold

RiskKind = Literal["dr", "rr"]  # DRPerm PO-risk vs RRPerm R-risk
PropensityMode = Literal["known", "fit"]
# ``known`` (default for stream L0): e ≡ n1/n — batch W is by design, so
# fitting e(X) absorbs covariate-shift and collapses PO → p≈1 on adjacent packs.
# ``fit``: classical cross-fit propensity (size/power studies, feature-rich e).


@dataclass
class CFPermGateResult:
    """One stream-step CFPerm decision."""

    statistic: float
    p_value: float
    reject: bool
    alpha: float
    risk_kind: str
    n_perm: int
    # L1 ingredients
    intensity: float
    p_strength: float
    t_strength: float
    po_gap: float
    # optional diagnostics
    mu_hat: Optional[np.ndarray] = None
    e_hat: Optional[np.ndarray] = None
    pseudo_outcome: Optional[np.ndarray] = None
    tau_score: Optional[np.ndarray] = None
    perm_stats: Optional[np.ndarray] = None


def _as_1d(a: np.ndarray) -> np.ndarray:
    return np.asarray(a, float).ravel()


def _make_folds(n: int, n_splits: int, seed: int) -> List[np.ndarray]:
    n_splits = int(max(2, min(n_splits, n // 2 if n >= 4 else 2)))
    if n < 4:
        # tiny batch: single fold = all rows (in-sample nuisances)
        return [np.arange(n)]
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    return [te for _, te in kf.split(np.arange(n))]


def _fit_mu(X: np.ndarray, y: np.ndarray, seed: int) -> RandomForestRegressor:
    m = RandomForestRegressor(
        n_estimators=40,
        max_depth=6,
        min_samples_leaf=2,
        random_state=seed,
        n_jobs=1,
    )
    m.fit(X, y)
    return m


def _fit_e(X: np.ndarray, w: np.ndarray, seed: int) -> RandomForestClassifier:
    clf = RandomForestClassifier(
        n_estimators=40,
        max_depth=6,
        min_samples_leaf=2,
        random_state=seed,
        n_jobs=1,
    )
    # need both classes for predict_proba
    if len(np.unique(w)) < 2:
        # degenerate: constant propensity
        class _Const:
            def predict_proba(self, X_):
                p = float(np.mean(w))
                out = np.zeros((len(X_), 2), float)
                out[:, 0] = 1.0 - p
                out[:, 1] = p
                return out

        return _Const()  # type: ignore[return-value]
    clf.fit(X, w)
    return clf


def _e_prob(model, X: np.ndarray) -> np.ndarray:
    proba = model.predict_proba(X)
    # column for class 1
    classes = getattr(model, "classes_", np.array([0, 1]))
    if 1 in list(classes):
        j = int(np.where(classes == 1)[0][0])
        return np.asarray(proba[:, j], float)
    return np.asarray(proba[:, -1], float)


def crossfit_nuisances(
    X: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    *,
    n_splits: int = 3,
    seed: int = 0,
    clip_e: float = 0.02,
) -> Tuple[np.ndarray, np.ndarray]:
    """Cross-fit μ(x)=E[Y|X] and e(x)=P(W=1|X)."""
    X = np.asarray(X, float)
    y = _as_1d(y)
    w = _as_1d(w).astype(int)
    n = len(y)
    mu = np.zeros(n, float)
    e = np.zeros(n, float)
    folds = _make_folds(n, n_splits, seed)
    if len(folds) == 1:
        mu_m = _fit_mu(X, y, seed)
        e_m = _fit_e(X, w, seed + 1)
        mu[:] = mu_m.predict(X)
        e[:] = _e_prob(e_m, X)
    else:
        idx_all = np.arange(n)
        for k, te in enumerate(folds):
            tr = np.setdiff1d(idx_all, te)
            mu_m = _fit_mu(X[tr], y[tr], seed + k)
            e_m = _fit_e(X[tr], w[tr], seed + 100 + k)
            mu[te] = mu_m.predict(X[te])
            e[te] = _e_prob(e_m, X[te])
    e = np.clip(e, clip_e, 1.0 - clip_e)
    return mu, e


def po_dr_scores(
    y: np.ndarray,
    w: np.ndarray,
    mu: np.ndarray,
    e: np.ndarray,
) -> np.ndarray:
    """Doubly-robust style pseudo-outcome residual product: (Y−μ)(W−e)."""
    return (y - mu) * (w - e)


def po_risk_statistic(
    X: np.ndarray,
    po: np.ndarray,
    *,
    seed: int = 0,
) -> Tuple[float, np.ndarray]:
    """Observed DRPerm statistic = mean(τ̂(X)^2), τ̂ := regress PO on X."""
    tau_m = _fit_mu(X, po, seed)
    tau = tau_m.predict(X)
    return float(np.mean(tau**2)), np.asarray(tau, float)


def rr_risk_statistic(
    y: np.ndarray,
    w: np.ndarray,
    mu: np.ndarray,
    e: np.ndarray,
    X: np.ndarray,
    *,
    seed: int = 0,
    clip: float = 0.02,
) -> Tuple[float, np.ndarray]:
    """R-risk: mean((Y−μ − τ̂(X)·(W−e))^2) with τ̂ from Ỹ/(W−e)."""
    y_t = y - mu
    w_t = w - e
    den = np.where(np.abs(w_t) < clip, np.sign(w_t) * clip + (w_t == 0) * clip, w_t)
    pseudo_tau = y_t / den
    tau_m = _fit_mu(X, pseudo_tau, seed)
    tau = tau_m.predict(X)
    risk = float(np.mean((y_t - tau * w_t) ** 2))
    return risk, np.asarray(tau, float)


def _split_idx(n: int, seed: int, frac: float = 0.5) -> Tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    n_tr = max(4, int(round(n * frac)))
    n_tr = min(n_tr, n - 4) if n >= 8 else max(1, n // 2)
    return idx[:n_tr], idx[n_tr:]


def _stat_on_split(
    X: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    tr: np.ndarray,
    te: np.ndarray,
    *,
    risk: RiskKind,
    seed: int,
    clip_e: float,
    e_mode: PropensityMode = "known",
) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Fit nuisances on ``tr``, score risk on ``te`` (honest for both obs & perm)."""
    mu_m = _fit_mu(X[tr], y[tr], seed)
    mu_te = mu_m.predict(X[te])
    if e_mode == "known":
        # Design propensity: W is batch assignment, not a treatment to estimate.
        e_const = float(np.clip(np.mean(w[tr]), clip_e, 1.0 - clip_e))
        e_te = np.full(len(te), e_const, float)
    else:
        e_m = _fit_e(X[tr], w[tr], seed + 1)
        e_te = np.clip(_e_prob(e_m, X[te]), clip_e, 1.0 - clip_e)
    y_te, w_te = y[te], w[te]
    X_te = X[te]
    if risk == "dr":
        po_te = po_dr_scores(y_te, w_te, mu_te, e_te)
        # fit τ on te (small-sample); statistic = mean(τ^2)
        if len(te) >= 8:
            tau_m = _fit_mu(X_te, po_te, seed + 2)
            tau = tau_m.predict(X_te)
            stat = float(np.mean(tau**2))
        else:
            tau = po_te
            stat = float(np.mean(po_te**2))
    else:
        po_te = po_dr_scores(y_te, w_te, mu_te, e_te)
        if len(te) >= 8:
            stat, tau = rr_risk_statistic(
                y_te, w_te, mu_te, e_te, X_te, seed=seed + 2, clip=clip_e
            )
        else:
            tau = po_te
            stat = float(np.mean((y_te - mu_te) ** 2))
    return stat, mu_te, e_te, po_te, tau


def cfperm_batch_test(
    X_recent: np.ndarray,
    y_recent: np.ndarray,
    X_ood: np.ndarray,
    y_ood: np.ndarray,
    *,
    risk: RiskKind = "dr",
    n_perm: int = 39,
    n_splits: int = 3,  # kept for API compat; split protocol uses 50/50
    alpha: float = 0.05,
    seed: int = 0,
    clip_e: float = 0.02,
    e_mode: PropensityMode = "known",
    return_detail: bool = False,
) -> CFPermGateResult:
    """CFPerm L0: permute batch labels W between recent (0) and current (1).

    Uses a **shared train/test split** for observed and permutations so the
    null is honest. Default ``e_mode='known'`` (constant design propensity)
    — required for stream packs where X predicts batch and fitted e collapses PO.

    Statistic large under shift → one-sided p = P(perm ≥ obs).
    """
    del n_splits  # split protocol supersedes K-fold here
    Xr = np.asarray(X_recent, float)
    Xo = np.asarray(X_ood, float)
    yr = _as_1d(y_recent)
    yo = _as_1d(y_ood)
    X = np.vstack([Xr, Xo])
    y = np.concatenate([yr, yo])
    w = np.concatenate([np.zeros(len(yr), int), np.ones(len(yo), int)])
    tr, te = _split_idx(len(y), seed)

    stat, mu_te, e_te, po_te, tau = _stat_on_split(
        X, y, w, tr, te, risk=risk, seed=seed + 200, clip_e=clip_e, e_mode=e_mode
    )

    rng = np.random.default_rng(seed + 999)
    perm_stats = np.zeros(n_perm, float)
    for b in range(n_perm):
        wb = rng.permutation(w)
        perm_stats[b], _, _, _, _ = _stat_on_split(
            X,
            y,
            wb,
            tr,
            te,
            risk=risk,
            seed=seed + 400 + b,
            clip_e=clip_e,
            e_mode=e_mode,
        )

    # Monte-Carlo p uses (1+#)/ (B+1); requires 1/(B+1) < alpha for any reject.
    p = float((1.0 + np.sum(perm_stats >= stat)) / (n_perm + 1))
    reject = bool(p < alpha)
    if n_perm + 1 <= 1.0 / max(alpha, 1e-12):
        # Caller should bump n_perm; keep reject=False rather than lie.
        reject = False

    # L1 intensity from CFPerm signals (not RFPerm MSE-gap)
    p_strength = float(np.clip((alpha - p) / max(alpha, 1e-8), 0.0, 1.0))
    null_scale = float(np.std(perm_stats) + 1e-8)
    t_strength = float(np.tanh(max(stat - float(np.mean(perm_stats)), 0.0) / null_scale))
    # PO |score| gap on held-out recent vs ood rows
    te_w = w[te]
    abs_po = np.abs(po_te)
    m0 = float(np.mean(abs_po[te_w == 0])) + 1e-8 if np.any(te_w == 0) else 1e-8
    m1 = float(np.mean(abs_po[te_w == 1])) if np.any(te_w == 1) else m0
    po_gap = float(np.clip(m1 / m0 - 1.0, 0.0, 2.0) / 2.0)
    intensity = float(
        np.clip(0.55 * po_gap + 0.30 * p_strength + 0.15 * t_strength, 0.0, 1.0)
    )

    out = CFPermGateResult(
        statistic=float(stat),
        p_value=p,
        reject=reject,
        alpha=float(alpha),
        risk_kind=risk,
        n_perm=int(n_perm),
        intensity=intensity,
        p_strength=p_strength,
        t_strength=t_strength,
        po_gap=po_gap,
    )
    if return_detail:
        out.mu_hat = mu_te
        out.e_hat = e_te
        out.pseudo_outcome = po_te
        out.tau_score = tau
        out.perm_stats = perm_stats
    return out


def cfperm_intensity_to_temper(
    intensity: float,
    *,
    lam_max: float = 0.75,
    gate: float = 0.20,
    beijing_gate: float = 0.45,
) -> Tuple[float, bool]:
    """Map CFPerm intensity → (λ, is_beijing).

    Same dual schedule as before, but driven by CFPerm L1 not RFPerm drift.
    """
    d = float(np.clip(intensity, 0.0, 1.0))
    is_bj = d > float(beijing_gate)
    if d <= gate:
        return 0.0, is_bj
    lam = float(lam_max * (d - gate) / max(1.0 - gate, 1e-8))
    return lam, is_bj


# ---------------------------------------------------------------------------
# Evaluation helpers (gate quality, not downstream MSE)
# ---------------------------------------------------------------------------


@dataclass
class CFPermEvalSummary:
    """Size / power / agreement diagnostics for the CFPerm gate."""

    size: float  # false-reject rate under null packs
    power: float  # true-reject rate under shift packs
    mean_p_null: float
    mean_p_alt: float
    mean_stat_null: float
    mean_stat_alt: float
    n_null: int
    n_alt: int


def eval_cfperm_size_power(
    trials: Sequence[Dict[str, float]],
) -> CFPermEvalSummary:
    """Aggregate per-trial dicts with keys: label∈{null,alt}, reject, p, statistic."""
    null = [t for t in trials if t.get("label") == "null"]
    alt = [t for t in trials if t.get("label") == "alt"]

    def _mean(xs, key, default=float("nan")):
        vals = [float(x[key]) for x in xs if key in x]
        return float(np.mean(vals)) if vals else default

    return CFPermEvalSummary(
        size=_mean(null, "reject", 0.0),
        power=_mean(alt, "reject", 0.0),
        mean_p_null=_mean(null, "p"),
        mean_p_alt=_mean(alt, "p"),
        mean_stat_null=_mean(null, "statistic"),
        mean_stat_alt=_mean(alt, "statistic"),
        n_null=len(null),
        n_alt=len(alt),
    )


def synthetic_shift_trial(
    *,
    n0: int = 128,
    n1: int = 128,
    p: int = 8,
    shift: float = 0.0,
    seed: int = 0,
    risk: RiskKind = "dr",
    n_perm: int = 39,
    alpha: float = 0.05,
    e_mode: PropensityMode = "known",
) -> Dict[str, float]:
    """Gaussian DGP for gate eval.

    Null (``shift=0``): same β on both batches.
    Alt (``shift>0``): **concept drift** — batch-1 uses β·(1+shift) and an
    additive label offset. Covariate law stays the same (avoids pure
    covariate-shift cases where fitted e(X) absorbs W and PO-risk collapses).
    """
    rng = np.random.default_rng(seed)
    beta = rng.normal(size=p)
    X0 = rng.normal(size=(n0, p))
    X1 = rng.normal(size=(n1, p))
    y0 = X0 @ beta + 0.5 * rng.normal(size=n0)
    beta1 = beta * (1.0 + float(shift))
    y1 = X1 @ beta1 + 1.0 * float(shift) + 0.5 * rng.normal(size=n1)
    res = cfperm_batch_test(
        X0,
        y0,
        X1,
        y1,
        risk=risk,
        n_perm=n_perm,
        alpha=alpha,
        seed=seed,
        e_mode=e_mode,
    )
    return {
        "label": "null" if abs(shift) < 1e-12 else "alt",
        "reject": float(res.reject),
        "p": res.p_value,
        "statistic": res.statistic,
        "intensity": res.intensity,
        "shift": float(shift),
    }
