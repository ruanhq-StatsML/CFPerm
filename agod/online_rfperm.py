"""OnlineRFPerm (OnlinePermOOB) — PDF Algorithm 1.

Gate distribution-shift by predictive degradation of a fixed f_ref,
NOT by permute-W Monte Carlo.

  T_j = E[MSE_j] - E_ref
  p_k = #{i < b+k : T_{b+k} <= T_i} / (b+k)     # rank / conformal-style
  optional EWMA-weighted p
  online FDR (fixed-α or simple Alpha-Investing-style wealth) → reject

When reject → enable √PO reweighting; else keep uniform.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np
from sklearn.ensemble import RandomForestRegressor


@dataclass
class OnlineRFPermState:
    e_ref: float
    f_ref: RandomForestRegressor
    T_hist: List[float] = field(default_factory=list)
    p_hist: List[float] = field(default_factory=list)
    reject_hist: List[int] = field(default_factory=list)
    wealth: float = 1.0  # alpha-investing style wealth for online FDR
    n_burn: int = 0


def fit_online_rfperm(
    X_ref: np.ndarray,
    y_ref: np.ndarray,
    *,
    seed: int = 0,
    n_estimators: int = 40,
) -> OnlineRFPermState:
    """S1: fit f_ref once; E_ref ≈ OOB / holdout MSE on ref."""
    rf = RandomForestRegressor(
        n_estimators=n_estimators,
        max_depth=8,
        min_samples_leaf=2,
        random_state=seed,
        n_jobs=1,
        oob_score=True,
        bootstrap=True,
    )
    rf.fit(X_ref, y_ref)
    # OOB MSE when available; else in-sample fallback
    if getattr(rf, "oob_prediction_", None) is not None:
        e_ref = float(np.mean((y_ref - rf.oob_prediction_) ** 2))
    else:
        e_ref = float(np.mean((y_ref - rf.predict(X_ref)) ** 2))
    return OnlineRFPermState(e_ref=e_ref, f_ref=rf)


def batch_T(state: OnlineRFPermState, X: np.ndarray, y: np.ndarray) -> float:
    """T = E[MSE_batch] - E_ref under fixed f_ref."""
    pred = state.f_ref.predict(X)
    mse = float(np.mean((np.asarray(y, float) - pred) ** 2))
    return mse - state.e_ref


def rank_pvalue(T_cur: float, T_hist: List[float], *, ewma: bool = True, lam: float = 1.0) -> float:
    """Empirical / EWMA p-value vs historical T (PDF §2).

    Uniform rank:
      p = #{i: T_cur <= T_i} / n_hist
    EWMA (recent-focused), same orientation as rank:
      p = sum_i w_i 1{T_cur <= T_i} / sum_i w_i
      (large degradation → small p).
    """
    if not T_hist:
        return 1.0
    hist = np.asarray(T_hist, float)
    n = len(hist)
    if not ewma:
        # PDF: p_k = #{i: T_{b+k} <= T_i} / (b+k)  with denom = n_hist(+current slot)
        return float((1.0 + np.sum(T_cur <= hist)) / (n + 1))
    # EWMA weights: most recent hist gets weight 1, older decay.
    # Same orientation as rank p: large T_cur (worse) → fewer past T_i >= T_cur → small p.
    ages = np.arange(n - 1, -1, -1, dtype=float)  # 0 for newest
    w = np.exp(-lam * ages)
    ind = (T_cur <= hist).astype(float)
    return float((1e-12 + np.sum(w * ind)) / (np.sum(w) + 1e-12))


def online_fdr_step(
    state: OnlineRFPermState,
    p: float,
    *,
    alpha: float = 0.05,
    procedure: str = "alpha_investing",
) -> bool:
    """Decide reject at this step with online FDR control.

    - fixed: reject if p < alpha
    - alpha_investing: simple wealth rule (Foster & Stine style)
      spend α_t = wealth * alpha / (1+alpha); if reject, wealth += alpha
    """
    if procedure == "fixed":
        rej = bool(p < alpha)
        state.reject_hist.append(int(rej))
        return rej
    # alpha-investing
    wealth = max(state.wealth, 1e-8)
    alpha_t = wealth * alpha / (1.0 + alpha)
    rej = bool(p < alpha_t)
    if rej:
        state.wealth = wealth + alpha
    else:
        state.wealth = max(wealth - alpha_t, 1e-8)
    state.reject_hist.append(int(rej))
    return rej


def update_online_rfperm(
    state: OnlineRFPermState,
    X: np.ndarray,
    y: np.ndarray,
    *,
    burn_in: bool = False,
    alpha: float = 0.05,
    ewma: bool = True,
    fdr: str = "alpha_investing",
) -> dict:
    """Process one batch: append T; if not burn-in, form p and FDR reject."""
    T = batch_T(state, X, y)
    out = {"T": T, "p": 1.0, "reject": False, "burn_in": burn_in}
    if burn_in:
        state.T_hist.append(T)
        state.n_burn += 1
        state.p_hist.append(1.0)
        state.reject_hist.append(0)
        return out
    # monitoring: p vs history (burn-in + past monitoring), then append
    p = rank_pvalue(T, state.T_hist, ewma=ewma)
    rej = online_fdr_step(state, p, alpha=alpha, procedure=fdr)
    state.T_hist.append(T)
    state.p_hist.append(p)
    out.update({"p": p, "reject": rej})
    return out


def gate_blend(reject: bool, po_w: np.ndarray, *, soft: bool = False, p: float = 1.0, alpha: float = 0.05) -> np.ndarray:
    """Mix uniform and √PO weights by RFPerm gate."""
    ones = np.ones_like(po_w)
    if soft:
        # a → 1 when p << alpha
        a = float(1.0 / (1.0 + np.exp(8.0 * (p / max(alpha, 1e-8) - 1.0))))
        return (1.0 - a) * ones + a * po_w
    if reject:
        return po_w
    return ones
