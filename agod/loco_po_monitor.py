"""LOCO PO-risk continuous-batch monitor for attributed tip timing.

After localization / FSDS delivers tip features ``J*``, stream consecutive
batches and ask:

  *When did the attributed tips' PO-risk contribution change?*

Per batch ``t`` (control = ``t-1``, probe = ``t``):

1. Tip-group LOCO
   ``L_tip = MAE(μ^{-J*}) − MAE(μ)`` on the probe (mask tips → control means).
2. Freeze burn-in center ``L̄ = mean(L_tip on burn)``.
3. Monitor score
   ``S_t = |L_tip(t) − L̄| + PO_t``.
4. **Hard timing gate** (primary):
   ``S_t > mean(S_burn) + k·std(S_burn)`` (default ``k=5``).
   First post-burn hit = tip change-point ``t*_tip``.
5. Peer: OnlineRFPerm on ``S_t`` and on ``PO_t`` (AR / FDR diagnostics).

At a tip mean/coef shift, ``L_tip`` often *collapses* on the transition
batch while ``PO_t`` spikes — ``S_t`` pins that accurate time point.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from sklearn.ensemble import RandomForestRegressor

from agod.online_rfperm import ScalarOnlineRFPerm, first_reject_index
from agod.po_vimp_weights import loco_po_vimp


def _fit_rf(X: np.ndarray, y: np.ndarray, *, seed: int, n_estimators: int) -> RandomForestRegressor:
    rf = RandomForestRegressor(
        n_estimators=n_estimators,
        max_depth=6,
        min_samples_leaf=2,
        random_state=seed,
        n_jobs=1,
    )
    rf.fit(X, y)
    return rf


def _mask_columns(
    X: np.ndarray,
    cols: Sequence[int],
    fill: np.ndarray,
) -> np.ndarray:
    Xm = np.asarray(X, float).copy()
    for j in cols:
        Xm[:, int(j)] = float(fill[int(j)])
    return Xm


def batch_po_risk_mae(
    X0: np.ndarray,
    y0: np.ndarray,
    X1: np.ndarray,
    y1: np.ndarray,
    *,
    seed: int = 0,
    n_estimators: int = 40,
) -> float:
    """Window PO-risk proxy: MAE gap of μ0(control) on probe vs control."""
    X0 = np.asarray(X0, float)
    X1 = np.asarray(X1, float)
    y0 = np.asarray(y0, float).ravel()
    y1 = np.asarray(y1, float).ravel()
    if len(X0) < 8 or len(X1) < 8:
        return 0.0
    mu = _fit_rf(X0, y0, seed=seed, n_estimators=n_estimators)
    e0 = float(np.mean(np.abs(y0 - mu.predict(X0))))
    e1 = float(np.mean(np.abs(y1 - mu.predict(X1))))
    return max(e1 - e0, 0.0)


def loco_po_risks(
    X0: np.ndarray,
    y0: np.ndarray,
    X1: np.ndarray,
    y1: np.ndarray,
    *,
    tip_idx: Sequence[int],
    seed: int = 0,
    n_estimators: int = 40,
    per_tip: bool = True,
) -> Dict[str, object]:
    """Raw LOCO PO-risk on tip set (and optional per-tip scores).

    Returns
    -------
    dict with
      base_mae : MAE of full μ0 on probe
      tip_group_loco : MAE(mask tips) − base_mae  (≥0 ⇒ tips help)
      po_full : max(base_mae − MAE_control, 0) style gap vs control fit
      per_tip_loco : dict j → LOCO_j (only if per_tip)
      tip_share : tip_group / (sum tip LOCO + eps)
    """
    tip = sorted({int(j) for j in tip_idx})
    X0 = np.asarray(X0, float)
    X1 = np.asarray(X1, float)
    y0 = np.asarray(y0, float).ravel()
    y1 = np.asarray(y1, float).ravel()
    n0, d = X0.shape
    if not tip or n0 < 8 or len(X1) < 8:
        return {
            "base_mae": 0.0,
            "tip_group_loco": 0.0,
            "po_full": 0.0,
            "per_tip_loco": {},
            "tip_share": 0.0,
        }

    fill = X0.mean(axis=0)
    mu = _fit_rf(X0, y0, seed=seed, n_estimators=n_estimators)
    base = float(np.mean(np.abs(y1 - mu.predict(X1))))
    e0 = float(np.mean(np.abs(y0 - mu.predict(X0))))
    po_full = max(base - e0, 0.0)

    X1_tip = _mask_columns(X1, tip, fill)
    X0_tip = _mask_columns(X0, tip, fill)
    mu_tip = _fit_rf(X0_tip, y0, seed=seed + 91, n_estimators=n_estimators)
    tip_mae = float(np.mean(np.abs(y1 - mu_tip.predict(X1_tip))))
    tip_group = max(tip_mae - base, 0.0)

    per: Dict[int, float] = {}
    tip_share = 0.0
    if per_tip:
        for j in tip:
            X0j = _mask_columns(X0, [j], fill)
            X1j = _mask_columns(X1, [j], fill)
            muj = _fit_rf(X0j, y0, seed=seed + 17 * (j + 1), n_estimators=n_estimators)
            per[j] = max(float(np.mean(np.abs(y1 - muj.predict(X1j)))) - base, 0.0)
        s = float(sum(per.values()))
        tip_share = float(tip_group / (s + 1e-12)) if s > 0 else (1.0 if tip_group > 0 else 0.0)
    else:
        tip_share = 1.0 if tip_group > 0 else 0.0

    return {
        "base_mae": base,
        "tip_group_loco": float(tip_group),
        "po_full": float(po_full),
        "per_tip_loco": per,
        "tip_share": float(tip_share),
    }


@dataclass
class LocoPOMonitorState:
    """Streaming tip-change score + OnlineRFPerm peers.

    Monitor scalar (frozen burn center ``L̄``):
      ``S_t = |L_tip(t) − L̄| + PO_t``

    Timing uses a burn-calibrated hard gate
      ``S_t > mean(S_burn) + k · std(S_burn)`` (default ``k=5``)
    so stably high tip-LOCO does not early-alarm; OnlineRFPerm on ``S_t``
    remains as a peer FDR stream.
    """

    tip_idx: List[int]
    tip_stream: ScalarOnlineRFPerm = field(default_factory=ScalarOnlineRFPerm)
    po_stream: ScalarOnlineRFPerm = field(default_factory=ScalarOnlineRFPerm)
    tip_loco_hist: List[float] = field(default_factory=list)
    tip_change_hist: List[float] = field(default_factory=list)
    score_hist: List[float] = field(default_factory=list)
    po_full_hist: List[float] = field(default_factory=list)
    per_tip_hist: List[Dict[int, float]] = field(default_factory=list)
    tip_share_hist: List[float] = field(default_factory=list)
    hard_reject_hist: List[int] = field(default_factory=list)
    loco_center: Optional[float] = None
    score_burn_mean: Optional[float] = None
    score_burn_std: Optional[float] = None
    score_threshold: Optional[float] = None
    hard_k: float = 5.0
    n_estimators: int = 30
    seed: int = 0


def init_loco_po_monitor(
    tip_idx: Sequence[int],
    *,
    seed: int = 0,
    n_estimators: int = 30,
    hard_k: float = 5.0,
) -> LocoPOMonitorState:
    tips = sorted({int(j) for j in tip_idx})
    if not tips:
        raise ValueError("tip_idx must be non-empty (attributed J*)")
    return LocoPOMonitorState(
        tip_idx=tips,
        tip_stream=ScalarOnlineRFPerm(),
        po_stream=ScalarOnlineRFPerm(),
        n_estimators=n_estimators,
        seed=seed,
        hard_k=float(hard_k),
    )


def update_loco_po_monitor(
    state: LocoPOMonitorState,
    X0: np.ndarray,
    y0: np.ndarray,
    X1: np.ndarray,
    y1: np.ndarray,
    *,
    burn_in: bool = False,
    alpha: float = 0.05,
    per_tip: bool = True,
    t: int = 0,
) -> dict:
    """One batch: tip-change score S=|L−L̄|+PO into hard gate + OnlineRFPerm."""
    risks = loco_po_risks(
        X0,
        y0,
        X1,
        y1,
        tip_idx=state.tip_idx,
        seed=state.seed + 1009 * int(t),
        n_estimators=state.n_estimators,
        per_tip=per_tip,
    )
    tip_v = float(risks["tip_group_loco"])  # type: ignore[arg-type]
    po_v = float(risks["po_full"])  # type: ignore[arg-type]
    state.tip_loco_hist.append(tip_v)
    state.po_full_hist.append(po_v)
    state.per_tip_hist.append(dict(risks["per_tip_loco"]))  # type: ignore[arg-type]
    state.tip_share_hist.append(float(risks["tip_share"]))  # type: ignore[arg-type]

    if burn_in:
        state.loco_center = float(np.mean(state.tip_loco_hist))
        tip_change = abs(tip_v - float(state.loco_center))
        tip_sign = 0.0
        score = tip_change + po_v
        state.tip_change_hist.append(tip_change)
        state.score_hist.append(score)
        sb = np.asarray(state.score_hist, float)
        state.score_burn_mean = float(sb.mean())
        state.score_burn_std = float(sb.std() + 1e-6)
        state.score_threshold = float(
            state.score_burn_mean + state.hard_k * state.score_burn_std
        )
        tip_out = state.tip_stream.step(score, burn_in=True, alpha=alpha)
        po_out = state.po_stream.step(po_v, burn_in=True, alpha=alpha)
        hard = False
    else:
        center = float(state.loco_center if state.loco_center is not None else 0.0)
        tip_change = abs(tip_v - center)
        tip_sign = float(np.sign(tip_v - center))
        score = tip_change + po_v
        state.tip_change_hist.append(tip_change)
        state.score_hist.append(score)
        thr = float(
            state.score_threshold
            if state.score_threshold is not None
            else (state.score_burn_mean or 0.0)
            + state.hard_k * (state.score_burn_std or 1e-3)
        )
        hard = bool(score > thr)
        tip_out = state.tip_stream.step(score, burn_in=False, alpha=alpha)
        po_out = state.po_stream.step(po_v, burn_in=False, alpha=alpha)

    state.hard_reject_hist.append(int(hard))
    return {
        "tip_group_loco": tip_v,
        "tip_change": float(state.tip_change_hist[-1]),
        "score": float(state.score_hist[-1]),
        "tip_loco_sign_vs_burn": tip_sign if not burn_in else 0.0,
        "po_full": po_v,
        "tip_share": float(risks["tip_share"]),  # type: ignore[arg-type]
        "per_tip_loco": risks["per_tip_loco"],
        "loco_center": state.loco_center,
        "score_threshold": state.score_threshold,
        "hard_reject": hard,
        "tip_rfperm": tip_out,
        "po_rfperm": po_out,
        "burn_in": burn_in,
    }


def run_loco_po_stream(
    stream: Sequence[Tuple[np.ndarray, np.ndarray]],
    tip_idx: Sequence[int],
    *,
    burn_in: int = 8,
    alpha: float = 0.05,
    seed: int = 0,
    n_estimators: int = 30,
    per_tip: bool = True,
    hard_k: float = 5.0,
    known_shift: Optional[int] = None,
) -> dict:
    """Walk a batch stream; return tip change-point vs full-PO peer.

    Primary timing = first hard-gate reject on ``S=|L−L̄|+PO``.
    OnlineRFPerm first-rejects are reported as peers.
    """
    if len(stream) < burn_in + 2:
        raise ValueError("stream too short for burn_in + one monitor batch")
    state = init_loco_po_monitor(
        tip_idx, seed=seed, n_estimators=n_estimators, hard_k=hard_k
    )
    rows: List[dict] = []
    for t in range(1, len(stream)):
        X0, y0 = stream[t - 1]
        X1, y1 = stream[t]
        burn = t <= burn_in
        out = update_loco_po_monitor(
            state,
            X0,
            y0,
            X1,
            y1,
            burn_in=burn,
            alpha=alpha,
            per_tip=per_tip,
            t=t,
        )
        rows.append(
            {
                "t": t,
                **{k: out[k] for k in out if k != "per_tip_loco"},
                "per_tip_loco": out["per_tip_loco"],
            }
        )

    hard_first = first_reject_index(state.hard_reject_hist, after=burn_in)
    tip_rf_first = first_reject_index(state.tip_stream.reject_hist, after=burn_in)
    po_first = first_reject_index(state.po_stream.reject_hist, after=burn_in)
    tip_t = None if hard_first is None else int(hard_first + 1)
    tip_rf_t = None if tip_rf_first is None else int(tip_rf_first + 1)
    po_t = None if po_first is None else int(po_first + 1)

    delay_tip = None
    delay_po = None
    if known_shift is not None:
        if tip_t is not None:
            delay_tip = int(tip_t - known_shift)
        if po_t is not None:
            delay_po = int(po_t - known_shift)

    return {
        "tip_idx": list(state.tip_idx),
        "burn_in": burn_in,
        "alpha": alpha,
        "hard_k": hard_k,
        "loco_center": state.loco_center,
        "score_threshold": state.score_threshold,
        "tip_first_reject_t": tip_t,
        "tip_rfperm_first_reject_t": tip_rf_t,
        "po_first_reject_t": po_t,
        "detection_delay_tip": delay_tip,
        "detection_delay_po": delay_po,
        "known_shift": known_shift,
        "lead_tip_vs_po": (
            None if tip_t is None or po_t is None else int(tip_t - po_t)
        ),
        "tip_loco_hist": list(state.tip_loco_hist),
        "tip_change_hist": list(state.tip_change_hist),
        "score_hist": list(state.score_hist),
        "po_full_hist": list(state.po_full_hist),
        "tip_share_hist": list(state.tip_share_hist),
        "hard_reject_hist": list(state.hard_reject_hist),
        "tip_reject_hist": list(state.tip_stream.reject_hist),
        "po_reject_hist": list(state.po_stream.reject_hist),
        "rows": rows,
    }



def make_tip_shift_stream(
    *,
    n_batches: int = 40,
    batch_size: int = 128,
    d: int = 8,
    tip_idx: Sequence[int] = (0, 1),
    shift_at: int = 20,
    seed: int = 0,
    shift_mean: float = 2.5,
    noise: float = 0.35,
) -> Tuple[List[Tuple[np.ndarray, np.ndarray]], List[int]]:
    """Synthetic stream: tip features drive y; mean-shift tips at ``shift_at``.

    Pre-shift: X ~ N(0,I), y = X_tip · 1 + noise.
    At/after ``shift_at``: tip columns get +shift_mean and coef flips sign
    on tip 0 so both distribution and PO structure move on J*.
    """
    rng = np.random.default_rng(seed)
    tips = sorted({int(j) for j in tip_idx})
    stream: List[Tuple[np.ndarray, np.ndarray]] = []
    for t in range(n_batches):
        X = rng.normal(size=(batch_size, d)).astype(float)
        if t >= shift_at:
            for j in tips:
                X[:, j] = X[:, j] + shift_mean
        coef = np.zeros(d, float)
        for j in tips:
            coef[j] = 1.0
        if t >= shift_at and tips:
            coef[tips[0]] = -1.2
        y = X @ coef + rng.normal(0.0, noise, size=batch_size)
        stream.append((X, y))
    return stream, tips


def tip_probs_from_loco_vimp(
    X0: np.ndarray,
    y0: np.ndarray,
    X1: np.ndarray,
    y1: np.ndarray,
    *,
    seed: int = 0,
) -> np.ndarray:
    """Convenience: LOCO-PO simplex (existing ``loco_po_vimp``)."""
    return loco_po_vimp(X0, y0, X1, y1, seed=seed)
