"""LOCO PO-risk continuous-batch monitor for attributed tip timing.

After localization / FSDS delivers tip features ``J*``, stream consecutive
batches and ask:

  *When did the attributed tips' PO-risk contribution change?*

Per batch ``t`` (control = ``t-1``, probe = ``t``):

1. Tip-group LOCO
   ``L_tip = MAE(μ^{-J*}) − MAE(μ)`` on the probe (mask tips → control means).
2. Freeze burn-in center ``L̄ = mean(L_tip on burn)``.
3. Score modes (v2):
   ``sum`` = |L−L̄|+PO | ``abs_dev`` | ``po`` |
   ``collapse`` = max(0,L̄−L)+PO (default) | ``cusum`` on sum base.
4. **Hard / CUSUM timing gate** → tip change-point ``t*_tip``.
5. Peer: OnlineRFPerm on score / PO (AR / FDR diagnostics).

At a tip mean/coef shift, ``L_tip`` often *collapses* on the transition
batch while ``PO_t`` spikes — ``collapse`` / ``sum`` pin that time.
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
    """Streaming tip-change score + OnlineRFPerm peers (v2)."""

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
    tip_change_threshold: Optional[float] = None
    po_threshold: Optional[float] = None
    tip_change_burn: List[float] = field(default_factory=list)
    po_burn: List[float] = field(default_factory=list)
    hard_k: float = 5.0
    score_mode: str = "confirm"
    cusum: float = 0.0
    cusum_delta: float = 0.05
    cusum_h: float = 1.0
    n_estimators: int = 30
    seed: int = 0


def init_loco_po_monitor(
    tip_idx: Sequence[int],
    *,
    seed: int = 0,
    n_estimators: int = 30,
    hard_k: float = 5.0,
    score_mode: str = "confirm",
    cusum_delta: float = 0.05,
    cusum_h: float = 1.0,
) -> LocoPOMonitorState:
    tips = sorted({int(j) for j in tip_idx})
    if not tips:
        raise ValueError("tip_idx must be non-empty (attributed J*)")
    mode = str(score_mode).lower()
    if mode not in {"sum", "abs_dev", "po", "collapse", "cusum", "confirm"}:
        raise ValueError(f"unknown score_mode={score_mode!r}")
    return LocoPOMonitorState(
        tip_idx=tips,
        tip_stream=ScalarOnlineRFPerm(),
        po_stream=ScalarOnlineRFPerm(),
        n_estimators=n_estimators,
        seed=seed,
        hard_k=float(hard_k),
        score_mode=mode,
        cusum_delta=float(cusum_delta),
        cusum_h=float(cusum_h),
    )


def _compose_score(
    *,
    mode: str,
    tip_v: float,
    center: float,
    po_v: float,
    tip_change: float,
) -> float:
    if mode == "abs_dev":
        return float(tip_change)
    if mode == "po":
        return float(po_v)
    if mode == "collapse":
        return float(max(0.0, center - tip_v) + po_v)
    return float(tip_change + po_v)  # sum / cusum / confirm base


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
    """One batch: tip-change score into hard/CUSUM gate + OnlineRFPerm."""
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
        center = float(state.loco_center)
        tip_change = abs(tip_v - center)
        tip_sign = 0.0
        score = _compose_score(
            mode=state.score_mode, tip_v=tip_v, center=center, po_v=po_v, tip_change=tip_change
        )
        state.tip_change_hist.append(tip_change)
        state.score_hist.append(score)
        tip_local_burn = float(max(0.0, center - tip_v)) if state.score_mode in {"confirm", "collapse"} else tip_change
        state.tip_change_burn.append(tip_local_burn)
        state.po_burn.append(po_v)
        sb = np.asarray(state.score_hist, float)
        state.score_burn_mean = float(sb.mean())
        state.score_burn_std = float(sb.std() + 1e-6)
        state.score_threshold = float(state.score_burn_mean + state.hard_k * state.score_burn_std)
        tc = np.asarray(state.tip_change_burn, float)
        pb = np.asarray(state.po_burn, float)
        state.tip_change_threshold = float(tc.mean() + state.hard_k * (tc.std() + 1e-6))
        state.po_threshold = float(pb.mean() + state.hard_k * (pb.std() + 1e-6))
        state.cusum = 0.0
        tip_out = state.tip_stream.step(score, burn_in=True, alpha=alpha)
        po_out = state.po_stream.step(po_v, burn_in=True, alpha=alpha)
        hard = False
    else:
        center = float(state.loco_center if state.loco_center is not None else 0.0)
        tip_change = abs(tip_v - center)
        tip_sign = float(np.sign(tip_v - center))
        score = _compose_score(
            mode=state.score_mode, tip_v=tip_v, center=center, po_v=po_v, tip_change=tip_change
        )
        state.tip_change_hist.append(tip_change)
        state.score_hist.append(score)
        if state.score_mode == "cusum":
            base = float(state.score_burn_mean or 0.0)
            state.cusum = max(0.0, state.cusum + (score - base - state.cusum_delta))
            hard = bool(state.cusum > state.cusum_h)
        elif state.score_mode == "confirm":
            # tip-local collapse AND global PO both exceed burn gates
            tip_local = float(max(0.0, center - tip_v))
            hard = bool(
                tip_local > float(state.tip_change_threshold or 0.0)
                and po_v > float(state.po_threshold or 0.0)
            )
        else:
            thr = float(
                state.score_threshold
                if state.score_threshold is not None
                else (state.score_burn_mean or 0.0) + state.hard_k * (state.score_burn_std or 1e-3)
            )
            hard = bool(score > thr)
        tip_out = state.tip_stream.step(score, burn_in=False, alpha=alpha)
        po_out = state.po_stream.step(po_v, burn_in=False, alpha=alpha)

    state.hard_reject_hist.append(int(hard))
    return {
        "tip_group_loco": tip_v,
        "tip_change": float(state.tip_change_hist[-1]),
        "score": float(state.score_hist[-1]),
        "cusum": float(state.cusum),
        "score_mode": state.score_mode,
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
    score_mode: str = "confirm",
    cusum_delta: float = 0.05,
    cusum_h: float = 1.0,
    known_shift: Optional[int] = None,
) -> dict:
    """Walk a batch stream; primary timing = hard/CUSUM gate."""
    if len(stream) < burn_in + 2:
        raise ValueError("stream too short for burn_in + one monitor batch")
    state = init_loco_po_monitor(
        tip_idx,
        seed=seed,
        n_estimators=n_estimators,
        hard_k=hard_k,
        score_mode=score_mode,
        cusum_delta=cusum_delta,
        cusum_h=cusum_h,
    )
    rows: List[dict] = []
    for t in range(1, len(stream)):
        X0, y0 = stream[t - 1]
        X1, y1 = stream[t]
        burn = t <= burn_in
        out = update_loco_po_monitor(
            state, X0, y0, X1, y1, burn_in=burn, alpha=alpha, per_tip=per_tip, t=t
        )
        rows.append({"t": t, **{k: out[k] for k in out if k != "per_tip_loco"}, "per_tip_loco": out["per_tip_loco"]})

    hard_first = first_reject_index(state.hard_reject_hist, after=burn_in)
    tip_rf_first = first_reject_index(state.tip_stream.reject_hist, after=burn_in)
    po_first = first_reject_index(state.po_stream.reject_hist, after=burn_in)
    tip_t = None if hard_first is None else int(hard_first + 1)
    tip_rf_t = None if tip_rf_first is None else int(tip_rf_first + 1)
    po_t = None if po_first is None else int(po_first + 1)

    delay_tip = delay_po = None
    early_false = 0
    if known_shift is not None:
        if tip_t is not None:
            delay_tip = int(tip_t - known_shift)
        if po_t is not None:
            delay_po = int(po_t - known_shift)
        for i, r in enumerate(state.hard_reject_hist):
            t_i = i + 1
            if t_i >= burn_in + 1 and t_i < known_shift and int(r) == 1:
                early_false += 1

    return {
        "tip_idx": list(state.tip_idx),
        "burn_in": burn_in,
        "alpha": alpha,
        "hard_k": hard_k,
        "score_mode": state.score_mode,
        "loco_center": state.loco_center,
        "score_threshold": state.score_threshold,
        "tip_first_reject_t": tip_t,
        "tip_rfperm_first_reject_t": tip_rf_t,
        "po_first_reject_t": po_t,
        "detection_delay_tip": delay_tip,
        "detection_delay_po": delay_po,
        "early_false_alarms": early_false,
        "known_shift": known_shift,
        "lead_tip_vs_po": (None if tip_t is None or po_t is None else int(tip_t - po_t)),
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


def evaluate_score_modes(
    *,
    seeds: Sequence[int] = (0, 1, 2, 3, 4),
    modes: Sequence[str] = ("sum", "abs_dev", "po", "collapse", "cusum", "confirm"),
    shift_at: int = 20,
    burn_in: int = 8,
    true_tips: Sequence[int] = (0, 1),
    wrong_tips: Sequence[int] = (6, 7),
    n_batches: int = 40,
    batch_size: int = 128,
    d: int = 8,
    n_estimators: int = 20,
) -> dict:
    """Multi-seed sweep: true tips vs wrong tips × score modes."""
    rows = []
    for mode in modes:
        delays_true, delays_wrong = [], []
        hit_true = hit_wrong = early_true = early_wrong = 0
        miss_true = miss_wrong = 0
        for s in seeds:
            stream, _ = make_tip_shift_stream(
                n_batches=n_batches,
                batch_size=batch_size,
                d=d,
                tip_idx=true_tips,
                shift_at=shift_at,
                seed=int(s),
            )
            for label, tips in (("true", true_tips), ("wrong", wrong_tips)):
                out = run_loco_po_stream(
                    stream,
                    tips,
                    burn_in=burn_in,
                    seed=int(s),
                    n_estimators=n_estimators,
                    per_tip=False,
                    score_mode=mode,
                    known_shift=shift_at,
                )
                delay = out["detection_delay_tip"]
                early = int(out["early_false_alarms"])
                if label == "true":
                    early_true += early
                    if delay is None:
                        miss_true += 1
                    else:
                        delays_true.append(delay)
                        if abs(delay) <= 1:
                            hit_true += 1
                else:
                    early_wrong += early
                    if delay is None:
                        miss_wrong += 1
                    else:
                        delays_wrong.append(delay)
                        if abs(delay) <= 1:
                            hit_wrong += 1
        n = float(len(seeds))
        rows.append(
            {
                "mode": mode,
                "n_seeds": len(seeds),
                "true_delay_mean": float(np.mean(delays_true)) if delays_true else None,
                "true_delay_std": float(np.std(delays_true)) if delays_true else None,
                "true_hit_pm1": hit_true / n,
                "true_miss": miss_true / n,
                "true_early_alarms_mean": early_true / n,
                "wrong_delay_mean": float(np.mean(delays_wrong)) if delays_wrong else None,
                "wrong_hit_pm1": hit_wrong / n,
                "wrong_miss": miss_wrong / n,
                "wrong_early_alarms_mean": early_wrong / n,
            }
        )
    return {"shift_at": shift_at, "burn_in": burn_in, "rows": rows}



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
