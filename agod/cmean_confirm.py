r"""Confirm + cmean joint formulation (prototype).

Locked loop (no extra modes):

1. **Attribution (offline / period):** cmean \(δ, D, r\) → Drill → \(K^\star\),
   tip set \(J^\star\), signed direction \(\mathrm{sign}(Δ\bar y)\),
   \(\mathrm{sign}(δ_j)\).
2. **Timing (online batches):** freeze \(J^\star\); each batch \(t\)
   (control \(t-1\), probe \(t\)) compute
   - tip-cmean magnitude
     \(Δ_{\mathrm{tip}}(t)=\|μ_t^{J^\star}-μ_{t-1}^{J^\star}\|_2\)
   - confirm LOCO: tip-collapse \(\max(0,\bar L-L_{\mathrm{tip}})\)
     **and** window PO both exceed burn gates.
3. **Joint fire (default):**
   \(\mathrm{Fire}_t = \mathbf{1}\{Δ_{\mathrm{tip}}\geε_Δ\}
   \cdot\mathbf{1}\{\mathrm{confirm}_t\}\).
   First post-burn fire = change-point \(t^\star\).

Magnitude thresholds are **burn-calibrated engineering gates**
(mean + \(k\cdot\mathrm{sd}\)), not Type-I \(α\) tests. See
``calibrate_cmean_confirm_thresholds`` and the LaTeX note.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from agod.loco_po_monitor import (
    loco_po_risks,
    make_tip_shift_stream,
)
from agod.online_rfperm import first_reject_index


# ---------------------------------------------------------------------------
# Tip-cmean magnitude
# ---------------------------------------------------------------------------


def tip_cmean_delta(
    X0: np.ndarray,
    X1: np.ndarray,
    tip_idx: Sequence[int],
) -> float:
    """‖μ1 − μ0‖₂ restricted to tip columns (standardized or raw OK)."""
    tips = sorted({int(j) for j in tip_idx})
    if not tips:
        return 0.0
    X0 = np.asarray(X0, float)
    X1 = np.asarray(X1, float)
    d0 = X0[:, tips].mean(axis=0)
    d1 = X1[:, tips].mean(axis=0)
    return float(np.linalg.norm(d1 - d0))


def tip_cmean_signed(
    X0: np.ndarray,
    X1: np.ndarray,
    tip_idx: Sequence[int],
) -> Dict[int, float]:
    """Per-tip mean shift (W_probe − W_control); sign = 正向/负向 on feature."""
    tips = sorted({int(j) for j in tip_idx})
    X0 = np.asarray(X0, float)
    X1 = np.asarray(X1, float)
    return {j: float(X1[:, j].mean() - X0[:, j].mean()) for j in tips}


# ---------------------------------------------------------------------------
# Threshold calibration (how magnitude cutoffs are set)
# ---------------------------------------------------------------------------


@dataclass
class CmeanConfirmThresholds:
    """Burn-calibrated magnitude gates.

    Recipe (default ``k=5``):
      ε_Δ     = mean(Δ_tip on burn) + k · sd(Δ_tip on burn)
      thr_col = mean(collapse on burn) + k · sd(collapse)
      thr_PO  = mean(PO on burn) + k · sd(PO)

    Justification: under a near-stationary burn window the consecutive-batch
    tip-cmean / LOCO / PO fluctuate at a stable scale; a tip change-point is
    an excursion beyond that empirical scale. This is an **engineering
    gate**, not a Type-I α claim (stream is updating; no stationary null).
    """

    eps_delta: float
    thr_collapse: float
    thr_po: float
    hard_k: float = 5.0
    burn_n: int = 0
    delta_burn_mean: float = 0.0
    delta_burn_std: float = 0.0
    collapse_burn_mean: float = 0.0
    collapse_burn_std: float = 0.0
    po_burn_mean: float = 0.0
    po_burn_std: float = 0.0
    loco_center: float = 0.0
    recipe: str = "burn_mean_plus_k_sd"


def calibrate_cmean_confirm_thresholds(
    stream: Sequence[Tuple[np.ndarray, np.ndarray]],
    tip_idx: Sequence[int],
    *,
    burn_in: int = 8,
    hard_k: float = 5.0,
    seed: int = 0,
    n_estimators: int = 20,
) -> CmeanConfirmThresholds:
    """Fit magnitude thresholds from burn consecutive batches only."""
    tips = sorted({int(j) for j in tip_idx})
    if burn_in < 2 or len(stream) < burn_in + 1:
        raise ValueError("need burn_in>=2 and stream longer than burn_in")

    deltas: List[float] = []
    collapses: List[float] = []
    pos: List[float] = []
    L_vals: List[float] = []

    for t in range(1, burn_in + 1):
        X0, y0 = stream[t - 1]
        X1, y1 = stream[t]
        dlt = tip_cmean_delta(X0, X1, tips)
        risks = loco_po_risks(
            X0, y0, X1, y1, tip_idx=tips, seed=seed + 17 * t, n_estimators=n_estimators, per_tip=False
        )
        L = float(risks["tip_group_loco"])  # type: ignore[arg-type]
        L_vals.append(L)
        po = float(risks["po_full"])  # type: ignore[arg-type]
        # provisional center = running mean of L
        center = float(np.mean(L_vals))
        col = max(0.0, center - L)
        deltas.append(dlt)
        collapses.append(col)
        pos.append(po)

    def _ms(xs: List[float]) -> Tuple[float, float]:
        a = np.asarray(xs, float)
        return float(a.mean()), float(a.std() + 1e-6)

    dm, ds = _ms(deltas)
    cm, cs = _ms(collapses)
    pm, ps = _ms(pos)
    k = float(hard_k)
    return CmeanConfirmThresholds(
        eps_delta=dm + k * ds,
        thr_collapse=cm + k * cs,
        thr_po=pm + k * ps,
        hard_k=k,
        burn_n=burn_in,
        delta_burn_mean=dm,
        delta_burn_std=ds,
        collapse_burn_mean=cm,
        collapse_burn_std=cs,
        po_burn_mean=pm,
        po_burn_std=ps,
        loco_center=float(np.mean(L_vals)),
        recipe="burn_mean_plus_k_sd",
    )


# ---------------------------------------------------------------------------
# Streaming state
# ---------------------------------------------------------------------------


@dataclass
class CmeanConfirmState:
    tip_idx: List[int]
    thr: CmeanConfirmThresholds
    delta_hist: List[float] = field(default_factory=list)
    collapse_hist: List[float] = field(default_factory=list)
    po_hist: List[float] = field(default_factory=list)
    L_hist: List[float] = field(default_factory=list)
    fire_hist: List[int] = field(default_factory=list)
    confirm_hist: List[int] = field(default_factory=list)
    cmean_gate_hist: List[int] = field(default_factory=list)
    signed_hist: List[Dict[int, float]] = field(default_factory=list)
    n_estimators: int = 20
    seed: int = 0


def init_cmean_confirm(
    tip_idx: Sequence[int],
    thr: CmeanConfirmThresholds,
    *,
    seed: int = 0,
    n_estimators: int = 20,
) -> CmeanConfirmState:
    tips = sorted({int(j) for j in tip_idx})
    if not tips:
        raise ValueError("tip_idx empty")
    return CmeanConfirmState(
        tip_idx=tips, thr=thr, seed=seed, n_estimators=n_estimators
    )


def update_cmean_confirm(
    state: CmeanConfirmState,
    X0: np.ndarray,
    y0: np.ndarray,
    X1: np.ndarray,
    y1: np.ndarray,
    *,
    burn_in: bool = False,
    t: int = 0,
) -> dict:
    """One batch of the joint formulation.

    During burn_in: record series only (thresholds already frozen).
    After burn: Fire = 1{Δ_tip ≥ ε_Δ} · 1{confirm}.
    """
    tips = state.tip_idx
    dlt = tip_cmean_delta(X0, X1, tips)
    signed = tip_cmean_signed(X0, X1, tips)
    risks = loco_po_risks(
        X0,
        y0,
        X1,
        y1,
        tip_idx=tips,
        seed=state.seed + 1009 * int(t),
        n_estimators=state.n_estimators,
        per_tip=False,
    )
    L = float(risks["tip_group_loco"])  # type: ignore[arg-type]
    po = float(risks["po_full"])  # type: ignore[arg-type]
    center = float(state.thr.loco_center)
    collapse = max(0.0, center - L)

    state.delta_hist.append(dlt)
    state.collapse_hist.append(collapse)
    state.po_hist.append(po)
    state.L_hist.append(L)
    state.signed_hist.append(signed)

    if burn_in:
        state.confirm_hist.append(0)
        state.cmean_gate_hist.append(0)
        state.fire_hist.append(0)
        return {
            "delta_tip": dlt,
            "collapse": collapse,
            "po": po,
            "L_tip": L,
            "signed_tips": signed,
            "cmean_gate": False,
            "confirm": False,
            "fire": False,
            "burn_in": True,
        }

    cmean_gate = bool(dlt >= state.thr.eps_delta)
    confirm = bool(collapse >= state.thr.thr_collapse and po >= state.thr.thr_po)
    fire = bool(cmean_gate and confirm)
    state.cmean_gate_hist.append(int(cmean_gate))
    state.confirm_hist.append(int(confirm))
    state.fire_hist.append(int(fire))
    return {
        "delta_tip": dlt,
        "collapse": collapse,
        "po": po,
        "L_tip": L,
        "signed_tips": signed,
        "cmean_gate": cmean_gate,
        "confirm": confirm,
        "fire": fire,
        "burn_in": False,
        "eps_delta": state.thr.eps_delta,
        "thr_collapse": state.thr.thr_collapse,
        "thr_po": state.thr.thr_po,
    }


def run_cmean_confirm_prototype(
    stream: Sequence[Tuple[np.ndarray, np.ndarray]],
    tip_idx: Sequence[int],
    *,
    burn_in: int = 8,
    hard_k: float = 5.0,
    seed: int = 0,
    n_estimators: int = 20,
    known_shift: Optional[int] = None,
) -> dict:
    """End-to-end prototype: calibrate on burn → monitor → first fire t*."""
    thr = calibrate_cmean_confirm_thresholds(
        stream,
        tip_idx,
        burn_in=burn_in,
        hard_k=hard_k,
        seed=seed,
        n_estimators=n_estimators,
    )
    state = init_cmean_confirm(tip_idx, thr, seed=seed, n_estimators=n_estimators)
    rows: List[dict] = []
    for t in range(1, len(stream)):
        out = update_cmean_confirm(
            state,
            stream[t - 1][0],
            stream[t - 1][1],
            stream[t][0],
            stream[t][1],
            burn_in=(t <= burn_in),
            t=t,
        )
        rows.append({"t": t, **out})

    fire_t = first_reject_index(state.fire_hist, after=burn_in)
    t_star = None if fire_t is None else int(fire_t + 1)
    delay = None if known_shift is None or t_star is None else int(t_star - known_shift)

    # signed direction at fire (or last row)
    signed_at = {}
    if t_star is not None:
        signed_at = state.signed_hist[t_star - 1]
    elif state.signed_hist:
        signed_at = state.signed_hist[-1]

    return {
        "tip_idx": list(state.tip_idx),
        "burn_in": burn_in,
        "hard_k": hard_k,
        "thresholds": {
            "eps_delta": thr.eps_delta,
            "thr_collapse": thr.thr_collapse,
            "thr_po": thr.thr_po,
            "recipe": thr.recipe,
            "delta_burn_mean": thr.delta_burn_mean,
            "delta_burn_std": thr.delta_burn_std,
            "collapse_burn_mean": thr.collapse_burn_mean,
            "collapse_burn_std": thr.collapse_burn_std,
            "po_burn_mean": thr.po_burn_mean,
            "po_burn_std": thr.po_burn_std,
            "loco_center": thr.loco_center,
        },
        "t_star": t_star,
        "detection_delay": delay,
        "known_shift": known_shift,
        "signed_tips_at_fire": signed_at,
        "fire_hist": list(state.fire_hist),
        "confirm_hist": list(state.confirm_hist),
        "cmean_gate_hist": list(state.cmean_gate_hist),
        "delta_hist": list(state.delta_hist),
        "rows": rows,
    }
