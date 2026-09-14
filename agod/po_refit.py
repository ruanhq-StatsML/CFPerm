"""Streaming PO-learner re-fit after OnlineRFPerm reject.

Protocol (preferred)
--------------------
At stream index ``t``, OnlineRFPerm marks the **current** batch as OOD
(significant). Roles are intentionally separate:

  recent_control  R = batch_{t-k}, …, batch_{t-1}
      → re-fit μ0 (green outcome under the *recent* regime)

  ood_batch       O = batch_t
      → score PO_i = |Y_i − μ0(X_i)| (optional μ-gap blend)
      → IPTW weights w ∝ √PO / PO^{1/3} **only on O**
      → fit the downstream learner on O with those weights

Non-reject steps stay uniform (w=1) and are excluded from sig-only MSE.

Legacy pair mode
----------------
``window_mode="prev_cur"`` keeps the older cut:
  T=0 = batches immediately *before* the (prev, cur) pair
  T=1 = prev ∪ cur
  weights still taken only from the current slice of T=1.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Literal, Sequence, Tuple

import numpy as np
from sklearn.ensemble import RandomForestRegressor

from agod.po_iptw import po_iptw_weights

WindowMode = Literal["recent_ood", "prev_cur"]


@dataclass(frozen=True)
class RefitWindows:
    """Named recent-control vs OOD cut at time t."""

    X_recent: np.ndarray
    y_recent: np.ndarray
    X_ood: np.ndarray
    y_ood: np.ndarray
    recent_idx: Tuple[int, ...]
    ood_idx: int
    window_mode: WindowMode

    # Optional larger pool used only to fit μ1 / score (legacy prev∪cur).
    X_score: np.ndarray | None = None
    y_score: np.ndarray | None = None

    @property
    def X0(self) -> np.ndarray:
        return self.X_recent

    @property
    def y0(self) -> np.ndarray:
        return self.y_recent

    @property
    def X1(self) -> np.ndarray:
        return self.X_ood if self.X_score is None else self.X_score

    @property
    def y1(self) -> np.ndarray:
        return self.y_ood if self.y_score is None else self.y_score


def _stack_batches(stream: Sequence, indices: Sequence[int]) -> Tuple[np.ndarray, np.ndarray]:
    Xs, ys = [], []
    for i in indices:
        Xi, yi = stream[i]
        Xs.append(np.asarray(Xi, float))
        ys.append(np.asarray(yi, float).ravel())
    return np.vstack(Xs), np.concatenate(ys)


def _fit_mu(X: np.ndarray, y: np.ndarray, seed: int) -> RandomForestRegressor:
    rf = RandomForestRegressor(
        n_estimators=30,
        max_depth=6,
        min_samples_leaf=2,
        random_state=seed,
        n_jobs=1,
    )
    rf.fit(X, y)
    return rf


def build_recent_ood_windows(
    stream: list,
    t: int,
    *,
    n_recent: int = 1,
    window_mode: WindowMode = "recent_ood",
) -> RefitWindows:
    """Build recent-control + OOD windows at current index ``t``.

    Parameters
    ----------
    n_recent :
        How many *immediately preceding* batches form the control for μ0.
        ``recent_ood``: indices ``[t-n_recent, …, t-1]`` (includes prev).
        ``prev_cur``:   indices ``[t-1-n_recent, …, t-2]`` (before the pair).
    window_mode :
        ``recent_ood`` — OOD = current only; recent = last n_recent batches.
        ``prev_cur``   — score pool = prev∪current; control ends at t-2.
    """
    if t < 1:
        raise ValueError("need t>=1 so a previous batch exists")
    if n_recent < 1:
        raise ValueError("n_recent must be >= 1")

    X_ood, y_ood = stream[t]
    X_ood = np.asarray(X_ood, float)
    y_ood = np.asarray(y_ood, float).ravel()

    if window_mode == "recent_ood":
        recent_idx = tuple(range(max(0, t - n_recent), t))
        if not recent_idx:
            recent_idx = (0,)
        X_recent, y_recent = _stack_batches(stream, recent_idx)
        return RefitWindows(
            X_recent=X_recent,
            y_recent=y_recent,
            X_ood=X_ood,
            y_ood=y_ood,
            recent_idx=recent_idx,
            ood_idx=t,
            window_mode=window_mode,
        )

    if window_mode == "prev_cur":
        # Legacy: control before the (prev, cur) pair; score on prev∪cur.
        recent_idx = tuple(range(max(0, t - 1 - n_recent), t - 1))
        if not recent_idx:
            recent_idx = (0,)
        X_recent, y_recent = _stack_batches(stream, recent_idx)
        X_prev, y_prev = stream[t - 1]
        X_score = np.vstack([np.asarray(X_prev, float), X_ood])
        y_score = np.concatenate([np.asarray(y_prev, float).ravel(), y_ood])
        return RefitWindows(
            X_recent=X_recent,
            y_recent=y_recent,
            X_ood=X_ood,
            y_ood=y_ood,
            recent_idx=recent_idx,
            ood_idx=t,
            window_mode=window_mode,
            X_score=X_score,
            y_score=y_score,
        )

    raise ValueError(f"unknown window_mode={window_mode!r}")


def build_t01_windows(
    stream: list,
    t: int,
    *,
    n_control: int = 1,
    window_mode: WindowMode = "recent_ood",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Backward-compatible (X0,y0,X1,y1) view of :func:`build_recent_ood_windows`.

    For ``recent_ood``, ``(X1,y1)`` is the OOD batch only.
    For ``prev_cur``, ``(X1,y1)`` is prev∪cur (legacy).
    """
    w = build_recent_ood_windows(
        stream, t, n_recent=n_control, window_mode=window_mode
    )
    return w.X0, w.y0, w.X1, w.y1


def refit_po_risk_t01(
    X0: np.ndarray,
    y0: np.ndarray,
    X1: np.ndarray,
    y1: np.ndarray,
    *,
    seed: int = 0,
    blend_mu_gap: float = 0.25,
) -> np.ndarray:
    """Re-fit PO proxy: μ0 on recent control; score on the T=1 / OOD pool.

    Returns per-row PO-risk on ``X1`` (same length as ``y1``).
    """
    mu0 = _fit_mu(X0, y0, seed)
    resid = np.abs(np.asarray(y1, float) - mu0.predict(X1))
    if blend_mu_gap > 0 and len(X1) >= 16:
        mu1 = _fit_mu(X1, y1, seed + 1)
        gap = np.abs(mu1.predict(X1) - mu0.predict(X1))
        po = (1.0 - blend_mu_gap) * resid + blend_mu_gap * gap
    else:
        po = resid
    return np.asarray(po, float)


def refit_po_on_windows(
    windows: RefitWindows,
    *,
    seed: int = 0,
    blend_mu_gap: float = 0.25,
) -> np.ndarray:
    """PO-risk on the **OOD batch only** (length = len(y_ood))."""
    po_score = refit_po_risk_t01(
        windows.X0,
        windows.y0,
        windows.X1,
        windows.y1,
        seed=seed,
        blend_mu_gap=blend_mu_gap,
    )
    if windows.window_mode == "prev_cur":
        n_ood = len(windows.y_ood)
        return po_score[-n_ood:]
    return po_score


def current_batch_po_weights(
    stream: list,
    t: int,
    *,
    seed: int = 0,
    n_control: int = 1,
    n_recent: int | None = None,
    blend_mu_gap: float = 0.25,
    mode: str = "sqrt",
    power: float | None = None,
    window_mode: WindowMode = "recent_ood",
) -> np.ndarray:
    """Post-hoc PO weights for the OOD (current) batch after recent re-fit.

    ``n_recent`` aliases ``n_control`` (how many preceding batches re-fit μ0).
    """
    k = n_control if n_recent is None else int(n_recent)
    windows = build_recent_ood_windows(
        stream, t, n_recent=k, window_mode=window_mode
    )
    po_ood = refit_po_on_windows(windows, seed=seed, blend_mu_gap=blend_mu_gap)
    return po_iptw_weights(po_ood, mode=mode, power=power)  # type: ignore[arg-type]


def current_batch_sqrt_po_weights(
    stream: list,
    t: int,
    *,
    seed: int = 0,
    n_control: int = 1,
    n_recent: int | None = None,
    blend_mu_gap: float = 0.25,
    window_mode: WindowMode = "recent_ood",
) -> np.ndarray:
    """√PO weights on the OOD batch."""
    return current_batch_po_weights(
        stream,
        t,
        seed=seed,
        n_control=n_control,
        n_recent=n_recent,
        blend_mu_gap=blend_mu_gap,
        mode="sqrt",
        window_mode=window_mode,
    )


def current_batch_cbrt_po_weights(
    stream: list,
    t: int,
    *,
    seed: int = 0,
    n_control: int = 1,
    n_recent: int | None = None,
    blend_mu_gap: float = 0.25,
    window_mode: WindowMode = "recent_ood",
) -> np.ndarray:
    """PO^{1/3} weights on the OOD batch (softer than √PO)."""
    return current_batch_po_weights(
        stream,
        t,
        seed=seed,
        n_control=n_control,
        n_recent=n_recent,
        blend_mu_gap=blend_mu_gap,
        mode="cbrt",
        window_mode=window_mode,
    )


def describe_windows(windows: RefitWindows) -> str:
    """One-line human summary of a cut."""
    return (
        f"mode={windows.window_mode} recent={list(windows.recent_idx)} "
        f"ood={windows.ood_idx} "
        f"|R|={len(windows.y_recent)} |O|={len(windows.y_ood)}"
    )
