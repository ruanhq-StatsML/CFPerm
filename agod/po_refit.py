"""Streaming PO-learner re-fit for post-hoc IPTW (after OnlineRFPerm reject).

When a new batch t arrives and the shift test is significant:

  T=0  : recent control window  — batch(es) just before the pair
  T=1  : previous + current     — batch_{t-1} ∪ batch_t

Re-fit a green PO-risk proxy on this labeled pool, then use
  w_i ∝ √PO(X_i, Y_i | T=1)
only on the current batch for the downstream learner.

PO proxy (green, no DRE):
  fit μ0 on T=0;  PO_i = |Y_i − μ0(X_i)|  for i in T=1
  (optionally blend with |μ1(X_i)−μ0(X_i)| if μ1 is fit on T=1)

This matches: uniform most of the time; re-adjust only on clearly
different batches; PO-learner always refreshed from the latest T=0/1 cut.
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
from sklearn.ensemble import RandomForestRegressor

from agod.po_iptw import po_iptw_weights


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


def build_t01_windows(
    stream: list,
    t: int,
    *,
    n_control: int = 1,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Construct (X0,y0,X1,y1) at stream index t (current = stream[t]).

    T=0: up to ``n_control`` *recent* batches immediately before the T=1 pair
         (indices ending at t-2). Cold start (t==1): reuse stream[0] as T=0
         while T=1 still = stream[0]∪stream[1] (overlap is acceptable as a
         weak control; gate usually off early anyway).
    T=1: stream[t-1] ∪ stream[t]
    """
    if t < 1:
        raise ValueError("need t>=1 for prev+cur as T=1")
    Xp, yp = stream[t - 1]
    Xc, yc = stream[t]
    X1 = np.vstack([Xp, Xc])
    y1 = np.concatenate([np.asarray(yp, float), np.asarray(yc, float)])

    # recent control: batches just before the (prev, cur) pair
    ctrl_idx = [i for i in range(max(0, t - 1 - n_control), t - 1)]
    if not ctrl_idx:
        ctrl_idx = [0]

    Xs, ys = [], []
    for i in ctrl_idx:
        Xi, yi = stream[i]
        Xs.append(Xi)
        ys.append(np.asarray(yi, float))
    X0 = np.vstack(Xs)
    y0 = np.concatenate(ys)
    return X0, y0, X1, y1


def refit_po_risk_t01(
    X0: np.ndarray,
    y0: np.ndarray,
    X1: np.ndarray,
    y1: np.ndarray,
    *,
    seed: int = 0,
    blend_mu_gap: float = 0.25,
) -> np.ndarray:
    """Re-fit PO proxy with T=0 control vs T=1 (prev+cur).

    Returns per-row PO-risk on the T=1 pool (len = len(X1)).
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


def current_batch_sqrt_po_weights(
    stream: list,
    t: int,
    *,
    seed: int = 0,
    n_control: int = 1,
    blend_mu_gap: float = 0.25,
) -> np.ndarray:
    """Post-hoc √PO weights for *current* batch only, after T=0/T=1 re-fit.

    T=1 = prev∪cur; weights returned have length = len(current batch),
    taken from the second half of the T=1 PO vector.
    """
    X0, y0, X1, y1 = build_t01_windows(stream, t, n_control=n_control)
    po1 = refit_po_risk_t01(X0, y0, X1, y1, seed=seed, blend_mu_gap=blend_mu_gap)
    n_cur = len(stream[t][1])
    po_cur = po1[-n_cur:]
    return po_iptw_weights(po_cur, mode="sqrt")
