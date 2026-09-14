"""Online RF probe → gated √po_risk0, plus PO-tail localization.

This is the IPTW stream logic in ``run_agod_po_ood_stream_multi`` /
``run_agod_po_iptw_mse``, with one change: do **not** apply √PO every
batch. The shallow RF on 上一批 is T=0. Instance

  po_risk0_i = |Y_i − μ0(X_i)|   (mixed with batch gap, same as IPTW)

and ``w = np.sqrt(po_risk0)`` only when the *consecutive OOS* probe
error jumps — in-sample e1>e0 always holds for trees and would fire
every hop.

  e_now  = err(μ0 fitted on B_{t-1}, scored on B_t)
  e_prev = err(μ0 fitted on B_{t-2}, scored on B_{t-1})
  fire iff e_now / e_prev ≥ γ   (default 1.5; skip the first hop)

Quiet → last two batches, uniform. Fire → current batch with √po_risk0
or the high-po_risk0 tail (post-hoc subset localization).
"""
from __future__ import annotations

import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import accuracy_score

from agod.po_iptw import instance_po_risk, po_iptw_weights
from agod.po_refit import (
    Stream,
    _pack_stream,
    _stream_task,
    assign_hop,
    fit_predict,
    hop_score,
)


def fit_online_rf(X, y, *, seed=0, task="mse"):
    """Same shallow RF as the IPTW stream probe (not the next-step model)."""
    X = np.asarray(X, dtype=float)
    y = np.asarray(y).ravel()
    if str(task) == "acc":
        clf = RandomForestClassifier(
            n_estimators=20,
            max_depth=4,
            min_samples_leaf=2,
            random_state=int(seed),
            n_jobs=1,
        )
        clf.fit(X, y.astype(int))
        return clf
    rf = RandomForestRegressor(
        n_estimators=20,
        max_depth=4,
        min_samples_leaf=3,
        random_state=int(seed),
        n_jobs=1,
    )
    rf.fit(X, y.astype(float))
    return rf


fit_online_probe = fit_online_rf


def probe_err(model, X, y, task="mse"):
    X = np.asarray(X, dtype=float)
    y = np.asarray(y).ravel()
    if str(task) == "acc":
        pred = model.predict(X).astype(int)
        return 1.0 - float(accuracy_score(y.astype(int), pred))
    pred = np.asarray(model.predict(X), dtype=float).ravel()
    return float(np.mean(np.abs(y.astype(float) - pred)))


score_probe = probe_err


def hop_fires(e_now, e_prev, gate=1.5):
    """Consecutive OOS gate. First hop has no previous OOS → quiet."""
    if e_prev is None:
        return False
    ratio = shift_ratio(e_now, e_prev)
    return bool(np.isfinite(ratio) and ratio >= float(gate))


def po_risk0_rows(model, X, y, *, batch_po=0.0, task="mse"):
    """Instance PO vs T=0 probe. This is po_risk0 in the IPTW stream."""
    X = np.asarray(X, dtype=float)
    y = np.asarray(y).ravel()
    if str(task) == "acc":
        y = y.astype(int)
        if hasattr(model, "predict_proba"):
            proba = model.predict_proba(X)
            classes = list(model.classes_)
            p_true = np.zeros(len(y), dtype=float)
            for i, yi in enumerate(y):
                if yi in classes:
                    p_true[i] = float(proba[i, classes.index(yi)])
            po = 1.0 - p_true
        else:
            po = (y != model.predict(X).astype(int)).astype(float)
        return instance_po_risk(po, np.zeros_like(po), batch_po=batch_po, mix=0.5)
    pred = np.asarray(model.predict(X), dtype=float).ravel()
    return instance_po_risk(y, pred, batch_po=batch_po, mix=0.5)


def po_tail_mask(risk, q=0.30, min_n=8):
    """High-PO tail inside the current batch (post-hoc subset)."""
    r = np.asarray(risk, dtype=float).ravel()
    n = r.size
    if n == 0:
        return np.zeros(0, dtype=bool)
    keep = max(int(min_n), int(np.ceil(float(q) * n)))
    keep = min(keep, n)
    out = np.zeros(n, dtype=bool)
    if keep >= n:
        out[:] = True
        return out
    out[np.argpartition(r, -keep)[-keep:]] = True
    return out


def shift_ratio(e_now, e_prev):
    return float(e_now) / (float(e_prev) + 1e-8)


def run_rfperm_stream(
    stream: Stream,
    *,
    gate=1.5,
    localize=False,
    q=0.30,
    learner="rf",
    seed=0,
):
    """Gated online-RF √po_risk0. ``localize`` trains on the PO tail only."""
    task = _stream_task(stream)
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y).ravel()
    batch = np.asarray(stream.batch, dtype=int)
    k = int(batch.max()) + 1
    history = []
    e_prev = None
    for t in range(1, k - 1):
        t0, t1 = assign_hop(batch, t)
        te = batch == (t + 1)
        if not np.any(t0) or not np.any(t1) or not np.any(te):
            continue
        probe = fit_online_rf(X[t0], y[t0], seed=int(seed) + t, task=task)
        e_now = probe_err(probe, X[t1], y[t1], task=task)
        fired = hop_fires(e_now, e_prev, gate=gate)
        if e_prev is None:
            ratio, po_b = 1.0, 0.0
        else:
            ratio = shift_ratio(e_now, e_prev)
            po_b = max(float(e_now) - float(e_prev), 0.0)
        po0 = po_risk0_rows(probe, X[t1], y[t1], batch_po=po_b, task=task)
        if fired:
            if localize:
                tail = po_tail_mask(po0, q=q)
                tr = np.zeros(t1.shape[0], dtype=bool)
                tr[np.flatnonzero(t1)[tail]] = True
                if str(task) == "acc" and np.unique(y[tr]).size < 2:
                    # one-class tail: keep the hop (drop old batch) but don't subset
                    tr = t1
                w_tr = None
            else:
                tr = t1
                w_tr = po_iptw_weights(po0, mode="sqrt")
        else:
            tr = (batch == (t - 1)) | (batch == t)
            w_tr = None
        pred = fit_predict(
            X[tr], y[tr], X[te], w=w_tr, learner=learner, seed=int(seed) + t, task=task
        )
        mse = hop_score(y[te], pred, task=task)
        history.append(
            {
                "t": int(t),
                "next_mse": mse,
                "fired": bool(fired),
                "ratio": float(ratio),
                "lambda": 1.0 if fired else 0.0,
                "mean_r0": float("nan") if e_prev is None else float(e_prev),
                "mean_r1": float(e_now),
                "n_train": int(tr.sum()),
                "n_test": int(te.sum()),
                "q_tail": float(q),
                "localize": bool(localize),
            }
        )
        e_prev = e_now
    tag = "local" if localize else "rfperm"
    return _pack_stream(
        history,
        assign=tag,
        gate=gate,
        always=False,
        mode=tag,
        learner=learner,
        metric=task,
    )
