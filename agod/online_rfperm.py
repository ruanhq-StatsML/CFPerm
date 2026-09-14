"""Online RFPerm + PO-risk (the AGOD stream method).

Probe = the same shallow RF as ``run_agod_po_ood_stream_multi``
(``n_estimators=20``, ``max_depth=4``) on 上一批 as T=0.
Score = instance ``po_risk0`` from ``instance_po_risk``.
Weight = ``w = sqrt(po_risk0)`` on T=1, mean 1.

Not every batch. Consecutive OOS (in-sample e1>e0 fires every hop
on trees):

  e_now  = err(μ0 fitted on B_{t-1}, scored on B_t)
  e_prev = err(μ0 fitted on B_{t-2}, scored on B_{t-1})
  fire iff e_now / e_prev ≥ γ   (default 1.5; skip the first hop)

Quiet → last two batches, w=1.
Fire → same rows; T=0 stays 1, T=1 gets √po_risk0.
Then the next hop cools (anneal). No subset localization.
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


def shift_ratio(e_now, e_prev):
    return float(e_now) / (float(e_prev) + 1e-8)


def _quantiles(x):
    """Per-vector quantification: mean / spread / tails."""
    a = np.asarray(x, dtype=float).ravel()
    a = a[np.isfinite(a)]
    if a.size == 0:
        return {
            "n": 0,
            "mean": float("nan"),
            "std": float("nan"),
            "min": float("nan"),
            "p10": float("nan"),
            "p25": float("nan"),
            "p50": float("nan"),
            "p75": float("nan"),
            "p90": float("nan"),
            "max": float("nan"),
        }
    qs = np.quantile(a, [0.10, 0.25, 0.50, 0.75, 0.90])
    return {
        "n": int(a.size),
        "mean": float(a.mean()),
        "std": float(a.std()),
        "min": float(a.min()),
        "p10": float(qs[0]),
        "p25": float(qs[1]),
        "p50": float(qs[2]),
        "p75": float(qs[3]),
        "p90": float(qs[4]),
        "max": float(a.max()),
    }


def quantify_last_two(
    probe,
    X,
    y,
    batch,
    t,
    *,
    task="mse",
    batch_po=0.0,
    fired=False,
):
    """Per-observation PO-risk and IPTW weight on B_{t-1} ∪ B_t.

    Every row gets ``po_risk0_i`` vs the T=0 probe. Quiet → ``w_i=1``.
    Fire → T=0 stays 1, T=1 gets ``w_i = √po_risk0`` (mean 1).
    """
    batch = np.asarray(batch, dtype=int)
    tr = (batch == int(t) - 1) | (batch == int(t))
    idx = np.flatnonzero(tr)
    treated = (batch[tr] == int(t)).astype(float)
    po = po_risk0_rows(probe, X[tr], y[tr], batch_po=batch_po, task=task)
    if fired:
        w = po_iptw_weights(po, mode="sqrt", treated=treated)
    else:
        w = np.ones(idx.size, dtype=float)
    t1 = treated > 0.5
    return {
        "index": idx,
        "batch": batch[tr],
        "treated": t1.astype(int),
        "po_risk0": np.asarray(po, dtype=float),
        "w": np.asarray(w, dtype=float),
        "fired": bool(fired),
        "po_t1": _quantiles(po[t1]),
        "po_t0": _quantiles(po[~t1]),
        "w_t1": _quantiles(w[t1]),
        "w_t0": _quantiles(w[~t1]),
        "w_all": _quantiles(w),
    }


def last_two_sqrt_weights(batch, t, po1):
    """Last two batches: T=0 w=1, T=1 w=√po_risk0, mean 1."""
    batch = np.asarray(batch, dtype=int)
    tr = (batch == int(t) - 1) | (batch == int(t))
    treated = (batch[tr] == int(t)).astype(float)
    po = np.ones(int(tr.sum()), dtype=float)
    po[treated > 0.5] = np.asarray(po1, dtype=float).ravel()
    w = po_iptw_weights(po, mode="sqrt", treated=treated)
    return tr, w


def run_rfperm_stream(
    stream: Stream,
    *,
    gate=1.5,
    learner="rf",
    seed=0,
    detail=False,
):
    """Online RFPerm + PO-risk. Last two batches; √po_risk0 on T=1 iff gated.

    Every hop records per-observation PO/weight *quantiles*. ``detail=True``
    also keeps the full ``po_risk0`` and ``w`` vectors on last-two rows.
    """
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
        q = quantify_last_two(
            probe, X, y, batch, t, task=task, batch_po=po_b, fired=fired
        )
        tr = (batch == (t - 1)) | (batch == t)
        w_tr = q["w"] if fired else None
        pred = fit_predict(
            X[tr], y[tr], X[te], w=w_tr, learner=learner, seed=int(seed) + t, task=task
        )
        mse = hop_score(y[te], pred, task=task)
        rec = {
            "t": int(t),
            "next_mse": mse,
            "fired": bool(fired),
            "ratio": float(ratio),
            "lambda": 1.0 if fired else 0.0,
            "mean_r0": float("nan") if e_prev is None else float(e_prev),
            "mean_r1": float(e_now),
            "n_train": int(tr.sum()),
            "n_test": int(te.sum()),
            "mean_w": float(q["w_all"]["mean"]),
            "po_t1": q["po_t1"],
            "w_t1": q["w_t1"],
        }
        if detail:
            rec["po_risk0"] = q["po_risk0"].tolist()
            rec["w"] = q["w"].tolist()
            rec["treated"] = q["treated"].tolist()
            rec["index"] = q["index"].tolist()
        history.append(rec)
        e_prev = e_now
    return _pack_stream(
        history,
        assign="rfperm",
        gate=gate,
        always=False,
        mode="sqrt",
        learner=learner,
        metric=task,
    )


run_online_rfperm = run_rfperm_stream
