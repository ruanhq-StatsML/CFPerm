"""Rolling PO-learner refit → gated sample weights.

Uniform is the default. Re-adjust only when the new batch is clearly
different. When batch t arrives, *refit* the PO-learner (do not reuse a
frozen probe):

  T = 0  most recent control  (batch t-2, or batch t-1 on the first hop)
  T = 1  上一批 ∪ 这一批     (batch t-1 and batch t)

Instance PO-risk r_i = φ_i² on T=1 rows becomes √PO weights iff the batch
contrast clears a gate; otherwise w = 1.

Next-step model trains on the T=1 rows and is scored on batch t+1.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from sklearn.linear_model import Ridge

from agod.po_iptw import dre_weights, po_iptw_weights

Assign = Literal["hop", "pair"]


def dr_pseudo_outcome(X, y, t, clip=0.05, ridge_alpha=3.0):
    """DR φ = (μ1−μ0) + (T−π)(Y−μ_T)/(π(1−π)). π = clip(mean T)."""
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float).ravel()
    t = np.asarray(t, dtype=int).ravel()
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    n0 = int((t == 0).sum())
    n1 = int((t == 1).sum())
    if n0 < 3 or n1 < 3:
        z = np.zeros(y.shape[0], dtype=float)
        return z, {"ok": False, "pi": 0.5, "mu0": z, "mu1": z}
    mu0 = Ridge(alpha=float(ridge_alpha)).fit(X[t == 0], y[t == 0]).predict(X)
    mu1 = Ridge(alpha=float(ridge_alpha)).fit(X[t == 1], y[t == 1]).predict(X)
    pi = float(np.clip(n1 / float(n0 + n1), clip, 1.0 - clip))
    mu_t = np.where(t == 1, mu1, mu0)
    phi = (mu1 - mu0) + (t - pi) * (y - mu_t) / (pi * (1.0 - pi))
    return np.asarray(phi, dtype=float), {
        "ok": True,
        "pi": pi,
        "mu0": np.asarray(mu0, dtype=float),
        "mu1": np.asarray(mu1, dtype=float),
    }


def instance_phi2(phi):
    return np.square(np.asarray(phi, dtype=float).ravel())


def batch_contrast(risk, t):
    """How different T=1 looks vs T=0: mean r1 / (mean r0 + eps)."""
    t = np.asarray(t, dtype=int).ravel()
    r = np.asarray(risk, dtype=float).ravel()
    r0 = r[t == 0]
    r1 = r[t == 1]
    m0 = float(r0.mean()) if r0.size else 0.0
    m1 = float(r1.mean()) if r1.size else 0.0
    return m1 / (m0 + 1e-8), m0, m1


def should_readjust(ratio, gate=1.25):
    """Re-adjust iff T=1 PO-risk is clearly larger than T=0."""
    return bool(np.isfinite(ratio) and ratio >= float(gate))


def assign_hop(batch, t):
    """T=0: 上一批 (t-1). T=1: 这一批 (t)."""
    batch = np.asarray(batch, dtype=int)
    t0 = batch == int(t) - 1
    t1 = batch == int(t)
    return t0, t1


def assign_pair(batch, t):
    """T=0: most recent leftover (t-2). T=1: 上一批 ∪ 这一批."""
    batch = np.asarray(batch, dtype=int)
    if int(t) < 2:
        return assign_hop(batch, t)
    t0 = batch == int(t) - 2
    t1 = (batch == int(t) - 1) | (batch == int(t))
    return t0, t1


def mix_lambda(ratio, gate=1.25, soft_scale=2.0):
    """Soft gate: 0 below γ, then ramp to 1 over ``soft_scale * γ`` extra ratio."""
    g = float(gate)
    den = float(soft_scale) * g
    if den <= 0:
        return 1.0 if float(ratio) >= g else 0.0
    return float(np.clip((float(ratio) - g) / den, 0.0, 1.0))


def residual_hop_ratio(X, y, t0, t1, ridge_alpha=3.0):
    """How badly a model fit on 上一批 predicts 这一批.

    Global concept flip raises this a lot; a global flip does *not*
    raise mean φ²_{T=1} / mean φ²_{T=0} (both sides jump together).
    """
    t0 = np.asarray(t0, dtype=bool).ravel()
    t1 = np.asarray(t1, dtype=bool).ravel()
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float).ravel()
    if int(t0.sum()) < 3 or int(t1.sum()) < 3:
        return 1.0, 0.0, 0.0
    clf = Ridge(alpha=float(ridge_alpha)).fit(X[t0], y[t0])
    mse0 = float(np.mean((clf.predict(X[t0]) - y[t0]) ** 2))
    mse1 = float(np.mean((clf.predict(X[t1]) - y[t1]) ** 2))
    return mse1 / (mse0 + 1e-8), mse0, mse1


def refit_po_weights(
    X,
    y,
    t0,
    t1,
    *,
    gate=1.25,
    always=False,
    mode="sqrt",
    ridge_alpha=3.0,
    soft=False,
    soft_scale=2.0,
):
    """Refit PO-learner on T0 vs T1; return weights on T1 (and gate flag).

    Hard gate: uniform unless ρ ≥ γ, then full √PO.
    Soft gate: w = (1-λ) + λ w_PO with λ = clip((ρ-γ)/(2γ), 0, 1).
    """
    t0 = np.asarray(t0, dtype=bool).ravel()
    t1 = np.asarray(t1, dtype=bool).ravel()
    on = t0 | t1
    T = np.zeros(int(on.sum()), dtype=int)
    T[t1[on]] = 1
    phi, nuis = dr_pseudo_outcome(X[on], y[on], T, ridge_alpha=ridge_alpha)
    risk_on = instance_phi2(phi)
    ratio, m0, m1 = batch_contrast(risk_on, T)
    w_all = np.ones(t1.shape[0], dtype=float)
    lam = 0.0
    if nuis["ok"]:
        if always:
            lam = 1.0
        elif soft:
            lam = mix_lambda(ratio, gate=gate, soft_scale=soft_scale)
        elif should_readjust(ratio, gate=gate):
            lam = 1.0
        if lam > 0:
            r_full = np.zeros(t1.shape[0], dtype=float)
            r_full[on] = risk_on
            w_po = po_iptw_weights(r_full[t1], mode=mode)  # type: ignore[arg-type]
            mixed = (1.0 - lam) + lam * w_po
            w_all[t1] = mixed / (mixed.mean() + 1e-8)
    fired = bool(lam > 0)
    return {
        "weights": w_all,
        "fired": fired,
        "ok": nuis["ok"],
        "ratio": float(ratio),
        "lambda": float(lam),
        "mean_r0": float(m0),
        "mean_r1": float(m1),
        "pi": float(nuis["pi"]),
    }


@dataclass
class Stream:
    X: np.ndarray
    y: np.ndarray
    batch: np.ndarray
    name: str = ""
    meta: dict | None = None


def make_batch_stream(
    *,
    n_batches=8,
    n_per=100,
    p=12,
    seed=0,
    cov=0.0,
    concept_at=None,
    concept=0.0,
    noise=0.45,
    rank=4,
):
    """Same P(Y|X) until ``concept_at``, then flip the rating map.

    cov>0 adds a per-batch mean hop (covariate). concept_at=None → similar.
    """
    rng = np.random.default_rng(seed)
    beta = np.zeros(p)
    beta[: int(rank)] = 0.85
    rows, ys, batches = [], [], []
    dirs = rng.normal(size=(int(n_batches), p))
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True) + 1e-12
    for b in range(int(n_batches)):
        X = rng.normal(size=(int(n_per), p))
        if cov:
            X = X + float(cov) * dirs[b] * np.sqrt(p)
        use = -beta if (concept_at is not None and b >= int(concept_at) and concept) else beta
        y = 3.0 + X @ use + rng.normal(scale=float(noise), size=int(n_per))
        y = np.clip(y, 1.0, 5.0)
        rows.append(X)
        ys.append(y)
        batches.append(np.full(int(n_per), b, dtype=int))
    return Stream(
        X=np.vstack(rows),
        y=np.concatenate(ys),
        batch=np.concatenate(batches),
        name="synth",
        meta={
            "cov": float(cov),
            "concept": float(concept),
            "concept_at": None if concept_at is None else int(concept_at),
            "seed": int(seed),
            "n_batches": int(n_batches),
            "n_per": int(n_per),
        },
    )


def _ridge_predict(Xtr, ytr, Xte, w=None, alpha=3.0):
    clf = Ridge(alpha=float(alpha))
    if w is None:
        clf.fit(Xtr, ytr)
    else:
        clf.fit(Xtr, ytr, sample_weight=w)
    return clf.predict(Xte)


def _pack_stream(history, *, assign, gate, always, mode):
    mses = np.array([h["next_mse"] for h in history], dtype=float)
    fires = np.array([h["fired"] for h in history], dtype=float)
    return {
        "assign": assign,
        "gate": float(gate),
        "always": bool(always),
        "mode": mode,
        "online_mse": float(mses.mean()) if mses.size else float("nan"),
        "mse_std": float(mses.std()) if mses.size else float("nan"),
        "fire_rate": float(fires.mean()) if fires.size else 0.0,
        "path": mses.tolist(),
        "history": history,
        "n_hops": int(mses.size),
    }


def run_refit_stream(
    stream: Stream,
    *,
    assign: Assign = "pair",
    gate=1.25,
    always=False,
    mode="sqrt",
    ridge_alpha=3.0,
    soft=False,
    soft_scale=2.0,
):
    """Online next-batch MSE. Train on T=1 rows, score batch t+1."""
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    k = int(batch.max()) + 1
    assign_fn = assign_pair if assign == "pair" else assign_hop
    history = []
    for t in range(1, k - 1):
        t0, t1 = assign_fn(batch, t)
        rec = refit_po_weights(
            X,
            y,
            t0,
            t1,
            gate=gate,
            always=always,
            mode=mode,
            ridge_alpha=ridge_alpha,
            soft=soft,
            soft_scale=soft_scale,
        )
        tr = t1
        te = batch == (t + 1)
        if not np.any(tr) or not np.any(te):
            continue
        pred = _ridge_predict(
            X[tr],
            y[tr],
            X[te],
            w=rec["weights"][tr],
            alpha=ridge_alpha,
        )
        mse = float(np.mean((pred - y[te]) ** 2))
        history.append(
            {
                "t": int(t),
                "next_mse": mse,
                "fired": rec["fired"],
                "ratio": rec["ratio"],
                "lambda": rec["lambda"],
                "mean_r0": rec["mean_r0"],
                "mean_r1": rec["mean_r1"],
                "n_train": int(tr.sum()),
                "n_test": int(te.sum()),
            }
        )
    tag = "uniform" if (not always and not soft and gate >= 1e6) else mode
    if soft and not always:
        tag = f"soft_{mode}"
    return _pack_stream(history, assign=assign, gate=gate, always=always, mode=tag)


def run_uniform_on_same_rows(stream: Stream, assign: Assign = "pair", ridge_alpha=3.0):
    """Same train rows as refit (T=1), but w=1. Fair uniform baseline."""
    return run_refit_stream(
        stream,
        assign=assign,
        gate=1e9,
        always=False,
        mode="sqrt",
        ridge_alpha=ridge_alpha,
    )


def run_adaptive_stream(
    stream: Stream,
    *,
    gate=1.25,
    mode="sqrt",
    ridge_alpha=3.0,
):
    """Similar hops: pool last two, uniform. Different hops: train on new batch + √PO.

    Gate is the hop contrast (T=0=上一批, T=1=这一批). That is the
    re-adjustment the quiet stream was missing.
    """
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    k = int(batch.max()) + 1
    history = []
    for t in range(1, k - 1):
        t0, t1 = assign_hop(batch, t)
        rec = refit_po_weights(
            X,
            y,
            t0,
            t1,
            gate=gate,
            always=False,
            mode=mode,
            ridge_alpha=ridge_alpha,
        )
        te = batch == (t + 1)
        if rec["fired"]:
            tr = t1
            w_tr = rec["weights"][tr]
        else:
            tr = (batch == (t - 1)) | (batch == t)
            w_tr = None
        if not np.any(tr) or not np.any(te):
            continue
        pred = _ridge_predict(X[tr], y[tr], X[te], w=w_tr, alpha=ridge_alpha)
        mse = float(np.mean((pred - y[te]) ** 2))
        history.append(
            {
                "t": int(t),
                "next_mse": mse,
                "fired": rec["fired"],
                "ratio": rec["ratio"],
                "lambda": rec.get("lambda", 1.0 if rec["fired"] else 0.0),
                "mean_r0": rec["mean_r0"],
                "mean_r1": rec["mean_r1"],
                "n_train": int(tr.sum()),
                "n_test": int(te.sum()),
            }
        )
    return _pack_stream(
        history, assign="adaptive", gate=gate, always=False, mode="adaptive"
    )


def run_switch_stream(stream: Stream, *, gate=1.25, ridge_alpha=3.0):
    """Same train-set switch as adaptive, but always uniform (no √PO).

    Isolates 'drop the old batch' from 'reweight the new batch'.
    """
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    k = int(batch.max()) + 1
    history = []
    for t in range(1, k - 1):
        t0, t1 = assign_hop(batch, t)
        rec = refit_po_weights(
            X, y, t0, t1, gate=gate, always=False, ridge_alpha=ridge_alpha
        )
        te = batch == (t + 1)
        tr = t1 if rec["fired"] else ((batch == (t - 1)) | (batch == t))
        if not np.any(tr) or not np.any(te):
            continue
        pred = _ridge_predict(X[tr], y[tr], X[te], w=None, alpha=ridge_alpha)
        mse = float(np.mean((pred - y[te]) ** 2))
        history.append(
            {
                "t": int(t),
                "next_mse": mse,
                "fired": rec["fired"],
                "ratio": rec["ratio"],
                "lambda": rec["lambda"],
                "mean_r0": rec["mean_r0"],
                "mean_r1": rec["mean_r1"],
                "n_train": int(tr.sum()),
                "n_test": int(te.sum()),
            }
        )
    return _pack_stream(
        history, assign="switch", gate=gate, always=False, mode="switch"
    )


def run_dre_hop(stream: Stream, *, ridge_alpha=3.0):
    """X-only density-ratio weights on the new batch. Ignores Y-shift."""
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    k = int(batch.max()) + 1
    history = []
    for t in range(1, k - 1):
        prev = batch == (t - 1)
        tr = batch == t
        te = batch == (t + 1)
        if not np.any(prev) or not np.any(tr) or not np.any(te):
            continue
        w = dre_weights(X[prev], X[tr], seed=int(t))
        pred = _ridge_predict(X[tr], y[tr], X[te], w=w, alpha=ridge_alpha)
        mse = float(np.mean((pred - y[te]) ** 2))
        history.append(
            {
                "t": int(t),
                "next_mse": mse,
                "fired": True,
                "ratio": float("nan"),
                "lambda": 1.0,
                "mean_r0": float("nan"),
                "mean_r1": float("nan"),
                "n_train": int(tr.sum()),
                "n_test": int(te.sum()),
            }
        )
    return _pack_stream(
        history, assign="hop", gate=0.0, always=True, mode="dre"
    )


def run_oracle_switch(stream: Stream, *, ridge_alpha=3.0):
    """Knows concept_at: after the cut, train on the new batch only.

    Diagnostic upper bound for *train-set* choice, not a deployable method.
    Quiet / covariate streams have no cut → always pair-uniform.
    """
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    k = int(batch.max()) + 1
    cut = None if stream.meta is None else stream.meta.get("concept_at")
    history = []
    for t in range(1, k - 1):
        te = batch == (t + 1)
        fired = cut is not None and int(t) >= int(cut)
        tr = (batch == t) if fired else ((batch == (t - 1)) | (batch == t))
        if not np.any(tr) or not np.any(te):
            continue
        pred = _ridge_predict(X[tr], y[tr], X[te], w=None, alpha=ridge_alpha)
        mse = float(np.mean((pred - y[te]) ** 2))
        history.append(
            {
                "t": int(t),
                "next_mse": mse,
                "fired": bool(fired),
                "ratio": float("nan"),
                "lambda": 1.0 if fired else 0.0,
                "mean_r0": float("nan"),
                "mean_r1": float("nan"),
                "n_train": int(tr.sum()),
                "n_test": int(te.sum()),
            }
        )
    return _pack_stream(
        history, assign="oracle", gate=0.0, always=False, mode="oracle"
    )


def run_resid_stream(
    stream: Stream,
    *,
    gate=2.0,
    po_on_fire=False,
    po_gate=1.25,
    mode="sqrt",
    ridge_alpha=3.0,
):
    """Re-adjust when 上一批→这一批 residual MSE jumps (not PO-ratio).

    Quiet: train pair, uniform. Fire: drop the old batch.
    ``po_on_fire`` then puts √PO weights on the new batch (refit PO-learner).
    """
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    k = int(batch.max()) + 1
    history = []
    for t in range(1, k - 1):
        t0, t1 = assign_hop(batch, t)
        rho, mse0, mse1 = residual_hop_ratio(X, y, t0, t1, ridge_alpha=ridge_alpha)
        fired = bool(np.isfinite(rho) and rho >= float(gate))
        te = batch == (t + 1)
        if fired:
            tr = t1
            w_tr = None
            if po_on_fire:
                rec = refit_po_weights(
                    X,
                    y,
                    t0,
                    t1,
                    gate=po_gate,
                    always=True,
                    mode=mode,
                    ridge_alpha=ridge_alpha,
                )
                w_tr = rec["weights"][tr]
        else:
            tr = (batch == (t - 1)) | (batch == t)
            w_tr = None
        if not np.any(tr) or not np.any(te):
            continue
        pred = _ridge_predict(X[tr], y[tr], X[te], w=w_tr, alpha=ridge_alpha)
        mse = float(np.mean((pred - y[te]) ** 2))
        history.append(
            {
                "t": int(t),
                "next_mse": mse,
                "fired": fired,
                "ratio": float(rho),
                "lambda": 1.0 if fired else 0.0,
                "mean_r0": float(mse0),
                "mean_r1": float(mse1),
                "n_train": int(tr.sum()),
                "n_test": int(te.sum()),
            }
        )
    tag = "resid_po" if po_on_fire else "resid"
    return _pack_stream(history, assign=tag, gate=gate, always=False, mode=tag)
