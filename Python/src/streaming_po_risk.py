"""Streaming PO-risk on a tabular table.

New batch is T=1, reference is T=0. Outcome model μ(Y|X) and
propensity e(T|X) are maintained separately — not the serving MLP,
not a single lstsq. φ = (Y − μ)(T − e), τ̂(X) ≈ φ, risk = mean(τ̂²).

Read the number. Raw PO-risk on a small n_new jitters; a causal
moving average is the stability readout. No online-bootstrap —
repeated MLP inference cannot be afforded.
"""
from __future__ import annotations

import numpy as np
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.preprocessing import StandardScaler

CLIP = 1e-3
REF_N = 10_000
MIN_STREAM_N = 5_000
# Stream vs ref-split baseline. Below this, all layers stay trainable.
DEVIATION_RATIO = 2.0
ACTION_KEEP = "keep_training"
ACTION_WATCH = "watch"
ACTION_FREEZE = "freeze"


def zscore(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    sd = X.std(axis=0)
    sd = np.where(sd < 1e-8, 1.0, sd)
    return (X - X.mean(axis=0)) / sd


def pack_ref_new(X_ref, Y_ref, X_new, Y_new):
    """Concatenate reference (T=0) and the incoming batch (T=1)."""
    X_ref = np.asarray(X_ref, dtype=float)
    X_new = np.asarray(X_new, dtype=float)
    Y = np.concatenate(
        [np.asarray(Y_ref, dtype=float).ravel(), np.asarray(Y_new, dtype=float).ravel()]
    )
    T = np.concatenate([np.zeros(len(X_ref)), np.ones(len(X_new))])
    X = np.vstack([X_ref, X_new])
    return X, Y, T


def _is_binary(Y: np.ndarray) -> bool:
    u = np.unique(np.asarray(Y, dtype=float).ravel())
    return u.size <= 2 and set(np.round(u, 8)).issubset({0.0, 1.0})


class TabularPORisk:
    """Own outcome model + own propensity model."""

    def __init__(self, clip: float = CLIP, seed: int = 2026):
        self.clip = float(clip)
        self.seed = int(seed)
        self.outcome = None
        self.propensity = None
        self.scaler = StandardScaler()

    def fit_nuisance(self, X, Y, T, X_outcome=None):
        """Fit μ(Y|X_outcome) and e(T|X). X_outcome defaults to X."""
        X = np.asarray(X, dtype=float)
        if X.ndim == 1:
            X = X.reshape(-1, 1)
        Y = np.asarray(Y, dtype=float).ravel()
        T = np.asarray(T, dtype=float).ravel().astype(int)
        Xo = X if X_outcome is None else np.asarray(X_outcome, dtype=float)
        if Xo.ndim == 1:
            Xo = Xo.reshape(-1, 1)
        Xs = self.scaler.fit_transform(X)
        self.propensity = LogisticRegression(max_iter=300, random_state=self.seed)
        self.propensity.fit(Xs, T)
        e = np.clip(self.propensity.predict_proba(Xs)[:, 1], self.clip, 1.0 - self.clip)
        scale_y = StandardScaler()
        Xos = scale_y.fit_transform(Xo)
        if _is_binary(Y):
            self.outcome = LogisticRegression(max_iter=300, random_state=self.seed)
            self.outcome.fit(Xos, Y.astype(int))
            mu = self.outcome.predict_proba(Xos)[:, 1]
        else:
            self.outcome = Ridge(alpha=1.0)
            self.outcome.fit(Xos, Y)
            mu = self.outcome.predict(Xos)
        self._outcome_scaler = scale_y
        return mu, e

    def tau_risk(self, X, Y, T, mu, e) -> float:
        X = np.asarray(X, dtype=float)
        if X.ndim == 1:
            X = X.reshape(-1, 1)
        Y = np.asarray(Y, dtype=float).ravel()
        T = np.asarray(T, dtype=float).ravel()
        mu = np.asarray(mu, dtype=float).ravel()
        e = np.clip(np.asarray(e, dtype=float).ravel(), self.clip, 1.0 - self.clip)
        phi = (Y - mu) * (T - e)
        Xd = np.column_stack([np.ones(len(X)), zscore(X)])
        tau, *_ = np.linalg.lstsq(Xd, phi, rcond=None)
        tau_hat = Xd @ tau
        return float(np.mean(tau_hat**2))

    def risk(self, X, Y, T, X_outcome=None) -> dict:
        mu, e = self.fit_nuisance(X, Y, T, X_outcome=X_outcome)
        value = self.tau_risk(X, Y, T, mu, e)
        return {"po_risk": value, "mu": mu, "e": e}


def po_risk(X, Y, T, mu=None, clip: float = CLIP, seed: int = 2026) -> float:
    """Tabular PO-risk. μ may be passed (conditional on a serving model);
    e is always a separate propensity. If μ is None, a separate outcome
    model is fit too.
    """
    est = TabularPORisk(clip=clip, seed=seed)
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    T = np.asarray(T, dtype=float).ravel()
    if mu is None:
        return float(est.risk(X, Y, T)["po_risk"])
    mu = np.asarray(mu, dtype=float).ravel()
    X_outcome = np.column_stack([X if X.ndim == 2 else X.reshape(-1, 1), mu.reshape(-1, 1)])
    return float(est.risk(X, Y, T, X_outcome=X_outcome)["po_risk"])


def batch_mse(Y, mu) -> float:
    """mean((Y − μ)²). Layer-wise attribution companion to PO-risk."""
    Y = np.asarray(Y, dtype=float).ravel()
    mu = np.asarray(mu, dtype=float).ravel()
    return float(np.mean((Y - mu) ** 2))


def streaming_po_risk(X_ref, Y_ref, X_new, Y_new, mu_fn=None, clip: float = CLIP, seed: int = 2026) -> float:
    """PO-risk on (ref ∪ new) with T=1 on the new batch."""
    X, Y, T = pack_ref_new(X_ref, Y_ref, X_new, Y_new)
    mu = None if mu_fn is None else np.asarray(mu_fn(X), dtype=float).ravel()
    return po_risk(X, Y, T, mu=mu, clip=clip, seed=seed)


def streaming_po_and_mse(X_ref, Y_ref, X_new, Y_new, mu_fn=None, clip: float = CLIP, seed: int = 2026):
    """PO-risk on (ref ∪ new) and MSE on the new batch, one μ forward.

    Freeze-depth clones share this call so PO-risk and MSE stay paired.
    No extra bootstrap inference.
    """
    X, Y, T = pack_ref_new(X_ref, Y_ref, X_new, Y_new)
    n_ref = len(np.asarray(Y_ref, dtype=float).ravel())
    if mu_fn is None:
        est = TabularPORisk(clip=clip, seed=seed)
        out = est.risk(X, Y, T)
        mu = np.asarray(out["mu"], dtype=float).ravel()
        po = float(out["po_risk"])
    else:
        mu = np.asarray(mu_fn(X), dtype=float).ravel()
        po = po_risk(X, Y, T, mu=mu, clip=clip, seed=seed)
    mse = batch_mse(Y_new, mu[n_ref : n_ref + len(np.asarray(Y_new).ravel())])
    return float(po), float(mse)


def ref_split_baseline(X_ref, Y_ref, *, seed: int = 2026, clip: float = CLIP) -> float:
    """PO-risk on a fake T split of D_ref. The quiet level to read against."""
    X_ref = np.asarray(X_ref, dtype=float)
    Y_ref = np.asarray(Y_ref, dtype=float).ravel()
    n = len(Y_ref)
    half = n // 2
    return streaming_po_risk(X_ref[:half], Y_ref[:half], X_ref[half:], Y_ref[half:], clip=clip, seed=seed)


def moving_average(y, window: int) -> np.ndarray:
    """Causal moving average. No bootstrap, no extra inference."""
    y = np.asarray(y, dtype=float).ravel()
    if y.size == 0:
        return y
    w = max(1, int(window))
    out = np.empty(len(y), dtype=float)
    c = np.cumsum(y)
    for i in range(len(y)):
        lo = i - w + 1
        if lo <= 0:
            out[i] = c[i] / (i + 1)
        else:
            out[i] = (c[i] - c[lo - 1]) / w
    return out


def ma_window(n_new: int) -> int:
    """Cover ~1000 stream rows. n_new=20 → window 50. No extra inference."""
    return max(5, int(round(1000 / max(int(n_new), 1))))


def annotate_moving_average(rows, baseline, n_new=None, ratio: float = DEVIATION_RATIO) -> dict:
    """Attach causal MA to each row. MA below 2× baseline → all-layer backprop."""
    rows = list(rows)
    if n_new is None:
        n_new = int(rows[0]["n_new"]) if rows else 1
    w = ma_window(n_new)
    ys = np.asarray([float(r["po_stream"]) for r in rows], dtype=float) if rows else np.array([])
    ma = moving_average(ys, w) if ys.size else np.array([])
    for r, m in zip(rows, ma):
        r["po_ma"] = float(m)
        r["ma_large"] = bool(large_deviation(float(m), baseline, ratio=ratio))
    thresh = float(ratio) * max(float(baseline), 1e-12)
    stable = bool(ma.size == 0 or float(np.nanmax(ma)) < thresh)
    return {
        "ma_window": int(w),
        "ma_max": float(np.nanmax(ma)) if ma.size else 0.0,
        "ma_stable": stable,
        "all_layer_backprop": stable,
        "frac_ma_large": float(np.mean([r["ma_large"] for r in rows])) if rows else 0.0,
    }


def large_deviation(stream_po: float, baseline_po: float, ratio: float = DEVIATION_RATIO) -> bool:
    """Read the two numbers. Large iff stream is clearly above the ref-split baseline."""
    denom = max(float(baseline_po), 1e-12)
    return float(stream_po) >= float(ratio) * denom


def po_mse_action(po_broken: bool, mse_broken: bool) -> str:
    """PO × MSE contrast on the board.

    PO-risk: did P(Y|X) hop (T=1 new vs T=0 ref)?
    Serving MSE: is the current MLP still within its D_ref error?
    Freeze only if both flags fire — a hop with a quiet MSE is not a
    failed all-layer update. MSE-only bumps stay keep-training.
    """
    if po_broken and mse_broken:
        return ACTION_FREEZE
    if po_broken:
        return ACTION_WATCH
    return ACTION_KEEP


def annotate_po_mse_contrast(rows, po_base, mse_base=None, n_new=None, ratio: float = DEVIATION_RATIO) -> dict:
    """Causal MA of PO-risk and of serving MSE, then the three-way action."""
    po_info = annotate_moving_average(rows, po_base, n_new=n_new, ratio=ratio)
    w = int(po_info["ma_window"])
    has_mse = bool(rows) and all(r.get("mse_stream") is not None for r in rows)
    if has_mse and mse_base is not None:
        ys = np.asarray([float(r["mse_stream"]) for r in rows], dtype=float)
        ma = moving_average(ys, w) if ys.size else np.array([])
        for r, m in zip(rows, ma):
            r["mse_ma"] = float(m)
            r["mse_large"] = bool(large_deviation(float(m), mse_base, ratio=ratio))
            r["po_broken"] = bool(r.get("ma_large"))
            r["mse_broken"] = bool(r["mse_large"])
            r["action"] = po_mse_action(r["po_broken"], r["mse_broken"])
        mse_max = float(np.nanmax(ma)) if ma.size else 0.0
        mse_thresh = float(ratio) * max(float(mse_base), 1e-12)
        mse_stable = bool(ma.size == 0 or mse_max < mse_thresh)
    else:
        for r in rows:
            r["po_broken"] = bool(r.get("ma_large", r.get("large_deviation")))
            r["mse_broken"] = False
            r["action"] = po_mse_action(bool(r["po_broken"]), False)
        mse_max, mse_stable = 0.0, True
    counts = {
        ACTION_KEEP: sum(1 for r in rows if r.get("action") == ACTION_KEEP),
        ACTION_WATCH: sum(1 for r in rows if r.get("action") == ACTION_WATCH),
        ACTION_FREEZE: sum(1 for r in rows if r.get("action") == ACTION_FREEZE),
    }
    if counts[ACTION_FREEZE]:
        board_action = ACTION_FREEZE
    elif counts[ACTION_WATCH]:
        board_action = ACTION_WATCH
    else:
        board_action = ACTION_KEEP
    return {
        **po_info,
        "mse_ma_max": mse_max,
        "mse_stable": mse_stable,
        "n_keep_training": int(counts[ACTION_KEEP]),
        "n_watch": int(counts[ACTION_WATCH]),
        "n_freeze": int(counts[ACTION_FREEZE]),
        "board_action": board_action,
        "all_layer_backprop": bool(board_action == ACTION_KEEP),
    }
