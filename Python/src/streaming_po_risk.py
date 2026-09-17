"""Streaming PO-risk on a tabular table.

New batch is T=1, reference is T=0. Outcome model μ(Y|X) and
propensity e(T|X) are Random Forests, maintained separately — not
the serving MLP, not a single lstsq.

φ = (Y − μ)(T − e), τ̂(X) ≈ φ, risk = mean(τ̂²).

RF PO-risk is not expected to collapse suddenly. Serving MSE of the
MLP is the series more likely to break first.

MMD口径 is fixed: RBF MMD²(X_new, X_ref). Same T=0 reference as
PO-risk. Not the mean pairwise MMD against all previous batches, and
not MMD on a layer representation vs the last batch. Those answer
different questions.

No online-bootstrap.
"""
from __future__ import annotations

import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor

CLIP = 1e-3
REF_N = 10_000
MIN_STREAM_N = 5_000
# Stream vs ref-split baseline. Below this, all layers stay trainable.
DEVIATION_RATIO = 2.0
ACTION_KEEP = "keep_training"
ACTION_WATCH = "watch"
ACTION_FREEZE = "freeze"
ACTION_XSHIFT = "x_shift"
ACTION_TRICKY = "tricky"
MMD_MAX_N = 512


def _sqdist(A, B):
    aa = np.sum(np.asarray(A, dtype=float) ** 2, axis=1, keepdims=True)
    bb = np.sum(np.asarray(B, dtype=float) ** 2, axis=1, keepdims=True).T
    return np.maximum(aa + bb - 2.0 * (A @ B.T), 0.0)


def _subsample_rows(X, n, rng):
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    if len(X) <= n:
        return X
    return X[rng.choice(len(X), size=int(n), replace=False)]


def rbf_bandwidth(X, max_n: int = MMD_MAX_N, seed: int = 2026) -> float:
    """Median heuristic on a subsample of X. Fixed for the stream so MMD is comparable."""
    rng = np.random.default_rng(seed)
    Xs = _subsample_rows(X, max_n, rng)
    D = _sqdist(Xs, Xs)
    tri = D[np.triu_indices_from(D, k=1)]
    med = float(np.median(tri)) if tri.size else 0.0
    if med <= 0:
        return 1.0
    return float(np.sqrt(med / 2.0))


def rbf_mmd2(X, Y, sigma: float, max_n: int = MMD_MAX_N, seed: int = 2026) -> float:
    """Unbiased RBF MMD²(X, Y). Call as MMD²(X_new, X_ref) only."""
    rng = np.random.default_rng(seed)
    Xs = _subsample_rows(X, max_n, rng)
    Ys = _subsample_rows(Y, max_n, rng)
    sig = max(float(sigma), 1e-8)
    gamma = 1.0 / (2.0 * sig * sig)
    Kxx = np.exp(-gamma * _sqdist(Xs, Xs))
    Kyy = np.exp(-gamma * _sqdist(Ys, Ys))
    Kxy = np.exp(-gamma * _sqdist(Xs, Ys))
    n, m = len(Xs), len(Ys)
    xx = 0.0 if n < 2 else (Kxx.sum() - np.trace(Kxx)) / (n * (n - 1))
    yy = 0.0 if m < 2 else (Kyy.sum() - np.trace(Kyy)) / (m * (m - 1))
    return float(xx + yy - 2.0 * float(Kxy.mean()))


def mmd_vs_reference(X_ref, X_new, sigma: float, max_n: int = MMD_MAX_N, seed: int = 2026) -> float:
    """Covariate-shift readout: MMD²(X_new, X_ref).

    Same reference batch as PO-risk's T=0. Not pairwise-vs-history, not
    last-batch layer representations.
    """
    return rbf_mmd2(X_new, X_ref, sigma=sigma, max_n=max_n, seed=seed)


def ref_split_mmd(X_ref, *, sigma: float, seed: int = 2026, max_n: int = MMD_MAX_N) -> float:
    """MMD of a fake split of D_ref. Quiet level for P(X) vs itself."""
    X_ref = np.asarray(X_ref, dtype=float)
    n = len(X_ref)
    half = n // 2
    return rbf_mmd2(X_ref[:half], X_ref[half:], sigma=sigma, max_n=max_n, seed=seed)


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
    """Own RF outcome model + own RF propensity. Not the serving MLP."""

    def __init__(self, clip: float = CLIP, seed: int = 2026):
        self.clip = float(clip)
        self.seed = int(seed)
        self.outcome = None
        self.propensity = None

    def _rf_cls(self):
        return RandomForestClassifier(
            n_estimators=60,
            max_depth=8,
            min_samples_leaf=10,
            n_jobs=1,
            random_state=self.seed,
        )

    def _rf_reg(self):
        return RandomForestRegressor(
            n_estimators=60,
            max_depth=8,
            min_samples_leaf=10,
            n_jobs=1,
            random_state=self.seed,
        )

    def fit_nuisance(self, X, Y, T, X_outcome=None):
        """Fit μ(Y|X_outcome) and e(T|X) with RF. X_outcome defaults to X."""
        X = np.asarray(X, dtype=float)
        if X.ndim == 1:
            X = X.reshape(-1, 1)
        Y = np.asarray(Y, dtype=float).ravel()
        T = np.asarray(T, dtype=float).ravel().astype(int)
        Xo = X if X_outcome is None else np.asarray(X_outcome, dtype=float)
        if Xo.ndim == 1:
            Xo = Xo.reshape(-1, 1)
        self.propensity = self._rf_cls()
        self.propensity.fit(X, T)
        e = np.clip(self.propensity.predict_proba(X)[:, 1], self.clip, 1.0 - self.clip)
        if _is_binary(Y):
            self.outcome = self._rf_cls()
            self.outcome.fit(Xo, Y.astype(int))
            mu = self.outcome.predict_proba(Xo)[:, 1]
        else:
            self.outcome = self._rf_reg()
            self.outcome.fit(Xo, Y)
            mu = self.outcome.predict(Xo)
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


def po_mse_action(po_broken: bool, mse_broken: bool, mmd_broken: bool | None = None) -> str:
    """PO × MSE contrast. Freeze-from-layer only when both PO and MSE break.

    RF PO-risk is not expected to collapse first. Serving MSE breaking
    while PO holds is the main cell — then read MMD²(X_new, X_ref).
    PO broken and MSE holds → watch, do not freeze yet.
    """
    if po_broken and mse_broken:
        return ACTION_FREEZE
    if po_broken:
        return ACTION_WATCH
    if mse_broken:
        if mmd_broken:
            return ACTION_XSHIFT
        return ACTION_TRICKY
    return ACTION_KEEP


def annotate_po_mse_contrast(rows, po_base, mse_base=None, mmd_base=None, n_new=None, ratio: float = DEVIATION_RATIO) -> dict:
    """Causal MA of PO-risk, serving MSE, and MMD of X, then the action."""
    po_info = annotate_moving_average(rows, po_base, n_new=n_new, ratio=ratio)
    w = int(po_info["ma_window"])

    def _ma_flag(key, baseline):
        has = bool(rows) and all(r.get(key) is not None for r in rows) and baseline is not None
        if not has:
            return 0.0, True, None
        ys = np.asarray([float(r[key]) for r in rows], dtype=float)
        ma = moving_average(ys, w) if ys.size else np.array([])
        ma_key = "mse_ma" if key == "mse_stream" else "mmd_ma"
        large_key = "mse_large" if key == "mse_stream" else "mmd_large"
        for r, m in zip(rows, ma):
            r[ma_key] = float(m)
            r[large_key] = bool(large_deviation(float(m), baseline, ratio=ratio))
        mx = float(np.nanmax(ma)) if ma.size else 0.0
        thresh = float(ratio) * max(float(baseline), 1e-12)
        return mx, bool(ma.size == 0 or mx < thresh), True

    def _mmd_key(rows):
        if rows and all(r.get("mmd_vs_ref") is not None for r in rows):
            return "mmd_vs_ref"
        if rows and all(r.get("mmd_stream") is not None for r in rows):
            return "mmd_stream"
        return None

    mse_max, mse_stable, has_mse = _ma_flag("mse_stream", mse_base)
    mmd_key = _mmd_key(rows)
    mmd_max, mmd_stable, has_mmd = _ma_flag(mmd_key, mmd_base) if mmd_key else (0.0, True, None)
    for r in rows:
        r["po_broken"] = bool(r.get("ma_large", r.get("large_deviation")))
        r["mse_broken"] = bool(r.get("mse_large")) if has_mse else False
        r["mmd_broken"] = bool(r.get("mmd_large")) if has_mmd else False
        r["action"] = po_mse_action(r["po_broken"], r["mse_broken"], r["mmd_broken"])
    counts = {
        ACTION_KEEP: sum(1 for r in rows if r.get("action") == ACTION_KEEP),
        ACTION_WATCH: sum(1 for r in rows if r.get("action") == ACTION_WATCH),
        ACTION_XSHIFT: sum(1 for r in rows if r.get("action") == ACTION_XSHIFT),
        ACTION_TRICKY: sum(1 for r in rows if r.get("action") == ACTION_TRICKY),
        ACTION_FREEZE: sum(1 for r in rows if r.get("action") == ACTION_FREEZE),
    }
    if counts[ACTION_FREEZE]:
        board_action = ACTION_FREEZE
    elif counts[ACTION_XSHIFT]:
        board_action = ACTION_XSHIFT
    elif counts[ACTION_TRICKY]:
        board_action = ACTION_TRICKY
    elif counts[ACTION_WATCH]:
        board_action = ACTION_WATCH
    else:
        board_action = ACTION_KEEP
    return {
        **po_info,
        "mse_ma_max": mse_max,
        "mse_stable": mse_stable,
        "mmd_ma_max": mmd_max,
        "mmd_stable": mmd_stable,
        "n_keep_training": int(counts[ACTION_KEEP]),
        "n_watch": int(counts[ACTION_WATCH]),
        "n_x_shift": int(counts[ACTION_XSHIFT]),
        "n_tricky": int(counts[ACTION_TRICKY]),
        "n_freeze": int(counts[ACTION_FREEZE]),
        "board_action": board_action,
        "all_layer_backprop": bool(board_action == ACTION_KEEP),
    }
