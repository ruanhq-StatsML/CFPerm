"""Amazon MSE prototype: SDC star-prototypes and GPM-Ridge, nothing else.

The only score is online rating MSE (predict the arriving batch, then update).
Plateau RidgeSGD is the current winner on this stream; these two existing
methods are plugged into the same predict-then-update loop so you can read
the shapes and the MSE in one place.

Shapes
------
Amazon (this file):
    X      (N, P)   TF-IDF of title+body, vocab frozen on batch 0
                    P ≤ 128 (Gift-card head often yields P ≈ 109)
    y      (N,)     stars in [1, 5]
    batch  (N,)     category index 0..K-1, K=9
    N      = K * n_per   (n_per=240 in the full run, 80 in --quick)
    mu     (5, P)   one prototype per rounded star
    M      (P, k)   GPM core-gradient bases of the current X batch

MSR-VTT (not scored here):
    s      (n, 2049) = [768 video | 512 audio | 768 text | 1 label]
    X      (n, 2048), W early/late inside each video
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from amazon_continuous_batches import (
    CACHE_DIR,
    CATEGORIES,
    ETA0,
    RidgeSGD,
    batch_mean_cosine,
    concept_intensity_mse,
    featurize_reviews,
    heatmap_hop_c,
    load_amazon_reviews,
    make_amazon_like_stream,
    run_amazon_method,
)
from gpm_fsds import gpm_extend, gpm_gate, project_grad
from instance_discrimination import l2_normalize
from msrvtt_multimodal_attribution import D_AUDIO, D_TEXT, D_VIDEO, P_X, SEED
from prototype_drift import location_shift_field, sdc_compensate

STARS = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
METHODS = ("plateau", "bank", "sdc", "gpm_typed")


def describe_shapes(stream=None):
    """The arrays this prototype actually touches."""
    out = {
        "msrvtt": {
            "s": "(n, 2049) = 768 video + 512 audio + 768 text + 1 label",
            "X": "(n, %d)" % P_X,
            "blocks": {"video": D_VIDEO, "audio": D_AUDIO, "text": D_TEXT},
        },
        "amazon": {
            "X": "(N, P) TF-IDF, P<=128, vocab frozen on batch 0",
            "y": "(N,) stars in [1, 5]",
            "batch": "(N,) category id 0..K-1, K=9",
            "mu": "(5, P) star prototypes",
            "M": "(P, k) GPM bases",
        },
    }
    if stream is not None:
        X = np.asarray(stream.X)
        y = np.asarray(stream.y)
        batch = np.asarray(stream.batch)
        k = int(batch.max()) + 1 if batch.size else 0
        out["amazon_live"] = {
            "X": tuple(int(z) for z in X.shape),
            "y": tuple(int(z) for z in y.shape),
            "batch": tuple(int(z) for z in batch.shape),
            "K": k,
            "P": int(X.shape[1]) if X.ndim == 2 else 0,
            "N": int(X.shape[0]),
            "n_per": {int(t): int((batch == t).sum()) for t in range(k)},
            "y_range": [float(y.min()) if y.size else None, float(y.max()) if y.size else None],
            "categories": list(stream.categories or ()),
        }
    return out


def _softmax(z):
    z = np.asarray(z, dtype=float)
    z = z - z.max(axis=1, keepdims=True)
    e = np.exp(np.clip(z, -50.0, 50.0))
    return e / (e.sum(axis=1, keepdims=True) + 1e-12)


def star_prototypes(X, y):
    """mu[r-1] = mean TF-IDF of reviews with rounded star r."""
    X = np.asarray(X, dtype=float)
    yb = np.clip(np.round(np.asarray(y, dtype=float)), 1, 5).astype(int)
    mu = np.zeros((5, X.shape[1]))
    present = np.zeros(5, dtype=bool)
    for r in range(1, 6):
        sl = yb == r
        if sl.any():
            mu[r - 1] = X[sl].mean(axis=0)
            present[r - 1] = True
        else:
            mu[r - 1] = X.mean(axis=0)
    return mu, present


def sdc_predict(X, mu, tau=0.20):
    """Cosine-weighted star rating. ŷ ∈ [1, 5]."""
    Xn, _ = l2_normalize(X)
    Mn, _ = l2_normalize(mu)
    w = _softmax(Xn @ Mn.T / float(tau))
    pred = w @ STARS
    return np.clip(pred, 1.0, 5.0)


def run_bank_mse(stream, tau=0.20, cap=800):
    """Wu-style instance bank as a rating regressor.

    Store (x, y) from seen batches. On the arriving batch, cosine-weight the
    stored ratings (Nadaraya–Watson). Then enqueue the new rows. This is the
    memory bank used for MSE, not for NCE.
    """
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    n_batches = int(batch.max()) + 1
    i0 = np.flatnonzero(batch == 0)
    Z, yz = X[i0].copy(), y[i0].copy()
    history = []
    for t in range(1, n_batches):
        ic = np.flatnonzero(batch == t)
        Xt, yt = X[ic], y[ic]
        pred = _nw_predict(Xt, Z, yz, tau=tau)
        online = _mse(pred, yt)
        Z = np.vstack([Z, Xt])
        yz = np.concatenate([yz, yt])
        if len(yz) > int(cap):
            Z, yz = Z[-int(cap) :], yz[-int(cap) :]
        history.append(
            {
                "round": int(t),
                "phase": "adapt",
                "method": "bank",
                "online_mse": online,
                "bank": int(len(yz)),
                "n": int(ic.size),
            }
        )
    online = np.array([h["online_mse"] for h in history], dtype=float)
    return {
        "method": "bank",
        "n_batches": n_batches,
        "online_mse": float(online.mean()) if online.size else float("nan"),
        "cum_mse": float(online.sum()) if online.size else float("nan"),
        "last_online_mse": float(online[-1]) if online.size else float("nan"),
        "bwt": float("nan"),
        "online_path": online.tolist(),
        "history": history,
        "meta": dict(stream.meta),
    }


def _nw_predict(X, Z, yz, tau=0.20):
    Xn, _ = l2_normalize(X)
    Zn, _ = l2_normalize(Z)
    w = _softmax(Xn @ Zn.T / float(tau))
    return np.clip(w @ np.asarray(yz, dtype=float), 1.0, 5.0)


def _mse(pred, y):
    return float(np.mean((np.asarray(pred) - np.asarray(y, dtype=float)) ** 2))


def run_sdc_mse(stream, tau=0.20):
    """Yu SDC on star prototypes. Predict batch t, then replace μ with t's means."""
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    n_batches = int(batch.max()) + 1
    i0 = np.flatnonzero(batch == 0)
    mu, _ = star_prototypes(X[i0], y[i0])
    X_prev = X[i0]
    history = []
    for t in range(1, n_batches):
        ic = np.flatnonzero(batch == t)
        Xt, yt = X[ic], y[ic]
        field = location_shift_field(X_prev, Xt)
        mu_use = sdc_compensate(mu, Xt - field, Xt)
        pred = sdc_predict(Xt, mu_use, tau=tau)
        online = _mse(pred, yt)
        mu, _ = star_prototypes(Xt, yt)
        X_prev = Xt
        history.append(
            {
                "round": int(t),
                "phase": "adapt",
                "method": "sdc",
                "online_mse": online,
                "n": int(ic.size),
            }
        )
    online = np.array([h["online_mse"] for h in history], dtype=float)
    return {
        "method": "sdc",
        "n_batches": n_batches,
        "online_mse": float(online.mean()) if online.size else float("nan"),
        "cum_mse": float(online.sum()) if online.size else float("nan"),
        "last_online_mse": float(online[-1]) if online.size else float("nan"),
        "bwt": float("nan"),
        "online_path": online.tolist(),
        "history": history,
        "meta": dict(stream.meta),
    }


def _ridge_gpm_step(probe, X, y, eta, M, gate, ridge=1e-3):
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    n = max(len(y), 1)
    err = probe.predict(X) - y
    dw = (X.T @ err) / n + float(ridge) * probe.w
    if gate:
        dw = project_grad(dw, M)
    probe.w = probe.w - float(eta) * dw
    probe.b = probe.b - float(eta) * float(err.mean())


def run_gpm_mse(
    stream,
    typed=True,
    eta0=ETA0,
    steps_per_batch=8,
    warmup_steps=6,
    seed=SEED,
    ridge=1e-3,
    batch_size=32,
    thresh=0.90,
    max_k=12,
):
    """Plateau RidgeSGD with GPM projection on dw. typed uses (ĉ, δ̂) as the gate."""
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    n_batches = int(batch.max()) + 1
    probe = RidgeSGD(X.shape[1], seed=seed)
    rng = np.random.default_rng(seed)
    i0 = np.flatnonzero(batch == 0)
    probe.b = float(y[i0].mean())
    eta = float(eta0)
    chunk = min(int(batch_size), max(4, i0.size))
    for _ in range(int(max(1, warmup_steps))):
        sl = rng.choice(i0.size, size=chunk, replace=False)
        probe.step(X[i0][sl], y[i0][sl], eta, ridge=ridge)
    M = gpm_extend(None, X[i0], thresh=thresh, max_k=max_k)
    R = batch_mean_cosine(X, batch, n_batches=n_batches)
    plateau_eta = float(eta0)
    best = np.inf
    stall = 0
    history = []
    for t in range(1, n_batches):
        ip, ic = np.flatnonzero(batch == t - 1), np.flatnonzero(batch == t)
        c_hat = heatmap_hop_c(R, t)
        d_hat, _ = concept_intensity_mse(X[ip], y[ip], X[ic], y[ic])
        gate = bool(gpm_gate(c_hat, d_hat)) if typed else True
        online = probe.mse(X[ic], y[ic])
        chunk_t = min(int(batch_size), max(4, ic.size))
        for _ in range(max(1, int(steps_per_batch))):
            sl = rng.choice(ic.size, size=chunk_t, replace=False)
            _ridge_gpm_step(probe, X[ic][sl], y[ic][sl], plateau_eta, M, gate, ridge=ridge)
        history.append(
            {
                "round": int(t),
                "phase": "adapt",
                "method": "gpm_typed" if typed else "gpm",
                "online_mse": float(online),
                "c": float(c_hat),
                "delta": float(d_hat),
                "gate": gate,
                "rank": int(0 if M is None else M.shape[1]),
                "n": int(ic.size),
            }
        )
        if online < best - 1e-4:
            best = online
            stall = 0
        else:
            stall += 1
            if stall >= 2:
                plateau_eta *= 0.5
                stall = 0
        M = gpm_extend(M, X[ic], thresh=thresh, max_k=max_k)
    online = np.array([h["online_mse"] for h in history], dtype=float)
    name = "gpm_typed" if typed else "gpm"
    return {
        "method": name,
        "n_batches": n_batches,
        "online_mse": float(online.mean()) if online.size else float("nan"),
        "cum_mse": float(online.sum()) if online.size else float("nan"),
        "last_online_mse": float(online[-1]) if online.size else float("nan"),
        "bwt": probe.mse(X[i0], y[i0]),
        "mean_gate": float(np.mean([h["gate"] for h in history])) if history else 0.0,
        "online_path": online.tolist(),
        "history": history,
        "meta": dict(stream.meta),
    }


def run_mse_method(stream, method, seed=SEED, **kw):
    if method == "plateau":
        rec, _ = run_amazon_method(stream, method="plateau", seed=seed, **kw)
        return rec
    if method == "sdc":
        return run_sdc_mse(stream)
    if method == "bank":
        return run_bank_mse(stream)
    if method == "gpm":
        return run_gpm_mse(stream, typed=False, seed=seed, **kw)
    if method == "gpm_typed":
        return run_gpm_mse(stream, typed=True, seed=seed, **kw)
    raise ValueError("unknown method %s" % method)


def load_amazon_stream(n_per=240, n_batches=9, seed=SEED, max_features=128):
    raw = load_amazon_reviews(
        categories=CATEGORIES[: int(n_batches)], n_per=n_per, seed=seed, cache_dir=CACHE_DIR
    )
    return featurize_reviews(raw, max_features=max_features)


def run_mse_suite(seeds=None, n_per=240, n_batches=9, source="amazon", steps_per_batch=8):
    seeds = list(seeds if seeds is not None else range(SEED, SEED + 4))
    rows = []
    traces = {m: [] for m in METHODS}
    shapes = None
    for s in seeds:
        if source == "synthetic":
            stream = make_amazon_like_stream(n_batches=n_batches, n_per=n_per, seed=int(s), cov=0.35)
        else:
            stream = load_amazon_stream(n_per=n_per, n_batches=n_batches, seed=int(s))
        if shapes is None:
            shapes = describe_shapes(stream)
        for method in METHODS:
            rec = run_mse_method(stream, method, seed=int(s), steps_per_batch=steps_per_batch)
            traces[method].append(rec)
            rows.append(
                {
                    "seed": int(s),
                    "method": method,
                    "online_mse": rec["online_mse"],
                    "cum_mse": rec["cum_mse"],
                    "last_online_mse": rec["last_online_mse"],
                }
            )
    table = {}
    for method in METHODS:
        v = np.array([r["online_mse"] for r in rows if r["method"] == method], dtype=float)
        table[method] = {
            "online_mse": {"mean": float(v.mean()), "sd": float(v.std(ddof=1) if len(v) > 1 else 0.0)},
            "n": int(len(v)),
        }
    return {
        "table": table,
        "rows": rows,
        "traces": traces,
        "shapes": shapes,
        "seeds": [int(s) for s in seeds],
        "source": source,
    }


def plot_mse_suite(suite, path):
    import matplotlib.pyplot as plt
    from msrvtt_attribution_plots import GRID, INK, MUTED, _save, _style

    _style()
    colors = {"plateau": "#C4892A", "bank": "#2F6B4F", "sdc": "#2C4A6E", "gpm_typed": "#C45C26"}
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    traces = suite["traces"]
    t = None
    for method in METHODS:
        recs = traces[method]
        paths = np.array([r["online_path"] for r in recs], dtype=float)
        t = np.arange(1, paths.shape[1] + 1)
        mu, sd = paths.mean(0), paths.std(0)
        ax.plot(t, mu, color=colors[method], lw=2.1, label=method)
        ax.fill_between(t, mu - sd, mu + sd, color=colors[method], alpha=0.12, lw=0)
    ax.set_xlabel("batch")
    ax.set_ylabel("online MSE")
    ax.set_title("Amazon rating MSE prototype", loc="left", fontsize=13, fontweight="bold")
    ax.grid(True, color=GRID)
    ax.legend(frameon=False, fontsize=8.5)
    fig.text(
        0.04,
        0.01,
        "Predict the arriving category, then update. Plateau RidgeSGD is the baseline. "
        "bank = instance memory (cosine-weighted ratings); sdc = five star means; gpm_typed = GPM on dw.",
        fontsize=8.0,
        color=MUTED,
    )
    return _save(fig, Path(path))
