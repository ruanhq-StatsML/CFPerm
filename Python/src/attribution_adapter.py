"""Attribution as next-batch training adapter (Amazon rating MSE).

FSDS / TSS already estimates the hop (ĉ from the batch-relationship heatmap,
δ̂ from excess risk after mean-align). This module turns that pair into the
*training weights* for the next batch — not a new loss.

Existing pieces only:
  - heatmap hop cosine → sample weights for closed-form Ridge (domain proximity IW)
  - δ̂ / ĉ → gate between instance bank (Wu NW) and hop-weighted Ridge (TSS signs)
  - rolling EWMA of per-predictor MSE → online stacking weight
  - river.PARegressor as an off-the-shelf online baseline

Locked signs (unchanged):
  ĉ large  → prefer bank / shrink SGD step (do not chase P(X|W))
  δ̂ large  → prefer Ridge refit / raise step (relearn P(Y|X))
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from sklearn.linear_model import Ridge

from amazon_continuous_batches import (
    CACHE_DIR,
    CATEGORIES,
    ETA0,
    batch_mean_cosine,
    concept_intensity_mse,
    featurize_reviews,
    heatmap_hop_c,
    load_amazon_reviews,
    make_amazon_like_stream,
    run_amazon_method,
)
from amazon_mse_prototype import _mse, _nw_predict, describe_shapes, run_bank_mse
from msrvtt_multimodal_attribution import SEED

METHODS = (
    "plateau",
    "bank",
    "ridge_past",
    "hop_ridge",
    "attr_adapter",
    "river_pa",
)


def load_amazon_stream(n_per=240, n_batches=9, seed=SEED, max_features=128):
    raw = load_amazon_reviews(
        categories=CATEGORIES[: int(n_batches)],
        n_per=n_per,
        seed=seed,
        cache_dir=CACHE_DIR,
    )
    return featurize_reviews(raw, max_features=max_features)


def rolling_hop_stats(stream):
    """Online-causal (ĉ_t, δ̂_t) on each hop t-1 → t, plus batch means."""
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    k = int(batch.max()) + 1
    R = batch_mean_cosine(X, batch, n_batches=k)
    mus = np.stack([X[batch == t].mean(axis=0) for t in range(k)])
    c = np.zeros(k)
    d = np.zeros(k)
    for t in range(1, k):
        ip, ic = np.flatnonzero(batch == t - 1), np.flatnonzero(batch == t)
        c[t] = heatmap_hop_c(R, t)
        d[t], _ = concept_intensity_mse(X[ip], y[ip], X[ic], y[ic])
    return {"R": R, "mus": mus, "c": c, "delta": d, "n_batches": k}


def hop_sample_weights(mus, batch, t, gamma=4.0):
    """Importance weights for past rows: exp(γ (cos(μ_s, μ_{t-1}) − 1)).

    Causal: the arriving batch mean is unknown, so the reference is μ_{t-1}
    (last observed category). Same cosine board as the FSDS heatmap.
    """
    batch = np.asarray(batch, dtype=int)
    w = np.zeros(batch.shape[0], dtype=float)
    ref = mus[t - 1]
    ref = ref / (np.linalg.norm(ref) + 1e-12)
    for s in range(t):
        mu = mus[s] / (np.linalg.norm(mus[s]) + 1e-12)
        cos = float(np.dot(mu, ref))
        w[batch == s] = float(np.exp(float(gamma) * (cos - 1.0)))
    return np.maximum(w, 1e-3)


def _ridge_predict(Xtr, ytr, Xte, sample_weight=None, alpha=3.0):
    clf = Ridge(alpha=float(alpha))
    if sample_weight is None:
        clf.fit(Xtr, ytr)
    else:
        clf.fit(Xtr, ytr, sample_weight=sample_weight)
    return np.clip(clf.predict(Xte), 1.0, 5.0)


def run_ridge_past(stream, alpha=3.0):
    """Closed-form Ridge on all past batches (uniform weights)."""
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    k = int(batch.max()) + 1
    history = []
    for t in range(1, k):
        tr, ic = batch < t, batch == t
        pred = _ridge_predict(X[tr], y[tr], X[ic], alpha=alpha)
        history.append(
            {
                "round": int(t),
                "online_mse": _mse(pred, y[ic]),
                "n": int(ic.sum()),
            }
        )
    online = np.array([h["online_mse"] for h in history], dtype=float)
    return _summary("ridge_past", online, history, stream)


def run_hop_ridge(stream, alpha=3.0, gamma=4.0):
    """Heatmap-hop importance-weighted Ridge — attribution as sample weights."""
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    stats = rolling_hop_stats(stream)
    history = []
    for t in range(1, stats["n_batches"]):
        tr, ic = batch < t, batch == t
        w = hop_sample_weights(stats["mus"], batch, t, gamma=gamma)
        pred = _ridge_predict(X[tr], y[tr], X[ic], sample_weight=w[tr], alpha=alpha)
        history.append(
            {
                "round": int(t),
                "online_mse": _mse(pred, y[ic]),
                "c": float(stats["c"][t]),
                "delta": float(stats["delta"][t]),
                "w_mean": float(w[tr].mean()),
                "w_max": float(w[tr].max()),
                "n": int(ic.sum()),
            }
        )
    online = np.array([h["online_mse"] for h in history], dtype=float)
    return _summary("hop_ridge", online, history, stream, extras={"gamma": gamma, "alpha": alpha})


def typed_bank_weight(c, delta, ewma_bank=None, ewma_ridge=None):
    """TSS signs → mixture weight on the instance bank.

    ĉ large / δ̂ quiet → bank up. δ̂ large → bank down (Ridge refit).
    Optional EWMA stack blends with inverse recent error.
    """
    c = float(c)
    delta = float(delta)
    w_cd = (1.0 + 2.0 * c) / (1.0 + 2.0 * c + 4.0 * delta)
    if ewma_bank is None or ewma_ridge is None:
        return float(np.clip(w_cd, 0.05, 0.70))
    inv_b = 1.0 / (float(ewma_bank) + 1e-3)
    inv_r = 1.0 / (float(ewma_ridge) + 1e-3)
    w_err = inv_b / (inv_b + inv_r)
    return float(np.clip(0.5 * w_cd + 0.5 * w_err, 0.05, 0.70))


def run_attr_adapter(stream, alpha=3.0, gamma=4.0, tau=0.20, cap=800):
    """Next-batch adapter: typed mix of Wu bank and hop-weighted Ridge.

    Predict the arriving batch with
        ŷ = w · bank_NW + (1−w) · Ridge_hop(past),
    w from (ĉ, δ̂) and rolling predictor EWMA, then enqueue the batch.
    """
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    stats = rolling_hop_stats(stream)
    i0 = np.flatnonzero(batch == 0)
    Z, yz = X[i0].copy(), y[i0].copy()
    ewma_bank = ewma_ridge = None
    history = []
    for t in range(1, stats["n_batches"]):
        ic = np.flatnonzero(batch == t)
        tr = batch < t
        Xt, yt = X[ic], y[ic]
        c, d = float(stats["c"][t]), float(stats["delta"][t])
        w = typed_bank_weight(c, d, ewma_bank, ewma_ridge)
        pb = _nw_predict(Xt, Z, yz, tau=tau)
        sw = hop_sample_weights(stats["mus"], batch, t, gamma=gamma)
        pr = _ridge_predict(X[tr], y[tr], Xt, sample_weight=sw[tr], alpha=alpha)
        pred = w * pb + (1.0 - w) * pr
        mb, mr = _mse(pb, yt), _mse(pr, yt)
        ewma_bank = mb if ewma_bank is None else 0.6 * ewma_bank + 0.4 * mb
        ewma_ridge = mr if ewma_ridge is None else 0.6 * ewma_ridge + 0.4 * mr
        history.append(
            {
                "round": int(t),
                "online_mse": _mse(pred, yt),
                "mse_bank": mb,
                "mse_ridge": mr,
                "w_bank": w,
                "c": c,
                "delta": d,
                "n": int(ic.size),
            }
        )
        Z = np.vstack([Z, Xt])
        yz = np.concatenate([yz, yt])
        if len(yz) > int(cap):
            Z, yz = Z[-int(cap) :], yz[-int(cap) :]
    online = np.array([h["online_mse"] for h in history], dtype=float)
    return _summary(
        "attr_adapter",
        online,
        history,
        stream,
        extras={
            "alpha": alpha,
            "gamma": gamma,
            "mean_w_bank": float(np.mean([h["w_bank"] for h in history])),
        },
    )


def run_river_pa(stream, C=0.01):
    """river PassiveAwareRegressor — package baseline on the same stream."""
    from river import compose, linear_model, preprocessing

    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    k = int(batch.max()) + 1
    p = X.shape[1]
    model = compose.Pipeline(
        preprocessing.StandardScaler(),
        linear_model.PARegressor(C=float(C), mode=2),
    )

    def row(i):
        return {j: float(X[i, j]) for j in range(p)}

    for i in np.flatnonzero(batch == 0):
        model.learn_one(row(i), float(y[i]))
    history = []
    for t in range(1, k):
        ic = np.flatnonzero(batch == t)
        preds = []
        for i in ic:
            yhat = model.predict_one(row(i))
            preds.append(float(np.clip(3.0 if yhat is None else yhat, 1.0, 5.0)))
        history.append({"round": int(t), "online_mse": _mse(preds, y[ic]), "n": int(ic.size)})
        for i in ic:
            model.learn_one(row(i), float(y[i]))
    online = np.array([h["online_mse"] for h in history], dtype=float)
    return _summary("river_pa", online, history, stream, extras={"C": C, "package": "river"})


def _summary(method, online, history, stream, extras=None):
    out = {
        "method": method,
        "n_batches": int(np.asarray(stream.batch).max()) + 1,
        "online_mse": float(online.mean()) if online.size else float("nan"),
        "cum_mse": float(online.sum()) if online.size else float("nan"),
        "last_online_mse": float(online[-1]) if online.size else float("nan"),
        "online_path": online.tolist(),
        "history": history,
        "meta": dict(getattr(stream, "meta", {}) or {}),
    }
    if extras:
        out.update(extras)
    return out


def run_adapter_method(stream, method, seed=SEED, **kw):
    if method == "plateau":
        rec, _ = run_amazon_method(stream, method="plateau", seed=seed, **kw)
        return {
            "method": "plateau",
            "n_batches": rec["n_batches"],
            "online_mse": rec["online_mse"],
            "cum_mse": rec["cum_mse"],
            "last_online_mse": rec["last_online_mse"],
            "online_path": rec.get("online_path", []),
            "history": rec.get("history", []),
            "meta": dict(getattr(stream, "meta", {}) or {}),
        }
    if method == "bank":
        return run_bank_mse(stream)
    if method == "ridge_past":
        return run_ridge_past(stream)
    if method == "hop_ridge":
        return run_hop_ridge(stream)
    if method == "attr_adapter":
        return run_attr_adapter(stream)
    if method == "river_pa":
        return run_river_pa(stream)
    raise ValueError("unknown method %s" % method)


def run_adapter_suite(
    seeds=None,
    n_per=240,
    n_batches=9,
    source="amazon",
    methods=None,
    steps_per_batch=8,
):
    seeds = list(seeds if seeds is not None else range(SEED, SEED + 4))
    methods = list(methods or METHODS)
    rows = []
    traces = {m: [] for m in methods}
    shapes = None
    hop_logs = []
    for s in seeds:
        if source == "synthetic":
            stream = make_amazon_like_stream(
                n_batches=n_batches, n_per=n_per, seed=int(s), cov=0.35
            )
        else:
            stream = load_amazon_stream(n_per=n_per, n_batches=n_batches, seed=int(s))
        if shapes is None:
            shapes = describe_shapes(stream)
        hops = rolling_hop_stats(stream)
        hop_logs.append(
            {
                "seed": int(s),
                "c": hops["c"].tolist(),
                "delta": hops["delta"].tolist(),
            }
        )
        for method in methods:
            kw = {}
            if method == "plateau":
                kw["steps_per_batch"] = steps_per_batch
            rec = run_adapter_method(stream, method, seed=int(s), **kw)
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
    for method in methods:
        v = np.array([r["online_mse"] for r in rows if r["method"] == method], dtype=float)
        table[method] = {
            "online_mse": {
                "mean": float(v.mean()),
                "sd": float(v.std(ddof=1) if len(v) > 1 else 0.0),
            },
            "n": int(len(v)),
        }
    # paired lift vs bank / plateau
    lifts = {}
    for baseline in ("bank", "plateau"):
        if baseline not in methods:
            continue
        lifts[baseline] = {}
        for method in methods:
            if method == baseline:
                continue
            diffs = []
            for s in seeds:
                a = next(r["online_mse"] for r in rows if r["seed"] == s and r["method"] == method)
                b = next(r["online_mse"] for r in rows if r["seed"] == s and r["method"] == baseline)
                diffs.append(b - a)  # positive => method better (lower MSE)
            d = np.asarray(diffs, dtype=float)
            lifts[baseline][method] = {
                "mean_mse_drop": float(d.mean()),
                "sd": float(d.std(ddof=1) if len(d) > 1 else 0.0),
            }
    return {
        "table": table,
        "rows": rows,
        "traces": traces,
        "shapes": shapes,
        "lifts": lifts,
        "hop_logs": hop_logs,
        "seeds": [int(s) for s in seeds],
        "source": source,
        "methods": methods,
    }


def plot_adapter_suite(suite, path):
    import matplotlib.pyplot as plt
    from msrvtt_attribution_plots import GRID, INK, MUTED, _save, _style

    _style()
    colors = {
        "plateau": "#C4892A",
        "bank": "#2F6B4F",
        "ridge_past": "#6B7C8A",
        "hop_ridge": "#C45C26",
        "attr_adapter": "#2C4A6E",
        "river_pa": "#9AA3AE",
    }
    labels = {
        "plateau": "plateau RidgeSGD",
        "bank": "instance bank",
        "ridge_past": "Ridge on past",
        "hop_ridge": "hop-weighted Ridge",
        "attr_adapter": "attr adapter (bank⊕hop)",
        "river_pa": "river PA",
    }
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2), gridspec_kw={"width_ratios": [1.15, 1.0]})
    ax = axes[0]
    methods = [m for m in suite["methods"] if m in suite["traces"]]
    for method in methods:
        paths = np.array([r["online_path"] for r in suite["traces"][method]], dtype=float)
        if paths.size == 0 or paths.shape[1] == 0:
            continue
        t = np.arange(1, paths.shape[1] + 1)
        mu, sd = paths.mean(0), paths.std(0)
        ax.plot(t, mu, color=colors.get(method, INK), lw=2.1, label=labels.get(method, method))
        ax.fill_between(t, mu - sd, mu + sd, color=colors.get(method, INK), alpha=0.12, lw=0)
    ax.set_xlabel("batch")
    ax.set_ylabel("online MSE")
    ax.set_title("Attribution adapter on Amazon ratings", loc="left", fontsize=12.5, fontweight="bold")
    ax.grid(True, color=GRID)
    ax.legend(frameon=False, fontsize=7.8)

    ax = axes[1]
    means = [suite["table"][m]["online_mse"]["mean"] for m in methods]
    sds = [suite["table"][m]["online_mse"]["sd"] for m in methods]
    ypos = np.arange(len(methods))[::-1]
    ax.barh(
        ypos,
        means,
        xerr=sds,
        color=[colors.get(m, INK) for m in methods],
        height=0.62,
        error_kw={"ecolor": "#555", "lw": 1.0, "capsize": 2},
    )
    ax.set_yticks(ypos)
    ax.set_yticklabels([labels.get(m, m) for m in methods], fontsize=8.2)
    ax.set_xlabel("mean online MSE")
    ax.set_title("Lower is better", loc="left", fontsize=12.5, fontweight="bold")
    ax.grid(True, color=GRID, axis="x")
    fig.text(
        0.04,
        0.01,
        "hop-weighted Ridge uses heatmap cos(μ_s, μ_{t−1}) as sample weights. "
        "attr adapter mixes Wu bank and hop Ridge with TSS signs on (ĉ, δ̂).",
        fontsize=8.0,
        color=MUTED,
    )
    return _save(fig, Path(path))


def write_adapter_tex(suite, path):
    methods = suite["methods"]
    labels = {
        "plateau": "plateau RidgeSGD",
        "bank": "instance bank (Wu NW)",
        "ridge_past": "Ridge on all past",
        "hop_ridge": "hop-weighted Ridge",
        "attr_adapter": "attr.\\ adapter (bank$\\oplus$hop)",
        "river_pa": "river PARegressor",
    }
    lines = [
        r"% Attribution adapter on Amazon. Auto-generated.",
        r"\begin{table}[ht]\centering",
        r"\caption{Multimodal / batch attribution as a next-batch training adapter on Amazon Reviews 2023.",
        r"Online rating MSE (predict the arriving category, then update). Hop-weighted Ridge reuses the",
        r"batch-relationship heatmap as sample weights $w_s=\exp(\gamma(\cos(\bar x_s,\bar x_{t-1})-1))$.",
        r"The typed adapter mixes the Wu instance bank with hop Ridge using the locked TSS signs on $(\hat c,\hat\delta)$.}",
        r"\label{tab:attr-adapter-amazon}",
        r"\small",
        r"\begin{tabular}{@{}lcc@{}}\toprule",
        r"Method & online MSE & vs bank \\",
        r"\midrule",
    ]
    bank_mu = suite["table"].get("bank", {}).get("online_mse", {}).get("mean", float("nan"))
    for m in methods:
        cell = suite["table"][m]["online_mse"]
        drop = bank_mu - cell["mean"] if np.isfinite(bank_mu) else float("nan")
        lines.append(
            r"%s & $%.3f$ ($%.3f$) & $%+.3f$ \\"
            % (labels.get(m, m), cell["mean"], cell["sd"], drop)
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))
    return path
