"""PO-risk mask attribution + modality-π hop (on top of heatmap hop-ridge).

Locked recipe:
  Amazon          → hop_ridge (heatmap IW)
  Multimodal      → modality-π hop; π_m from ``benchmark_feature_selection``
                    via ``attribution_adapter.modality_pi_shares``
  PO-risk VIMP    → mask tokens/patches only (zero / noise / missing)

  w_s = exp(γ Σ_m π_m (cos(μ_s^m, μ_{t-1}^m) − 1))

Swap the body of ``benchmark_feature_selection`` to change the selector;
shares still come from ``modality_mass(vimp)``. No new loss.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from sklearn.linear_model import Ridge

from amazon_continuous_batches import make_amazon_like_stream
from amazon_mse_prototype import _mse, describe_shapes
from attribution_adapter import (
    attribution_ewma_rate,
    block_means,
    hop_sample_weights,
    load_amazon_stream,
    modality_hop_weights,
    modality_pi_shares,
    pi_l1_change,
    rolling_hop_stats,
    run_hop_ridge,
    run_ridge_past,
)
from msrvtt_multimodal_attribution import (
    GROUP_NAMES,
    GROUPS,
    SEED,
    assign_temporal_batches,
    crossfit_po,
    load_window_bundle,
    make_synthetic_bundle,
    po_tau_vimp,
    standardize_columns,
)

METHODS = (
    "ridge_past",
    "hop_ridge",
    "hop_po_feat",
    "hop_po_amp",
    "hop_po_mod",
    "hop_po_both",
)


def po_feature_vimp_safe(X0, X1, y0, y1, seed=SEED, n_estimators=60):
    """PO-risk feature VIMP on a two-batch hop (W∈{0,1}).

    Uses ``po_tau_vimp`` (RF importances of τ̂ on the PO score). Do not call
    ``logo_po_risk`` here: that helper hard-codes MSR-VTT GROUPS (2048-d).
    """
    X0 = np.asarray(X0, dtype=float)
    X1 = np.asarray(X1, dtype=float)
    if X0.shape[0] < 8 or X1.shape[0] < 8 or X0.shape[1] != X1.shape[1]:
        p = int(X0.shape[1] if X0.ndim == 2 else X1.shape[1])
        return np.ones(max(p, 1), dtype=float)
    X = np.vstack([X0, X1])
    y = np.concatenate([np.asarray(y0, dtype=float), np.asarray(y1, dtype=float)])
    W = np.concatenate([np.zeros(len(X0), dtype=int), np.ones(len(X1), dtype=int)])
    po, _, _ = crossfit_po(X, y, W, seed=seed, n_splits=3, n_estimators=n_estimators)
    _, vimp, _ = po_tau_vimp(X, po, seed=seed + 4, n_estimators=n_estimators)
    return np.maximum(np.asarray(vimp, dtype=float), 0.0)


def feature_weights_from_vimp(vimp, floor=0.05, power=1.0, invert=True):
    """Map PO-risk VIMP → feature weights.

    Default ``invert=True`` (TSS-aligned): high PO / domain VIMP → *lower*
    weight. ``invert=False`` amplifies high-VIMP coords (ablation).
    """
    v = np.asarray(vimp, dtype=float)
    v = np.maximum(v, 0.0)
    if float(v.sum()) <= 0:
        return np.ones_like(v)
    z = (v / (v.mean() + 1e-12)) ** float(power)
    z = np.maximum(z, float(floor))
    if invert:
        z = 1.0 / z
    return z / (z.mean() + 1e-12)


def modality_shares_from_vimp(vimp, groups=None):
    groups = groups or GROUPS
    v = np.asarray(vimp, dtype=float)
    mass = {g: float(np.maximum(v[sl], 0.0).sum()) for g, sl in groups.items()}
    tot = sum(mass.values()) + 1e-12
    return {g: mass[g] / tot for g in mass}


def apply_mask(X, cols, mode="zero", seed=SEED, missing_fill="mean"):
    X = np.asarray(X, dtype=float).copy()
    cols = np.asarray(cols, dtype=int)
    if cols.size == 0:
        return X
    rng = np.random.default_rng(seed)
    if mode == "zero":
        X[:, cols] = 0.0
    elif mode == "noise":
        scale = X[:, cols].std(axis=0, ddof=0) + 1e-6
        X[:, cols] = rng.normal(0.0, scale, size=(X.shape[0], cols.size))
    elif mode == "missing":
        if missing_fill == "mean":
            X[:, cols] = X[:, cols].mean(axis=0)
        else:
            X[:, cols] = 0.0
    else:
        raise ValueError("unknown mask mode %s" % mode)
    return X


def topk_from_vimp(vimp, k=16):
    v = np.asarray(vimp, dtype=float)
    k = int(max(1, min(k, v.size)))
    return np.argpartition(v, -k)[-k:]


def _scale_features(X, feat_w):
    return np.asarray(X, dtype=float) * np.asarray(feat_w, dtype=float)


def run_hop_po_ridge(
    stream,
    use_feat=True,
    use_mod=False,
    groups=None,
    alpha=3.0,
    gamma=4.0,
    n_estimators=40,
    seed=SEED,
    invert_feat=True,
):
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    k = int(batch.max()) + 1
    groups = groups or {"text": slice(0, X.shape[1])}
    stats = rolling_hop_stats(stream)
    history = []
    feat_w = np.ones(X.shape[1], dtype=float)
    shares = {g: 1.0 / max(len(groups), 1) for g in groups}
    block_mus = block_means(X, batch, groups)

    for t in range(1, k):
        tr, ic = batch < t, batch == t
        if t >= 2:
            X0, y0 = X[batch == t - 2], y[batch == t - 2]
            X1, y1 = X[batch == t - 1], y[batch == t - 1]
            vimp = po_feature_vimp_safe(X0, X1, y0, y1, seed=seed + t, n_estimators=n_estimators)
        else:
            X0 = X1 = X[batch == 0]
            vimp = np.ones(X.shape[1], dtype=float)

        if use_feat:
            feat_w = feature_weights_from_vimp(vimp, invert=invert_feat)

        if use_mod and len(groups) > 1 and t >= 2 and len(X0) >= 8 and len(X1) >= 8:
            # π_m from benchmark_feature_selection (drop-in selector hook)
            shares, fs_vimp, _ = modality_pi_shares(
                X0,
                X1,
                seed=seed + 17 + t,
                n_estimators=n_estimators,
                prev_shares=shares,
                ewma=0.3,
            )
            vimp = fs_vimp
            w = modality_hop_weights(block_mus, batch, t, shares, gamma=gamma)
        else:
            w = hop_sample_weights(stats["mus"], batch, t, gamma=gamma)

        Xs = _scale_features(X, feat_w) if use_feat else X
        clf = Ridge(alpha=float(alpha))
        clf.fit(Xs[tr], y[tr], sample_weight=w[tr])
        pred = np.clip(clf.predict(Xs[ic]), 1.0, 5.0)
        history.append(
            {
                "round": int(t),
                "online_mse": _mse(pred, y[ic]),
                "c": float(stats["c"][t]),
                "delta": float(stats["delta"][t]),
                "shares": {g: float(shares.get(g, 0.0)) for g in groups},
                "topk": [int(i) for i in topk_from_vimp(vimp, k=8)],
                "n": int(ic.sum()),
            }
        )

    online = np.array([h["online_mse"] for h in history], dtype=float)
    if use_feat and use_mod:
        name = "hop_po_both"
    elif use_feat and not invert_feat:
        name = "hop_po_amp"
    elif use_feat:
        name = "hop_po_feat"
    else:
        name = "hop_po_mod"
    return {
        "method": name,
        "online_mse": float(online.mean()) if online.size else float("nan"),
        "cum_mse": float(online.sum()) if online.size else float("nan"),
        "last_online_mse": float(online[-1]) if online.size else float("nan"),
        "online_path": online.tolist(),
        "history": history,
        "meta": dict(getattr(stream, "meta", {}) or {}),
    }


def run_mask_ablation(
    stream,
    k_tokens=16,
    modes=("none", "zero", "noise", "missing"),
    alpha=3.0,
    gamma=4.0,
    n_estimators=40,
    seed=SEED,
):
    """Refit hop-ridge after masking top-k PO tokens; report ΔMSE vs none."""
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    k = int(batch.max()) + 1
    stats = rolling_hop_stats(stream)
    vimp = po_feature_vimp_safe(
        X[batch == 0],
        X[batch == min(1, k - 1)],
        y[batch == 0],
        y[batch == min(1, k - 1)],
        seed=seed,
        n_estimators=n_estimators,
    )
    cols = topk_from_vimp(vimp, k=k_tokens)
    out = {
        "topk": [int(i) for i in cols],
        "vimp_topk": [float(vimp[i]) for i in cols],
        "modes": {},
    }
    for mode in modes:
        history = []
        for t in range(1, k):
            tr, ic = batch < t, batch == t
            w = hop_sample_weights(stats["mus"], batch, t, gamma=gamma)
            Xt = X if mode == "none" else apply_mask(X, cols, mode=mode, seed=seed + t)
            clf = Ridge(alpha=float(alpha))
            clf.fit(Xt[tr], y[tr], sample_weight=w[tr])
            pred = np.clip(clf.predict(Xt[ic]), 1.0, 5.0)
            history.append(_mse(pred, y[ic]))
        online = float(np.mean(history)) if history else float("nan")
        out["modes"][mode] = {"online_mse": online, "path": history}
    base = out["modes"]["none"]["online_mse"]
    for mode, cell in out["modes"].items():
        cell["delta_mse"] = float(cell["online_mse"] - base) if np.isfinite(base) else float("nan")
    return out


def msrvtt_as_stream(n_batches=5, seed=SEED, synthetic=False):
    if synthetic:
        bundle = make_synthetic_bundle(n_videos=8, n_windows=20, seed=seed)
    else:
        bundle = load_window_bundle()
    X = standardize_columns(bundle.X)
    _, y_id = np.unique(np.asarray(bundle.video_id), return_inverse=True)
    batch = assign_temporal_batches(bundle.window_idx, n_batches=int(n_batches))
    return {
        "X": X,
        "y_id": y_id.astype(int),
        "batch": batch,
        "groups": GROUPS,
        "meta": {
            "source": "synthetic" if synthetic else "msrvtt",
            "n_batches": int(n_batches),
            "n": int(X.shape[0]),
            "p": int(X.shape[1]),
            "n_classes": int(y_id.max()) + 1,
        },
    }


def _multi_mse(pred, Y):
    return float(np.mean((np.asarray(pred) - np.asarray(Y)) ** 2))


def _group_input_index(groups, target):
    idx = []
    for g, sl in groups.items():
        if g == target:
            continue
        idx.extend(range(sl.start, sl.stop))
    return np.asarray(idx, dtype=int)


def run_msrvtt_hop_po(
    stream,
    use_mod=True,
    alpha=3.0,
    gamma=4.0,
    seed=SEED,
    n_estimators=40,
    target="text",
    ewma_pi="auto",
):
    """Online cross-modal Ridge: other modalities → target block.

    π_m from ``modality_pi_shares`` → ``benchmark_feature_selection``.
    ``ewma_pi``:
      - float in [0,1]: fixed persistence on previous π
      - ``"auto"`` (default): attribution sets λ_new via
        ``attribution_ewma_rate(ĉ, δ̂, Δπ)``; persistence = 1−λ_new
        so concept / π jumps → forget fast; stable hops → keep ensemble π
    """
    from amazon_continuous_batches import (
        batch_mean_cosine,
        concept_intensity_mse,
        heatmap_hop_c,
    )

    X, batch = stream["X"], stream["batch"]
    groups = stream["groups"]
    if target not in groups:
        raise ValueError("unknown target modality %s" % target)
    in_idx = _group_input_index(groups, target)
    Y = X[:, groups[target]]
    Xin = X[:, in_idx]
    k = int(batch.max()) + 1
    block_mus = block_means(X, batch, groups)
    mus = np.stack([X[batch == t].mean(0) for t in range(k)])
    R = batch_mean_cosine(X, batch, n_batches=k)
    history = []
    shares = {g: 1.0 / 3.0 for g in GROUP_NAMES}
    auto_pi = isinstance(ewma_pi, str) and ewma_pi.lower() == "auto"

    for t in range(1, k):
        tr, ic = batch < t, batch == t
        lam = None
        c_t = d_t = 0.0
        if t >= 2:
            c_t = float(heatmap_hop_c(R, t - 1))
            d_t, _ = concept_intensity_mse(
                X[batch == t - 2],
                stream["y_id"][batch == t - 2].astype(float),
                X[batch == t - 1],
                stream["y_id"][batch == t - 1].astype(float),
            )
            raw, vimp, _meta = modality_pi_shares(
                X[batch == t - 2],
                X[batch == t - 1],
                seed=seed + t,
                n_estimators=n_estimators,
                ewma=0.0,
            )
            dpi = pi_l1_change(shares, raw)
            if auto_pi:
                lam = attribution_ewma_rate(c_t, d_t, pi_change=dpi)
                persist = 1.0 - lam
            else:
                persist = float(ewma_pi)
                lam = 1.0 - persist
            shares, vimp, _meta = modality_pi_shares(
                X[batch == t - 2],
                X[batch == t - 1],
                seed=seed + t,
                n_estimators=n_estimators,
                prev_shares=shares,
                ewma=persist,
            )
        else:
            vimp = None
        if use_mod:
            w = modality_hop_weights(block_mus, batch, t, shares, gamma=gamma)
        else:
            w = hop_sample_weights(mus, batch, t, gamma=gamma)
        clf = Ridge(alpha=float(alpha))
        clf.fit(Xin[tr], Y[tr], sample_weight=w[tr])
        pred = clf.predict(Xin[ic])
        history.append(
            {
                "round": int(t),
                "online_mse": _multi_mse(pred, Y[ic]),
                "shares": {g: float(shares[g]) for g in GROUP_NAMES},
                "topk": [] if vimp is None else [int(i) for i in topk_from_vimp(vimp, k=8)],
                "c": c_t,
                "delta": float(d_t),
                "ewma_rate": None if lam is None else float(lam),
                "ewma_pi_persist": None if lam is None else float(1.0 - lam),
            }
        )

    online = np.array([h["online_mse"] for h in history], dtype=float)
    rates = [h["ewma_rate"] for h in history if h["ewma_rate"] is not None]
    return {
        "method": "msrvtt_hop_pi" if use_mod else "msrvtt_hop",
        "target": target,
        "online_mse": float(online.mean()) if online.size else float("nan"),
        "online_path": online.tolist(),
        "history": history,
        "meta": stream["meta"],
        "ewma_pi": "auto" if auto_pi else float(ewma_pi),
        "mean_ewma_rate": float(np.mean(rates)) if rates else float("nan"),
        "gamma": float(gamma),
    }



def run_msrvtt_patch_mask(
    stream,
    modes=("none", "zero", "noise", "missing"),
    k_per_mod=32,
    seed=SEED,
    target="text",
    n_estimators=40,
):
    """Mask top-k PO coords inside each input modality; ΔMSE on cross-modal hop Ridge."""
    X, batch = stream["X"], stream["batch"]
    groups = stream["groups"]
    in_idx = _group_input_index(groups, target)
    Y = X[:, groups[target]]
    k = int(batch.max()) + 1
    mus = np.stack([X[batch == t].mean(0) for t in range(k)])
    vimp = po_feature_vimp_safe(
        X[batch == 0],
        X[batch == min(1, k - 1)],
        stream["y_id"][batch == 0],
        stream["y_id"][batch == min(1, k - 1)],
        seed=seed,
        n_estimators=n_estimators,
    )
    cols = {}
    for g, sl in groups.items():
        if g == target:
            continue
        local = vimp[sl]
        top = topk_from_vimp(local, k=k_per_mod)
        cols[g] = (sl.start + top).astype(int)

    out = {
        "target": target,
        "cols": {g: [int(i) for i in cols[g]] for g in cols},
        "by_modality": {},
    }
    for g in cols:
        out["by_modality"][g] = {}
        for mode in modes:
            path = []
            for t in range(1, k):
                tr, ic = batch < t, batch == t
                w = hop_sample_weights(mus, batch, t, gamma=4.0)
                Xt = X if mode == "none" else apply_mask(X, cols[g], mode=mode, seed=seed + t)
                Xin = Xt[:, in_idx]
                clf = Ridge(alpha=3.0)
                clf.fit(Xin[tr], Y[tr], sample_weight=w[tr])
                path.append(_multi_mse(clf.predict(Xin[ic]), Y[ic]))
            mse = float(np.mean(path)) if path else float("nan")
            out["by_modality"][g][mode] = {"online_mse": mse, "path": path}
        base = out["by_modality"][g]["none"]["online_mse"]
        for mode, cell in out["by_modality"][g].items():
            cell["delta_mse"] = (
                float(cell["online_mse"] - base) if np.isfinite(base) else float("nan")
            )
    return out


def run_amazon_po_suite(seeds=None, n_per=240, n_batches=9, source="amazon", n_estimators=40):
    seeds = list(seeds if seeds is not None else range(SEED, SEED + 4))
    rows = []
    traces = {m: [] for m in METHODS}
    masks = []
    shapes = None
    for s in seeds:
        if source == "synthetic":
            stream = make_amazon_like_stream(
                n_batches=n_batches, n_per=n_per, seed=int(s), cov=0.35
            )
        else:
            stream = load_amazon_stream(n_per=n_per, n_batches=n_batches, seed=int(s))
        if shapes is None:
            shapes = describe_shapes(stream)
        text_g = {"text": slice(0, int(stream.X.shape[1]))}
        recs = {
            "ridge_past": run_ridge_past(stream),
            "hop_ridge": run_hop_ridge(stream),
            "hop_po_feat": run_hop_po_ridge(
                stream,
                use_feat=True,
                use_mod=False,
                seed=int(s),
                n_estimators=n_estimators,
                invert_feat=True,
            ),
            "hop_po_amp": run_hop_po_ridge(
                stream,
                use_feat=True,
                use_mod=False,
                seed=int(s),
                n_estimators=n_estimators,
                invert_feat=False,
            ),
            "hop_po_mod": run_hop_po_ridge(
                stream,
                use_feat=False,
                use_mod=True,
                groups=text_g,
                seed=int(s),
                n_estimators=n_estimators,
            ),
            "hop_po_both": run_hop_po_ridge(
                stream,
                use_feat=True,
                use_mod=True,
                groups=text_g,
                seed=int(s),
                n_estimators=n_estimators,
                invert_feat=True,
            ),
        }
        for m, rec in recs.items():
            traces[m].append(rec)
            rows.append({"seed": int(s), "method": m, "online_mse": rec["online_mse"]})
        masks.append(
            {"seed": int(s), **run_mask_ablation(stream, seed=int(s), n_estimators=n_estimators)}
        )

    table = {}
    for m in METHODS:
        v = np.array([r["online_mse"] for r in rows if r["method"] == m], dtype=float)
        table[m] = {
            "online_mse": {
                "mean": float(v.mean()),
                "sd": float(v.std(ddof=1) if len(v) > 1 else 0.0),
            },
            "n": int(len(v)),
        }
    mask_summary = {}
    for mode in masks[0]["modes"]:
        deltas = np.array([m["modes"][mode]["delta_mse"] for m in masks], dtype=float)
        mses = np.array([m["modes"][mode]["online_mse"] for m in masks], dtype=float)
        mask_summary[mode] = {
            "online_mse": {
                "mean": float(mses.mean()),
                "sd": float(mses.std(ddof=1) if len(mses) > 1 else 0.0),
            },
            "delta_mse": {
                "mean": float(deltas.mean()),
                "sd": float(deltas.std(ddof=1) if len(deltas) > 1 else 0.0),
            },
        }
    return {
        "table": table,
        "rows": rows,
        "traces": traces,
        "mask": mask_summary,
        "mask_runs": masks,
        "shapes": shapes,
        "seeds": [int(s) for s in seeds],
        "source": source,
        "methods": list(METHODS),
    }


def run_msrvtt_po_suite(n_batches=5, seed=SEED, synthetic=False, target="text", n_estimators=40):
    stream = msrvtt_as_stream(n_batches=n_batches, seed=seed, synthetic=synthetic)
    base = run_msrvtt_hop_po(
        stream, use_mod=False, seed=seed, target=target, n_estimators=n_estimators
    )
    mod = run_msrvtt_hop_po(
        stream, use_mod=True, seed=seed, target=target, n_estimators=n_estimators
    )
    mask = run_msrvtt_patch_mask(
        stream, seed=seed, target=target, n_estimators=n_estimators
    )
    return {
        "hop": base,
        "hop_po_mod": mod,
        "mask": mask,
        "lift_mod": float(base["online_mse"] - mod["online_mse"]),
        "meta": {**stream["meta"], "probe": "crossmodal→%s" % target},
        "shares_last": mod["history"][-1]["shares"] if mod["history"] else {},
    }


def plot_po_suite(amazon, msrvtt, path):
    import matplotlib.pyplot as plt
    from msrvtt_attribution_plots import COLORS, GRID, INK, MUTED, _save, _style

    _style()
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.0))

    ax = axes[0]
    methods = amazon["methods"]
    means = [amazon["table"][m]["online_mse"]["mean"] for m in methods]
    sds = [amazon["table"][m]["online_mse"]["sd"] for m in methods]
    y = np.arange(len(methods))[::-1]
    ax.barh(
        y,
        means,
        xerr=sds,
        color="#C45C26",
        height=0.62,
        error_kw={"ecolor": "#555", "capsize": 2},
    )
    ax.set_yticks(y)
    ax.set_yticklabels(methods, fontsize=8)
    ax.set_xlabel("online MSE")
    ax.set_title("Amazon: PO×hop-ridge", loc="left", fontweight="bold")
    ax.grid(True, color=GRID, axis="x")

    ax = axes[1]
    modes = [m for m in ("zero", "noise", "missing") if m in amazon["mask"]]
    vals = [amazon["mask"][m]["delta_mse"]["mean"] for m in modes]
    ax.bar(np.arange(len(modes)), vals, color="#2C4A6E")
    ax.set_xticks(np.arange(len(modes)))
    ax.set_xticklabels(modes)
    ax.set_ylabel("ΔMSE vs unmasked")
    ax.set_title("Amazon token mask (top-PO)", loc="left", fontweight="bold")
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.grid(True, color=GRID, axis="y")

    ax = axes[2]
    mods = [g for g in GROUP_NAMES if g in msrvtt["mask"]["by_modality"]]
    x = np.arange(len(mods))
    width = 0.25
    for i, mode in enumerate(("zero", "noise", "missing")):
        vals = [msrvtt["mask"]["by_modality"][g][mode]["delta_mse"] for g in mods]
        ax.bar(
            x + (i - 1) * width,
            vals,
            width,
            label=mode,
            color=["#9AA3AE", "#2C4A6E", "#C45C26"][i],
        )
    ax.set_xticks(x)
    ax.set_xticklabels(mods)
    ax.set_ylabel("ΔMSE")
    ax.set_title("MSR-VTT patch mask (→text)", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.5)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.grid(True, color=GRID, axis="y")
    _ = COLORS

    fig.suptitle(
        "PO-risk feature selection × heatmap hop-ridge",
        fontsize=13,
        fontweight="bold",
        color=INK,
        x=0.04,
        ha="left",
    )
    fig.text(
        0.04,
        0.01,
        "Feature weights and modality π_m from PO-risk VIMP on the causal hop. "
        "Mask = zero / noise / column-mean on top-PO tokens or patches; ΔMSE is the attribution.",
        fontsize=8.0,
        color=MUTED,
    )
    return _save(fig, Path(path))


def write_po_tex(amazon, msrvtt, path):
    lines = [
        r"% PO-risk × hop-ridge. Auto-generated.",
        r"\begin{table}[ht]\centering",
        r"\caption{PO-risk feature selection on the heatmap hop. Amazon online rating MSE;",
        r"token mask $\Delta$MSE for top-PO coordinates; MSR-VTT modality-specific hop weights",
        r"and per-block patch mask $\Delta$MSE (cross-modal Ridge: video+audio $\to$ text).}",
        r"\label{tab:po-hop-ridge}",
        r"\small\begin{tabular}{@{}lc@{}}\toprule",
        r"Amazon method & online MSE \\",
        r"\midrule",
    ]
    for m in amazon["methods"]:
        cell = amazon["table"][m]["online_mse"]
        lines.append(
            r"%s & $%.3f$ ($%.3f$) \\" % (m.replace("_", r"\_"), cell["mean"], cell["sd"])
        )
    lines.append(r"\midrule")
    lines.append(r"\multicolumn{2}{@{}l}{Amazon token mask $\Delta$MSE (top-PO)} \\")
    for mode, cell in amazon["mask"].items():
        if mode == "none":
            continue
        d = cell["delta_mse"]
        lines.append(r"mask %s & $%+.3f$ ($%.3f$) \\" % (mode, d["mean"], d["sd"]))
    lines.append(r"\midrule")
    lines.append(
        r"MSR-VTT hop / hop+$\pi_m$ (cross-modal) & $%.4f$ / $%.4f$ \\"
        % (msrvtt["hop"]["online_mse"], msrvtt["hop_po_mod"]["online_mse"])
    )
    for g in msrvtt["mask"]["by_modality"]:
        z = msrvtt["mask"]["by_modality"][g]["zero"]["delta_mse"]
        lines.append(r"mask %s patches (zero) $\Delta$MSE & $%+.4f$ \\" % (g, z))
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))
    return path
