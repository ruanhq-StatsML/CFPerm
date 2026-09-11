"""Typed shift stepsize (TSS): a statistical LR map, not a product board.

Two identified channels, opposite signs:

    covariate intensity  c_m  = |Cohen's d| of the modality-mean under W
    concept intensity    δ_m  = two-fold excess 0-1 risk of a unimodal
                                ridge head after aligning the block mean

    η_m = η0 · (1 + β δ_m) / (1 + λ c_m)
    η_m = 0  if both channels are below threshold (quiet freeze)

Large c_m lowers the step; large δ_m raises it. The map is not η ∝ δ,
which would freeze under pure covariate shift. π_m = c_m / ∑ c is a
which-head budget only; η ∝ π has the wrong sign for pure covariate shift.

Comparators: constant, step, cosine, inverse-time, plateau, and η ∝ π.
Oracle TSS substitutes the known (c, δ) of the simulation DGP.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy import stats

from msrvtt_continuous_trainer import SeparateHeadProbe, _minibatch, _softmax, _split_modalities
from msrvtt_multimodal_attribution import (
    GROUP_NAMES,
    GROUPS,
    P_X,
    SEED,
    _cohens_d,
    assign_temporal_batches,
    standardize_columns,
)

BETA = 2.5
LAMBDA = 1.5
TAU_C = 0.25
TAU_D = 0.05
ETA0 = 0.10
ETA_MAX_MULT = 2.5
METHODS = (
    "constant",
    "step",
    "cosine",
    "invtime",
    "plateau",
    "fsds_pi",
    "tss",
    "oracle_tss",
)


def _ce(logits, y):
    y = np.asarray(y, dtype=int)
    P = _softmax(logits)
    p = np.clip(P[np.arange(len(y)), y], 1e-8, 1.0)
    return float(-np.mean(np.log(p)))


def _acc(logits, y):
    y = np.asarray(y, dtype=int)
    if y.size == 0:
        return float("nan")
    return float(np.mean(logits.argmax(axis=1) == y))


def modality_mean(X):
    X = np.asarray(X, dtype=float)
    return X.mean(axis=1)


def covariate_intensity(X0, X1):
    """Absolute covariate intensity per modality: |d| on the block mean.

    Returns ``c``, ``share`` (π from c), and the raw signed d.
    """
    X0 = np.asarray(X0, dtype=float)
    X1 = np.asarray(X1, dtype=float)
    signed = {}
    c = {}
    for name, sl in GROUPS.items():
        d = _cohens_d(modality_mean(X0[:, sl]), modality_mean(X1[:, sl]))
        signed[name] = float(d)
        c[name] = float(abs(d))
    tot = sum(c.values()) + 1e-12
    share = {k: c[k] / tot for k in c}
    return c, share, signed


def _gaussian_iw(z0, z1):
    """Density ratio p_1(z)/p_0(z) under two 1-d Gaussians, mean-normalized."""
    z0 = np.asarray(z0, dtype=float).reshape(-1)
    z1 = np.asarray(z1, dtype=float).reshape(-1)
    mu0, mu1 = float(z0.mean()), float(z1.mean())
    v0 = float(z0.var()) + 1e-6
    v1 = float(z1.var()) + 1e-6
    logw = -0.5 * ((z0 - mu1) ** 2 / v1 - (z0 - mu0) ** 2 / v0)
    logw = logw - np.max(logw)
    w = np.exp(np.clip(logw, -8.0, 8.0))
    w = w / (w.mean() + 1e-12)
    return np.clip(w, 0.05, 20.0)


def _anova_f(Xm, y):
    Xm = np.asarray(Xm, dtype=float)
    y = np.asarray(y, dtype=int)
    classes = np.unique(y)
    n = Xm.shape[0]
    if len(classes) < 2 or n <= len(classes):
        return np.zeros(Xm.shape[1])
    grand = Xm.mean(axis=0)
    ssb = np.zeros(Xm.shape[1])
    ssw = np.zeros(Xm.shape[1])
    for c in classes:
        xc = Xm[y == c]
        nc = xc.shape[0]
        if nc == 0:
            continue
        mc = xc.mean(axis=0)
        ssb += nc * (mc - grand) ** 2
        ssw += ((xc - mc) ** 2).sum(axis=0)
    dfb = len(classes) - 1
    dfw = max(n - len(classes), 1)
    return (ssb / dfb) / (ssw / dfw + 1e-12)


def _screen_coords(Xm, y, k=12):
    """Supervised coordinate screen (ANOVA F) for a unimodal plug-in."""
    k = int(max(1, min(k, Xm.shape[1], max(Xm.shape[0] - 1, 1))))
    fvals = _anova_f(Xm, y)
    return np.argpartition(fvals, -k)[-k:]


def _ridge_head(Xm, y, cols, ridge=0.8, n_classes=None):
    X = np.asarray(Xm, dtype=float)[:, cols]
    y = np.asarray(y, dtype=int)
    mu = X.mean(axis=0)
    Z = np.c_[X - mu, np.ones(len(X))]
    if n_classes is None:
        n_classes = int(y.max()) + 1
    yoh = np.eye(n_classes)[np.clip(y, 0, n_classes - 1)]
    gram = Z.T @ Z + float(ridge) * np.eye(Z.shape[1])
    b = np.linalg.solve(gram, Z.T @ yoh)
    return mu, b, n_classes


def _ridge_logits(Xm, cols, mu, b):
    X = np.asarray(Xm, dtype=float)[:, cols] - mu
    Z = np.c_[X, np.ones(len(X))]
    return Z @ b


def concept_intensity(probe, Xs0, y0, Xs1, y1, k=12, ridge=0.8):
    """Two-fold excess 0-1 risk after mean-aligning X.

    Screen coordinates by ANOVA F on batch t-1, fit a ridge head, evaluate
    OOS on the other fold (R_{t-1}) and on mean-aligned batch t (R_t^{al}).
    """
    y0 = np.asarray(y0, dtype=int)
    y1 = np.asarray(y1, dtype=int)
    out = {}
    extras = {}
    n0 = len(y0)
    n_classes = int(max(int(y0.max() if n0 else 0), int(y1.max() if len(y1) else 0))) + 1
    rng = np.random.default_rng(n0 * 17 + 3 * len(y1) + int(y0[:1].sum() if n0 else 0))
    fold = (rng.random(n0) >= 0.5).astype(int)
    if fold.min() == fold.max() and n0 > 1:
        fold[0] = 1 - fold[0]
    for name in GROUP_NAMES:
        mu1 = Xs1[name].mean(axis=0)
        mu0 = Xs0[name].mean(axis=0)
        x1_al = Xs1[name] - (mu1 - mu0)
        a0s, aals = [], []
        for f in (0, 1):
            tr, te = np.flatnonzero(fold == f), np.flatnonzero(fold != f)
            if tr.size < 4 or te.size < 2:
                continue
            cols = _screen_coords(Xs0[name][tr], y0[tr], k=k)
            mu, b, _ = _ridge_head(Xs0[name][tr], y0[tr], cols, ridge=ridge, n_classes=n_classes)
            a0s.append(_acc(_ridge_logits(Xs0[name][te], cols, mu, b), y0[te]))
            aals.append(_acc(_ridge_logits(x1_al, cols, mu, b), y1))
        acc0 = float(np.mean(a0s)) if a0s else float("nan")
        acc_al = float(np.mean(aals)) if aals else float("nan")
        delta = max(0.0, acc0 - acc_al) if np.isfinite(acc0) and np.isfinite(acc_al) else 0.0
        out[name] = float(delta)
        extras[name] = {"acc0": acc0, "acc_al": acc_al}
    return out, extras


def tss_lr(
    c,
    delta,
    eta0=ETA0,
    lam=LAMBDA,
    beta=BETA,
    tau_c=TAU_C,
    tau_d=TAU_D,
    eta_max_mult=ETA_MAX_MULT,
):
    """Signed per-head stepsize. Quiet freeze when both channels are off.

    η ∝ (1+βδ)/(1+λc): covariate intensity lowers the step, concept intensity
    raises it. A numerator of δ alone would freeze under pure covariate shift.
    """
    lrs = {}
    cap = float(eta0) * float(eta_max_mult)
    for name in GROUP_NAMES:
        cm = float(c.get(name, 0.0))
        dm = float(delta.get(name, 0.0))
        if cm < tau_c and dm < tau_d:
            lrs[name] = 0.0
            continue
        eta = float(eta0) * (1.0 + float(beta) * dm) / (1.0 + float(lam) * cm)
        lrs[name] = float(min(max(eta, 0.0), cap))
    return lrs


def scheduler_lrs(
    method,
    t,
    n_batches,
    eta0,
    c,
    share,
    delta,
    true_c=None,
    true_delta=None,
    plateau_eta=None,
):
    """Return a dict of per-head learning rates for one round."""
    if method == "tss":
        return tss_lr(c, delta, eta0=eta0)
    if method == "oracle_tss":
        return tss_lr(true_c or c, true_delta or delta, eta0=eta0)
    if method == "fsds_pi":
        return {g: float(eta0) * float(share[g]) for g in GROUP_NAMES}
    if method == "constant":
        eta = float(eta0)
    elif method == "step":
        eta = float(eta0) * (0.3 if t >= max(2, n_batches // 2) else 1.0)
    elif method == "cosine":
        eta = float(eta0) * 0.5 * (1.0 + np.cos(np.pi * t / max(n_batches - 1, 1)))
    elif method == "invtime":
        eta = float(eta0) / (1.0 + 0.45 * t)
    elif method == "plateau":
        eta = float(plateau_eta if plateau_eta is not None else eta0)
    else:
        raise ValueError("unknown method %s" % method)
    return {g: float(eta) for g in GROUP_NAMES}


@dataclass
class TypedStream:
    X: np.ndarray
    y: np.ndarray
    batch: np.ndarray
    true_c: dict = field(default_factory=dict)
    true_delta: dict = field(default_factory=dict)
    meta: dict = field(default_factory=dict)


def make_typed_stream(
    n_batches=12,
    n_per=40,
    n_classes=4,
    seed=SEED,
    cov=None,
    concept=None,
    concept_at=6,
    noise=1.0,
    signal=1.15,
    rank=12,
):
    """Oracle DGP with known per-modality covariate and concept schedules.

    Class means live in the first ``rank`` coordinates of each block.
    Covariate shift is a class-independent mean drift on the *whole* block
    (so |d| of the coordinate-mean is an absolute intensity, not a diluted
    subspace statistic). Concept drift is a cyclic permutation of the
    class-mean assignment after ``concept_at``. Remaining coordinates are
    high-d noise, so an oversized step overfits batch-specific directions.
    Oracle ``true_c`` / ``true_delta`` are on the same scale as the TSS map,
    not the raw drift coefficients.
    """
    cov = {g: float((cov or {}).get(g, 0.0)) for g in GROUP_NAMES}
    concept = {g: float((concept or {}).get(g, 0.0)) for g in GROUP_NAMES}
    rng = np.random.default_rng(seed)
    rank = int(rank)
    means = {g: rng.normal(scale=float(signal), size=(n_classes, rank)) for g in GROUP_NAMES}
    rows, labs, batches = [], [], []
    true_c = {g: [] for g in GROUP_NAMES}
    true_delta = {g: [] for g in GROUP_NAMES}
    y_per = n_per // n_classes
    n_per = y_per * n_classes
    for t in range(n_batches):
        y = np.repeat(np.arange(n_classes), y_per)
        rng.shuffle(y)
        X = rng.normal(scale=noise, size=(n_per, P_X))
        for name, sl in GROUPS.items():
            use = y.copy()
            drifted = t >= int(concept_at) and concept[name] > 0
            if drifted:
                use = (y + int(round(concept[name]))) % n_classes
            X[:, sl.start : sl.start + rank] += means[name][use]
            if cov[name] > 0:
                X[:, sl] += cov[name] * float(t)
            true_c[name].append(1.0 if cov[name] > 0 else 0.0)
            true_delta[name].append(0.40 if drifted else 0.0)
        rows.append(X)
        labs.append(y)
        batches.append(np.full(n_per, t, dtype=int))
    return TypedStream(
        X=np.vstack(rows),
        y=np.concatenate(labs),
        batch=np.concatenate(batches),
        true_c={g: np.asarray(true_c[g], dtype=float) for g in GROUP_NAMES},
        true_delta={g: np.asarray(true_delta[g], dtype=float) for g in GROUP_NAMES},
        meta={
            "n_batches": n_batches,
            "n_per": n_per,
            "n_classes": n_classes,
            "cov": cov,
            "concept": concept,
            "concept_at": int(concept_at),
            "seed": int(seed),
        },
    )


def graft_midclip_concept(bundle):
    """Controlled concept drift on real X: flip the video-id map at mid-clip."""
    vid = np.asarray(bundle.video_id)
    w = np.asarray(bundle.window_idx, dtype=float)
    _, y0 = np.unique(vid, return_inverse=True)
    n = int(y0.max()) + 1
    mid = float(np.median(w))
    return (y0 + (w >= mid).astype(int)) % n


def stream_from_bundle(bundle, n_batches=10, y=None):
    """MSR-VTT (or any window bundle) as a typed stream. Oracle channels unknown."""
    X = standardize_columns(bundle.X)
    if y is None:
        _, y = np.unique(np.asarray(bundle.video_id), return_inverse=True)
    else:
        _, y = np.unique(np.asarray(y), return_inverse=True)
    batch = assign_temporal_batches(bundle.window_idx, n_batches=n_batches)
    n_batches = int(batch.max()) + 1 if batch.size else int(n_batches)
    return TypedStream(
        X=X,
        y=np.asarray(y, dtype=int),
        batch=batch,
        meta={
            "n_batches": n_batches,
            "n_classes": int(len(np.unique(y))),
            "source": "bundle",
            "n": int(len(y)),
        },
    )


def _true_at(stream, t):
    if not stream.true_c:
        return None, None
    tc = {g: float(stream.true_c[g][t]) for g in GROUP_NAMES}
    td = {g: float(stream.true_delta[g][t]) for g in GROUP_NAMES}
    return tc, td


def run_method(
    stream,
    method="tss",
    eta0=ETA0,
    steps_per_batch=6,
    warmup_steps=4,
    seed=SEED,
    lam=LAMBDA,
    beta=BETA,
):
    """Warmup on B0, then one SGD round per later batch under ``method``."""
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=int)
    batch = np.asarray(stream.batch, dtype=int)
    n_batches = int(batch.max()) + 1
    n_classes = int(y.max()) + 1
    probe = SeparateHeadProbe(n_classes=n_classes, seed=seed)
    rng = np.random.default_rng(seed)
    plateau_eta = float(eta0)
    best_ce = np.inf
    stall = 0
    i0 = np.flatnonzero(batch == 0)
    Xs0, y0 = _split_modalities(X[i0]), y[i0]
    chunk = max(4, i0.size // 2)
    warm = {g: eta0 / 3.0 for g in GROUP_NAMES}
    for head in GROUP_NAMES:
        for _ in range(int(max(1, warmup_steps))):
            sl = _minibatch(rng, i0.size, chunk)
            probe.step({g: Xs0[g][sl] for g in GROUP_NAMES}, y0[sl], warm, active=head)

    history = []
    c_hat = {g: 0.0 for g in GROUP_NAMES}
    d_hat = {g: 0.0 for g in GROUP_NAMES}
    share = {g: 1.0 / 3.0 for g in GROUP_NAMES}

    def snapshot(round_id, phase, idx, lrs, auc_w=float("nan")):
        Xs = _split_modalities(X[idx])
        yy = y[idx]
        logits = probe.logits(Xs)
        bwt = _acc(probe.logits(Xs0), y0)
        return {
            "round": int(round_id),
            "phase": phase,
            "method": method,
            "lr": {g: float(lrs[g]) for g in GROUP_NAMES},
            "c": {g: float(c_hat[g]) for g in GROUP_NAMES},
            "delta": {g: float(d_hat[g]) for g in GROUP_NAMES},
            "pi": {g: float(share[g]) for g in GROUP_NAMES},
            "acc": _acc(logits, yy),
            "ce": _ce(logits, yy),
            "bwt": float(bwt),
            "n": int(idx.size),
        }

    history.append(snapshot(0, "warmup", i0, warm))

    for t in range(1, n_batches):
        ip, ic = np.flatnonzero(batch == t - 1), np.flatnonzero(batch == t)
        Xp, yp = X[ip], y[ip]
        Xc, yc = X[ic], y[ic]
        Xsp, Xsc = _split_modalities(Xp), _split_modalities(Xc)
        c_hat, share, _ = covariate_intensity(Xp, Xc)
        d_hat, _ = concept_intensity(probe, Xsp, yp, Xsc, yc)
        tc, td = _true_at(stream, t)
        lrs = scheduler_lrs(
            method,
            t,
            n_batches,
            eta0,
            c_hat,
            share,
            d_hat,
            true_c=tc,
            true_delta=td,
            plateau_eta=plateau_eta,
        )
        if method == "tss":
            lrs = tss_lr(c_hat, d_hat, eta0=eta0, lam=lam, beta=beta)
        n_steps = max(1, int(steps_per_batch) // 3)
        chunk_t = max(4, ic.size // 2)
        for head in GROUP_NAMES:
            for _ in range(n_steps):
                sl = _minibatch(rng, ic.size, chunk_t)
                probe.step({g: Xsc[g][sl] for g in GROUP_NAMES}, yc[sl], lrs, active=head)
        row = snapshot(t, "adapt", ic, lrs)
        if tc is not None:
            row["true_c"] = tc
            row["true_delta"] = td
        history.append(row)
        if row["ce"] < best_ce - 1e-4:
            best_ce = row["ce"]
            stall = 0
        else:
            stall += 1
            if stall >= 2:
                plateau_eta *= 0.5
                stall = 0

    last = np.flatnonzero(batch == n_batches - 1)
    adapt = [h for h in history if h["phase"] == "adapt"]
    split = int(stream.meta.get("concept_at") or max(n_batches // 2, 1))
    post = [h for h in adapt if h["round"] >= split]
    summary = {
        "method": method,
        "n_batches": n_batches,
        "n_classes": n_classes,
        "eta0": float(eta0),
        "online_acc": float(np.mean([h["acc"] for h in adapt])) if adapt else float("nan"),
        "online_ce": float(np.mean([h["ce"] for h in adapt])) if adapt else float("nan"),
        "post_acc": float(np.mean([h["acc"] for h in post])) if post else float("nan"),
        "last_acc": _acc(probe.logits(_split_modalities(X[last])), y[last]) if last.size else float("nan"),
        "bwt": history[-1]["bwt"] if history else float("nan"),
        "mean_c": {g: float(np.mean([h["c"][g] for h in adapt])) for g in GROUP_NAMES} if adapt else {},
        "mean_delta": {g: float(np.mean([h["delta"][g] for h in adapt])) for g in GROUP_NAMES} if adapt else {},
        "mean_lr": {g: float(np.mean([h["lr"][g] for h in adapt])) for g in GROUP_NAMES} if adapt else {},
        "history": history,
        "meta": dict(stream.meta),
    }
    return summary, probe


REGIMES = {
    "cov_only": dict(
        cov={"video": 0.12, "audio": 0.04, "text": 0.0},
        concept={},
        n_classes=6,
        signal=0.80,
        noise=1.15,
        rank=10,
    ),
    "concept_only": dict(
        cov={},
        concept={"video": 1.0, "audio": 0.0, "text": 0.0},
        concept_at=6,
        n_classes=6,
        signal=0.80,
        noise=1.15,
        rank=10,
    ),
    "both": dict(
        cov={"video": 0.12, "audio": 0.04, "text": 0.0},
        concept={"video": 1.0},
        concept_at=6,
        n_classes=6,
        signal=0.80,
        noise=1.15,
        rank=10,
    ),
}


def run_suite(
    seeds=None,
    n_batches=12,
    n_per=64,
    methods=None,
    regimes=None,
    eta0=ETA0,
    steps_per_batch=6,
):
    """Monte Carlo comparison: regimes × methods × seeds."""
    seeds = list(seeds if seeds is not None else range(SEED, SEED + 8))
    methods = list(methods or METHODS)
    regimes = regimes or REGIMES
    rows = []
    traces = {}
    for rname, spec in regimes.items():
        traces[rname] = {}
        for method in methods:
            traces[rname][method] = []
            for s in seeds:
                stream = make_typed_stream(n_batches=n_batches, n_per=n_per, seed=int(s), **spec)
                summary, _ = run_method(
                    stream,
                    method=method,
                    eta0=eta0,
                    steps_per_batch=steps_per_batch,
                    seed=int(s),
                )
                traces[rname][method].append(summary)
                rows.append(
                    {
                        "regime": rname,
                        "method": method,
                        "seed": int(s),
                        "online_acc": summary["online_acc"],
                        "online_ce": summary["online_ce"],
                        "post_acc": summary["post_acc"],
                        "bwt": summary["bwt"],
                        "last_acc": summary["last_acc"],
                        "mean_lr_video": summary["mean_lr"]["video"],
                        "mean_c_video": summary["mean_c"]["video"],
                        "mean_delta_video": summary["mean_delta"]["video"],
                    }
                )
    table = _summarize_rows(rows)
    ident = identification_report(traces)
    tests = paired_tests(rows)
    return {"rows": rows, "table": table, "identification": ident, "tests": tests, "traces": traces}


def _summarize_rows(rows):
    table = {}
    regimes = sorted({r["regime"] for r in rows})
    methods = []
    for r in rows:
        if r["method"] not in methods:
            methods.append(r["method"])
    for regime in regimes:
        table[regime] = {}
        for method in methods:
            sub = [r for r in rows if r["regime"] == regime and r["method"] == method]
            cell = {}
            for key in ("online_acc", "online_ce", "post_acc", "bwt", "last_acc", "mean_lr_video"):
                v = np.array([r[key] for r in sub], dtype=float)
                cell[key] = {"mean": float(v.mean()), "sd": float(v.std(ddof=1)) if len(v) > 1 else 0.0, "n": int(len(v))}
            table[regime][method] = cell
    return table


def paired_tests(rows, baseline="tss"):
    """Wilcoxon signed-rank of TSS vs each comparator, by regime, on BWT and online acc."""
    out = {}
    regimes = sorted({r["regime"] for r in rows})
    methods = sorted({r["method"] for r in rows if r["method"] != baseline})
    for regime in regimes:
        out[regime] = {}
        base = {(r["seed"], r["method"]): r for r in rows if r["regime"] == regime}
        seeds = sorted({r["seed"] for r in rows if r["regime"] == regime})
        for method in methods:
            rec = {}
            for metric in ("bwt", "online_acc", "post_acc"):
                a = np.array([base[(s, baseline)][metric] for s in seeds if (s, baseline) in base and (s, method) in base])
                b = np.array([base[(s, method)][metric] for s in seeds if (s, baseline) in base and (s, method) in base])
                if len(a) < 4 or np.allclose(a, b):
                    rec[metric] = {"n": int(len(a)), "mean_diff": float(np.mean(a - b)) if len(a) else float("nan"), "p": float("nan")}
                    continue
                try:
                    stat, p = stats.wilcoxon(a, b, zero_method="wilcox", alternative="two-sided")
                except ValueError:
                    stat, p = float("nan"), float("nan")
                rec[metric] = {"n": int(len(a)), "mean_diff": float(np.mean(a - b)), "stat": float(stat), "p": float(p)}
            out[regime][method] = rec
    return out


def identification_report(traces):
    """Does ĉ track true covariate, and δ̂ jump at concept_at?"""
    out = {}
    for rname, by_m in traces.items():
        # any method's history has c/delta; use tss if present else first
        method = "tss" if "tss" in by_m else next(iter(by_m))
        cs, ds, tcs, tds = [], [], [], []
        pre_d, post_d = [], []
        concept_at = None
        for summary in by_m[method]:
            concept_at = summary["meta"].get("concept_at")
            for h in summary["history"]:
                if h["phase"] != "adapt":
                    continue
                cs.append(h["c"]["video"])
                ds.append(h["delta"]["video"])
                if "true_c" in h:
                    tcs.append(h["true_c"]["video"])
                    tds.append(h["true_delta"]["video"])
                    if concept_at is not None and h["round"] < concept_at:
                        pre_d.append(h["delta"]["video"])
                    elif concept_at is not None:
                        post_d.append(h["delta"]["video"])
        rec = {
            "mean_c_video": float(np.mean(cs)) if cs else float("nan"),
            "mean_delta_video": float(np.mean(ds)) if ds else float("nan"),
        }
        if tcs:
            rec["corr_c"] = float(np.corrcoef(cs, tcs)[0, 1]) if len(set(tcs)) > 1 else float("nan")
        if pre_d and post_d:
            rec["delta_pre"] = float(np.mean(pre_d))
            rec["delta_post"] = float(np.mean(post_d))
            rec["delta_jump"] = rec["delta_post"] - rec["delta_pre"]
        out[rname] = rec
    return out


def run_bundle_suite(bundle, n_batches=10, methods=None, eta0=ETA0, seed=SEED, y=None, label="bundle"):
    stream = stream_from_bundle(bundle, n_batches=n_batches, y=y)
    methods = [m for m in (methods or METHODS) if m != "oracle_tss"]
    out = {}
    for method in methods:
        summary, _ = run_method(stream, method=method, eta0=eta0, seed=seed)
        out[method] = {k: v for k, v in summary.items() if k != "history"}
        out[method]["history"] = summary["history"]
        out[method]["label"] = label
    return out


def plot_comparison(suite, path):
    import matplotlib.pyplot as plt
    from msrvtt_attribution_plots import GRID, INK, MUTED, _save, _style

    _style()
    table = suite["table"]
    regimes = [r for r in ("cov_only", "concept_only", "both") if r in table]
    methods = [m for m in METHODS if m in table[regimes[0]]]
    labels = {
        "constant": "constant",
        "step": "step",
        "cosine": "cosine",
        "invtime": "inv-time",
        "plateau": "plateau",
        "fsds_pi": r"$\eta\propto\pi$ (wrong sign)",
        "tss": "TSS",
        "oracle_tss": "oracle TSS",
    }
    colors = {
        "constant": "#9AA3AE",
        "step": "#6B7C8A",
        "cosine": "#2C4A6E",
        "invtime": "#7A8B9A",
        "plateau": "#8A6A3D",
        "fsds_pi": "#A33B24",
        "tss": "#C45C26",
        "oracle_tss": "#2F6B4F",
    }
    fig, axes = plt.subplots(2, len(regimes), figsize=(4.1 * len(regimes), 7.2), sharey="row")
    if len(regimes) == 1:
        axes = np.array(axes).reshape(2, 1)
    titles = {"cov_only": "Covariate only", "concept_only": "Concept only", "both": "Both"}
    x = np.arange(len(methods))
    for j, regime in enumerate(regimes):
        for i, metric in enumerate(("bwt", "post_acc")):
            ax = axes[i, j]
            means = [table[regime][m][metric]["mean"] for m in methods]
            sds = [table[regime][m][metric]["sd"] for m in methods]
            ax.bar(x, means, yerr=sds, color=[colors[m] for m in methods], ecolor=MUTED, capsize=2.5, width=0.78)
            ax.set_xticks(x)
            ax.set_xticklabels([labels[m] for m in methods], rotation=55, ha="right", fontsize=8)
            ax.grid(True, axis="y", color=GRID)
            if j == 0:
                ax.set_ylabel("BWT (acc. on batch 0)" if metric == "bwt" else "post-change accuracy")
            if i == 0:
                ax.set_title(titles[regime], loc="left", fontsize=12, fontweight="bold")
            ax.set_ylim(0.0, 1.05)
    fig.suptitle("Typed shift stepsize vs LR schedulers", fontsize=13.5, fontweight="bold", color=INK, x=0.04, ha="left")
    fig.text(
        0.04,
        0.01,
        "TSS: η_m = η0 (1+β δ_m) / (1+λ c_m).  Covariate c_m lowers the step; concept δ_m raises it.  "
        "η∝π uses the same c_m as a simplex share with the opposite sign.  Error bars are seed s.d.",
        fontsize=8.2,
        color=MUTED,
    )
    return _save(fig, path)


def plot_eta_paths(suite, path):
    import matplotlib.pyplot as plt
    from msrvtt_attribution_plots import GRID, INK, MUTED, _save, _style

    _style()
    show = ["cosine", "plateau", "fsds_pi", "tss"]
    colors = {"cosine": "#2C4A6E", "plateau": "#8A6A3D", "fsds_pi": "#A33B24", "tss": "#C45C26"}
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.0), sharey=True)
    for ax, regime, title in zip(axes, ("cov_only", "concept_only"), ("Covariate only", "Concept only")):
        for method in show:
            recs = suite["traces"][regime][method]
            rounds = recs[0]["history"]
            t = [h["round"] for h in rounds]
            lr = np.array([[h["lr"]["video"] for h in r["history"]] for r in recs], dtype=float)
            mu, sd = lr.mean(axis=0), lr.std(axis=0)
            ax.plot(t, mu, color=colors[method], lw=2.1, label=method if method != "fsds_pi" else r"$\eta\propto\pi$")
            ax.fill_between(t, mu - sd, mu + sd, color=colors[method], alpha=0.14, lw=0)
        if regime == "concept_only":
            at = recs[0]["meta"].get("concept_at", 7)
            ax.axvline(at, color=MUTED, ls="--", lw=0.9)
        ax.set_title(title, loc="left", fontsize=12, fontweight="bold")
        ax.set_xlabel("round")
        ax.grid(True, color=GRID)
        ax.legend(frameon=False, fontsize=8.5)
    axes[0].set_ylabel("video-head  η")
    fig.suptitle("Video-head stepsize paths", fontsize=13.2, fontweight="bold", color=INK, x=0.06, ha="left")
    return _save(fig, path)


def write_tex_table(suite, path):
    table = suite["table"]
    tests = suite["tests"]
    ident = suite["identification"]
    methods = [m for m in METHODS if m in table[next(iter(table))]]
    labels = {
        "constant": "constant",
        "step": "step",
        "cosine": "cosine",
        "invtime": "inv-time",
        "plateau": "plateau",
        "fsds_pi": r"$\eta\propto\pi$",
        "tss": "TSS",
        "oracle_tss": "oracle TSS",
    }
    regimes = [r for r in ("cov_only", "concept_only", "both") if r in table]
    rt = {"cov_only": "covariate only", "concept_only": "concept only", "both": "both"}
    lines = [
        r"% Typed shift stepsize vs LR schedulers. Auto-generated.",
        r"\begin{table}[ht]\centering",
        r"\caption{Typed shift stepsize (TSS) versus standard LR schedulers on an oracle multimodal stream.",
        r"Covariate intensity $c_m$ is $|d|$ of the modality mean; concept intensity $\delta_m$ is two-fold excess 0-1 risk after mean-alignment of $X$.",
        r"TSS uses $\eta_m=\eta_0(1+\beta\delta_m)/(1+\lambda c_m)$ with a quiet freeze.",
        r"$\eta\propto\pi$ uses the same $c_m$ as a simplex share (wrong sign for pure covariate shift).",
        r"Entries are mean (s.d.) over seeds. BWT is accuracy on batch 0 after the stream;",
        r"post-change accuracy is the mean on rounds at or after the known concept time.}",
        r"\label{tab:tss-vs-schedulers}",
        r"\small",
        r"\begin{tabular}{@{}ll cccc@{}}\toprule",
        r"Regime & Method & online acc. & post-change & BWT & video $\bar\eta$ \\",
        r"\midrule",
    ]
    for regime in regimes:
        first = True
        for method in methods:
            cell = table[regime][method]
            def fmt(k):
                return r"$%.3f$ ($%.3f$)" % (cell[k]["mean"], cell[k]["sd"])
            name = rt[regime] if first else ""
            first = False
            star = ""
            if method != "tss" and regime in tests and method in tests[regime]:
                key = "bwt" if regime == "cov_only" else "post_acc"
                rec = tests[regime][method].get(key, {})
                p = rec.get("p", float("nan"))
                diff = rec.get("mean_diff", float("nan"))
                if p == p and p < 0.05 and diff == diff and diff > 0:
                    star = r"$^{\ast}$"
            lines.append(
                r"%s & %s%s & %s & %s & %s & $%.4f$ \\"
                % (name, labels[method], star, fmt("online_acc"), fmt("post_acc"), fmt("bwt"), cell["mean_lr_video"]["mean"])
            )
        lines.append(r"\midrule")
    if lines[-1] == r"\midrule":
        lines[-1] = r"\bottomrule"
    lines.append(r"\end{tabular}\\[0.4em]")
    lines.append(r"{\footnotesize Wilcoxon signed-rank, TSS $-$ comparator: $^{\ast}$ $p<0.05$ and mean difference $>0$ on BWT (covariate-only) or post-change accuracy (concept/both).")
    if "concept_only" in ident and "delta_jump" in ident["concept_only"]:
        lines.append(
            r" Concept-only $\hat\delta_{\mathrm{v}}$ jump at the known change point: $%.3f \to %.3f$."
            % (ident["concept_only"]["delta_pre"], ident["concept_only"]["delta_post"])
        )
    if "cov_only" in ident:
        lines.append(r" Covariate-only mean $\hat c_{\mathrm{v}}=%.3f$." % ident["cov_only"]["mean_c_video"])
    lines.append(r"}")
    lines.append(r"\end{table}")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path
