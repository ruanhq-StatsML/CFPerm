"""Online quantification of OOD-risk-as-weights (MSE + regret + BWT).

Boards in this repo (Waymo / Interstate / raw stock feeds are *not* checked in —
use the drift surrogates or Amazon / Diabetes):

  amazon      — category stream (TF-IDF ratings)
  synthetic   — Amazon-like covariate/concept drift
  diabetes    — source→target Readmission from datasets/datasets.zip
  interstate  — synthetic multi-regime traffic-speed surrogate
  stock       — synthetic multi-regime return surrogate

Weight channels (existing methods only — no new loss):
  uniform  — ridge_past
  hop      — heatmap cosine IW (feature-space proximity)
  dga      — Fan–Grangier–Ablin gradient alignment
  mahal    — Mahalanobis proximity to last-batch mean (OOD score → weight)
  attr     — TSS-gated bank⊕hop adapter

Regret: cumulative online MSE minus uniform ridge_past on the same seed
(natural "no OOD weighting" comparator for closed-form Ridge).
"""
from __future__ import annotations

import csv
import io
import json
import zipfile
from pathlib import Path

import numpy as np

from amazon_continuous_batches import AmazonStream, make_amazon_like_stream
from attribution_adapter import (
    SEED,
    dga_alignments,
    dga_mirror_step,
    dga_sample_weights,
    hop_sample_weights,
    load_amazon_stream,
    rolling_hop_stats,
    run_attr_adapter,
    run_dga_ridge,
    run_hop_ridge,
    run_ridge_past,
    _mse,
    _ridge_fit,
    _ridge_predict,
)

ROOT = Path(__file__).resolve().parents[2]
DIABETES_ZIP = ROOT / "datasets" / "datasets.zip"
OUT = ROOT / "results" / "ood_weight_online_eval"

BOARDS = ("amazon", "synthetic", "diabetes", "interstate", "stock")
WEIGHT_METHODS = ("uniform", "hop", "dga", "mahal", "attr")


def _as_stream(X, y, batch, meta=None):
    return AmazonStream(
        X=np.asarray(X, dtype=float),
        y=np.asarray(y, dtype=float),
        batch=np.asarray(batch, dtype=int),
        meta=dict(meta or {}),
    )


def load_board(board, n_per=240, n_batches=9, seed=SEED, **kw):
    board = str(board).lower()
    if board == "amazon":
        return load_amazon_stream(n_per=n_per, n_batches=n_batches, seed=seed)
    if board == "synthetic":
        return make_amazon_like_stream(
            n_batches=n_batches,
            n_per=n_per,
            seed=seed,
            cov=float(kw.get("cov", 0.45)),
            concept=float(kw.get("concept", 0.15)),
        )
    if board == "diabetes":
        return load_diabetes_stream(
            n_per=n_per,
            n_batches=n_batches,
            seed=seed,
            target=str(kw.get("target", "num_medications")),
        )
    if board == "interstate":
        return make_interstate_like_stream(
            n_batches=n_batches, n_per=n_per, seed=seed, p=int(kw.get("p", 24))
        )
    if board == "stock":
        return make_stock_like_stream(
            n_batches=n_batches, n_per=n_per, seed=seed, p=int(kw.get("p", 24))
        )
    raise ValueError("unknown board %s; choose from %s" % (board, BOARDS))


def load_diabetes_stream(
    n_per=240,
    n_batches=9,
    seed=SEED,
    target="num_medications",
    zip_path=DIABETES_ZIP,
):
    """Source→target Diabetes Readmission as ordered batches (covariate shift)."""
    zip_path = Path(zip_path)
    if not zip_path.is_file():
        raise FileNotFoundError("missing %s" % zip_path)

    def _read(name):
        with zipfile.ZipFile(zip_path) as zf:
            raw = zf.read(name).decode("utf-8", errors="replace")
        reader = csv.DictReader(io.StringIO(raw))
        rows = list(reader)
        cols = [c for c in reader.fieldnames if c != target]
        X = np.zeros((len(rows), len(cols)), dtype=float)
        y = np.zeros(len(rows), dtype=float)
        for i, row in enumerate(rows):
            y[i] = float(row[target])
            for j, c in enumerate(cols):
                try:
                    X[i, j] = float(row[c])
                except (TypeError, ValueError):
                    X[i, j] = 0.0
        return X, y

    Xs, ys = _read("source_DiabetesReadmission.csv")
    Xt, yt = _read("target_DiabetesReadmission.csv")
    rng = np.random.default_rng(int(seed))
    n_src_b = max(n_batches // 2, 1)
    n_tgt_b = max(n_batches - n_src_b, 1)
    parts_X, parts_y, parts_b = [], [], []
    for b in range(n_src_b):
        idx = rng.choice(len(ys), size=min(n_per, len(ys)), replace=False)
        parts_X.append(Xs[idx])
        parts_y.append(ys[idx])
        parts_b.append(np.full(len(idx), b, dtype=int))
    for j in range(n_tgt_b):
        b = n_src_b + j
        idx = rng.choice(len(yt), size=min(n_per, len(yt)), replace=False)
        parts_X.append(Xt[idx])
        parts_y.append(yt[idx])
        parts_b.append(np.full(len(idx), b, dtype=int))
    X = np.vstack(parts_X)
    y = np.concatenate(parts_y)
    batch = np.concatenate(parts_b)
    mu, sd = X.mean(0), X.std(0) + 1e-6
    X = (X - mu) / sd
    # standardize Y so MSE is comparable across boards (raw count scale is ~200)
    y = (y - float(y.mean())) / (float(y.std()) + 1e-6)
    return _as_stream(
        X,
        y,
        batch,
        meta={
            "board": "diabetes",
            "target": target,
            "y_standardized": True,
            "n_source_batches": int(n_src_b),
            "n_target_batches": int(n_tgt_b),
            "p": int(X.shape[1]),
            "n": int(X.shape[0]),
        },
    )


def make_interstate_like_stream(n_batches=9, n_per=240, seed=SEED, p=24):
    """Synthetic traffic-speed board (Interstate / PeMS surrogate)."""
    rng = np.random.default_rng(int(seed))
    regimes = np.array([0.0, 0.35, -0.25, 0.55, -0.4, 0.2, -0.15, 0.45, -0.3])
    X_parts, y_parts, b_parts = [], [], []
    w_true = rng.normal(size=p)
    w_true /= np.linalg.norm(w_true) + 1e-12
    for t in range(int(n_batches)):
        shift = regimes[t % len(regimes)]
        base = rng.normal(size=(n_per, p)) + shift
        base[:, 0] += 0.15 * np.sin(0.7 * t)
        noise = 0.35 * rng.normal(size=n_per)
        ww = w_true * (1.0 + (0.2 if t % 3 == 2 else 0.0))
        y = base @ ww + noise + 1.2 * shift
        X_parts.append(base)
        y_parts.append(y)
        b_parts.append(np.full(n_per, t, dtype=int))
    return _as_stream(
        np.vstack(X_parts),
        np.concatenate(y_parts),
        np.concatenate(b_parts),
        meta={"board": "interstate", "surrogate": True, "p": p, "n_per": n_per},
    )


def make_stock_like_stream(n_batches=9, n_per=240, seed=SEED, p=24):
    """Synthetic return board (equity-panel surrogate)."""
    rng = np.random.default_rng(int(seed) + 17)
    vols = np.array([0.6, 0.6, 1.4, 1.4, 0.5, 1.8, 0.7, 1.2, 0.9])
    X_parts, y_parts, b_parts = [], [], []
    beta = rng.normal(size=p)
    beta /= np.linalg.norm(beta) + 1e-12
    for t in range(int(n_batches)):
        vol = float(vols[t % len(vols)])
        factors = rng.normal(scale=vol, size=(n_per, p))
        if t % 3 == 0:
            factors += rng.normal(scale=0.5, size=p)
        eps = rng.normal(scale=0.25 * vol, size=n_per)
        rot = np.cos(0.4 * t) * beta + np.sin(0.4 * t) * rng.normal(size=p)
        rot /= np.linalg.norm(rot) + 1e-12
        y = factors @ rot + eps
        X_parts.append(factors)
        y_parts.append(y)
        b_parts.append(np.full(n_per, t, dtype=int))
    return _as_stream(
        np.vstack(X_parts),
        np.concatenate(y_parts),
        np.concatenate(b_parts),
        meta={"board": "stock", "surrogate": True, "p": p, "n_per": n_per},
    )


def mahal_sample_weights(X, batch, t, gamma=1.0):
    """Near last-batch mean (Mahalanobis) → larger weight; far / OOD → smaller."""
    batch = np.asarray(batch, dtype=int)
    X = np.asarray(X, dtype=float)
    w = np.ones(batch.shape[0], dtype=float)
    ref_idx = batch == (t - 1)
    if not np.any(ref_idx) or t < 1:
        return w
    ref = X[ref_idx]
    mu = ref.mean(0)
    xc = ref - mu
    cov = (xc.T @ xc) / max(len(ref) - 1, 1) + 1e-3 * np.eye(X.shape[1])
    try:
        prec = np.linalg.inv(cov)
    except np.linalg.LinAlgError:
        prec = np.linalg.pinv(cov)
    for s in range(t):
        idx = batch == s
        if not np.any(idx):
            continue
        delta = X[idx].mean(0) - mu
        d2 = float(delta @ prec @ delta)
        w[idx] = float(np.exp(-float(gamma) * np.log1p(max(d2, 0.0))))
    return np.maximum(w, 1e-3)


def _pack(method, online, history, stream, extras=None):
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


def run_mahal_ridge(stream, alpha=3.0, gamma=1.0):
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    k = int(batch.max()) + 1
    history = []
    for t in range(1, k):
        tr, ic = batch < t, batch == t
        w = mahal_sample_weights(X, batch, t, gamma=gamma)
        pred = _ridge_predict(X[tr], y[tr], X[ic], sample_weight=w[tr], alpha=alpha)
        history.append(
            {
                "round": int(t),
                "online_mse": _mse(pred, y[ic]),
                "w_mean": float(w[tr].mean()),
                "n": int(ic.sum()),
            }
        )
    online = np.array([h["online_mse"] for h in history], dtype=float)
    return _pack("mahal", online, history, stream, extras={"gamma": gamma, "alpha": alpha})


def run_weight_method(stream, method, **kw):
    method = str(method).lower()
    if method in ("uniform", "ridge_past"):
        rec = run_ridge_past(stream, alpha=kw.get("alpha", 3.0))
        rec["method"] = "uniform"
        return rec
    if method == "hop":
        rec = run_hop_ridge(stream, alpha=kw.get("alpha", 3.0), gamma=kw.get("gamma", 4.0))
        rec["method"] = "hop"
        return rec
    if method == "dga":
        rec = run_dga_ridge(
            stream,
            alpha=kw.get("alpha", 3.0),
            eta=kw.get("eta", 1.0),
            ema_beta=kw.get("ema_beta", 0.35),
            align=kw.get("align", "cosine"),
        )
        rec["method"] = "dga"
        return rec
    if method == "mahal":
        return run_mahal_ridge(
            stream, alpha=kw.get("alpha", 3.0), gamma=kw.get("mahal_gamma", 1.0)
        )
    if method == "attr":
        rec = run_attr_adapter(stream)
        rec["method"] = "attr"
        return rec
    raise ValueError("unknown weight method %s" % method)


def backward_transfer_mse(stream, method="hop", alpha=3.0, **kw):
    """MSE on batch 0 after fitting on all past with the method's weight rule."""
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    k = int(batch.max()) + 1
    if k < 2:
        return float("nan")
    tr = batch < k
    i0 = batch == 0
    t_ref = max(k - 1, 1)
    if method in ("uniform", "ridge_past", "attr"):
        pred = _ridge_predict(X[tr], y[tr], X[i0], sample_weight=None, alpha=alpha)
    elif method == "hop":
        stats = rolling_hop_stats(stream)
        w = hop_sample_weights(stats["mus"], batch, t=t_ref, gamma=kw.get("gamma", 4.0))
        pred = _ridge_predict(X[tr], y[tr], X[i0], sample_weight=w[tr], alpha=alpha)
    elif method == "mahal":
        w = mahal_sample_weights(X, batch, t=t_ref, gamma=kw.get("mahal_gamma", 1.0))
        pred = _ridge_predict(X[tr], y[tr], X[i0], sample_weight=w[tr], alpha=alpha)
    elif method == "dga":
        domains = list(range(k))
        spe = k - 1
        clf = _ridge_fit(X[tr], y[tr], sample_weight=None, alpha=alpha)
        a = dga_alignments(
            clf, X[tr], y[tr], batch[tr], domains, spe, align=kw.get("align", "cosine")
        )
        alpha_w = dga_mirror_step(
            np.ones(k, dtype=float) / float(k), a, eta=kw.get("eta", 1.0)
        )
        w = dga_sample_weights(alpha_w, batch, domains)
        pred = _ridge_predict(X[tr], y[tr], X[i0], sample_weight=w[tr], alpha=alpha)
    else:
        pred = _ridge_predict(X[tr], y[tr], X[i0], sample_weight=None, alpha=alpha)
    return float(_mse(pred, y[i0]))


def attach_regret(rec, baseline_cum):
    out = dict(rec)
    out["regret"] = float(rec["cum_mse"] - baseline_cum)
    out["baseline_cum_mse"] = float(baseline_cum)
    return out


def run_board_suite(
    board="amazon",
    seeds=None,
    methods=None,
    n_per=240,
    n_batches=9,
    **board_kw,
):
    seeds = list(seeds if seeds is not None else range(SEED, SEED + 4))
    methods = list(methods or WEIGHT_METHODS)
    rows = []
    traces = {m: [] for m in methods}
    shapes = None
    for s in seeds:
        stream = load_board(board, n_per=n_per, n_batches=n_batches, seed=int(s), **board_kw)
        if shapes is None:
            shapes = {
                "board": board,
                "n": int(stream.X.shape[0]),
                "p": int(stream.X.shape[1]),
                "n_batches": int(stream.batch.max()) + 1,
                "meta": dict(getattr(stream, "meta", {}) or {}),
            }
        base = run_weight_method(stream, "uniform")
        base_cum = float(base["cum_mse"])
        per_seed = {}
        for method in methods:
            rec = run_weight_method(stream, method)
            rec = attach_regret(rec, base_cum)
            rec["bwt"] = backward_transfer_mse(stream, method=method)
            per_seed[method] = rec
            traces[method].append(rec)
        best_cum = min(per_seed[m]["cum_mse"] for m in methods)
        for method in methods:
            rec = per_seed[method]
            rec["oracle_gap"] = float(rec["cum_mse"] - best_cum)
            rows.append(
                {
                    "seed": int(s),
                    "board": board,
                    "method": method,
                    "online_mse": rec["online_mse"],
                    "cum_mse": rec["cum_mse"],
                    "regret": rec["regret"],
                    "bwt": rec["bwt"],
                    "oracle_gap": rec["oracle_gap"],
                    "last_online_mse": rec["last_online_mse"],
                }
            )

    table = {}
    for method in methods:
        def col(key, method=method):
            v = np.array([r[key] for r in rows if r["method"] == method], dtype=float)
            return {
                "mean": float(np.nanmean(v)),
                "sd": float(np.nanstd(v, ddof=1) if len(v) > 1 else 0.0),
            }

        table[method] = {
            "online_mse": col("online_mse"),
            "cum_mse": col("cum_mse"),
            "regret": col("regret"),
            "bwt": col("bwt"),
            "oracle_gap": col("oracle_gap"),
            "n": int(sum(1 for r in rows if r["method"] == method)),
        }
    return {
        "board": board,
        "table": table,
        "rows": rows,
        "traces": traces,
        "shapes": shapes,
        "seeds": [int(s) for s in seeds],
        "methods": methods,
    }


def run_multi_board(
    boards=None,
    seeds=None,
    methods=None,
    n_per=240,
    n_batches=9,
    quick=False,
):
    boards = list(boards or ("amazon", "synthetic", "diabetes", "interstate", "stock"))
    methods = list(methods or WEIGHT_METHODS)
    if quick:
        n_per = min(int(n_per), 80)
        seeds = list(seeds if seeds is not None else range(SEED, SEED + 2))
        n_batches = min(int(n_batches), 6)
    else:
        seeds = list(seeds if seeds is not None else range(SEED, SEED + 4))
    suites = {}
    for board in boards:
        print("board", board, flush=True)
        suites[board] = run_board_suite(
            board=board,
            seeds=seeds,
            methods=methods,
            n_per=n_per,
            n_batches=n_batches,
        )
    return {
        "boards": boards,
        "suites": suites,
        "seeds": [int(s) for s in seeds],
        "methods": methods,
    }


def plot_ood_eval(multi, path):
    import matplotlib.pyplot as plt

    suites = multi["suites"]
    boards = multi["boards"]
    methods = multi["methods"]
    colors = {
        "uniform": "#6B7C8A",
        "hop": "#C45C26",
        "dga": "#5B3A8C",
        "mahal": "#2F6B4F",
        "attr": "#2C4A6E",
    }
    labels = {
        "uniform": "uniform Ridge",
        "hop": "hop cosine IW",
        "dga": "DGA align",
        "mahal": "Mahalanobis IW",
        "attr": "attr adapter",
    }
    n_b = len(boards)
    fig, axes = plt.subplots(2, n_b, figsize=(3.15 * n_b, 6.0), squeeze=False)
    for j, board in enumerate(boards):
        suite = suites[board]
        ax = axes[0, j]
        means = [suite["table"][m]["online_mse"]["mean"] for m in methods]
        sds = [suite["table"][m]["online_mse"]["sd"] for m in methods]
        ypos = np.arange(len(methods))[::-1]
        ax.barh(
            ypos,
            means,
            xerr=sds,
            color=[colors.get(m, "#333") for m in methods],
            height=0.62,
            error_kw={"ecolor": "#555", "lw": 0.9, "capsize": 2},
        )
        ax.set_yticks(ypos)
        ax.set_yticklabels([labels.get(m, m) for m in methods], fontsize=8)
        ax.set_xlabel("online MSE")
        ax.set_title(board, loc="left", fontsize=11, fontweight="bold")
        ax.grid(True, axis="x", alpha=0.35)

        ax = axes[1, j]
        means = [suite["table"][m]["regret"]["mean"] for m in methods]
        sds = [suite["table"][m]["regret"]["sd"] for m in methods]
        ax.barh(
            ypos,
            means,
            xerr=sds,
            color=[colors.get(m, "#333") for m in methods],
            height=0.62,
            error_kw={"ecolor": "#555", "lw": 0.9, "capsize": 2},
        )
        ax.axvline(0.0, color="#444", lw=0.8, ls="--")
        ax.set_yticks(ypos)
        ax.set_yticklabels([labels.get(m, m) for m in methods], fontsize=8)
        ax.set_xlabel("regret vs uniform")
        ax.grid(True, axis="x", alpha=0.35)
    fig.suptitle(
        "OOD-risk as weights — online MSE (top) · regret vs uniform (bottom)",
        fontsize=12,
        fontweight="bold",
        y=0.995,
    )
    fig.tight_layout(rect=[0, 0.02, 1, 0.96])
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_paths(multi, path, board="amazon"):
    import matplotlib.pyplot as plt

    suite = multi["suites"][board]
    methods = multi["methods"]
    colors = {
        "uniform": "#6B7C8A",
        "hop": "#C45C26",
        "dga": "#5B3A8C",
        "mahal": "#2F6B4F",
        "attr": "#2C4A6E",
    }
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    for method in methods:
        paths = np.array(
            [r["online_path"] for r in suite["traces"][method]], dtype=float
        )
        if paths.size == 0:
            continue
        t = np.arange(1, paths.shape[1] + 1)
        mu, sd = paths.mean(0), paths.std(0)
        ax.plot(t, mu, color=colors.get(method, "#333"), lw=2.0, label=method)
        ax.fill_between(
            t, mu - sd, mu + sd, color=colors.get(method, "#333"), alpha=0.12, lw=0
        )
    ax.set_xlabel("batch hop")
    ax.set_ylabel("online MSE")
    ax.set_title("Online MSE path — %s" % board, loc="left", fontsize=12, fontweight="bold")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(True, alpha=0.35)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path


def write_readme(multi, path):
    lines = [
        "# OOD-risk as weights — online quantification",
        "",
        "Compare **OOD / proximity / gradient-alignment weights** as sample weights",
        "for closed-form Ridge on streaming boards.",
        "",
        "Metrics: **online MSE**, **cum regret vs uniform**, **BWT** (batch-0 MSE),",
        "**oracle gap** (vs best method on that seed).",
        "",
        "> Waymo / Interstate / raw stock feeds are **not** in this repo.",
        "> `interstate` and `stock` are **regime-shift surrogates**;",
        "> `diabetes` is the real source→target table from `datasets/datasets.zip`.",
        "",
        "## Methods",
        "",
        "| id | weight rule |",
        "| --- | --- |",
        "| uniform | equal weights on all past (`ridge_past`) |",
        "| hop | heatmap cosine IW |",
        "| dga | Fan–Grangier–Ablin gradient alignment + EMA |",
        "| mahal | Mahalanobis proximity to last-batch mean |",
        "| attr | TSS-gated bank⊕hop adapter |",
        "",
        "## Results",
        "",
    ]
    for board in multi["boards"]:
        suite = multi["suites"][board]
        lines.append("### `%s`" % board)
        lines.append("")
        lines.append("| method | online MSE | regret vs uniform | BWT | oracle gap |")
        lines.append("| --- | --- | --- | --- | --- |")
        for m in multi["methods"]:
            t = suite["table"][m]
            lines.append(
                "| %s | %.4f (%.4f) | %+.4f (%.4f) | %.4f (%.4f) | %.4f (%.4f) |"
                % (
                    m,
                    t["online_mse"]["mean"],
                    t["online_mse"]["sd"],
                    t["regret"]["mean"],
                    t["regret"]["sd"],
                    t["bwt"]["mean"],
                    t["bwt"]["sd"],
                    t["oracle_gap"]["mean"],
                    t["oracle_gap"]["sd"],
                )
            )
        lines.append("")
    lines.extend(
        [
            "## Run",
            "",
            "```bash",
            "python3 scripts/run_ood_weight_online_eval.py",
            "python3 scripts/run_ood_weight_online_eval.py --quick",
            "python3 scripts/run_ood_weight_online_eval.py --boards amazon,diabetes,interstate",
            "```",
            "",
            "Code: `Python/src/ood_weight_online_eval.py`",
            "Skill: `.cursor/skills/ood-risk-weights-online/SKILL.md`",
            "",
        ]
    )
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))
    return path


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    def _c(o):
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(type(o))

    # drop heavy history from traces for the summary json
    slim = {
        "boards": obj["boards"],
        "seeds": obj["seeds"],
        "methods": obj["methods"],
        "suites": {},
    }
    for board, suite in obj["suites"].items():
        slim["suites"][board] = {
            "board": suite["board"],
            "table": suite["table"],
            "rows": suite["rows"],
            "shapes": suite["shapes"],
            "seeds": suite["seeds"],
            "methods": suite["methods"],
            "paths": {
                m: [r["online_path"] for r in suite["traces"][m]] for m in suite["methods"]
            },
        }
    path.write_text(json.dumps(slim, indent=2, default=_c))
    return path
