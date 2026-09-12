"""Amazon review continuous batches: scalar TSS vs LR clocks on rating MSE.

Nine product categories from Amazon Reviews 2023 are consecutive batches.
X is TF-IDF of (title + body), frozen on batch 0. Y is the 1-5 star rating.
The typed map is the same next-epoch TSS rule as the multimodal probe, on a
single text head:

    η* = η0 (1 + β δ / s) / (1 + λ c s),   s = sqrt(n_iter)

c is |d| of the TF-IDF coordinate-mean; δ is two-fold excess ridge MSE after
mean-aligning X. Quiet holds last η. Train batch t with the previous η, then
set η for t+1.

Comparators: constant, cosine, plateau, Polyak (loss-proportional, wrong
covariate sign), inv-c (covariate half of TSS).

Metrics: online MSE (predict the arriving batch, then update), cumulative
MSE, BWT (MSE on batch 0 after the stream), and regret vs the hindsight-best
constant η on the same stream.
"""
from __future__ import annotations

import json
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy import stats
from sklearn.feature_extraction.text import TfidfVectorizer

from msrvtt_multimodal_attribution import SEED, _cohens_d, write_json
from typed_shift_stepsize import (
    BETA,
    ETA0,
    ETA_MAX_MULT,
    LAMBDA,
    METHOD_COLORS,
    METHOD_LABELS,
    RHO_DOWN,
    RHO_UP,
    TAU_C,
    TAU_D,
    _cosine_eta,
    tss_lr,
)

CATEGORIES = (
    "Gift_Cards",
    "Digital_Music",
    "All_Beauty",
    "Software",
    "Subscription_Boxes",
    "Musical_Instruments",
    "Magazine_Subscriptions",
    "Handmade_Products",
    "Amazon_Fashion",
)
HF_JSONL = (
    "https://huggingface.co/datasets/McAuley-Lab/Amazon-Reviews-2023/"
    "resolve/main/raw/review_categories/{name}.jsonl"
)
METHODS = ("constant", "cosine", "plateau", "polyak", "inv_c", "tss")
TIME_METHODS = ("constant", "cosine", "plateau")
HINDSIGHT_ETAS = (0.04, 0.06, 0.08, 0.10, 0.16, 0.20, 0.25)
CACHE_DIR = Path(__file__).resolve().parents[2] / "data" / "amazon" / "heads"


@dataclass
class AmazonStream:
    X: np.ndarray
    y: np.ndarray
    batch: np.ndarray
    texts: list | None = None
    categories: tuple = ()
    meta: dict = field(default_factory=dict)


class RidgeSGD:
    """Linear MSE probe with an explicit stepsize."""

    def __init__(self, p, seed=SEED):
        rng = np.random.default_rng(seed)
        self.w = rng.normal(scale=0.01, size=int(p))
        self.b = 0.0

    def predict(self, X):
        X = np.asarray(X, dtype=float)
        return X @ self.w + self.b

    def mse(self, X, y):
        y = np.asarray(y, dtype=float)
        err = self.predict(X) - y
        return float(np.mean(err**2))

    def step(self, X, y, eta, ridge=1e-3):
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float)
        n = max(len(y), 1)
        err = self.predict(X) - y
        self.w -= float(eta) * ((X.T @ err) / n + float(ridge) * self.w)
        self.b -= float(eta) * float(err.mean())


def tss_eta_scalar(c, delta, eta0=ETA0, prev=None, n_iter=None, beta=BETA, lam=LAMBDA):
    """Same TSS map as the multimodal probe, read off the text head."""
    prev = float(eta0 if prev is None else prev)
    lrs = tss_lr(
        {"video": 0.0, "audio": 0.0, "text": float(c)},
        {"video": 0.0, "audio": 0.0, "text": float(delta)},
        eta0=float(eta0),
        prev={"video": prev, "audio": prev, "text": prev},
        n_iter=n_iter,
        beta=beta,
        lam=lam,
        rho_up=RHO_UP,
        rho_down=RHO_DOWN,
    )
    return float(lrs["text"])


def covariate_intensity_text(X0, X1):
    """|d| of the TF-IDF coordinate-mean, matching the TSS block estimator."""
    z0 = np.asarray(X0, dtype=float).mean(axis=1)
    z1 = np.asarray(X1, dtype=float).mean(axis=1)
    return float(abs(_cohens_d(z0, z1)))


def _corr_screen(X, y, k=16):
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    y = y - y.mean()
    Xc = X - X.mean(axis=0)
    num = Xc.T @ y
    den = np.sqrt((Xc**2).sum(axis=0) * float(y @ y) + 1e-12)
    score = np.abs(num / den)
    k = int(max(1, min(k, X.shape[1], max(X.shape[0] - 1, 1))))
    return np.argpartition(score, -k)[-k:]


def _ridge_predict(Xtr, ytr, Xte, cols, ridge=1.0):
    Z = np.c_[np.asarray(Xtr, dtype=float)[:, cols], np.ones(len(Xtr))]
    gram = Z.T @ Z + float(ridge) * np.eye(Z.shape[1])
    b = np.linalg.solve(gram, Z.T @ np.asarray(ytr, dtype=float))
    Zte = np.c_[np.asarray(Xte, dtype=float)[:, cols], np.ones(len(Xte))]
    return Zte @ b


def concept_intensity_mse(X0, y0, X1, y1, k=8, ridge=2.0):
    """Two-fold excess ridge MSE after mean-aligning X."""
    X0 = np.asarray(X0, dtype=float)
    X1 = np.asarray(X1, dtype=float)
    y0 = np.asarray(y0, dtype=float)
    y1 = np.asarray(y1, dtype=float)
    n0 = len(y0)
    if n0 < 8 or len(y1) < 4:
        return 0.0, {"mse0": float("nan"), "mse_al": float("nan")}
    rng = np.random.default_rng(n0 * 13 + 7 * len(y1) + int(np.round(y0[:1].sum() if n0 else 0)))
    fold = (rng.random(n0) >= 0.5).astype(int)
    if fold.min() == fold.max():
        fold[0] = 1 - fold[0]
    mu0, mu1 = X0.mean(axis=0), X1.mean(axis=0)
    X1_al = X1 - (mu1 - mu0)
    mse0s, mse_als = [], []
    for f in (0, 1):
        tr, te = np.flatnonzero(fold == f), np.flatnonzero(fold != f)
        if tr.size < 8 or te.size < 4:
            continue
        k_use = int(max(2, min(k, tr.size // 6, X0.shape[1])))
        cols = _corr_screen(X0[tr], y0[tr], k=k_use)
        pred_te = _ridge_predict(X0[tr], y0[tr], X0[te], cols, ridge=ridge)
        pred_al = _ridge_predict(X0[tr], y0[tr], X1_al, cols, ridge=ridge)
        mse0s.append(float(np.mean((pred_te - y0[te]) ** 2)))
        mse_als.append(float(np.mean((pred_al - y1) ** 2)))
    mse0 = float(np.mean(mse0s)) if mse0s else float("nan")
    mse_al = float(np.mean(mse_als)) if mse_als else float("nan")
    delta = max(0.0, mse_al - mse0) if np.isfinite(mse0) and np.isfinite(mse_al) else 0.0
    return float(delta), {"mse0": mse0, "mse_al": mse_al}


def make_amazon_like_stream(
    n_batches=9,
    n_per=80,
    p=48,
    seed=SEED,
    cov=0.35,
    concept=0.0,
    concept_at=5,
    noise=0.45,
    rank=6,
    signal=0.90,
):
    """Gaussian Amazon-like rating stream with known covariate / concept."""
    rng = np.random.default_rng(seed)
    beta = np.zeros(p)
    beta[: int(rank)] = rng.normal(scale=float(signal), size=int(rank))
    rows, ys, batches = [], [], []
    true_c, true_delta = [], []
    for t in range(int(n_batches)):
        X = rng.normal(size=(n_per, p))
        if cov:
            X += float(cov) * float(t)
        use = beta
        drifted = bool(concept) and t >= int(concept_at)
        if drifted:
            # Flip the rating map: P(Y|X) changes, P(X) does not.
            use = -beta
        y = 3.0 + X @ use + rng.normal(scale=noise, size=n_per)
        y = np.clip(y, 1.0, 5.0)
        rows.append(X)
        ys.append(y)
        batches.append(np.full(n_per, t, dtype=int))
        true_c.append(1.0 if cov else 0.0)
        true_delta.append(0.4 if drifted else 0.0)
    cats = tuple(f"synth_{t}" for t in range(n_batches))
    return AmazonStream(
        X=np.vstack(rows),
        y=np.concatenate(ys),
        batch=np.concatenate(batches),
        categories=cats,
        meta={
            "n_batches": int(n_batches),
            "n_per": int(n_per),
            "p": int(p),
            "cov": float(cov),
            "concept": float(concept),
            "concept_at": int(concept_at),
            "seed": int(seed),
            "source": "synthetic",
            "true_c": np.asarray(true_c, dtype=float),
            "true_delta": np.asarray(true_delta, dtype=float),
        },
    )


def _review_text(rec):
    title = str(rec.get("title") or "").strip()
    body = str(rec.get("text") or "").strip()
    return (title + " " + body).strip()


def parse_jsonl_records(text, max_rows=None):
    rows = []
    for line in str(text).splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue
        body = _review_text(rec)
        rating = rec.get("rating")
        if not body or rating is None:
            continue
        try:
            y = float(rating)
        except (TypeError, ValueError):
            continue
        if not np.isfinite(y):
            continue
        rows.append({"text": body, "y": y, "timestamp": rec.get("timestamp")})
        if max_rows is not None and len(rows) >= int(max_rows):
            break
    return rows


def fetch_category_head(
    name,
    cache_dir=CACHE_DIR,
    max_bytes=2_000_000,
    max_rows=700,
    timeout=120,
):
    """HTTP range-read the start of a McAuley category JSONL; cache locally."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    dest = cache_dir / f"{name}.jsonl"
    if dest.exists() and dest.stat().st_size > 200:
        rows = parse_jsonl_records(dest.read_text(encoding="utf-8"), max_rows=max_rows)
        if rows:
            return rows
    url = HF_JSONL.format(name=name)
    req = urllib.request.Request(
        url,
        headers={
            "Range": "bytes=0-%d" % (int(max_bytes) - 1),
            "User-Agent": "cfperm-amazon-batches/0.1",
        },
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        raw = resp.read()
    rows = parse_jsonl_records(raw.decode("utf-8", errors="replace"), max_rows=max_rows)
    slim = [
        json.dumps({"rating": r["y"], "title": "", "text": r["text"], "timestamp": r["timestamp"]}, ensure_ascii=False)
        for r in rows
    ]
    dest.write_text("\n".join(slim) + ("\n" if slim else ""), encoding="utf-8")
    return rows


def load_amazon_reviews(
    categories=CATEGORIES,
    n_per=240,
    seed=SEED,
    cache_dir=CACHE_DIR,
    max_rows=700,
):
    """Sample n_per reviews from each of 8-10 product categories."""
    rng = np.random.default_rng(seed)
    cats = tuple(categories)
    texts, ys, batches, used = [], [], [], []
    pool_sizes = {}
    for t, name in enumerate(cats):
        rows = fetch_category_head(name, cache_dir=cache_dir, max_rows=max_rows)
        if len(rows) < 20:
            raise RuntimeError("category %s has too few reviews (%d)" % (name, len(rows)))
        idx = np.arange(len(rows))
        rng.shuffle(idx)
        take = idx[: min(int(n_per), len(idx))]
        pool_sizes[name] = int(len(rows))
        used.append(name)
        for i in take:
            texts.append(rows[i]["text"])
            ys.append(rows[i]["y"])
            batches.append(t)
    return {
        "texts": texts,
        "y": np.asarray(ys, dtype=float),
        "batch": np.asarray(batches, dtype=int),
        "categories": tuple(used),
        "meta": {
            "source": "amazon-reviews-2023",
            "n_batches": len(used),
            "n_per": int(n_per),
            "pool_sizes": pool_sizes,
            "seed": int(seed),
        },
    }


def featurize_reviews(payload, max_features=128, min_df=2):
    """TF-IDF frozen on batch 0. Rows stay l2-normalized; do not column-scale.

    Column-standardizing rare words inflates the gradient so η0=0.10 (the
    multimodal default) diverges. The TSS comparison is then just who decayed
    fast enough. Frozen TF-IDF keeps η0 on the same scale as the probe.
    """
    texts = list(payload["texts"])
    y = np.asarray(payload["y"], dtype=float)
    batch = np.asarray(payload["batch"], dtype=int)
    i0 = np.flatnonzero(batch == 0)
    vec = TfidfVectorizer(
        max_features=int(max_features),
        min_df=min_df,
        stop_words="english",
        sublinear_tf=True,
        ngram_range=(1, 1),
    )
    vec.fit([texts[i] for i in i0])
    X = np.asarray(vec.transform(texts).toarray(), dtype=float)
    meta = dict(payload.get("meta") or {})
    meta.update(
        {
            "n_features": int(X.shape[1]),
            "vocab": int(len(vec.get_feature_names_out())),
            "y_mean": float(y[i0].mean()),
            "y_std": float(y[i0].std()),
        }
    )
    return AmazonStream(
        X=X,
        y=y,
        batch=batch,
        texts=texts,
        categories=tuple(payload.get("categories") or ()),
        meta=meta,
    )


def _scheduler_eta(method, t, n_batches, eta0, c, delta, plateau_eta, prev, n_iter, online_mse=None):
    if method == "tss":
        return tss_eta_scalar(c, delta, eta0=eta0, prev=prev, n_iter=n_iter)
    if method == "inv_c":
        return tss_eta_scalar(c, 0.0, eta0=eta0, prev=prev, n_iter=n_iter, beta=0.0)
    if method == "polyak":
        mse = float(1.0 if online_mse is None else online_mse)
        if not np.isfinite(mse) or mse > 50.0:
            return float(eta0) * 0.25
        cap = float(eta0) * float(ETA_MAX_MULT)
        return float(min(cap, eta0 * (mse / 1.5)))
    if method == "constant":
        return float(eta0)
    if method == "cosine":
        return _cosine_eta(eta0, t, max(n_batches - 1, 1))
    if method == "plateau":
        return float(plateau_eta if plateau_eta is not None else eta0)
    raise ValueError("unknown method %s" % method)


def run_amazon_method(
    stream,
    method="tss",
    eta0=ETA0,
    steps_per_batch=8,
    warmup_steps=6,
    seed=SEED,
    ridge=1e-3,
    batch_size=32,
):
    """Warmup on B0. Typed maps set η for the next epoch; clocks use t."""
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    n_batches = int(batch.max()) + 1
    probe = RidgeSGD(X.shape[1], seed=seed)
    rng = np.random.default_rng(seed)
    i0 = np.flatnonzero(batch == 0)
    probe.b = float(y[i0].mean())
    eta = float(eta0)
    n_iter = 0
    chunk = min(int(batch_size), max(4, i0.size))
    for _ in range(int(max(1, warmup_steps))):
        sl = rng.choice(i0.size, size=chunk, replace=False)
        probe.step(X[i0][sl], y[i0][sl], eta, ridge=ridge)
        n_iter += 1

    history = [
        {
            "round": 0,
            "phase": "warmup",
            "method": method,
            "eta": float(eta),
            "c": 0.0,
            "delta": 0.0,
            "online_mse": float("nan"),
            "post_mse": probe.mse(X[i0], y[i0]),
            "bwt": probe.mse(X[i0], y[i0]),
            "n": int(i0.size),
        }
    ]
    plateau_eta = float(eta0)
    best_mse = np.inf
    stall = 0
    next_eta = float(eta0)
    c_hat, d_hat = 0.0, 0.0

    for t in range(1, n_batches):
        ip, ic = np.flatnonzero(batch == t - 1), np.flatnonzero(batch == t)
        c_hat = covariate_intensity_text(X[ip], X[ic])
        d_hat, _ = concept_intensity_mse(X[ip], y[ip], X[ic], y[ic])
        if method in TIME_METHODS:
            eta = _scheduler_eta(method, t, n_batches, eta0, c_hat, d_hat, plateau_eta, next_eta, n_iter)
        else:
            eta = float(next_eta)
        online = probe.mse(X[ic], y[ic])
        chunk_t = min(int(batch_size), max(4, ic.size))
        n_steps = max(1, int(steps_per_batch))
        for _ in range(n_steps):
            sl = rng.choice(ic.size, size=chunk_t, replace=False)
            probe.step(X[ic][sl], y[ic][sl], eta, ridge=ridge)
            n_iter += 1
        post = probe.mse(X[ic], y[ic])
        bwt = probe.mse(X[i0], y[i0])
        history.append(
            {
                "round": int(t),
                "phase": "adapt",
                "method": method,
                "eta": float(eta),
                "c": float(c_hat),
                "delta": float(d_hat),
                "online_mse": float(online),
                "null_mse": float(np.mean((y[ic] - y[i0].mean()) ** 2)),
                "post_mse": float(post),
                "bwt": float(bwt),
                "n": int(ic.size),
            }
        )
        if online < best_mse - 1e-4:
            best_mse = online
            stall = 0
        else:
            stall += 1
            if stall >= 2:
                plateau_eta *= 0.5
                stall = 0
        next_eta = _scheduler_eta(
            method,
            t + 1,
            n_batches,
            eta0,
            c_hat,
            d_hat,
            plateau_eta,
            eta,
            n_iter,
            online_mse=online,
        )

    adapt = [h for h in history if h["phase"] == "adapt"]
    online = np.array([h["online_mse"] for h in adapt], dtype=float)
    null = np.array([h.get("null_mse", np.nan) for h in adapt], dtype=float)
    cum = np.cumsum(online)
    summary = {
        "method": method,
        "n_batches": n_batches,
        "eta0": float(eta0),
        "online_mse": float(online.mean()) if online.size else float("nan"),
        "null_mse": float(null.mean()) if null.size else float("nan"),
        "cum_mse": float(cum[-1]) if cum.size else float("nan"),
        "last_online_mse": float(online[-1]) if online.size else float("nan"),
        "post_mse": float(np.mean([h["post_mse"] for h in adapt])) if adapt else float("nan"),
        "bwt": float(history[-1]["bwt"]) if history else float("nan"),
        "mean_eta": float(np.mean([h["eta"] for h in adapt])) if adapt else float("nan"),
        "mean_c": float(np.mean([h["c"] for h in adapt])) if adapt else float("nan"),
        "mean_delta": float(np.mean([h["delta"] for h in adapt])) if adapt else float("nan"),
        "online_path": online.tolist(),
        "null_path": null.tolist(),
        "cum_path": cum.tolist(),
        "history": history,
        "meta": dict(stream.meta),
    }
    return summary, probe


def hindsight_constant(stream, etas=HINDSIGHT_ETAS, **kwargs):
    """Best constant η in the TSS-feasible box [min etas, η0·ETA_MAX_MULT].

    Unconstrained larger constants can still cut online MSE on this weak
    linear probe; that is underfit, not a scheduler ranking. Cap the grid
    at the same 2.5×η0 used by TSS.
    """
    best = None
    traces = []
    for eta in etas:
        rec, _ = run_amazon_method(stream, method="constant", eta0=float(eta), **kwargs)
        rec = dict(rec)
        rec["hindsight_eta"] = float(eta)
        traces.append(rec)
        if best is None or rec["cum_mse"] < best["cum_mse"]:
            best = rec
    best = dict(best)
    best["grid"] = [{"eta": r["hindsight_eta"], "cum_mse": r["cum_mse"], "online_mse": r["online_mse"]} for r in traces]
    return best


def _attach_regret(summary, hindsight):
    out = dict(summary)
    out["regret"] = float(summary["cum_mse"] - hindsight["cum_mse"])
    out["hindsight_eta"] = float(hindsight.get("hindsight_eta", hindsight.get("eta0", float("nan"))))
    out["hindsight_cum_mse"] = float(hindsight["cum_mse"])
    return out


def run_amazon_suite(
    seeds=None,
    n_batches=9,
    n_per=240,
    methods=None,
    eta0=ETA0,
    source="amazon",
    categories=CATEGORIES,
    cache_dir=CACHE_DIR,
    steps_per_batch=8,
    cov=0.35,
    concept=0.0,
):
    """Monte Carlo over seeds. ``source='amazon'`` or ``'synthetic'``."""
    seeds = list(seeds if seeds is not None else range(SEED, SEED + 6))
    methods = list(methods or METHODS)
    rows = []
    traces = {m: [] for m in methods}
    hindsight_rows = []
    ident = []
    for s in seeds:
        if source == "synthetic":
            stream = make_amazon_like_stream(
                n_batches=n_batches, n_per=n_per, seed=int(s), cov=cov, concept=concept
            )
        else:
            raw = load_amazon_reviews(
                categories=categories[: int(n_batches)],
                n_per=n_per,
                seed=int(s),
                cache_dir=cache_dir,
            )
            stream = featurize_reviews(raw)
        hind = hindsight_constant(
            stream, steps_per_batch=steps_per_batch, seed=int(s)
        )
        hindsight_rows.append(
            {
                "seed": int(s),
                "hindsight_eta": hind["hindsight_eta"],
                "cum_mse": hind["cum_mse"],
                "online_mse": hind["online_mse"],
            }
        )
        for method in methods:
            summary, _ = run_amazon_method(
                stream,
                method=method,
                eta0=eta0,
                steps_per_batch=steps_per_batch,
                seed=int(s),
            )
            summary = _attach_regret(summary, hind)
            traces[method].append(summary)
            rows.append(
                {
                    "seed": int(s),
                    "method": method,
                    "online_mse": summary["online_mse"],
                    "cum_mse": summary["cum_mse"],
                    "regret": summary["regret"],
                    "bwt": summary["bwt"],
                    "last_online_mse": summary["last_online_mse"],
                    "mean_eta": summary["mean_eta"],
                    "mean_c": summary["mean_c"],
                    "mean_delta": summary["mean_delta"],
                    "hindsight_eta": summary["hindsight_eta"],
                }
            )
        ident.append(
            {
                "seed": int(s),
                "mean_c": traces[methods[0]][-1]["mean_c"],
                "mean_delta": traces[methods[0]][-1]["mean_delta"],
                "n_batches": traces[methods[0]][-1]["n_batches"],
                "n_features": int(stream.X.shape[1]),
                "categories": list(stream.categories),
                "source": stream.meta.get("source"),
            }
        )
    table = _summarize_rows(rows)
    tests = paired_tests(rows)
    return {
        "rows": rows,
        "table": table,
        "tests": tests,
        "hindsight": hindsight_rows,
        "identification": ident,
        "traces": traces,
        "methods": methods,
        "meta": {
            "source": source,
            "n_batches": int(n_batches),
            "n_per": int(n_per),
            "seeds": [int(s) for s in seeds],
            "eta0": float(eta0),
            "categories": list(categories[: int(n_batches)]),
        },
    }


def _summarize_rows(rows):
    methods = []
    for r in rows:
        if r["method"] not in methods:
            methods.append(r["method"])
    table = {}
    for method in methods:
        sub = [r for r in rows if r["method"] == method]
        cell = {}
        for key in ("online_mse", "cum_mse", "regret", "bwt", "last_online_mse", "mean_eta", "mean_c", "mean_delta"):
            v = np.array([r[key] for r in sub], dtype=float)
            cell[key] = {
                "mean": float(v.mean()),
                "sd": float(v.std(ddof=1)) if len(v) > 1 else 0.0,
                "n": int(len(v)),
            }
        table[method] = cell
    return table


def paired_tests(rows, baseline="tss"):
    out = {}
    methods = [m for m in dict.fromkeys(r["method"] for r in rows) if m != baseline]
    seeds = sorted({r["seed"] for r in rows})
    by = {(r["seed"], r["method"]): r for r in rows}
    for method in methods:
        rec = {}
        for metric in ("online_mse", "regret", "bwt"):
            a = np.array([by[(s, baseline)][metric] for s in seeds if (s, baseline) in by and (s, method) in by])
            b = np.array([by[(s, method)][metric] for s in seeds if (s, baseline) in by and (s, method) in by])
            if len(a) < 4 or np.allclose(a, b):
                rec[metric] = {
                    "n": int(len(a)),
                    "mean_diff": float(np.mean(a - b)) if len(a) else float("nan"),
                    "p": float("nan"),
                }
                continue
            try:
                stat, p = stats.wilcoxon(a, b, zero_method="wilcox", alternative="two-sided")
            except ValueError:
                stat, p = float("nan"), float("nan")
            rec[metric] = {
                "n": int(len(a)),
                "mean_diff": float(np.mean(a - b)),
                "stat": float(stat),
                "p": float(p),
            }
        out[method] = rec
    return out


DIVERGE_MSE = 50.0


def _stable_methods(table, methods):
    out = []
    for m in methods:
        mu = table[m]["online_mse"]["mean"]
        if np.isfinite(mu) and abs(mu) < DIVERGE_MSE:
            out.append(m)
    return out


def plot_amazon_suite(suite, path):
    import matplotlib.pyplot as plt
    from msrvtt_attribution_plots import GRID, INK, MUTED, _save, _style

    _style()
    methods = [m for m in METHODS if m in suite["table"]]
    bars = _stable_methods(suite["table"], methods)
    colors = {m: METHOD_COLORS.get(m, METHOD_COLORS.get("polyak_m" if m == "polyak" else m, "#9AA3AE")) for m in methods}
    labels = {m: METHOD_LABELS.get(m, m) for m in methods}
    labels["polyak"] = r"Polyak"
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.6))

    ax = axes[0, 0]
    x = np.arange(len(bars))
    means = [suite["table"][m]["online_mse"]["mean"] for m in bars]
    sds = [suite["table"][m]["online_mse"]["sd"] for m in bars]
    ax.bar(x, means, yerr=sds, color=[colors[m] for m in bars], ecolor=MUTED, capsize=2.5, width=0.78)
    ax.set_xticks(x)
    ax.set_xticklabels([labels[m] for m in bars], rotation=35, ha="right")
    ax.set_ylabel("online MSE")
    ax.set_title("Mean online MSE", loc="left", fontsize=12, fontweight="bold")
    ax.grid(True, axis="y", color=GRID)

    ax = axes[0, 1]
    means = [suite["table"][m]["regret"]["mean"] for m in bars]
    sds = [suite["table"][m]["regret"]["sd"] for m in bars]
    ax.bar(x, means, yerr=sds, color=[colors[m] for m in bars], ecolor=MUTED, capsize=2.5, width=0.78)
    ax.set_xticks(x)
    ax.set_xticklabels([labels[m] for m in bars], rotation=35, ha="right")
    ax.set_ylabel("cumulative regret")
    ax.set_title("Regret vs hindsight-best constant η", loc="left", fontsize=12, fontweight="bold")
    ax.grid(True, axis="y", color=GRID)

    ax = axes[1, 0]
    for method in bars:
        recs = suite["traces"][method]
        path_mse = np.array([r["online_path"] for r in recs], dtype=float)
        t = np.arange(1, path_mse.shape[1] + 1)
        mu, sd = path_mse.mean(axis=0), path_mse.std(axis=0)
        ax.plot(t, mu, color=colors[method], lw=2.0, label=labels[method])
        ax.fill_between(t, mu - sd, mu + sd, color=colors[method], alpha=0.12, lw=0)
    recs0 = suite["traces"][bars[0]] if bars else None
    if recs0 and recs0[0].get("null_path"):
        null = np.array([r["null_path"] for r in recs0], dtype=float)
        t = np.arange(1, null.shape[1] + 1)
        ax.plot(t, null.mean(axis=0), color=MUTED, ls="--", lw=1.2, label="intercept")
    ax.set_xlabel("batch")
    ax.set_ylabel("online MSE")
    ax.set_title("Per-batch online MSE", loc="left", fontsize=12, fontweight="bold")
    ax.grid(True, color=GRID)
    ax.legend(frameon=False, fontsize=7.5, ncol=2)

    ax = axes[1, 1]
    for method in methods:
        recs = suite["traces"][method]
        t = [h["round"] for h in recs[0]["history"]]
        eta = np.array([[h["eta"] for h in r["history"]] for r in recs], dtype=float)
        mu, sd = eta.mean(axis=0), eta.std(axis=0)
        ax.plot(t, mu, color=colors[method], lw=2.0, label=labels[method])
        ax.fill_between(t, mu - sd, mu + sd, color=colors[method], alpha=0.12, lw=0)
    ax.set_xlabel("round")
    ax.set_ylabel("η")
    ax.set_title("Next-epoch stepsize", loc="left", fontsize=12, fontweight="bold")
    ax.grid(True, color=GRID)

    n_b = suite["meta"].get("n_batches", "")
    fig.suptitle(
        "Amazon reviews, %s batches — TSS vs LR clocks" % n_b,
        fontsize=13.2,
        fontweight="bold",
        color=INK,
        x=0.04,
        ha="left",
    )
    cats = suite["meta"].get("categories") or []
    omitted = [m for m in methods if m not in bars]
    note = ""
    if omitted:
        note = "  Omitted from MSE panels (diverged): %s." % ", ".join(omitted)
    fig.text(
        0.04,
        0.01,
        "Amazon Reviews 2023, %d category batches.  Online MSE is predict-then-update.  "
        "Dashed line is the batch-0 intercept.  Regret vs hindsight-best constant η.%s  Error bars are seed s.d."
        % (len(cats) or int(n_b or 0), note),
        fontsize=7.8,
        color=MUTED,
    )
    fig.tight_layout(rect=(0, 0.06, 1, 0.94))
    return _save(fig, path)


def write_tex_table(suite, path):
    table = suite["table"]
    tests = suite["tests"]
    methods = [m for m in METHODS if m in table]
    labels = {m: METHOD_LABELS.get(m, m) for m in methods}
    labels["polyak"] = r"Polyak"
    meta = suite["meta"]
    cats = meta.get("categories") or []
    ident = suite.get("identification") or []
    mean_c = float(np.mean([r["mean_c"] for r in ident])) if ident else float("nan")
    mean_d = float(np.mean([r["mean_delta"] for r in ident])) if ident else float("nan")
    hind = suite.get("hindsight") or []
    hind_eta = float(np.mean([r["hindsight_eta"] for r in hind])) if hind else float("nan")
    lines = [
        r"% Amazon continuous batches: online MSE and regret. Auto-generated.",
        r"\begin{table}[ht]\centering",
        r"\caption{Amazon Reviews 2023, %d product-category batches. Linear TF-IDF MSE, predict-then-update."
        % int(meta.get("n_batches") or len(cats) or 0),
        r"TSS is the same next-epoch map as the multimodal probe, on a scalar text head.",
        r"Regret is cumulative online MSE minus the hindsight-best constant $\eta$ inside the TSS box $[0.04,2.5\eta_0]$.",
        r"Entries are mean (s.d.) over seeds.}",
        r"\label{tab:amazon-continuous-batches}",
        r"\small",
        r"\begin{tabular}{@{}l cccc c@{}}\toprule",
        r"Method & online MSE & regret & BWT MSE & $\bar\eta$ & $\bar c$ \\",
        r"\midrule",
    ]
    for method in methods:
        cell = table[method]

        def fmt(k):
            mu, sd = cell[k]["mean"], cell[k]["sd"]
            if not np.isfinite(mu) or abs(mu) >= DIVERGE_MSE:
                return r"diverged"
            return r"$%.3f$ ($%.3f$)" % (mu, sd)

        star = ""
        rec = tests.get(method, {}).get("regret", {})
        p = rec.get("p", float("nan"))
        diff = rec.get("mean_diff", float("nan"))
        # mean_diff is TSS − comparator; negative regret-diff means TSS has lower regret
        if method != "tss" and p == p and p < 0.05 and diff == diff and diff < 0:
            star = r"$^{\ast}$"
        lines.append(
            r"%s%s & %s & %s & %s & $%.4f$ & $%.3f$ \\"
            % (
                labels[method],
                star,
                fmt("online_mse"),
                fmt("regret"),
                fmt("bwt"),
                cell["mean_eta"]["mean"],
                cell["mean_c"]["mean"],
            )
        )
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}\\[0.4em]")
    cat_tex = ", ".join(c.replace("_", r"\_") for c in cats)
    lines.append(
        r"{\footnotesize Wilcoxon signed-rank on regret, TSS $-$ comparator: $^{\ast}$ $p<0.05$ and TSS lower."
        r" Mean $\hat c=%.3f$, $\hat\delta=%.3f$. Hindsight-best constant $\bar\eta=%.3f$. Categories: %s.}"
        % (mean_c, mean_d, hind_eta, cat_tex)
    )
    lines.append(r"\end{table}")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


def strip_traces(suite):
    out = {k: v for k, v in suite.items() if k != "traces"}
    slim = {}
    for method, recs in (suite.get("traces") or {}).items():
        slim[method] = [
            {k: v for k, v in r.items() if k != "history"}
            | {
                "last_eta": r["history"][-1]["eta"],
                "last_c": r["history"][-1]["c"],
                "last_delta": r["history"][-1]["delta"],
            }
            for r in recs
        ]
    out["traces_slim"] = slim
    return out
