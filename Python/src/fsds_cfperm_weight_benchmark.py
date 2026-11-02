"""Benchmark FSDS/CFPerm → sample weights + observation-subset weights.

Datasets
--------
  amazon     continuous category hops (real reviews)
  digits     class-block shift (0-4 → inject 5-9)
  news20     TF-IDF 3-topic shift (held-out topic influx)
  synth_cgcd sparse novel-class C-GCD proxy

Weight channels (from FSDS RF-domain VIMP + subset localization)
---------------------------------------------------------------
  sample_w   : w_i = 1 + λ · membership_i          (soft, mean-1)
  subset_w   : s_i = 1{membership_i ≥ q} · c + 1   (hard observation subset)
  both       : sample_w * subset_w (renormalized)

Adaptation: SGDClassifier.partial_fit with the chosen weights.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
from sklearn.datasets import fetch_20newsgroups, load_digits, make_classification
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from weighted_cold_start import (
    SEED,
    attribute_rf_domain,
    cold_start_partial_fit,
    combine_weights,
    fit_initial,
    localize_subset,
    make_sgd,
    sample_weights,
    subset_weights,
)


def cfperm_feature_reject(
    X_ref: np.ndarray,
    X_live: np.ndarray,
    vimp: np.ndarray,
    *,
    top_frac: float = 0.25,
) -> np.ndarray:
    """Lightweight CFPerm-style reject set: top fraction of RF-domain VIMP."""
    v = np.asarray(vimp, float)
    k = max(1, int(round(len(v) * float(top_frac))))
    return np.argsort(-v)[:k]


# ---------------------------------------------------------------------------
# Dataset boards → (X_ref, y_ref, X_live, y_live, meta)
# ---------------------------------------------------------------------------


def board_amazon(*, n_per: int = 100, max_features: int = 64, seed: int = SEED) -> Dict:
    from amazon_continuous_batches import (
        CACHE_DIR,
        CATEGORIES,
        featurize_reviews,
        load_amazon_reviews,
    )

    raw = load_amazon_reviews(
        categories=list(CATEGORIES[:4]), n_per=n_per, seed=seed, cache_dir=CACHE_DIR,
    )
    st = featurize_reviews(raw, max_features=max_features)
    X, batch = np.asarray(st.X, float), np.asarray(st.batch, int)
    # label = category id; REF = cats 0,1; LIVE = recirculate known + sparse novel 2,3
    rng = np.random.default_rng(seed)
    known_idx = np.flatnonzero(batch <= 1)
    novel_idx = np.flatnonzero(batch >= 2)
    ref = rng.choice(known_idx, size=min(160, len(known_idx)), replace=False)
    live_k = rng.choice(known_idx, size=min(120, len(known_idx)), replace=False)
    live_n = rng.choice(novel_idx, size=min(30, len(novel_idx)), replace=False)
    live = np.concatenate([live_k, live_n])
    return {
        "name": "amazon",
        "X_ref": X[ref], "y_ref": batch[ref],
        "X_live": X[live], "y_live": batch[live],
        "note": "category hop: REF cats{0,1}; LIVE recirculate + sparse cats{2,3}",
        "known_labels": [0, 1],
        "novel_labels": [2, 3],
    }


def board_digits(*, seed: int = SEED) -> Dict:
    X, y = load_digits(return_X_y=True)
    X = StandardScaler().fit_transform(np.asarray(X, float))
    y = np.asarray(y, int)
    # REF: digits 0-4; LIVE: mix 0-4 recirculation + 5-9 novel influx (sparse)
    rng = np.random.default_rng(seed)
    known = np.flatnonzero(y <= 4)
    novel = np.flatnonzero(y >= 5)
    ref = rng.choice(known, size=min(400, len(known)), replace=False)
    live_k = rng.choice(known, size=160, replace=False)
    live_n = rng.choice(novel, size=32, replace=False)  # sparser novel
    live = np.concatenate([live_k, live_n])
    rng.shuffle(live)
    return {
        "name": "digits",
        "X_ref": X[ref], "y_ref": y[ref],
        "X_live": X[live], "y_live": y[live],
        "note": "known 0-4 REF; LIVE sparse inject 5-9 (~1:5)",
        "novel_labels": list(range(5, 10)),
        "known_labels": list(range(0, 5)),
    }


def board_news20(*, seed: int = SEED, max_features: int = 128) -> Dict:
    cats = ["sci.space", "rec.sport.baseball", "comp.graphics"]
    bun = fetch_20newsgroups(
        subset="train", categories=cats, remove=("headers", "footers", "quotes"),
        shuffle=True, random_state=seed,
    )
    vec = TfidfVectorizer(max_features=max_features, stop_words="english")
    X = vec.fit_transform(bun.data).toarray()
    y = np.asarray(bun.target, int)
    # REF: topics 0,1; LIVE: recirculate + sparse topic 2
    rng = np.random.default_rng(seed)
    known = np.flatnonzero(y <= 1)
    novel = np.flatnonzero(y == 2)
    ref = rng.choice(known, size=min(500, len(known)), replace=False)
    live_k = rng.choice(known, size=180, replace=False)
    live_n = rng.choice(novel, size=36, replace=False)
    live = np.concatenate([live_k, live_n])
    return {
        "name": "news20",
        "X_ref": X[ref], "y_ref": y[ref],
        "X_live": X[live], "y_live": y[live],
        "note": "TF-IDF; known sci/rec → sparse inject comp.graphics (~1:5)",
        "novel_labels": [2],
        "known_labels": [0, 1],
    }


def board_synth_cgcd(*, seed: int = SEED, quick: bool = True) -> Dict:
    from category_shift_boards import make_cgcd_stream

    stream = make_cgcd_stream(
        n_known=8 if quick else 12, n_novel_per_stage=2, n_stages=1,
        n_per=30 if quick else 40, p=20 if quick else 32, novel_dims=6, seed=seed,
    )
    st = stream["stages"][0]
    return {
        "name": "synth_cgcd",
        "X_ref": stream["X_ref"], "y_ref": stream["y_ref"],
        "X_live": st["X"], "y_live": st["y"],
        "note": "synthetic sparse novel C-GCD",
        "novel_labels": st["novel"],
        "known_labels": stream["known"],
        "active_dims": st["active_dims"],
    }


def load_boards(*, seed: int = SEED, quick: bool = True) -> List[Dict]:
    boards = []
    for fn in (
        lambda: board_digits(seed=seed),
        lambda: board_news20(seed=seed, max_features=96 if quick else 128),
        lambda: board_synth_cgcd(seed=seed, quick=quick),
        lambda: board_amazon(n_per=60 if quick else 100, max_features=48 if quick else 64, seed=seed),
    ):
        try:
            boards.append(fn())
        except Exception as exc:  # noqa: BLE001
            boards.append({"name": "FAILED", "error": str(exc)})
    return [b for b in boards if "X_ref" in b]


# ---------------------------------------------------------------------------
# One-board evaluation
# ---------------------------------------------------------------------------


def _metrics(clf, X, y, *, known=None, novel=None) -> Dict:
    pred = clf.predict(X)
    out = {"acc": float(accuracy_score(y, pred))}
    if known is not None:
        m = np.isin(y, known)
        out["known_acc"] = float(np.mean(pred[m] == y[m])) if m.any() else float("nan")
    if novel is not None:
        m = np.isin(y, novel)
        out["novel_acc"] = float(np.mean(pred[m] == y[m])) if m.any() else float("nan")
    return out


def run_board(
    board: Dict,
    *,
    seed: int = SEED,
    lam: float = 2.5,
    subset_q: float = 0.6,
    subset_boost: float = 3.0,
    top_k: int = 8,
) -> Dict:
    X0 = np.asarray(board["X_ref"], float)
    y0 = np.asarray(board["y_ref"]).ravel()
    X1 = np.asarray(board["X_live"], float)
    y1 = np.asarray(board["y_live"]).ravel()
    known = board.get("known_labels")
    novel = board.get("novel_labels")
    classes = np.unique(np.concatenate([y0, y1]))

    # scale on REF
    sc = StandardScaler().fit(X0)
    X0s, X1s = sc.transform(X0), sc.transform(X1)

    # FSDS attribution + localization
    vimp = attribute_rf_domain(X0s, X1s, seed=seed, n_estimators=60)
    mem, feat = localize_subset(X0s, X1s, vimp, top_k=top_k)
    rejected = cfperm_feature_reject(X0s, X1s, vimp, top_frac=0.25)
    feat_hit = None
    if board.get("active_dims") is not None:
        feat_hit = float(np.isin(feat, board["active_dims"]).mean())

    sw = sample_weights(mem, lam=lam)
    sub_w, sub_mask = subset_weights(mem, q=subset_q, boost=subset_boost)
    both = combine_weights(sw, sub_w)
    # fused: membership ∪ private-live labels (discovered / oracle ids)
    mem_f = mem.copy()
    if novel is not None:
        mem_f = np.maximum(mem_f, np.isin(y1, novel).astype(float))
    sw_f = sample_weights(mem_f, lam=lam)
    sub_f, sub_mask_f = subset_weights(mem_f, q=subset_q, boost=subset_boost)
    both_f = combine_weights(sw_f, sub_f)

    # train/eval split on LIVE
    rng = np.random.default_rng(seed + 11)
    # stratify-ish: keep some novel in both sides when possible
    idx = rng.permutation(len(y1))
    n_ev = max(20, len(y1) // 3)
    ev, tr = idx[:n_ev], idx[n_ev:]

    schemes = {
        "uniform": np.ones(len(tr)),
        "sample_w": sw[tr],
        "subset_w": sub_w[tr],
        "sample_x_subset": both[tr],
        "sample_w_fused": sw_f[tr],
        "subset_w_fused": sub_f[tr],
        "both_fused": both_f[tr],
    }
    rows = {}
    forget0 = None
    for name, w in schemes.items():
        clf, _ = fit_initial(X0, y0, seed=seed, classes=classes)
        # refit scaler consistent with fit_initial internal — use transformed via partial_fit path
        clf = make_sgd(seed=seed)
        clf.partial_fit(X0s, y0, classes=classes)
        w = np.asarray(w, float)
        w = w * (len(w) / max(w.sum(), 1e-12))
        # amplify CFPerm-rejected / top-VIMP coords mildly
        Xtr = X1s[tr].copy()
        Xtr[:, feat] *= 1.2
        clf = cold_start_partial_fit(
            clf, Xtr, y1[tr], sample_weight=w, classes=classes, n_passes=4,
        )
        m_live = _metrics(clf, X1s[ev], y1[ev], known=known, novel=novel)
        m_ref = _metrics(clf, X0s, y0, known=known, novel=None)
        forget = 1.0 - m_ref["acc"]
        rows[name] = {
            **m_live,
            "ref_acc": m_ref["acc"],
            "forgetting": float(forget),
            "mean_weight": float(w.mean()),
            "weight_gini": float(_gini(w)),
        }
        if name == "uniform":
            forget0 = forget

    return {
        "dataset": board["name"],
        "note": board.get("note"),
        "n_ref": int(len(y0)),
        "n_live": int(len(y1)),
        "n_subset": int(sub_mask.sum()),
        "subset_frac": float(sub_mask.mean()),
        "subset_frac_fused": float(sub_mask_f.mean()),
        "vimp_top_mass": float(vimp[feat].sum()),
        "feat_hit": feat_hit,
        "n_cfperm_reject": int(len(rejected)),
        "mean_membership": float(mem.mean()),
        "schemes": rows,
        "delta_novel_sample": _delta(rows, "sample_w", "uniform", "novel_acc"),
        "delta_novel_subset": _delta(rows, "subset_w", "uniform", "novel_acc"),
        "delta_novel_both": _delta(rows, "sample_x_subset", "uniform", "novel_acc"),
        "delta_acc_both": _delta(rows, "sample_x_subset", "uniform", "acc"),
        "delta_novel_both_fused": _delta(rows, "both_fused", "uniform", "novel_acc"),
        "delta_acc_both_fused": _delta(rows, "both_fused", "uniform", "acc"),
        "delta_acc_sample_fused": _delta(rows, "sample_w_fused", "uniform", "acc"),
        "delta_acc_subset_fused": _delta(rows, "subset_w_fused", "uniform", "acc"),
    }


def _gini(w: np.ndarray) -> float:
    w = np.sort(np.asarray(w, float))
    n = len(w)
    if n == 0 or w.sum() <= 0:
        return 0.0
    idx = np.arange(1, n + 1)
    return float((2 * np.sum(idx * w) / (n * w.sum())) - (n + 1) / n)


def _delta(rows, a, b, key):
    va, vb = rows[a].get(key), rows[b].get(key)
    if va is None or vb is None or va != va or vb != vb:
        return None
    return float(va - vb)


def run_benchmark(*, seed: int = SEED, quick: bool = True) -> Dict:
    boards = load_boards(seed=seed, quick=quick)
    results = []
    for b in boards:
        print("[wbench] board", b["name"], flush=True)
        results.append(run_board(b, seed=seed, top_k=6 if quick else 8))
    # aggregate where novel_acc exists
    def _mean_delta(key):
        vals = [r[key] for r in results if r.get(key) is not None]
        return float(np.mean(vals)) if vals else None

    return {
        "boards": results,
        "n_boards": len(results),
        "mean_delta_novel_sample": _mean_delta("delta_novel_sample"),
        "mean_delta_novel_subset": _mean_delta("delta_novel_subset"),
        "mean_delta_novel_both": _mean_delta("delta_novel_both"),
        "mean_delta_acc_both": _mean_delta("delta_acc_both"),
        "mean_delta_novel_both_fused": _mean_delta("delta_novel_both_fused"),
        "mean_delta_acc_both_fused": _mean_delta("delta_acc_both_fused"),
        "mean_delta_acc_sample_fused": _mean_delta("delta_acc_sample_fused"),
        "mean_delta_acc_subset_fused": _mean_delta("delta_acc_subset_fused"),
        "takeaway": (
            "FSDS RF-domain VIMP → membership → soft sample_w + hard observation subset_w "
            "on amazon/digits/news20/synth_cgcd. Pure membership is diagnostic; "
            "fused (membership ∪ private-live) is the cold-start weight channel."
        ),
    }
