"""Benchmark FSDS/CFPerm → sample weights + observation-subset weights.

Datasets (v5 swap)
------------------
  wine       UCI wine; known {0,1} → sparse inject {2}
  covtype    covertype subsample; known {1,2,3} → sparse {4,5}
  news20b    TF-IDF 4-topic board (atheism/religion/med/hardware)
  cancer     breast-cancer; REF class0 + LIVE recirculate + sparse class1

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
from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import StandardScaler

from weighted_cold_start import (
    SEED,
    attribute_rf_domain,
    cold_start_partial_fit,
    fit_initial,
    localize_subset,
    make_sgd,
    weights_from_attribution,
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


def _sparse_novel_split(X, y, known, novel, *, seed, n_ref, n_live_k, n_live_n):
    rng = np.random.default_rng(seed)
    X, y = np.asarray(X, float), np.asarray(y).ravel()
    k_idx = np.flatnonzero(np.isin(y, known))
    n_idx = np.flatnonzero(np.isin(y, novel))
    ref = rng.choice(k_idx, size=min(n_ref, len(k_idx)), replace=False)
    live_k = rng.choice(k_idx, size=min(n_live_k, len(k_idx)), replace=False)
    live_n = rng.choice(n_idx, size=min(n_live_n, len(n_idx)), replace=False)
    live = np.concatenate([live_k, live_n])
    rng.shuffle(live)
    return X[ref], y[ref], X[live], y[live]


def board_wine(*, seed: int = SEED) -> Dict:
    from sklearn.datasets import load_wine
    X, y = load_wine(return_X_y=True)
    X = StandardScaler().fit_transform(np.asarray(X, float))
    Xr, yr, Xl, yl = _sparse_novel_split(
        X, y, known=[0, 1], novel=[2], seed=seed, n_ref=90, n_live_k=50, n_live_n=12,
    )
    return {
        "name": "wine", "X_ref": Xr, "y_ref": yr, "X_live": Xl, "y_live": yl,
        "note": "UCI wine; known {0,1} → sparse inject {2}",
        "known_labels": [0, 1], "novel_labels": [2],
    }


def board_covtype(*, seed: int = SEED, quick: bool = True) -> Dict:
    from sklearn.datasets import fetch_covtype
    data = fetch_covtype()
    X, y = np.asarray(data.data, float), np.asarray(data.target, int)
    # subsample for speed
    rng = np.random.default_rng(seed)
    keep = np.isin(y, [1, 2, 3, 4, 5])
    idx = np.flatnonzero(keep)
    idx = rng.choice(idx, size=min(4000 if quick else 8000, len(idx)), replace=False)
    X, y = StandardScaler().fit_transform(X[idx]), y[idx]
    Xr, yr, Xl, yl = _sparse_novel_split(
        X, y, known=[1, 2, 3], novel=[4, 5], seed=seed + 1,
        n_ref=800 if quick else 1600, n_live_k=300, n_live_n=60,
    )
    return {
        "name": "covtype", "X_ref": Xr, "y_ref": yr, "X_live": Xl, "y_live": yl,
        "note": "covertype subsample; known {1,2,3} → sparse {4,5}",
        "known_labels": [1, 2, 3], "novel_labels": [4, 5],
    }


def board_news20b(*, seed: int = SEED, max_features: int = 96) -> Dict:
    cats = ["alt.atheism", "talk.religion.misc", "sci.med", "comp.sys.ibm.pc.hardware"]
    bun = fetch_20newsgroups(
        subset="train", categories=cats, remove=("headers", "footers", "quotes"),
        shuffle=True, random_state=seed,
    )
    vec = TfidfVectorizer(max_features=max_features, stop_words="english")
    X = vec.fit_transform(bun.data).toarray()
    y = np.asarray(bun.target, int)
    Xr, yr, Xl, yl = _sparse_novel_split(
        X, y, known=[0, 1], novel=[2, 3], seed=seed,
        n_ref=500, n_live_k=180, n_live_n=40,
    )
    return {
        "name": "news20b", "X_ref": Xr, "y_ref": yr, "X_live": Xl, "y_live": yl,
        "note": "TF-IDF 4-topic; known atheism/religion → sparse med/hardware",
        "known_labels": [0, 1], "novel_labels": [2, 3],
    }


def board_cancer(*, seed: int = SEED) -> Dict:
    from sklearn.datasets import load_breast_cancer
    X, y = load_breast_cancer(return_X_y=True)
    X = StandardScaler().fit_transform(np.asarray(X, float))
    # REF on benign(0); LIVE recirculate benign + sparse malignant(1) as "novel"
    Xr, yr, Xl, yl = _sparse_novel_split(
        X, y, known=[0], novel=[1], seed=seed, n_ref=220, n_live_k=100, n_live_n=24,
    )
    return {
        "name": "cancer", "X_ref": Xr, "y_ref": yr, "X_live": Xl, "y_live": yl,
        "note": "breast-cancer; REF benign → LIVE + sparse malignant",
        "known_labels": [0], "novel_labels": [1],
    }


def load_boards(*, seed: int = SEED, quick: bool = True) -> List[Dict]:
    boards = []
    for fn in (
        lambda: board_wine(seed=seed),
        lambda: board_covtype(seed=seed, quick=quick),
        lambda: board_news20b(seed=seed, max_features=64 if quick else 96),
        lambda: board_cancer(seed=seed),
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

    # membership-only channels (diagnostic)
    wt_s = weights_from_attribution(mem, feat, X1s.shape[1], lam=lam, subset_q=subset_q, subset_boost=subset_boost, mode="sample", temper=False)
    wt_u = weights_from_attribution(mem, feat, X1s.shape[1], lam=lam, subset_q=subset_q, subset_boost=subset_boost, mode="subset", temper=False)
    wt_b = weights_from_attribution(mem, feat, X1s.shape[1], lam=lam, subset_q=subset_q, subset_boost=subset_boost, mode="fused", temper=True)
    # private-live fused channels (cold-start default)
    mem_f = mem.copy()
    if novel is not None:
        mem_f = np.maximum(mem_f, np.isin(y1, novel).astype(float))
    wt_sf = weights_from_attribution(mem_f, feat, X1s.shape[1], lam=lam, subset_q=subset_q, subset_boost=subset_boost, mode="sample", temper=False)
    wt_uf = weights_from_attribution(mem_f, feat, X1s.shape[1], lam=lam, subset_q=subset_q, subset_boost=subset_boost, mode="subset", temper=False)
    wt_bf = weights_from_attribution(mem_f, feat, X1s.shape[1], lam=lam, subset_q=subset_q, subset_boost=subset_boost, mode="fused", temper=True)

    # train/eval split on LIVE
    rng = np.random.default_rng(seed + 11)
    idx = rng.permutation(len(y1))
    n_ev = max(20, len(y1) // 3)
    ev, tr = idx[:n_ev], idx[n_ev:]

    schemes = {
        "uniform": np.ones(len(tr)),
        "sample_w": wt_s["sample_weight"][tr],
        "subset_w": wt_u["sample_weight"][tr],
        "sample_x_subset": wt_b["sample_weight"][tr],
        "sample_w_fused": wt_sf["sample_weight"][tr],
        "subset_w_fused": wt_uf["sample_weight"][tr],
        "both_fused": wt_bf["sample_weight"][tr],
    }
    sub_mask = wt_b["subset_mask"]
    sub_mask_f = wt_bf["subset_mask"]
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
        Xtr[:, feat] *= 1.08
        clf = cold_start_partial_fit(
            clf, Xtr, y1[tr], sample_weight=w, classes=classes, n_passes=3,
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
