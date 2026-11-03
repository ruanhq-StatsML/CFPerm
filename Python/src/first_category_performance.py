"""Per-category *first* indices → record performance.

Idea
----
For each category ``c``, take the earliest slice of its stream indices
(the cold-start window of that class), then score accuracy on the union
(or per-class) of those indices.

This is the right KPI for ``initial stage → directed weights``:
not whole-LIVE novel_acc, but performance on the *first part* of each
new category.

API
---
  first_indices_per_category(y, n_first=..., frac=...)
  record_first_category_performance(y_true, y_pred, first_idx)
  evaluate_board_first_category(...)   # per-dataset, weights vs no-weights
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from sklearn.linear_model import SGDClassifier
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


Label = Union[int, np.integer]


# ---------------------------------------------------------------------------
# Core: find first indices + record performance
# ---------------------------------------------------------------------------


def first_indices_per_category(
    y: Sequence,
    *,
    n_first: Optional[int] = None,
    frac: Optional[float] = 0.25,
    categories: Optional[Sequence] = None,
    min_n: int = 1,
) -> Dict[int, np.ndarray]:
    """For each category, return the earliest indices in stream order.

    Parameters
    ----------
    y : label stream (order = time / arrival)
    n_first : keep this many first hits per class (overrides frac if set)
    frac : keep this fraction of each class's indices from the start
    categories : restrict to these labels (default = all in y)
    min_n : lower bound on kept count per class

    Returns
    -------
    dict: category -> 1d int array of indices into ``y`` (ascending).
    """
    y = np.asarray(y).ravel()
    if categories is None:
        cats = [int(c) for c in np.unique(y)]
    else:
        cats = [int(c) for c in categories]

    out: Dict[int, np.ndarray] = {}
    for c in cats:
        idx = np.flatnonzero(y == c)  # already stream-ordered
        if idx.size == 0:
            out[c] = np.asarray([], dtype=int)
            continue
        if n_first is not None:
            k = int(n_first)
        else:
            f = 0.25 if frac is None else float(frac)
            k = int(max(min_n, round(len(idx) * f)))
        k = int(np.clip(k, min_n, len(idx)))
        out[c] = idx[:k]
    return out


def union_indices(first_idx: Dict[int, np.ndarray]) -> np.ndarray:
    """Flatten + sort unique indices across categories."""
    if not first_idx:
        return np.asarray([], dtype=int)
    parts = [np.asarray(v, int) for v in first_idx.values() if len(v)]
    if not parts:
        return np.asarray([], dtype=int)
    return np.unique(np.concatenate(parts))


def record_first_category_performance(
    y_true: Sequence,
    y_pred: Sequence,
    first_idx: Dict[int, np.ndarray],
    *,
    categories: Optional[Sequence] = None,
) -> Dict:
    """Record accuracy on each category's first-index slice + micro average.

    Returns
    -------
    dict with:
      per_category: {c: {n, n_correct, acc, indices}}
      micro_acc: accuracy on union of first indices
      n_total: |union|
      categories: list of cats scored
    """
    y_true = np.asarray(y_true).ravel()
    y_pred = np.asarray(y_pred).ravel()
    if len(y_true) != len(y_pred):
        raise ValueError("y_true/y_pred length mismatch")

    cats = list(categories) if categories is not None else list(first_idx.keys())
    per: Dict[int, Dict] = {}
    for c in cats:
        c = int(c)
        idx = np.asarray(first_idx.get(c, []), int)
        if idx.size == 0:
            per[c] = {"n": 0, "n_correct": 0, "acc": float("nan"), "indices": []}
            continue
        correct = int(np.sum(y_pred[idx] == y_true[idx]))
        per[c] = {
            "n": int(idx.size),
            "n_correct": correct,
            "acc": float(correct / idx.size),
            "indices": idx.tolist(),
        }

    uni = union_indices({c: first_idx[c] for c in cats if c in first_idx})
    if uni.size == 0:
        micro = float("nan")
    else:
        micro = float(np.mean(y_pred[uni] == y_true[uni]))

    return {
        "per_category": per,
        "micro_acc": micro,
        "n_total": int(uni.size),
        "categories": [int(c) for c in cats],
        "union_indices": uni.tolist(),
    }


# ---------------------------------------------------------------------------
# Dataset helpers: build an arrival-ordered LIVE stream
# ---------------------------------------------------------------------------


def arrival_order_live(
    X_live: np.ndarray,
    y_live: np.ndarray,
    *,
    known: Sequence,
    novel: Sequence,
    seed: int = SEED,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Re-order LIVE so known recirculation comes first, then each novel class
    arrives as a contiguous block (stable cold-start windows).

    Returns (X_ord, y_ord, order_idx) where order_idx maps into original LIVE.
    """
    X = np.asarray(X_live, float)
    y = np.asarray(y_live).ravel()
    rng = np.random.default_rng(int(seed))
    known = [int(c) for c in known]
    novel = [int(c) for c in novel]

    parts: List[np.ndarray] = []
    # known first (shuffled)
    k_idx = np.flatnonzero(np.isin(y, known))
    if k_idx.size:
        parts.append(rng.permutation(k_idx))
    # each novel category as an arrival block (internal shuffle, blocks in novel order)
    for c in novel:
        n_idx = np.flatnonzero(y == c)
        if n_idx.size:
            parts.append(rng.permutation(n_idx))
    if not parts:
        order = np.arange(len(y))
    else:
        order = np.concatenate(parts)
    return X[order], y[order], order


# ---------------------------------------------------------------------------
# Per-dataset evaluation: first-category performance, weights vs no-weights
# ---------------------------------------------------------------------------


def _fit_scheme(
    X_ref: np.ndarray,
    y_ref: np.ndarray,
    X_tr: np.ndarray,
    y_tr: np.ndarray,
    *,
    scheme: str,
    seed: int,
    classes: np.ndarray,
    top_k: int = 8,
    lam: float = 2.5,
) -> SGDClassifier:
    sc = StandardScaler().fit(X_ref)
    X0s, Xtrs = sc.transform(X_ref), sc.transform(X_tr)
    classes = np.asarray(classes, int)
    clf = make_sgd(seed=seed)
    clf.partial_fit(X0s, y_ref, classes=classes)

    if scheme == "uniform":
        w = np.ones(len(y_tr), dtype=float)
        feat = np.arange(min(top_k, Xtrs.shape[1]))
    else:
        vimp = attribute_rf_domain(X0s, Xtrs, seed=seed, n_estimators=60)
        mem, feat = localize_subset(X0s, Xtrs, vimp, top_k=top_k)
        # fuse private-live if labels reveal novel in train (practical cold-start)
        # caller may already have only known+novel; boost high-mem
        mode = "fused" if scheme in ("both_fused", "repr_sep", "weighted") else "sample"
        if scheme == "sample_w":
            mode = "sample"
        wt = weights_from_attribution(
            mem, feat, Xtrs.shape[1], lam=lam, mode=mode, temper=True,
        )
        w = wt["sample_weight"]
        Xtrs = Xtrs.copy()
        Xtrs[:, feat] *= 1.08

    w = np.asarray(w, float)
    w = w * (len(w) / max(w.sum(), 1e-12))
    clf = cold_start_partial_fit(
        clf, Xtrs, y_tr, sample_weight=w, classes=classes, n_passes=3,
    )
    # stash scaler on clf for predict helper
    clf._first_cat_scaler_ = sc  # type: ignore[attr-defined]
    return clf


def predict_scaled(clf: SGDClassifier, X: np.ndarray) -> np.ndarray:
    sc = getattr(clf, "_first_cat_scaler_", None)
    Xs = sc.transform(X) if sc is not None else np.asarray(X, float)
    return clf.predict(Xs)


def split_first_indices(
    first_idx: Dict[int, np.ndarray],
    *,
    seed: int = SEED,
    eval_frac: float = 0.35,
    min_eval: int = 1,
) -> Tuple[Dict[int, np.ndarray], Dict[int, np.ndarray]]:
    """Within each category's *first* window, hold out a slice for eval.

    Never substitutes later-than-first samples — cold-start KPI stays honest.
    """
    rng = np.random.default_rng(int(seed))
    tr_map: Dict[int, np.ndarray] = {}
    ev_map: Dict[int, np.ndarray] = {}
    for c, idx in first_idx.items():
        idx = np.asarray(idx, int)
        if idx.size == 0:
            tr_map[int(c)] = idx
            ev_map[int(c)] = idx
            continue
        if idx.size == 1:
            # too small to split: adapt on it, also score on it (optimistic; flagged n=1)
            tr_map[int(c)] = idx
            ev_map[int(c)] = idx
            continue
        perm = rng.permutation(idx)
        n_ev = int(max(min_eval, round(len(perm) * float(eval_frac))))
        n_ev = min(n_ev, len(perm) - 1)  # keep ≥1 for adapt
        ev_map[int(c)] = np.sort(perm[:n_ev])
        tr_map[int(c)] = np.sort(perm[n_ev:])
    return tr_map, ev_map


def evaluate_board_first_category(
    board: Dict,
    *,
    seed: int = SEED,
    n_first: Optional[int] = None,
    frac: float = 0.25,
    schemes: Sequence[str] = ("uniform", "weighted", "both_fused"),
    top_k: int = 8,
    eval_frac: float = 0.35,
) -> Dict:
    """Per-dataset: arrival-order LIVE → first-indices per novel class → score.

    1. Re-order LIVE: known block, then each novel class arrival block.
    2. ``first_indices_per_category`` on that stream (earliest frac/n_first).
    3. Split those first-indices into adapt_tr / first_eval (no later substitutes).
    4. Fit schemes on known + adapt_tr; record performance on first_eval.
    """
    known = list(board.get("known_labels") or [])
    novel = list(board.get("novel_labels") or [])
    X_ref = np.asarray(board["X_ref"], float)
    y_ref = np.asarray(board["y_ref"]).ravel()
    X_live = np.asarray(board["X_live"], float)
    y_live = np.asarray(board["y_live"]).ravel()

    X_ord, y_ord, _order = arrival_order_live(
        X_live, y_live, known=known, novel=novel, seed=seed,
    )
    first_all = first_indices_per_category(
        y_ord, n_first=n_first, frac=frac, categories=novel if novel else None,
    )
    first_tr, first_eval = split_first_indices(
        first_all, seed=seed + 7, eval_frac=eval_frac,
    )

    # adapt set = all known in ordered LIVE + first-window train split
    known_idx = np.flatnonzero(np.isin(y_ord, known)) if known else np.asarray([], dtype=int)
    tr_first_u = union_indices(first_tr)
    tr = np.unique(np.concatenate([known_idx, tr_first_u])) if tr_first_u.size else known_idx
    if tr.size == 0:
        tr = np.arange(min(16, len(y_ord)))

    classes = np.unique(np.concatenate([y_ref, y_ord]))
    scheme_rows = {}
    for sch in schemes:
        sch_key = {
            "uniform": "uniform",
            "weighted": "weighted",
            "both_fused": "both_fused",
            "sample_w": "sample_w",
            "repr_sep": "repr_sep",
        }.get(sch, sch)
        clf = _fit_scheme(
            X_ref, y_ref, X_ord[tr], y_ord[tr],
            scheme=sch_key, seed=seed, classes=classes, top_k=top_k,
        )
        pred = predict_scaled(clf, X_ord)
        rec = record_first_category_performance(
            y_ord, pred, first_eval, categories=novel if novel else None,
        )
        # reference: all novel in ordered stream
        if novel:
            m = np.isin(y_ord, novel)
            novel_acc_all = float(np.mean(pred[m] == y_ord[m])) if m.any() else float("nan")
        else:
            novel_acc_all = float("nan")
        scheme_rows[sch] = {
            "first_micro_acc": rec["micro_acc"],
            "first_n_total": rec["n_total"],
            "per_category": {
                str(c): {"n": v["n"], "acc": v["acc"]} for c, v in rec["per_category"].items()
            },
            "novel_acc_all": novel_acc_all,
        }

    base = scheme_rows.get("uniform", {}).get("first_micro_acc")
    deltas = {}
    for sch, row in scheme_rows.items():
        if sch == "uniform" or base is None or base != base:
            continue
        v = row["first_micro_acc"]
        deltas[sch] = float(v - base) if v == v else None

    return {
        "dataset": board.get("name"),
        "note": board.get("note"),
        "known_labels": known,
        "novel_labels": novel,
        "n_live": int(len(y_ord)),
        "n_adapt": int(len(tr)),
        "n_known_adapt": int(len(known_idx)),
        "n_first_adapt": int(tr_first_u.size),
        "first_spec": {"n_first": n_first, "frac": frac, "eval_frac": eval_frac},
        "first_indices_stream": {str(c): v.tolist() for c, v in first_all.items()},
        "first_indices_adapt": {str(c): v.tolist() for c, v in first_tr.items()},
        "first_indices_eval": {str(c): v.tolist() for c, v in first_eval.items()},
        "schemes": scheme_rows,
        "delta_first_vs_uniform": deltas,
    }


def run_all_boards_first_category(
    *,
    seed: int = SEED,
    quick: bool = True,
    frac: float = 0.25,
    n_first: Optional[int] = None,
    include_legacy: bool = True,
) -> Dict:
    """Run first-category performance across datasets."""
    from fsds_cfperm_weight_benchmark import (
        board_cancer,
        board_cnae9,
        board_covtype,
        board_har,
        board_jvowel,
        board_mice,
        board_news20b,
        board_wine,
        load_boards,
    )

    boards = load_boards(seed=seed, quick=quick)
    if include_legacy:
        for fn in (
            lambda: board_wine(seed=seed),
            lambda: board_cancer(seed=seed),
            lambda: board_news20b(seed=seed),
            lambda: board_covtype(seed=seed, quick=quick),
        ):
            try:
                boards.append(fn())
            except Exception as exc:  # noqa: BLE001
                print("[first-cat] legacy board failed:", exc, flush=True)

    rows = []
    for b in boards:
        print("[first-cat] dataset", b.get("name"), flush=True)
        rows.append(
            evaluate_board_first_category(
                b, seed=seed, frac=frac, n_first=n_first,
                schemes=("uniform", "weighted", "both_fused"),
            )
        )
    # aggregate
    def _mean_delta(sch):
        vals = [
            r["delta_first_vs_uniform"].get(sch)
            for r in rows
            if r.get("delta_first_vs_uniform", {}).get(sch) is not None
        ]
        return float(np.mean(vals)) if vals else None

    return {
        "seed": seed,
        "frac": frac,
        "n_first": n_first,
        "boards": rows,
        "mean_delta_weighted": _mean_delta("weighted"),
        "mean_delta_both_fused": _mean_delta("both_fused"),
        "takeaway": (
            "Per-dataset first-category indices (earliest frac/n_first per novel class); "
            "KPI = micro_acc on that cold-start window; compare weights vs uniform."
        ),
    }
