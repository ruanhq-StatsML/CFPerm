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
from sklearn.metrics import roc_auc_score
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

# Head is still CE / logistic: SGDClassifier(loss="log_loss").
LOSS_NAME = "CE (log_loss)"


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


def predict_proba_scaled(clf: SGDClassifier, X: np.ndarray) -> np.ndarray:
    sc = getattr(clf, "_first_cat_scaler_", None)
    Xs = sc.transform(X) if sc is not None else np.asarray(X, float)
    return clf.predict_proba(Xs)


def _class_column(clf: SGDClassifier, c: int) -> Optional[int]:
    classes = list(np.asarray(clf.classes_).ravel())
    try:
        return int(classes.index(int(c)))
    except ValueError:
        return None


def record_ova_first_category(
    y_true: Sequence,
    proba: np.ndarray,
    y_pred: Sequence,
    clf: SGDClassifier,
    first_idx: Dict[int, np.ndarray],
    *,
    neg_idx: np.ndarray,
    categories: Optional[Sequence] = None,
) -> Dict:
    """One-vs-all AUC + acc for each new category on its first-index window.

    For category ``c``:
      - positives = first_idx[c]
      - negatives = ``neg_idx`` (typically known / other-class rows)
      - score = P(y=c) from the CE multiclass head
      - ova_auc = ROC-AUC(1{y==c}, score)
      - ova_acc = accuracy of argmax==c on the positive first-window
        (also report threshold-0.5 OVA acc on pos∪neg)
    """
    y_true = np.asarray(y_true).ravel()
    y_pred = np.asarray(y_pred).ravel()
    proba = np.asarray(proba, float)
    neg_idx = np.asarray(neg_idx, int)
    cats = list(categories) if categories is not None else list(first_idx.keys())

    per: Dict[int, Dict] = {}
    aucs, accs_pos, accs_ova = [], [], []
    for c in cats:
        c = int(c)
        pos = np.asarray(first_idx.get(c, []), int)
        col = _class_column(clf, c)
        if pos.size == 0 or col is None:
            per[c] = {
                "n_pos": int(pos.size), "n_neg": int(neg_idx.size),
                "ova_auc": float("nan"), "acc_on_first": float("nan"),
                "ova_acc_thr05": float("nan"),
            }
            continue
        # OVA pool
        pool = np.unique(np.concatenate([pos, neg_idx])) if neg_idx.size else pos
        binary = (y_true[pool] == c).astype(int)
        score = proba[pool, col]
        # AUC needs both classes present
        if binary.min() == binary.max():
            auc = float("nan")
        else:
            auc = float(roc_auc_score(binary, score))
        acc_first = float(np.mean(y_pred[pos] == c))  # multiclass hit on first window
        pred_ova = (score >= 0.5).astype(int)
        ova_acc = float(np.mean(pred_ova == binary))
        per[c] = {
            "n_pos": int(pos.size),
            "n_neg": int(neg_idx.size),
            "ova_auc": auc,
            "acc_on_first": acc_first,
            "ova_acc_thr05": ova_acc,
        }
        if auc == auc:
            aucs.append(auc)
        if acc_first == acc_first:
            accs_pos.append(acc_first)
        if ova_acc == ova_acc:
            accs_ova.append(ova_acc)

    return {
        "loss": LOSS_NAME,
        "per_category": per,
        "mean_ova_auc": float(np.mean(aucs)) if aucs else float("nan"),
        "mean_acc_on_first": float(np.mean(accs_pos)) if accs_pos else float("nan"),
        "mean_ova_acc_thr05": float(np.mean(accs_ova)) if accs_ova else float("nan"),
        "categories": [int(c) for c in cats],
    }


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


# ---------------------------------------------------------------------------
# OVA AUC/acc on first-n (default 100) — where the method can empower
# ---------------------------------------------------------------------------


def _sample_board_n_per_novel(
    X: np.ndarray,
    y: np.ndarray,
    *,
    known: Sequence[int],
    novel: Sequence[int],
    n_first: int,
    n_ref: int,
    n_live_k: int,
    seed: int,
    name: str,
    note: str,
) -> Optional[Dict]:
    """Build a board with ≥ n_first rows per novel class (else skip)."""
    rng = np.random.default_rng(int(seed))
    X = np.asarray(X, float)
    y = np.asarray(y).ravel()
    known = [int(c) for c in known]
    novel = [int(c) for c in novel]

    # need enough per novel class
    per_n = {}
    for c in novel:
        idx = np.flatnonzero(y == c)
        per_n[c] = len(idx)
        if len(idx) < int(n_first):
            return None

    k_idx = np.flatnonzero(np.isin(y, known))
    if len(k_idx) < max(n_ref, n_live_k):
        return None
    ref = rng.choice(k_idx, size=min(n_ref, len(k_idx)), replace=False)
    live_k = rng.choice(k_idx, size=min(n_live_k, len(k_idx)), replace=False)

    live_n_parts = []
    for c in novel:
        idx = np.flatnonzero(y == c)
        # take n_first in a shuffled order (= arrival order within class)
        take = rng.choice(idx, size=int(n_first), replace=False)
        live_n_parts.append(take)
    live_n = np.concatenate(live_n_parts)
    live = np.concatenate([live_k, live_n])
    # do NOT shuffle away class blocks yet — arrival_order_live will structure it
    return {
        "name": name,
        "X_ref": X[ref],
        "y_ref": y[ref],
        "X_live": X[live],
        "y_live": y[live],
        "known_labels": known,
        "novel_labels": novel,
        "note": note,
        "n_first_target": int(n_first),
        "n_per_novel_available": {str(c): int(per_n[c]) for c in novel},
    }


def load_boards_n_first(*, seed: int = SEED, n_first: int = 100) -> List[Dict]:
    """Datasets that can supply n_first samples per novel category."""
    from fsds_cfperm_weight_benchmark import (
        _openml_xy,
        board_cancer,
        board_covtype,
        board_news20b,
        board_wine,
    )

    boards: List[Dict] = []

    # HAR
    try:
        X, y = _openml_xy("har", version=1)
        X = StandardScaler().fit_transform(X)
        b = _sample_board_n_per_novel(
            X, y, known=[0, 1, 2], novel=[3, 4, 5], n_first=n_first,
            n_ref=1400, n_live_k=400, seed=seed,
            name="har", note="HAR; locomotion→posture; n_first per novel",
        )
        if b:
            boards.append(b)
        else:
            print("[ova] har skip: not enough per novel", flush=True)
    except Exception as exc:  # noqa: BLE001
        print("[ova] har failed:", exc, flush=True)

    # MiceProtein
    try:
        X, y = _openml_xy("MiceProtein", version=4)
        X = StandardScaler().fit_transform(X)
        b = _sample_board_n_per_novel(
            X, y, known=[0, 1, 2, 3], novel=[4, 5, 6, 7], n_first=n_first,
            n_ref=320, n_live_k=160, seed=seed,
            name="mice", note="MiceProtein; n_first per trisomic class",
        )
        if b:
            boards.append(b)
        else:
            print("[ova] mice skip: not enough per novel for n_first=%d" % n_first, flush=True)
    except Exception as exc:  # noqa: BLE001
        print("[ova] mice failed:", exc, flush=True)

    # CNAE-9
    try:
        X, y = _openml_xy("cnae-9", version=1)
        if X.shape[1] > 128:
            keep = np.argsort(-X.var(0))[:128]
            X = X[:, keep]
        X = StandardScaler().fit_transform(X)
        b = _sample_board_n_per_novel(
            X, y, known=[0, 1, 2, 3, 4], novel=[5, 6, 7, 8], n_first=n_first,
            n_ref=400, n_live_k=200, seed=seed,
            name="cnae9", note="CNAE-9; n_first per novel doc class",
        )
        if b:
            boards.append(b)
        else:
            print("[ova] cnae9 skip: not enough per novel", flush=True)
    except Exception as exc:  # noqa: BLE001
        print("[ova] cnae9 failed:", exc, flush=True)

    # JapaneseVowels
    try:
        X, y = _openml_xy("JapaneseVowels", version=1)
        X = StandardScaler().fit_transform(X)
        b = _sample_board_n_per_novel(
            X, y, known=list(range(5)), novel=list(range(5, 9)), n_first=n_first,
            n_ref=1000, n_live_k=300, seed=seed,
            name="jvowel", note="JapaneseVowels; n_first per novel speaker",
        )
        if b:
            boards.append(b)
        else:
            print("[ova] jvowel skip: not enough per novel", flush=True)
    except Exception as exc:  # noqa: BLE001
        print("[ova] jvowel failed:", exc, flush=True)

    # Covtype (often has enough for novel 4,5)
    try:
        from sklearn.datasets import fetch_covtype
        data = fetch_covtype()
        X = np.asarray(data.data, float)
        y = np.asarray(data.target, int)
        rng = np.random.default_rng(seed)
        keep = np.isin(y, [1, 2, 3, 4, 5])
        idx = np.flatnonzero(keep)
        idx = rng.choice(idx, size=min(12000, len(idx)), replace=False)
        X = StandardScaler().fit_transform(X[idx])
        y = y[idx]
        b = _sample_board_n_per_novel(
            X, y, known=[1, 2, 3], novel=[4, 5], n_first=n_first,
            n_ref=1600, n_live_k=400, seed=seed,
            name="covtype", note="covertype; n_first per novel cover",
        )
        if b:
            boards.append(b)
        else:
            print("[ova] covtype skip: not enough per novel", flush=True)
    except Exception as exc:  # noqa: BLE001
        print("[ova] covtype failed:", exc, flush=True)

    # news20b — may not have 100/class in 4-topic subset; try
    try:
        b0 = board_news20b(seed=seed, max_features=96)
        # rebuild with larger novel if possible from same fetch path
        from sklearn.datasets import fetch_20newsgroups
        from sklearn.feature_extraction.text import TfidfVectorizer
        cats = ["alt.atheism", "talk.religion.misc", "sci.med", "comp.sys.ibm.pc.hardware"]
        bun = fetch_20newsgroups(
            subset="train", categories=cats, remove=("headers", "footers", "quotes"),
            shuffle=True, random_state=seed,
        )
        vec = TfidfVectorizer(max_features=96, stop_words="english")
        X = vec.fit_transform(bun.data).toarray()
        y = np.asarray(bun.target, int)
        b = _sample_board_n_per_novel(
            X, y, known=[0, 1], novel=[2, 3], n_first=n_first,
            n_ref=500, n_live_k=200, seed=seed,
            name="news20b", note="20newsgroups; n_first per novel topic",
        )
        if b:
            boards.append(b)
        else:
            print("[ova] news20b skip: not enough per novel (need %d)" % n_first, flush=True)
    except Exception as exc:  # noqa: BLE001
        print("[ova] news20b failed:", exc, flush=True)

    # wine / cancer usually too small for n_first=100 — skip silently unless they pass
    for fn, label in (
        (lambda: board_wine(seed=seed), "wine"),
        (lambda: board_cancer(seed=seed), "cancer"),
    ):
        try:
            b0 = fn()
            # check counts
            y = np.asarray(b0["y_live"]).ravel()
            ok = all((y == c).sum() >= n_first for c in b0["novel_labels"])
            if ok:
                boards.append(b0)
            else:
                print("[ova] %s skip: live novel < n_first=%d" % (label, n_first), flush=True)
        except Exception as exc:  # noqa: BLE001
            print("[ova] %s failed:" % label, exc, flush=True)

    return boards


def evaluate_board_ova_n_first(
    board: Dict,
    *,
    seed: int = SEED,
    n_first: int = 100,
    schemes: Sequence[str] = ("uniform", "weighted", "both_fused"),
    top_k: int = 8,
    eval_frac: float = 0.35,
    n_neg: int = 200,
) -> Dict:
    """CE head + OVA AUC/acc on first ``n_first`` of each novel category."""
    known = list(board.get("known_labels") or [])
    novel = list(board.get("novel_labels") or [])
    X_ref = np.asarray(board["X_ref"], float)
    y_ref = np.asarray(board["y_ref"]).ravel()
    X_live = np.asarray(board["X_live"], float)
    y_live = np.asarray(board["y_live"]).ravel()

    X_ord, y_ord, _ = arrival_order_live(
        X_live, y_live, known=known, novel=novel, seed=seed,
    )
    first_all = first_indices_per_category(
        y_ord, n_first=int(n_first), frac=None, categories=novel,
    )
    # actual n used (may be < n_first if board undersized — should not happen)
    actual_n = {int(c): int(len(idx)) for c, idx in first_all.items()}
    first_tr, first_eval = split_first_indices(
        first_all, seed=seed + 7, eval_frac=eval_frac, min_eval=max(5, int(0.2 * n_first)),
    )

    known_idx = np.flatnonzero(np.isin(y_ord, known)) if known else np.asarray([], dtype=int)
    tr = np.unique(np.concatenate([known_idx, union_indices(first_tr)]))
    rng = np.random.default_rng(seed + 13)
    if known_idx.size > n_neg:
        neg_idx = rng.choice(known_idx, size=int(n_neg), replace=False)
    else:
        neg_idx = known_idx

    classes = np.unique(np.concatenate([y_ref, y_ord]))
    scheme_rows = {}
    for sch in schemes:
        sch_key = sch
        clf = _fit_scheme(
            X_ref, y_ref, X_ord[tr], y_ord[tr],
            scheme=sch_key, seed=seed, classes=classes, top_k=top_k,
        )
        pred = predict_scaled(clf, X_ord)
        proba = predict_proba_scaled(clf, X_ord)
        ova = record_ova_first_category(
            y_ord, proba, pred, clf, first_eval,
            neg_idx=neg_idx, categories=novel,
        )
        scheme_rows[sch] = ova

    def _delta(metric):
        out = {}
        base = scheme_rows.get("uniform", {}).get(metric)
        for sch, row in scheme_rows.items():
            if sch == "uniform" or base is None or base != base:
                continue
            v = row.get(metric)
            out[sch] = float(v - base) if v is not None and v == v else None
        return out

    return {
        "dataset": board.get("name"),
        "note": board.get("note"),
        "loss": LOSS_NAME,
        "n_first": int(n_first),
        "actual_n_first": actual_n,
        "known_labels": known,
        "novel_labels": novel,
        "n_adapt": int(len(tr)),
        "n_neg_ova": int(len(neg_idx)),
        "schemes": {
            sch: {
                "mean_ova_auc": row["mean_ova_auc"],
                "mean_acc_on_first": row["mean_acc_on_first"],
                "mean_ova_acc_thr05": row["mean_ova_acc_thr05"],
                "per_category": {
                    str(c): v for c, v in row["per_category"].items()
                },
            }
            for sch, row in scheme_rows.items()
        },
        "delta_ova_auc_vs_uniform": _delta("mean_ova_auc"),
        "delta_acc_on_first_vs_uniform": _delta("mean_acc_on_first"),
    }


def run_ova_n_first_benchmark(
    *,
    seed: int = SEED,
    n_first: int = 100,
    schemes: Sequence[str] = ("uniform", "weighted", "both_fused"),
) -> Dict:
    boards = load_boards_n_first(seed=seed, n_first=n_first)
    rows = []
    for b in boards:
        print("[ova-nfirst] dataset", b.get("name"), "n_first", n_first, flush=True)
        rows.append(
            evaluate_board_ova_n_first(
                b, seed=seed, n_first=n_first, schemes=schemes,
            )
        )

    def _mean_delta(key, sch):
        vals = [
            r[key].get(sch)
            for r in rows
            if r.get(key, {}).get(sch) is not None
        ]
        return float(np.mean(vals)) if vals else None

    return {
        "seed": seed,
        "n_first": n_first,
        "loss": LOSS_NAME,
        "boards": rows,
        "n_boards": len(rows),
        "mean_delta_ova_auc_weighted": _mean_delta("delta_ova_auc_vs_uniform", "weighted"),
        "mean_delta_ova_auc_fused": _mean_delta("delta_ova_auc_vs_uniform", "both_fused"),
        "mean_delta_acc_first_weighted": _mean_delta("delta_acc_on_first_vs_uniform", "weighted"),
        "mean_delta_acc_first_fused": _mean_delta("delta_acc_on_first_vs_uniform", "both_fused"),
        "takeaway": (
            "Still CE (log_loss). Empowerment KPI = one-vs-all AUC / acc on the "
            "first n_first=%d indices of each new category; weights vs uniform."
            % int(n_first)
        ),
    }
