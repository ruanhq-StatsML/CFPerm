"""Text multiclass: feature-level vs sample-level weights — FULL accuracy.

Not first-n subset. KPI = overall / known / novel accuracy on the full
held-out LIVE split.

Schemes
-------
  uniform     : CE, no extra weights
  sample_w    : sample-level w_i from FSDS membership (partial_fit sample_weight)
  feature_w   : feature-level scale x_j ← x_j · (1 + α v_j) from RF-domain VIMP
  both        : feature_w + sample_w
  ewma        : novel EWMA boost w_i = 1 + β · exp(-λ · n_idx_c)
  ewma_pure   : literal w_i = exp(-λ · n_idx_c) (down-weights late novel)
  ewma_feat   : ewma + feature_w
  ewma_both   : ewma × sample_w (+ optional feature_w)

Head remains CE (log_loss). Tabular boards are out of scope here.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import StandardScaler

from first_category_performance import arrival_order_live
from text_ova_nfirst import (
    board_amazon_category,
    board_cnae9_text,
    board_news20_comp,
    board_news20_mix,
    board_news20_rec,
)
from weighted_cold_start import (
    SEED,
    attribute_rf_domain,
    cold_start_partial_fit,
    localize_subset,
    make_sgd,
    sample_weights,
)

LOSS_NAME = "CE (log_loss)"

# Full-acc boards: take as many novel samples as available (cap), not first-100.
_FULL_N_PER_NOVEL = 400
# EWMA half-life ≈ ln(2)/λ ≈ 35 samples → λ ≈ 0.02 (intuitive cold-start focus)
_DEFAULT_EWMA_LAM = 0.02
_DEFAULT_EWMA_BETA = 3.0


def feature_level_weights(
    vimp: np.ndarray,
    *,
    alpha: float = 3.0,
    floor: float = 0.25,
    cap: float = 4.0,
) -> np.ndarray:
    """Turn RF-domain VIMP into per-feature multiplicative weights.

    ``x_j ← x_j * w^feat_j``, with
        w^feat_j = clip(1 + α · ṽ_j, floor, cap),
    ``ṽ`` = vimp / max(vimp) so the top attributed coord gets ~1+α.
    """
    v = np.asarray(vimp, float).ravel()
    vmax = float(v.max()) if v.size else 1.0
    vn = v / max(vmax, 1e-12)
    w = 1.0 + float(alpha) * vn
    return np.clip(w, float(floor), float(cap))


def apply_feature_weights(X: np.ndarray, feat_w: np.ndarray) -> np.ndarray:
    X = np.asarray(X, float)
    fw = np.asarray(feat_w, float).ravel()
    if fw.shape[0] != X.shape[1]:
        raise ValueError("feat_w dim mismatch")
    return X * fw.reshape(1, -1)


def novel_category_arrival_index(
    y: np.ndarray,
    novel: Sequence,
) -> np.ndarray:
    """Within-category arrival index n_idx for novel labels (0,1,2,…); -1 if known."""
    y = np.asarray(y).ravel()
    novel_set = {int(c) for c in novel}
    idx = np.full(len(y), -1, dtype=int)
    counts: Dict[int, int] = {}
    for i, yi in enumerate(y):
        c = int(yi)
        if c not in novel_set:
            continue
        k = counts.get(c, 0)
        idx[i] = k
        counts[c] = k + 1
    return idx


def ewma_novel_weights(
    y: np.ndarray,
    novel: Sequence,
    *,
    lam: float = _DEFAULT_EWMA_LAM,
    beta: float = _DEFAULT_EWMA_BETA,
    known_weight: float = 1.0,
    arrival_idx: Optional[np.ndarray] = None,
    mode: str = "boost",
) -> np.ndarray:
    """EWMA decay on novel-category arrival index ``n_idx``.

    Modes
    -----
    boost : ``w = 1 + β · exp(-λ · n_idx)`` for novel (early↑, late→1); known=1.
            Keeps late novel in the CE mass — the practical cold-start form.
    pure  : ``w = exp(-λ · n_idx)`` literally (late novel → 0). Useful ablation.
    """
    y = np.asarray(y).ravel()
    if arrival_idx is None:
        arrival_idx = novel_category_arrival_index(y, novel)
    arrival_idx = np.asarray(arrival_idx, int).ravel()
    w = np.full(len(y), float(known_weight), dtype=float)
    m = arrival_idx >= 0
    decay = np.exp(-float(lam) * arrival_idx[m].astype(float))
    if str(mode).lower() in ("pure", "raw", "exp"):
        w[m] = decay
    else:
        w[m] = 1.0 + float(beta) * decay
    return w


def _metrics_full(y_true: np.ndarray, y_pred: np.ndarray, *, known, novel) -> Dict:
    y_true = np.asarray(y_true).ravel()
    y_pred = np.asarray(y_pred).ravel()
    out = {"acc": float(accuracy_score(y_true, y_pred))}
    if known is not None:
        m = np.isin(y_true, known)
        out["known_acc"] = float(np.mean(y_pred[m] == y_true[m])) if m.any() else float("nan")
    if novel is not None:
        m = np.isin(y_true, novel)
        out["novel_acc"] = float(np.mean(y_pred[m] == y_true[m])) if m.any() else float("nan")
        out["n_novel_eval"] = int(m.sum())
    out["n_eval"] = int(len(y_true))
    return out


def fit_text_scheme(
    X_ref: np.ndarray,
    y_ref: np.ndarray,
    X_tr: np.ndarray,
    y_tr: np.ndarray,
    *,
    scheme: str,
    seed: int,
    classes: np.ndarray,
    known: Sequence,
    novel: Sequence,
    top_k: int = 16,
    lam_sample: float = 2.5,
    alpha_feat: float = 3.0,
    ewma_lam: float = _DEFAULT_EWMA_LAM,
    ewma_beta: float = _DEFAULT_EWMA_BETA,
    arrival_idx: Optional[np.ndarray] = None,
) -> Tuple[SGDClassifier, Dict]:
    """Fit one scheme; returns clf + diagnostics (feat_w, sample stats, vimp)."""
    scheme = str(scheme).lower()
    sc = StandardScaler().fit(X_ref)
    X0 = sc.transform(np.asarray(X_ref, float))
    X1 = sc.transform(np.asarray(X_tr, float))
    y0 = np.asarray(y_ref).ravel()
    y1 = np.asarray(y_tr).ravel()
    classes = np.asarray(classes, int)

    clf = make_sgd(seed=seed)
    clf.partial_fit(X0, y0, classes=classes)

    use_feat = scheme in (
        "feature_w", "both", "feature", "feat", "ewma_feat", "ewma_both",
    )
    use_sample = scheme in (
        "sample_w", "both", "weighted", "sample", "both_fused", "ewma_both",
    )
    use_ewma = scheme in (
        "ewma", "ewma_feat", "ewma_both", "ewma_w", "ewma_pure", "ewma_boost",
    )
    ewma_mode = "pure" if scheme == "ewma_pure" else "boost"
    if scheme == "uniform":
        use_feat = use_sample = use_ewma = False

    feat_w = np.ones(X1.shape[1], dtype=float)
    sample_w = np.ones(len(y1), dtype=float)
    vimp = None
    mem = None
    feat_idx = None

    if use_feat or use_sample:
        vimp = attribute_rf_domain(X0, X1, seed=seed, n_estimators=80)
        mem, feat_idx = localize_subset(X0, X1, vimp, top_k=top_k)
        # private-live fuse for membership (cold-start)
        if novel:
            mem = np.maximum(mem, np.isin(y1, list(novel)).astype(float))

    if use_feat:
        feat_w = feature_level_weights(vimp, alpha=alpha_feat)
        X0 = apply_feature_weights(X0, feat_w)
        X1 = apply_feature_weights(X1, feat_w)
        # re-init head on feature-weighted REF so coefs live in same geometry
        clf = make_sgd(seed=seed)
        clf.partial_fit(X0, y0, classes=classes)

    parts = []
    if use_sample:
        parts.append(sample_weights(mem, lam=lam_sample, cap=5.0))
    if use_ewma:
        parts.append(
            ewma_novel_weights(
                y1, novel, lam=ewma_lam, beta=ewma_beta,
                arrival_idx=arrival_idx, mode=ewma_mode,
            )
        )
    if parts:
        sample_w = np.ones(len(y1), dtype=float)
        for p in parts:
            sample_w = sample_w * np.asarray(p, float).ravel()
    else:
        sample_w = np.ones(len(y1), dtype=float)
    sample_w = sample_w * (len(sample_w) / max(sample_w.sum(), 1e-12))

    clf = cold_start_partial_fit(
        clf, X1, y1, sample_weight=sample_w, classes=classes, n_passes=4,
    )
    # stash transform for full-eval inference
    clf._tw_scaler_ = sc  # type: ignore[attr-defined]
    clf._tw_feat_w_ = feat_w  # type: ignore[attr-defined]
    clf._tw_use_feat_ = bool(use_feat)  # type: ignore[attr-defined]

    n_idx = arrival_idx if arrival_idx is not None else novel_category_arrival_index(y1, novel)
    novel_m = n_idx >= 0
    info = {
        "scheme": scheme,
        "use_feature_weights": bool(use_feat),
        "use_sample_weights": bool(use_sample),
        "use_ewma_weights": bool(use_ewma),
        "ewma_mode": ewma_mode if use_ewma else None,
        "ewma_lam": float(ewma_lam) if use_ewma else None,
        "ewma_beta": float(ewma_beta) if use_ewma and ewma_mode == "boost" else None,
        "feat_w_mean": float(feat_w.mean()),
        "feat_w_max": float(feat_w.max()),
        "feat_w_top_mass": float(np.sort(feat_w)[-max(1, len(feat_w)//10):].sum() / feat_w.sum()),
        "sample_w_mean": float(sample_w.mean()),
        "sample_w_max": float(sample_w.max()),
        "ewma_novel_w_mean": float(sample_w[novel_m].mean()) if novel_m.any() else None,
        "vimp_top_mass": float(vimp[feat_idx].sum()) if vimp is not None and feat_idx is not None else None,
        "mean_membership": float(np.mean(mem)) if mem is not None else None,
    }
    return clf, info


def predict_text_scheme(clf: SGDClassifier, X: np.ndarray) -> np.ndarray:
    sc = getattr(clf, "_tw_scaler_", None)
    Xs = sc.transform(np.asarray(X, float)) if sc is not None else np.asarray(X, float)
    if getattr(clf, "_tw_use_feat_", False):
        Xs = apply_feature_weights(Xs, clf._tw_feat_w_)  # type: ignore[attr-defined]
    return clf.predict(Xs)


def evaluate_board_full_acc(
    board: Dict,
    *,
    seed: int = SEED,
    schemes: Sequence[str] = (
        "uniform", "sample_w", "feature_w", "both", "ewma", "ewma_feat",
    ),
    eval_frac: float = 0.35,
    top_k: int = 16,
    alpha_feat: float = 3.0,
    ewma_lam: float = _DEFAULT_EWMA_LAM,
) -> Dict:
    """Train/eval split on LIVE; report FULL held-out accuracy (not first-n).

    LIVE is arrival-ordered so EWMA ``n_idx`` is the within-novel-category
    arrival index (0 = first hit of that new category).
    """
    known = list(board.get("known_labels") or [])
    novel = list(board.get("novel_labels") or [])
    X_ref = np.asarray(board["X_ref"], float)
    y_ref = np.asarray(board["y_ref"]).ravel()
    X_live_raw = np.asarray(board["X_live"], float)
    y_live_raw = np.asarray(board["y_live"]).ravel()
    X_live, y_live, _ = arrival_order_live(
        X_live_raw, y_live_raw, known=known, novel=novel, seed=seed,
    )
    # stream-level arrival index (preserved under train/eval subsetting)
    stream_idx = novel_category_arrival_index(y_live, novel)

    rng = np.random.default_rng(int(seed) + 11)
    idx = rng.permutation(len(y_live))
    n_ev = max(30, int(round(len(y_live) * float(eval_frac))))
    ev, tr = idx[:n_ev], idx[n_ev:]
    classes = np.unique(np.concatenate([y_ref, y_live]))

    rows = {}
    for sch in schemes:
        clf, info = fit_text_scheme(
            X_ref, y_ref, X_live[tr], y_live[tr],
            scheme=sch, seed=seed, classes=classes,
            known=known, novel=novel, top_k=top_k, alpha_feat=alpha_feat,
            ewma_lam=ewma_lam, arrival_idx=stream_idx[tr],
        )
        pred = predict_text_scheme(clf, X_live[ev])
        m = _metrics_full(y_live[ev], pred, known=known, novel=novel)
        # also REF retention
        pred_ref = predict_text_scheme(clf, X_ref)
        ref_acc = float(accuracy_score(y_ref, pred_ref))
        rows[sch] = {
            **m,
            "ref_acc": ref_acc,
            "forgetting": float(1.0 - ref_acc),
            **info,
        }

    base = rows.get("uniform", {})
    deltas = {}
    for sch, row in rows.items():
        if sch == "uniform":
            continue
        deltas[sch] = {
            "delta_acc": float(row["acc"] - base["acc"]) if base else None,
            "delta_novel_acc": (
                float(row["novel_acc"] - base["novel_acc"])
                if base and row.get("novel_acc") == row.get("novel_acc")
                and base.get("novel_acc") == base.get("novel_acc")
                else None
            ),
            "delta_known_acc": (
                float(row["known_acc"] - base["known_acc"])
                if base and row.get("known_acc") == row.get("known_acc")
                and base.get("known_acc") == base.get("known_acc")
                else None
            ),
        }

    return {
        "dataset": board.get("name"),
        "note": board.get("note"),
        "loss": LOSS_NAME,
        "ewma_lam": float(ewma_lam),
        "n_live": int(len(y_live)),
        "n_train": int(len(tr)),
        "n_eval": int(len(ev)),
        "n_novel_eval": int(np.isin(y_live[ev], novel).sum()) if novel else 0,
        "schemes": rows,
        "deltas_vs_uniform": deltas,
    }


def evaluate_board_incremental(
    board: Dict,
    *,
    seed: int = SEED,
    schemes: Sequence[str] = ("uniform", "ewma", "feature_w", "ewma_feat"),
    ewma_lam: float = _DEFAULT_EWMA_LAM,
    alpha_feat: float = 3.0,
    prefixes: Sequence[int] = (10, 25, 50, 100, 150, 200),
    eval_frac: float = 0.35,
) -> Dict:
    """Incremental cold-start: fit on first-k novel hits/category, score FULL eval.

    Shows EWMA's early-stage lift as ``k`` grows (增量效果).
    """
    known = list(board.get("known_labels") or [])
    novel = list(board.get("novel_labels") or [])
    X_ref = np.asarray(board["X_ref"], float)
    y_ref = np.asarray(board["y_ref"]).ravel()
    X_live_raw = np.asarray(board["X_live"], float)
    y_live_raw = np.asarray(board["y_live"]).ravel()
    X_live, y_live, _ = arrival_order_live(
        X_live_raw, y_live_raw, known=known, novel=novel, seed=seed,
    )
    stream_idx = novel_category_arrival_index(y_live, novel)
    classes = np.unique(np.concatenate([y_ref, y_live]))

    # reserve a late novel band for a fixed FULL eval pool (disjoint from train prefixes)
    rng = np.random.default_rng(int(seed) + 17)
    max_novel_idx = int(stream_idx.max()) if stream_idx.size and stream_idx.max() >= 0 else 0
    max_pref = int(max(prefixes))
    # keep at least ~25% of novel mass for eval when possible
    hold_from = min(max_pref, max(1, int(0.75 * (max_novel_idx + 1))))
    prefixes = [int(k) for k in prefixes if int(k) <= hold_from]
    late_novel = np.flatnonzero(stream_idx >= hold_from)
    known_idx = np.flatnonzero(stream_idx < 0)
    pool = np.concatenate([late_novel, known_idx]) if late_novel.size else np.arange(len(y_live))
    if pool.size < 40:
        pool = np.arange(len(y_live))
    n_ev = max(30, int(round(len(pool) * float(eval_frac))))
    n_ev = min(n_ev, len(pool))
    ev = rng.choice(pool, size=n_ev, replace=False)

    curve = {sch: [] for sch in schemes}
    for k in prefixes:
        # train = all known + novel with arrival_idx < k (and not in eval)
        tr_mask = (stream_idx < 0) | ((stream_idx >= 0) & (stream_idx < int(k)))
        tr = np.flatnonzero(tr_mask)
        tr = np.setdiff1d(tr, ev, assume_unique=False)
        if tr.size < 10:
            continue
        for sch in schemes:
            clf, info = fit_text_scheme(
                X_ref, y_ref, X_live[tr], y_live[tr],
                scheme=sch, seed=seed, classes=classes,
                known=known, novel=novel, alpha_feat=alpha_feat,
                ewma_lam=ewma_lam, arrival_idx=stream_idx[tr],
            )
            pred = predict_text_scheme(clf, X_live[ev])
            m = _metrics_full(y_live[ev], pred, known=known, novel=novel)
            curve[sch].append({
                "k": int(k),
                "n_train": int(tr.size),
                "acc": m["acc"],
                "novel_acc": m.get("novel_acc"),
                "known_acc": m.get("known_acc"),
                "sample_w_max": info.get("sample_w_max"),
            })

    # deltas vs uniform at each k
    deltas = {}
    for sch in schemes:
        if sch == "uniform":
            continue
        deltas[sch] = []
        u = {p["k"]: p for p in curve.get("uniform", [])}
        for p in curve[sch]:
            bu = u.get(p["k"])
            if not bu:
                continue
            deltas[sch].append({
                "k": p["k"],
                "delta_acc": float(p["acc"] - bu["acc"]),
                "delta_novel_acc": (
                    float(p["novel_acc"] - bu["novel_acc"])
                    if p.get("novel_acc") == p.get("novel_acc")
                    and bu.get("novel_acc") == bu.get("novel_acc")
                    else None
                ),
            })

    return {
        "dataset": board.get("name"),
        "ewma_lam": float(ewma_lam),
        "prefixes": [int(k) for k in prefixes],
        "n_eval": int(len(ev)),
        "curve": curve,
        "deltas_vs_uniform": deltas,
    }


def _expand_board_to_full_live(board: Dict, *, seed: int, n_per_novel: int) -> Dict:
    """Re-sample LIVE with up to ``n_per_novel`` per novel class (full mass).

    The n_first loaders only keep 100 novel rows; for overall accuracy we want
    the full held-out LIVE pool, not that cold-start subset.
    """
    # Prefer re-building via the board's available counts when present.
    avail = board.get("n_per_novel_available") or {}
    if not avail:
        return board
    # Cap by what the loader already put in X_live: if already large, keep.
    novel = list(board.get("novel_labels") or [])
    y_live = np.asarray(board["y_live"]).ravel()
    cur = min((int((y_live == c).sum()) for c in novel), default=0) if novel else 0
    if cur >= int(n_per_novel) * 0.9:
        return board
    # Cannot re-expand without raw X,y — callers should pass fuller boards.
    return board


def load_text_boards_full(
    *,
    seed: int = SEED,
    n_per_novel: int = _FULL_N_PER_NOVEL,
) -> List[Dict]:
    """Text boards with large LIVE novel mass for FULL accuracy eval."""
    # Reuse board builders but request n_first = n_per_novel so LIVE is full-sized.
    boards: List[Dict] = []
    for fn in (
        lambda: board_news20_comp(seed=seed, n_first=n_per_novel),
        lambda: board_news20_rec(seed=seed, n_first=n_per_novel),
        lambda: board_news20_mix(seed=seed, n_first=n_per_novel),
        lambda: board_cnae9_text(seed=seed, n_first=min(n_per_novel, 120)),
        lambda: board_amazon_category(seed=seed, n_first=min(n_per_novel, 300)),
    ):
        try:
            b = fn()
            if b is not None:
                boards.append(_expand_board_to_full_live(b, seed=seed, n_per_novel=n_per_novel))
                print(
                    "[text-full] loaded", b["name"],
                    "n_live=", len(b["y_live"]),
                    "novel_per≈", {
                        str(c): int((np.asarray(b["y_live"]) == c).sum())
                        for c in b.get("novel_labels") or []
                    },
                    flush=True,
                )
        except Exception as exc:  # noqa: BLE001
            print("[text-full] failed:", exc, flush=True)
    return boards


def run_text_full_acc_benchmark(
    *,
    seed: int = SEED,
    n_per_novel: int = _FULL_N_PER_NOVEL,
    schemes: Sequence[str] = (
        "uniform", "sample_w", "feature_w", "both", "ewma", "ewma_feat",
    ),
    alpha_feat: float = 3.0,
    ewma_lam: float = _DEFAULT_EWMA_LAM,
    run_incremental: bool = True,
) -> Dict:
    """Load FULL text LIVE boards → overall / known / novel accuracy compare."""
    boards = load_text_boards_full(seed=seed, n_per_novel=n_per_novel)
    rows = []
    incr = []
    for b in boards:
        print("[text-full] eval", b.get("name"), "schemes", list(schemes), flush=True)
        rows.append(
            evaluate_board_full_acc(
                b, seed=seed, schemes=schemes, alpha_feat=alpha_feat,
                ewma_lam=ewma_lam,
            )
        )
        if run_incremental:
            print("[text-full] incremental", b.get("name"), flush=True)
            incr.append(
                evaluate_board_incremental(
                    b, seed=seed, ewma_lam=ewma_lam, alpha_feat=alpha_feat,
                    schemes=("uniform", "ewma", "feature_w", "ewma_feat"),
                )
            )

    def _mean(sch, key):
        vals = []
        for r in rows:
            d = r.get("deltas_vs_uniform", {}).get(sch, {})
            v = d.get(key)
            if v is not None and v == v:
                vals.append(v)
        return float(np.mean(vals)) if vals else None

    def _mean_abs(sch, key):
        vals = []
        for r in rows:
            v = r.get("schemes", {}).get(sch, {}).get(key)
            if v is not None and v == v:
                vals.append(float(v))
        return float(np.mean(vals)) if vals else None

    sch_list = list(schemes)
    return {
        "domain": "text_multiclass_only",
        "tabular_excluded": True,
        "loss": LOSS_NAME,
        "kpi": "FULL held-out LIVE acc / novel_acc / known_acc (not first-n subset)",
        "seed": int(seed),
        "n_per_novel": int(n_per_novel),
        "alpha_feat": float(alpha_feat),
        "ewma_lam": float(ewma_lam),
        "ewma_formula": "w_i = exp(-lam * n_idx_novel_category)",
        "boards": rows,
        "incremental": incr,
        "n_boards": len(rows),
        "mean_acc": {sch: _mean_abs(sch, "acc") for sch in sch_list},
        "mean_novel_acc": {sch: _mean_abs(sch, "novel_acc") for sch in sch_list},
        "mean_delta_acc": {
            sch: _mean(sch, "delta_acc") for sch in sch_list if sch != "uniform"
        },
        "mean_delta_novel_acc": {
            sch: _mean(sch, "delta_novel_acc") for sch in sch_list if sch != "uniform"
        },
        "takeaway": (
            "EWMA w=exp(-λ·n_idx) focuses CE on earliest novel-category hits; "
            "compare incremental lift vs uniform / feature_w on FULL held-out LIVE acc."
        ),
    }
