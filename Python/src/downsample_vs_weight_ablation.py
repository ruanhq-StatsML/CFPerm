"""Downsample vs FSDS sample-weight ablation on category-shift boards.

Question
--------
Do we need label-based downsampling of known LIVE rows to cold-start novel
classes, or do FSDS-localized sample weights already deliver (and exceed)
that lift?

Schemes (same REF init + LIVE train/eval split)
----------------------------------------------
  uniform              full LIVE train, flat weights
  ds_known             keep all novel; downsample known to n_known=n_novel (1:1)
  ds_random            random downsample LIVE to same n as ds_known (n-matched control)
  ds_mem               membership hard downsample: keep top-|novel| by m (no labels)
  sample_w_fused       soft weights after private-live fuse (full n)
  subset_w_fused       hard subset weights after fuse (full n)
  both_fused           soft×hard fused (full n)
  ds_known+sample_f    downsample known, then soft-fused weights on remaining

Primary metric: novel_acc. Incremental value of weights = best_weight − ds_known.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from sklearn.preprocessing import StandardScaler

from fsds_cfperm_weight_benchmark import (
    SEED,
    _gini,
    _metrics,
    board_cancer,
    board_cnae9,
    board_covtype,
    board_har,
    board_isolet,
    board_jvowel,
    board_mice,
    board_news20b,
    board_wine,
    cfperm_feature_reject,
)
from weighted_cold_start import (
    attribute_rf_domain,
    cold_start_partial_fit,
    localize_subset,
    make_sgd,
    weights_from_attribution,
)


def load_ablation_boards(*, seed: int = SEED, quick: bool = True) -> List[Dict]:
    fns = [
        lambda: board_wine(seed=seed),
        lambda: board_cancer(seed=seed),
        lambda: board_har(seed=seed, quick=quick),
        lambda: board_mice(seed=seed),
        lambda: board_cnae9(seed=seed, quick=quick),
        lambda: board_jvowel(seed=seed, quick=quick),
        lambda: board_covtype(seed=seed, quick=quick),
        lambda: board_news20b(seed=seed, max_features=64 if quick else 96),
        lambda: board_isolet(seed=seed, quick=quick),
    ]
    out = []
    for fn in fns:
        try:
            out.append(fn())
        except Exception as exc:  # noqa: BLE001
            print("[ds-ablate] skip:", exc, flush=True)
    return out


def _fit_eval(
    X0s, y0, Xtr, ytr, Xev, yev, X0s_ref, y0_ref,
    *,
    known, novel, classes, feat, sample_weight=None, seed=SEED, n_passes=3,
) -> Dict:
    clf = make_sgd(seed=seed)
    clf.partial_fit(X0s, y0, classes=classes)
    Xuse = np.asarray(Xtr, float).copy()
    Xuse[:, feat] *= 1.08
    sw = None if sample_weight is None else np.asarray(sample_weight, float)
    if sw is not None:
        sw = sw * (len(sw) / max(sw.sum(), 1e-12))
    clf = cold_start_partial_fit(
        clf, Xuse, np.asarray(ytr).ravel(),
        sample_weight=sw, classes=classes, n_passes=n_passes,
    )
    m_live = _metrics(clf, Xev, yev, known=known, novel=novel)
    m_ref = _metrics(clf, X0s_ref, y0_ref, known=known, novel=None)
    return {
        **m_live,
        "ref_acc": m_ref["acc"],
        "forgetting": float(1.0 - m_ref["acc"]),
        "n_train": int(len(ytr)),
        "novel_frac_train": float(np.isin(ytr, novel).mean()) if novel is not None else float("nan"),
        "weight_gini": 0.0 if sw is None else float(_gini(sw)),
    }


def run_board_downsample_ablation(
    board: Dict,
    *,
    seed: int = SEED,
    top_k: int = 6,
    lam: float = 2.5,
    subset_q: float = 0.6,
    subset_boost: float = 3.0,
) -> Dict:
    X0 = np.asarray(board["X_ref"], float)
    y0 = np.asarray(board["y_ref"]).ravel()
    X1 = np.asarray(board["X_live"], float)
    y1 = np.asarray(board["y_live"]).ravel()
    known = board.get("known_labels")
    novel = board.get("novel_labels")
    classes = np.unique(np.concatenate([y0, y1]))

    sc = StandardScaler().fit(X0)
    X0s, X1s = sc.transform(X0), sc.transform(X1)

    vimp = attribute_rf_domain(X0s, X1s, seed=seed, n_estimators=60)
    mem, feat = localize_subset(X0s, X1s, vimp, top_k=top_k)
    mem_f = mem.copy()
    if novel is not None:
        mem_f = np.maximum(mem_f, np.isin(y1, novel).astype(float))

    wt_sf = weights_from_attribution(
        mem_f, feat, X1s.shape[1], lam=lam, subset_q=subset_q,
        subset_boost=subset_boost, mode="sample", temper=False,
    )
    wt_uf = weights_from_attribution(
        mem_f, feat, X1s.shape[1], lam=lam, subset_q=subset_q,
        subset_boost=subset_boost, mode="subset", temper=False,
    )
    wt_bf = weights_from_attribution(
        mem_f, feat, X1s.shape[1], lam=lam, subset_q=subset_q,
        subset_boost=subset_boost, mode="fused", temper=True,
    )

    rng = np.random.default_rng(seed + 11)
    idx = rng.permutation(len(y1))
    n_ev = max(20, len(y1) // 3)
    ev, tr = idx[:n_ev], idx[n_ev:]
    Xtr, ytr, mtr, mftr = X1s[tr], y1[tr], mem[tr], mem_f[tr]
    Xev, yev = X1s[ev], y1[ev]

    # --- downsample indices on TRAIN only ---
    is_nov = np.isin(ytr, novel) if novel is not None else np.zeros(len(ytr), bool)
    is_kn = np.isin(ytr, known) if known is not None else ~is_nov
    n_nov = int(is_nov.sum())
    n_kn = int(is_kn.sum())
    # 1:1 known:novel (keep all novel)
    n_kn_keep = min(n_kn, max(n_nov, 1))
    kn_idx = np.flatnonzero(is_kn)
    nov_idx = np.flatnonzero(is_nov)
    kn_keep = rng.choice(kn_idx, size=n_kn_keep, replace=False) if n_kn_keep else np.array([], dtype=int)
    ds_known_idx = np.concatenate([nov_idx, kn_keep])
    rng.shuffle(ds_known_idx)

    # n-matched random downsample
    n_ds = len(ds_known_idx)
    ds_rand_idx = rng.choice(len(ytr), size=min(n_ds, len(ytr)), replace=False)

    # membership downsample: keep top-n_nov by membership among all train,
    # plus enough to match n_ds total by filling with next-highest m
    # (smart downsample without using novel labels at selection time)
    order = np.argsort(-mtr)
    ds_mem_idx = order[:n_ds]

    # weights restricted to ds_known subset (for combo scheme)
    # map full-train fused sample weights onto ds_known rows
    w_sf_full = wt_sf["sample_weight"][tr]
    w_sf_on_ds = w_sf_full[ds_known_idx]
    w_sf_on_ds = w_sf_on_ds * (len(w_sf_on_ds) / max(w_sf_on_ds.sum(), 1e-12))

    schemes_idx_w = {
        "uniform": (np.arange(len(ytr)), None),
        "ds_known": (ds_known_idx, None),
        "ds_random": (ds_rand_idx, None),
        "ds_mem": (ds_mem_idx, None),
        "sample_w_fused": (np.arange(len(ytr)), w_sf_full),
        "subset_w_fused": (np.arange(len(ytr)), wt_uf["sample_weight"][tr]),
        "both_fused": (np.arange(len(ytr)), wt_bf["sample_weight"][tr]),
        "ds_known+sample_f": (ds_known_idx, w_sf_on_ds),
    }

    rows = {}
    for name, (sel, w) in schemes_idx_w.items():
        sel = np.asarray(sel, int)
        rows[name] = _fit_eval(
            X0s, y0, Xtr[sel], ytr[sel], Xev, yev, X0s, y0,
            known=known, novel=novel, classes=classes, feat=feat,
            sample_weight=w, seed=seed,
        )
        rows[name]["scheme"] = name

    u = rows["uniform"]["novel_acc"]
    ds = rows["ds_known"]["novel_acc"]

    def lift(a, b):
        return float(rows[a]["novel_acc"] - rows[b]["novel_acc"])

    best_w_name = max(
        ["sample_w_fused", "subset_w_fused", "both_fused"],
        key=lambda k: rows[k]["novel_acc"],
    )
    return {
        "dataset": board["name"],
        "note": board.get("note"),
        "n_ref": int(len(y0)),
        "n_live": int(len(y1)),
        "n_train_full": int(len(ytr)),
        "n_train_ds": int(n_ds),
        "n_novel_train": n_nov,
        "novel_frac_full": float(n_nov / max(len(ytr), 1)),
        "novel_frac_ds_known": float(rows["ds_known"]["novel_frac_train"]),
        "vimp_top_mass": float(vimp[feat].sum()),
        "mean_membership": float(mem.mean()),
        "schemes": rows,
        "delta_novel_ds_known_vs_uniform": lift("ds_known", "uniform"),
        "delta_novel_ds_random_vs_uniform": lift("ds_random", "uniform"),
        "delta_novel_ds_mem_vs_uniform": lift("ds_mem", "uniform"),
        "delta_novel_best_weight_vs_uniform": float(rows[best_w_name]["novel_acc"] - u),
        "delta_novel_best_weight_vs_ds_known": float(rows[best_w_name]["novel_acc"] - ds),
        "delta_novel_both_fused_vs_ds_known": lift("both_fused", "ds_known"),
        "delta_novel_sample_fused_vs_ds_known": lift("sample_w_fused", "ds_known"),
        "delta_novel_combo_vs_ds_known": lift("ds_known+sample_f", "ds_known"),
        "best_weight_scheme": best_w_name,
        "need_downsample": bool(lift("ds_known", "uniform") > 0.02 and lift("ds_known", "uniform") >= float(rows[best_w_name]["novel_acc"] - u) - 1e-9),
        # weights beat downsample if best weight > ds_known + margin
        "weights_beat_downsample": bool(rows[best_w_name]["novel_acc"] > ds + 1e-12),
    }


def run_downsample_ablation(*, seed: int = SEED, quick: bool = True) -> Dict:
    boards = load_ablation_boards(seed=seed, quick=quick)
    results = []
    for b in boards:
        print("[ds-ablate] board", b["name"], flush=True)
        results.append(run_board_downsample_ablation(b, seed=seed, top_k=6 if quick else 8))

    def mean_key(k):
        vals = [r[k] for r in results if r.get(k) is not None]
        return float(np.mean(vals)) if vals else None

    return {
        "boards": results,
        "n_boards": len(results),
        "mean_delta_ds_known_vs_uniform": mean_key("delta_novel_ds_known_vs_uniform"),
        "mean_delta_ds_random_vs_uniform": mean_key("delta_novel_ds_random_vs_uniform"),
        "mean_delta_ds_mem_vs_uniform": mean_key("delta_novel_ds_mem_vs_uniform"),
        "mean_delta_best_weight_vs_uniform": mean_key("delta_novel_best_weight_vs_uniform"),
        "mean_delta_best_weight_vs_ds_known": mean_key("delta_novel_best_weight_vs_ds_known"),
        "mean_delta_combo_vs_ds_known": mean_key("delta_novel_combo_vs_ds_known"),
        "n_weights_beat_ds": int(sum(1 for r in results if r.get("weights_beat_downsample"))),
        "takeaway": (
            "Label downsample (ds_known) is a blunt rebalance; FSDS sample weights "
            "are selective reweighting. Incremental value = best_weight − ds_known. "
            "Need downsample only when weights cannot recover novel mass and "
            "label balance itself is the bottleneck."
        ),
    }
