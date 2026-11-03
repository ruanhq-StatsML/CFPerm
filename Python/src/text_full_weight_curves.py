"""Text multiclass: feature / sample weights on FULL eval (v3 soft^γ gate).

No EWMA. OOD quantification (RF-domain VIMP + membership) → CE adapt via
partial_fit. Cold-start empowerment:

  soft^γ (default γ=0.5): boost *= ((1 - n_idx/N_cap)_+)^γ   ← linger
  soft^1 / soft^2        : linear / fast taper (ablations)
  hard / always          : cut / no-gate ablations

KPI = overall accuracy + mean novel OVA-AUC on full held-out LIVE.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

from text_feature_sample_weights import load_text_boards_full
from weighted_cold_start import (
    SEED,
    attribute_rf_domain,
    cold_start_partial_fit,
    localize_subset,
    make_sgd,
    weights_from_attribution,
)

LOSS_NAME = "CE (log_loss)"
# Primary: soft with linger γ=0.5. Ablations: γ=1, γ=2, hard, always.
SCHEMES = (
    "uniform", "feature_w", "sample_w", "both",
    "both_g1", "both_g2", "both_hard", "both_always",
)
COLD_START_CAP = 200
DEFAULT_SOFT_GAMMA = 0.5  # v3: linger — slower taper than linear
DEFAULT_BUDGETS = (20, 40, 60, 80, 100, 150, 200, 250, 300, 400)
GATE_MODES = ("soft", "hard", "always")


def feature_level_scale(
    X: np.ndarray,
    vimp: np.ndarray,
    *,
    alpha: float = 3.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Feature-level weights: scale each column by (1 + α · vimp_j)."""
    v = np.asarray(vimp, float).ravel()
    v = v / (v.sum() + 1e-12)
    scale = 1.0 + float(alpha) * (v / (v.max() + 1e-12))
    X = np.asarray(X, float) * scale.reshape(1, -1)
    return X, scale


def per_category_arrival_index(y: np.ndarray, cats: Sequence[int]) -> np.ndarray:
    """Within-category arrival index (0,1,2,…) for labels in ``cats``; -1 else."""
    y = np.asarray(y).ravel()
    cat_set = {int(c) for c in cats}
    idx = np.full(len(y), -1, dtype=int)
    counts: Dict[int, int] = {}
    for i, yi in enumerate(y):
        c = int(yi)
        if c not in cat_set:
            continue
        k = counts.get(c, 0)
        idx[i] = k
        counts[c] = k + 1
    return idx


def cold_start_decay(
    y: np.ndarray,
    novel: Sequence[int],
    *,
    cap: int = COLD_START_CAP,
    mode: str = "soft",
    gamma: float = DEFAULT_SOFT_GAMMA,
    arrival_idx: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Per-row cold-start intensity in [0, 1].

    soft   : ((1 - n_idx/cap)_+)^γ on novel rows (γ<1 lingers, γ>1 cuts faster)
    hard   : 1 if n_idx < cap else 0 on novel rows
    always : 1 on all novel rows (no gate; ablation)
    Known rows stay 0 (no novel empowerment).
    """
    y = np.asarray(y).ravel()
    if arrival_idx is None:
        arrival_idx = per_category_arrival_index(y, novel)
    arrival_idx = np.asarray(arrival_idx, int).ravel()
    decay = np.zeros(len(y), dtype=float)
    novel_m = arrival_idx >= 0
    if not novel_m.any():
        return decay
    mode = str(mode).lower()
    if mode == "always":
        decay[novel_m] = 1.0
    elif mode == "hard":
        decay[novel_m] = (arrival_idx[novel_m] < int(cap)).astype(float)
    else:  # soft^γ
        t = np.clip(
            1.0 - arrival_idx[novel_m].astype(float) / float(max(cap, 1)),
            0.0, 1.0,
        )
        decay[novel_m] = np.power(t, float(max(gamma, 1e-6)))
    return decay


def _parse_scheme(scheme: str) -> Tuple[str, str, float]:
    """Return (base_scheme, gate_mode, soft_gamma).

    base ∈ {uniform, feature_w, sample_w, both}
    gate ∈ {soft, hard, always}
    gamma: power on soft taper; default DEFAULT_SOFT_GAMMA (0.5 linger).
    Suffixes: _g05 / _g1 / _g2, _hard, _always, _soft.
    """
    s = str(scheme).lower().strip()
    gamma = float(DEFAULT_SOFT_GAMMA)
    if s == "uniform":
        return "uniform", "soft", gamma
    # gamma suffixes first (both_g1, both_g05, feature_w_g2, …)
    for tag, g in (("_g05", 0.5), ("_g0.5", 0.5), ("_g1", 1.0), ("_g2", 2.0)):
        if s.endswith(tag):
            return s[: -len(tag)], "soft", float(g)
    for gate in ("always", "hard", "soft"):
        suf = "_" + gate
        if s.endswith(suf):
            return s[: -len(suf)], gate, gamma
    # bare feature_w / sample_w / both → soft with linger γ
    return s, "soft", gamma


def _ova_auc_mean(
    y_true: np.ndarray,
    proba: np.ndarray,
    classes: np.ndarray,
    novel: Sequence[int],
) -> float:
    y_true = np.asarray(y_true).ravel()
    classes = list(np.asarray(classes).ravel())
    aucs = []
    for c in novel:
        c = int(c)
        if c not in classes:
            continue
        col = classes.index(c)
        binary = (y_true == c).astype(int)
        if binary.min() == binary.max():
            continue
        aucs.append(float(roc_auc_score(binary, proba[:, col])))
    return float(np.mean(aucs)) if aucs else float("nan")


def fit_full_scheme(
    X_ref: np.ndarray,
    y_ref: np.ndarray,
    X_tr: np.ndarray,
    y_tr: np.ndarray,
    *,
    scheme: str,
    seed: int,
    classes: np.ndarray,
    novel: Sequence[int],
    top_k: int = 16,
    lam: float = 2.5,
    feat_alpha: float = 3.0,
    cold_start_cap: int = COLD_START_CAP,
    arrival_idx: Optional[np.ndarray] = None,
    gate_mode: Optional[str] = None,
) -> Tuple[SGDClassifier, StandardScaler, Optional[np.ndarray], Dict]:
    """Fit CE head with OOD-driven weights under a cold-start gate."""
    base, gate, soft_gamma = _parse_scheme(scheme)
    if gate_mode is not None:
        gate = str(gate_mode).lower()
    sc = StandardScaler().fit(X_ref)
    X0 = sc.transform(X_ref)
    X1 = sc.transform(X_tr)
    y_tr = np.asarray(y_tr).ravel()
    classes = np.asarray(classes, int)
    clf = make_sgd(seed=seed)
    clf.partial_fit(X0, y_ref, classes=classes)

    if arrival_idx is None:
        arrival_idx = per_category_arrival_index(y_tr, novel)
    decay = cold_start_decay(
        y_tr, novel, cap=cold_start_cap, mode=gate, gamma=soft_gamma,
        arrival_idx=arrival_idx,
    )
    mean_decay = float(decay[decay > 0].mean()) if (decay > 0).any() else 0.0
    use_weights = base != "uniform" and mean_decay > 1e-8

    scale = None
    w = np.ones(len(y_tr), dtype=float)
    info = {
        "scheme_base": base,
        "gate_mode": gate,
        "soft_gamma": float(soft_gamma),
        "cold_start_cap": int(cold_start_cap),
        "mean_decay": mean_decay,
        "n_decay_pos": int((decay > 0).sum()),
        "n_train": int(len(y_tr)),
        "n_cold_start": int((decay > 0).sum()),
        "weights_active": False,
        "feature_w_active": False,
        "sample_w_active": False,
    }

    if use_weights:
        # attribute on rows that still carry cold-start mass when possible
        attr_m = decay > 0.05
        X_attr = X1[attr_m] if int(attr_m.sum()) >= 8 else X1
        vimp = attribute_rf_domain(
            X0, X_attr if len(X_attr) >= 8 else X1, seed=seed, n_estimators=80,
        )
        mem, feat = localize_subset(X0, X1, vimp, top_k=top_k)
        mem_f = mem.copy()
        # fuse membership with cold-start intensity (not a hard 0/1 only)
        mem_f = np.maximum(mem_f, decay)

        if base in ("feature_w", "both"):
            # taper feature stretch with mean cold-start intensity
            alpha_eff = float(feat_alpha) * max(mean_decay, 0.05)
            X1, scale = feature_level_scale(X1, vimp, alpha=alpha_eff)
            X0s = X0 * scale.reshape(1, -1)
            clf = make_sgd(seed=seed)
            clf.partial_fit(X0s, y_ref, classes=classes)
            info["feature_w_active"] = True
            info["feat_alpha_eff"] = alpha_eff

        if base in ("sample_w", "both"):
            wt = weights_from_attribution(
                mem_f, feat, X1.shape[1], lam=lam, mode="fused", temper=True,
            )
            w_raw = np.asarray(wt["sample_weight"], float).ravel()
            # excess over 1, tapered by per-row decay (warm novel → 1)
            excess = np.maximum(w_raw - 1.0, 0.0)
            w = 1.0 + excess * decay
            info["sample_w_active"] = True
        elif base == "feature_w":
            w = np.ones(len(y_tr), dtype=float)

        info["weights_active"] = True

    w = w * (len(w) / max(w.sum(), 1e-12))
    clf = cold_start_partial_fit(
        clf, X1, y_tr, sample_weight=w, classes=classes, n_passes=4,
    )
    return clf, sc, scale, info


def predict_full(
    clf: SGDClassifier,
    sc: StandardScaler,
    scale: Optional[np.ndarray],
    X: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    Xs = sc.transform(np.asarray(X, float))
    if scale is not None:
        Xs = Xs * scale.reshape(1, -1)
    return clf.predict(Xs), clf.predict_proba(Xs)


def eval_full_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    proba: np.ndarray,
    classes: np.ndarray,
    *,
    known: Sequence[int],
    novel: Sequence[int],
) -> Dict:
    y_true = np.asarray(y_true).ravel()
    out = {
        "acc": float(accuracy_score(y_true, y_pred)),
        "ova_auc_novel_mean": _ova_auc_mean(y_true, proba, classes, novel),
    }
    if known:
        m = np.isin(y_true, list(known))
        out["known_acc"] = float(np.mean(y_pred[m] == y_true[m])) if m.any() else float("nan")
    if novel:
        m = np.isin(y_true, list(novel))
        out["novel_acc"] = float(np.mean(y_pred[m] == y_true[m])) if m.any() else float("nan")
    return out


def _take_per_category(
    y_pool: np.ndarray,
    pool_idx: np.ndarray,
    cats: Sequence[int],
    n_per: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Take up to ``n_per`` pool rows per category (shuffled within cat)."""
    taken = []
    for c in cats:
        hits = pool_idx[np.asarray(y_pool[pool_idx] == int(c))]
        if hits.size == 0:
            continue
        hits = hits.copy()
        rng.shuffle(hits)
        taken.append(hits[: min(int(n_per), len(hits))])
    if not taken:
        return np.asarray([], dtype=int)
    return np.concatenate(taken)


def run_board_full_curves(
    board: Dict,
    *,
    seed: int = SEED,
    budgets: Optional[Sequence[int]] = None,
    schemes: Sequence[str] = SCHEMES,
    feat_alpha: float = 3.0,
    cold_start_cap: int = COLD_START_CAP,
) -> Dict:
    """Learning curves on FULL eval: budget = novel samples / category."""
    budgets = list(budgets or DEFAULT_BUDGETS)
    known = list(board["known_labels"])
    novel = list(board["novel_labels"])
    X_ref = np.asarray(board["X_ref"], float)
    y_ref = np.asarray(board["y_ref"]).ravel()
    X_live = np.asarray(board["X_live"], float)
    y_live = np.asarray(board["y_live"]).ravel()

    rng = np.random.default_rng(int(seed))
    idx = rng.permutation(len(y_live))
    n_ev = max(80, len(y_live) // 3)
    ev, pool = idx[:n_ev], idx[n_ev:]
    X_ev, y_ev = X_live[ev], y_live[ev]
    X_pool, y_pool = X_live[pool], y_live[pool]

    known_pool = np.flatnonzero(np.isin(y_pool, known))
    novel_pool = np.flatnonzero(np.isin(y_pool, novel))
    rng.shuffle(known_pool)

    classes = np.unique(np.concatenate([y_ref, y_live]))
    curves = {sch: [] for sch in schemes}

    for b in budgets:
        b = int(b)
        take_n = _take_per_category(y_pool, novel_pool, novel, b, rng)
        take_k = known_pool[: min(max(40, 2 * len(take_n)), len(known_pool))]
        tr = np.concatenate([take_k, take_n]) if len(take_n) else take_k
        if tr.size < 8:
            continue
        X_tr, y_tr = X_pool[tr], y_pool[tr]
        arrival = np.full(len(y_tr), -1, dtype=int)
        n_k = len(take_k)
        if len(take_n):
            arrival[n_k:] = per_category_arrival_index(y_tr[n_k:], novel)

        for sch in schemes:
            clf, sc, scale, info = fit_full_scheme(
                X_ref, y_ref, X_tr, y_tr,
                scheme=sch, seed=seed + b, classes=classes, novel=novel,
                feat_alpha=feat_alpha, cold_start_cap=cold_start_cap,
                arrival_idx=arrival,
            )
            pred, proba = predict_full(clf, sc, scale, X_ev)
            m = eval_full_metrics(
                y_ev, pred, proba, classes, known=known, novel=novel,
            )
            curves[sch].append({
                "budget": b,
                "budget_unit": "novel_per_category",
                "n_novel_train": int(len(take_n)),
                "n_known_train": int(len(take_k)),
                "gate_mode": info.get("gate_mode"),
                "soft_gamma": info.get("soft_gamma"),
                "mean_decay": info.get("mean_decay"),
                "n_cold_start": info["n_cold_start"],
                "weights_active": info["weights_active"],
                "feature_w_active": info["feature_w_active"],
                "sample_w_active": info["sample_w_active"],
                **m,
            })

    endpoint = {}
    for sch in schemes:
        endpoint[sch] = curves[sch][-1] if curves[sch] else {}

    at_cap = {}
    for sch in schemes:
        pts = [p for p in curves[sch] if int(p["budget"]) == int(cold_start_cap)]
        at_cap[sch] = pts[0] if pts else {}

    def _delta(metric, table):
        out = {}
        base = table.get("uniform", {}).get(metric)
        for sch in schemes:
            if sch == "uniform" or base is None or base != base:
                continue
            v = table.get(sch, {}).get(metric)
            out[sch] = float(v - base) if v is not None and v == v else None
        return out

    return {
        "dataset": board["name"],
        "note": board.get("note"),
        "loss": LOSS_NAME,
        "ewma": False,
        "cold_start_cap": int(cold_start_cap),
        "default_gate": "soft",
        "default_soft_gamma": float(DEFAULT_SOFT_GAMMA),
        "cold_start_rule": (
            f"soft^γ decay ((1 - n_idx/{cold_start_cap})_+)^γ, "
            f"default γ={DEFAULT_SOFT_GAMMA} (linger); "
            "both_g1 / both_g2 / both_hard / both_always = ablations"
        ),
        "eval": "full_holdout_LIVE",
        "n_eval": int(len(y_ev)),
        "budgets": list(budgets),
        "budget_unit": "novel_per_category",
        "schemes": list(schemes),
        "curves": curves,
        "endpoint": endpoint,
        "at_cold_start_cap": at_cap,
        "delta_acc_vs_uniform": _delta("acc", endpoint),
        "delta_ova_auc_vs_uniform": _delta("ova_auc_novel_mean", endpoint),
        "delta_novel_acc_vs_uniform": _delta("novel_acc", endpoint),
        "delta_acc_at_cap_vs_uniform": _delta("acc", at_cap),
        "delta_ova_at_cap_vs_uniform": _delta("ova_auc_novel_mean", at_cap),
    }


def run_text_full_weight_curves(
    *,
    seed: int = SEED,
    n_first_board: int = 400,
    budgets: Optional[Sequence[int]] = None,
    feat_alpha: float = 3.0,
    cold_start_cap: int = COLD_START_CAP,
    schemes: Sequence[str] = SCHEMES,
) -> Dict:
    """Text-only boards; full-eval curves with soft cold-start gate (v2)."""
    boards = load_text_boards_full(seed=seed, n_per_novel=n_first_board)
    rows = []
    for b in boards:
        print(
            "[full-curves]", b["name"], "n_live", len(b["y_live"]),
            "cold_start_cap", cold_start_cap, "schemes", list(schemes),
            flush=True,
        )
        rows.append(
            run_board_full_curves(
                b, seed=seed, budgets=budgets, feat_alpha=feat_alpha,
                cold_start_cap=cold_start_cap, schemes=schemes,
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
        "domain": "text_multiclass_only",
        "tabular_excluded": True,
        "loss": LOSS_NAME,
        "ewma": False,
        "cold_start_cap": int(cold_start_cap),
        "default_gate": "soft",
        "default_soft_gamma": float(DEFAULT_SOFT_GAMMA),
        "primary_kpi": "full_holdout overall acc + mean novel OVA-AUC",
        "seed": int(seed),
        "schemes": list(schemes),
        "boards": rows,
        "mean_delta_acc_feature_w": _mean_delta("delta_acc_vs_uniform", "feature_w"),
        "mean_delta_acc_sample_w": _mean_delta("delta_acc_vs_uniform", "sample_w"),
        "mean_delta_acc_both": _mean_delta("delta_acc_vs_uniform", "both"),
        "mean_delta_acc_both_g1": _mean_delta("delta_acc_vs_uniform", "both_g1"),
        "mean_delta_acc_both_g2": _mean_delta("delta_acc_vs_uniform", "both_g2"),
        "mean_delta_acc_both_hard": _mean_delta("delta_acc_vs_uniform", "both_hard"),
        "mean_delta_acc_both_always": _mean_delta("delta_acc_vs_uniform", "both_always"),
        "mean_delta_ova_feature_w": _mean_delta("delta_ova_auc_vs_uniform", "feature_w"),
        "mean_delta_ova_sample_w": _mean_delta("delta_ova_auc_vs_uniform", "sample_w"),
        "mean_delta_ova_both": _mean_delta("delta_ova_auc_vs_uniform", "both"),
        "takeaway": (
            f"v3 soft^γ gate: boost *= ((1-n_idx/{cold_start_cap})_+)^γ "
            f"with default γ={DEFAULT_SOFT_GAMMA} (linger). "
            "OOD VIMP → feature/sample weights via partial_fit. FULL holdout. No EWMA."
        ),
    }
