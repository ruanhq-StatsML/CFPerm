"""Representation-first cold-start under category shift.

Four steps (locked)
-------------------
  1. Detect + attribute
       RFPerm-style error rise → FSDS RF-domain VIMP → subspace S + membership.
       Role: *focus*, not classify.
  2. Representation separation  ← the "instant learn" step
       On S, margin-push high-membership rows away from known-class prototypes.
       Do **not** rush the classification head while novel is sparse.
  3. Dynamic head expand
       When novel buffer ≥ N* *or* the attributed footprint stabilizes,
       joint-fit new nodes with novel buffer + known REF replay.
  4. Optional TTA (cautious)
       Entropy-min only on coefs tied to S; freeze the rest.

Weights (sample_w / fused) still apply at step 3. Step 2 is geometry, not IW.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from sklearn.linear_model import SGDClassifier

from adaptation_opportunity import (
    adaptation_opportunity,
    opportunity_from_bound_labels,
)
from weighted_cold_start import (
    SEED,
    _membership_conc,
    attribute_rf_domain,
    boundary_check,
    choose_weight_mode,
    cold_start_partial_fit,
    detect_error_rise,
    localize_subset,
    mix_known_replay,
    weights_from_attribution,
)


# ---------------------------------------------------------------------------
# State
# ---------------------------------------------------------------------------


@dataclass
class ColdStartState:
    """Cross-batch memory for separation → expand."""

    novel_X: List[np.ndarray] = field(default_factory=list)
    novel_y: List[np.ndarray] = field(default_factory=list)
    novel_w: List[np.ndarray] = field(default_factory=list)
    feat_idx: Optional[np.ndarray] = None
    last_feat_idx: Optional[np.ndarray] = None
    # persistent mean displacement on S for high-mem rows (applied at infer)
    sep_offset: Optional[np.ndarray] = None  # shape (p,)
    n_sep_batches: int = 0
    last_vimp_top_mass: float = 0.0
    footprint_stable: bool = False
    expanded_labels: set = field(default_factory=set)
    total_novel_seen: int = 0
    last_opportunity: Optional[Dict] = None

    def buffer_size(self) -> int:
        return int(sum(len(y) for y in self.novel_y))

    def clear_buffer(self) -> None:
        self.novel_X.clear()
        self.novel_y.clear()
        self.novel_w.clear()


# ---------------------------------------------------------------------------
# Step 1 helpers already live in weighted_cold_start
# ---------------------------------------------------------------------------


def known_prototypes(
    X_ref: np.ndarray,
    y_ref: np.ndarray,
    feat_idx: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Class means on subspace S. Returns (proto[C, |S|], labels[C])."""
    X = np.asarray(X_ref, float)
    y = np.asarray(y_ref).ravel()
    feat = np.asarray(feat_idx, int)
    labs = np.unique(y)
    protos = np.zeros((len(labs), len(feat)), dtype=float)
    for i, c in enumerate(labs):
        sl = y == c
        protos[i] = X[sl][:, feat].mean(0) if sl.any() else 0.0
    return protos, labs


# ---------------------------------------------------------------------------
# Step 2 — representation separation (key)
# ---------------------------------------------------------------------------


def margin_separate(
    X_live: np.ndarray,
    membership: np.ndarray,
    prototypes_S: np.ndarray,
    feat_idx: np.ndarray,
    *,
    margin: float = 1.25,
    lr: float = 0.35,
    mem_thr: float = 0.55,
    n_steps: int = 3,
    max_step: float = 1.5,
    adaptive_margin: bool = True,
) -> Dict:
    """Margin push on S: high-mem rows move away from nearest known prototype.

    Goal is *not* "what is the new class" — only "stop collapsing into old
    centers". Sparse novel is enough because the loss is a geometric margin.

    Margin is adaptive by default: ``max(margin, 1.15 * median(d_highmem))``
    so already-separated novel still gets a soft repulsive lift on S.
    """
    X = np.asarray(X_live, float).copy()
    mem = np.asarray(membership, float).ravel()
    feat = np.asarray(feat_idx, int)
    proto = np.asarray(prototypes_S, float)
    if proto.size == 0 or len(feat) == 0:
        return {
            "X": X, "n_pushed": 0, "gap_before": 0.0, "gap_after": 0.0,
            "offset_S": np.zeros(len(feat)), "mask": np.zeros(len(mem), dtype=bool),
            "margin_eff": float(margin),
        }

    mask = mem >= float(mem_thr)
    if not mask.any():
        # fall back to top-quartile membership
        thr = float(np.quantile(mem, 0.75))
        mask = mem >= thr

    Xs = X[:, feat]
    # nearest known proto distance
    d0 = np.full(len(X), np.inf)
    nearest = np.zeros(len(X), dtype=int)
    for i in np.flatnonzero(mask):
        dif = Xs[i] - proto
        dist = np.linalg.norm(dif, axis=1)
        j = int(np.argmin(dist))
        d0[i] = float(dist[j])
        nearest[i] = j
    gap_before = float(np.mean(d0[mask])) if mask.any() else 0.0

    # adaptive target: push beyond the current high-mem median gap
    margin_eff = float(margin)
    if adaptive_margin and mask.any():
        med = float(np.median(d0[mask]))
        margin_eff = max(float(margin), 1.15 * med)

    offset_acc = np.zeros(len(feat), dtype=float)
    n_pushed = 0
    for _ in range(max(1, int(n_steps))):
        for i in np.flatnonzero(mask):
            j = int(nearest[i])
            v = Xs[i] - proto[j]
            d = float(np.linalg.norm(v) + 1e-8)
            # hard margin when inside; soft repulsion when outside but high-mem
            if d < margin_eff:
                strength = float(lr) * (margin_eff - d) * (0.5 + 0.5 * mem[i])
            else:
                # soft tail: still nudge so expand sees a clearer S gap
                strength = float(lr) * 0.15 * mem[i] * margin_eff / d
                if strength < 1e-4:
                    continue
            step = (v / d) * strength
            sn = float(np.linalg.norm(step))
            if sn > float(max_step):
                step *= float(max_step) / sn
            Xs[i] = Xs[i] + step
            offset_acc += step
            n_pushed += 1
            dif = Xs[i] - proto
            dist = np.linalg.norm(dif, axis=1)
            nearest[i] = int(np.argmin(dist))

    X[:, feat] = Xs
    d1 = np.full(len(X), np.inf)
    for i in np.flatnonzero(mask):
        dif = Xs[i] - proto
        d1[i] = float(np.linalg.norm(dif, axis=1).min())
    gap_after = float(np.mean(d1[mask])) if mask.any() else 0.0
    offset_S = offset_acc / max(n_pushed, 1)
    return {
        "X": X,
        "n_pushed": int(n_pushed),
        "gap_before": gap_before,
        "gap_after": gap_after,
        "gap_lift": float(gap_after - gap_before),
        "offset_S": offset_S,
        "mask": mask,
        "margin_eff": margin_eff,
        "mem_thr_eff": float(mem[mask].min()) if mask.any() else float(mem_thr),
    }


def apply_sep_offset(
    X: np.ndarray,
    membership: np.ndarray,
    state: ColdStartState,
    *,
    mem_thr: float = 0.55,
) -> np.ndarray:
    """Apply persistent S-offset to high-mem rows at inference."""
    if state.sep_offset is None or state.feat_idx is None:
        return np.asarray(X, float)
    X = np.asarray(X, float).copy()
    mem = np.asarray(membership, float).ravel()
    mask = mem >= float(mem_thr)
    if mask.any():
        X[mask] = X[mask] + state.sep_offset
    return X


# ---------------------------------------------------------------------------
# Step 3 — dynamic head expand
# ---------------------------------------------------------------------------


def footprint_stability(
    state: ColdStartState,
    vimp_top_mass: float,
    membership: np.ndarray,
    *,
    mass_tol: float = 0.08,
    min_conc: float = 0.50,
) -> bool:
    """True when attributed mass is stable and membership is concentrated."""
    conc = _membership_conc(membership)
    if state.n_sep_batches < 1:
        return False
    stable_mass = abs(float(vimp_top_mass) - float(state.last_vimp_top_mass)) <= float(mass_tol)
    return bool(stable_mass and conc >= float(min_conc))


def should_expand_head(
    state: ColdStartState,
    private_live: Sequence,
    membership: np.ndarray,
    vimp_top_mass: float,
    *,
    expand_min: int = 16,
    n_private_batch: int = 0,
    new_private: Optional[Sequence] = None,
) -> Tuple[bool, str]:
    """Expand when buffer ≥ N*, footprint stable, or new private-live is ready."""
    n_buf = state.buffer_size()
    if n_buf >= int(expand_min):
        return True, "buffer_ge_N"
    # new private labels in this batch with enough rows → expand now
    # (C-GCD stages introduce fresh novel ids; waiting forever starves the head)
    new_priv = [c for c in (new_private or []) if int(c) not in state.expanded_labels]
    if new_priv and int(n_private_batch) >= max(4, int(expand_min) // 2):
        return True, "new_private_batch"
    if private_live and footprint_stability(state, vimp_top_mass, membership):
        if n_buf >= max(4, int(expand_min) // 2):
            return True, "footprint_stable"
    return False, "wait"


def buffer_novel(
    state: ColdStartState,
    X: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    private_live: Sequence,
    membership: np.ndarray,
    *,
    mem_thr: float = 0.55,
) -> int:
    """Park novel / high-mem rows until head expand."""
    y = np.asarray(y).ravel()
    mem = np.asarray(membership, float).ravel()
    priv = set(int(c) for c in private_live)
    mask = np.isin(y, list(priv)) if priv else (mem >= float(mem_thr))
    if not mask.any():
        return 0
    state.novel_X.append(np.asarray(X, float)[mask])
    state.novel_y.append(y[mask])
    state.novel_w.append(np.asarray(w, float).ravel()[mask])
    state.total_novel_seen += int(mask.sum())
    return int(mask.sum())


def expand_and_fit(
    clf: SGDClassifier,
    state: ColdStartState,
    X_ref: np.ndarray,
    y_ref: np.ndarray,
    *,
    classes: Sequence[int],
    feature_mask: Optional[np.ndarray] = None,
    replay_frac: float = 0.35,
    seed: int = SEED,
    n_passes: int = 4,
) -> Tuple[SGDClassifier, Dict]:
    """Joint-fit buffered novel + known replay. Completes 'learn' after sep."""
    if state.buffer_size() == 0:
        return clf, {"fitted": False, "n_novel": 0, "n_replay": 0}
    Xn = np.vstack(state.novel_X)
    yn = np.concatenate(state.novel_y)
    wn = np.concatenate(state.novel_w)
    wn = wn * (len(wn) / max(wn.sum(), 1e-12))
    X_fit, y_fit, w_fit = mix_known_replay(
        Xn, yn, wn, X_ref, y_ref, replay_frac=replay_frac, seed=seed,
    )
    n_replay = int(len(y_fit) - len(yn))
    clf = cold_start_partial_fit(
        clf, X_fit, y_fit,
        sample_weight=w_fit,
        feature_mask=feature_mask,
        classes=classes,
        n_passes=n_passes,
    )
    state.expanded_labels |= set(int(c) for c in np.unique(yn))
    info = {
        "fitted": True,
        "n_novel": int(len(yn)),
        "n_replay": n_replay,
        "expanded_labels": sorted(state.expanded_labels),
    }
    state.clear_buffer()
    return clf, info


# ---------------------------------------------------------------------------
# Step 4 — optional TTA on subspace S
# ---------------------------------------------------------------------------


def tta_entropy_subspace(
    clf: SGDClassifier,
    X: np.ndarray,
    feat_idx: np.ndarray,
    *,
    steps: int = 3,
    lr: float = 5e-3,
    max_delta: float = 0.05,
) -> Dict:
    """TENT-style entropy min, but only on coef columns in S (rest frozen).

    Tabular SGD has no BatchNorm; we nudge ``coef_[:, S]`` + intercept lightly.
    """
    if not hasattr(clf, "coef_"):
        return {"applied": False, "reason": "no_coef"}
    X = np.asarray(X, float)
    feat = np.asarray(feat_idx, int)
    if len(feat) == 0 or len(X) == 0:
        return {"applied": False, "reason": "empty"}

    coef0 = clf.coef_.copy()
    intercept0 = clf.intercept_.copy()
    ent_before = _mean_entropy(clf, X)

    for _ in range(max(1, int(steps))):
        # finite-difference / score-function style nudge toward lower entropy
        proba = clf.predict_proba(X)
        # dH/dz ≈ (log p + 1) for each class; push logits of max-prob up a bit
        # on S features only via a small gradient on coef
        p = np.clip(proba, 1e-6, 1.0)
        # soft target = sharpen distribution
        sharp = p ** 2
        sharp = sharp / sharp.sum(1, keepdims=True)
        # residual (sharp - p): encourage confidence
        resid = sharp - p  # (n, C)
        # grad_coef[c, j] ≈ - mean_i resid[i,c] * X[i,j]
        grad = -(resid.T @ X) / max(len(X), 1)  # (C, p)
        step = float(lr) * grad
        # only touch S
        delta = np.zeros_like(clf.coef_)
        delta[:, feat] = step[:, feat]
        # clip
        nrm = float(np.linalg.norm(delta))
        if nrm > float(max_delta):
            delta *= float(max_delta) / nrm
        clf.coef_ = clf.coef_ + delta
        # tiny intercept nudge toward sharp
        clf.intercept_ = clf.intercept_ + float(lr) * resid.mean(0)

    ent_after = _mean_entropy(clf, X)
    if ent_after > ent_before * 1.05:
        # rolled worse — revert
        clf.coef_ = coef0
        clf.intercept_ = intercept0
        return {
            "applied": False, "reason": "entropy_rose",
            "ent_before": ent_before, "ent_after": ent_after,
        }
    return {
        "applied": True,
        "ent_before": ent_before,
        "ent_after": ent_after,
        "coef_delta_norm": float(np.linalg.norm(clf.coef_ - coef0)),
    }


def _mean_entropy(clf, X) -> float:
    p = np.clip(clf.predict_proba(X), 1e-6, 1.0)
    return float((-p * np.log(p)).sum(1).mean())


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------


def adapt_batch_repr_sep(
    clf: SGDClassifier,
    X_ref: np.ndarray,
    y_ref: np.ndarray,
    X_live: np.ndarray,
    y_live: np.ndarray,
    state: Optional[ColdStartState] = None,
    *,
    seed: int = SEED,
    top_k: int = 8,
    lam: float = 2.0,
    min_jaccard: float = 0.20,
    detect_thresh: float = 1.15,
    force_adapt: bool = False,
    use_private_labels: bool = True,
    weight_mode: str = "auto",
    subset_q: float = 0.65,
    subset_boost: float = 2.2,
    temper: bool = True,
    replay_frac: float = 0.35,
    expand_min: int = 16,
    margin: float = 1.25,
    sep_lr: float = 0.35,
    sep_steps: int = 3,
    mem_thr: float = 0.55,
    do_tta: bool = False,
    tta_steps: int = 2,
) -> Dict:
    """Full 4-step cold-start. Returns clf + state + diagnostics."""
    state = state or ColdStartState()
    y_ref = np.asarray(y_ref).ravel()
    y_live = np.asarray(y_live).ravel()
    X_ref = np.asarray(X_ref, float)
    X_live = np.asarray(X_live, float)

    # --- Step 1: detect + attribute ---
    gate = detect_error_rise(clf, X_ref, y_ref, X_live, y_live, thresh=detect_thresh)
    bound = boundary_check(y_ref, y_live, min_jaccard=min_jaccard)
    vimp = attribute_rf_domain(X_ref, X_live, seed=seed)
    mem_raw, feat = localize_subset(X_ref, X_live, vimp, top_k=top_k)
    mem = mem_raw.copy()
    if use_private_labels and bound["private_live"]:
        mem = np.maximum(mem, np.isin(y_live, bound["private_live"]).astype(float))
    vimp_top = float(vimp[feat].sum())

    mode_req = str(weight_mode).lower()
    mode_eff = choose_weight_mode(
        mode_req, mem, gate_trigger=bool(gate["trigger"] or force_adapt),
    )
    wt = weights_from_attribution(
        mem, feat, X_live.shape[1], lam=lam,
        subset_q=subset_q, subset_boost=subset_boost, mode=mode_eff,
        temper=temper,
    )

    phase = "skip"
    reason = "no_trigger"
    sep_info: Dict = {}
    expand_info: Dict = {"fitted": False}
    tta_info: Dict = {"applied": False}
    n_buffered = 0
    X_work = X_live
    opp: Dict = {}

    incomparable = bound["decision"] == "incomparable" and not force_adapt
    trigger = bool(gate["trigger"] or force_adapt)
    priv = list(bound.get("private_live") or [])
    n_priv_batch = opportunity_from_bound_labels(y_live, priv)
    protos, _labs = known_prototypes(X_ref, y_ref, feat)

    # Opportunity: can we adapt to new categories *this* batch?
    # O_persist discounts betting that S/offset still help next batch.
    opp = adaptation_opportunity(
        error_ratio=float(gate.get("ratio", 1.0)),
        vimp=vimp,
        feat_idx=feat,
        membership=mem,
        n_private=n_priv_batch,
        expand_min=expand_min,
        detect_thresh=detect_thresh,
        mem_thr=mem_thr,
        margin=margin,
        X_live=X_live,
        prototypes_S=protos,
        prev_vimp_top_mass=(state.last_vimp_top_mass if state.n_sep_batches else None),
        prev_feat=state.last_feat_idx,
        force=bool(force_adapt),
    )
    state.last_opportunity = {
        "O_local": opp["O_local"], "O_persist": opp["O_persist"],
        "action": opp["action"], "read": opp["read"],
    }

    if incomparable:
        reason = "incomparable"
        phase = "skip"
        # unknown / non-overlapping support: drop persisted geometry
        state.sep_offset = None
    elif trigger and opp["do_sep"]:
        # --- Step 2: representation separation (only if O_local says so) ---
        # Push only private-live / high-mem candidates — never move confirmed known
        # rows (that would fight retention).
        sep_mem = mem.copy()
        if priv:
            known_rows = ~np.isin(y_live, priv)
            sep_mem[known_rows] = 0.0  # do not push known recirculation
        sep_info = margin_separate(
            X_live, sep_mem, protos, feat,
            margin=margin, lr=sep_lr, mem_thr=mem_thr, n_steps=sep_steps,
        )
        X_work = sep_info["X"]
        off = np.zeros(X_live.shape[1], dtype=float)
        off[feat] = sep_info["offset_S"]
        # Persist offset only when O_persist is high; else short memory / replace.
        if opp["keep_offset"]:
            if state.sep_offset is None:
                state.sep_offset = off
            else:
                state.sep_offset = 0.7 * state.sep_offset + 0.3 * off
        else:
            state.sep_offset = off  # batch-local; next batch re-attributes
        state.last_feat_idx = np.asarray(state.feat_idx, int) if state.feat_idx is not None else None
        state.feat_idx = np.asarray(feat, int)
        state.footprint_stable = footprint_stability(state, vimp_top, mem)
        state.last_vimp_top_mass = vimp_top
        state.n_sep_batches += 1
        phase = "separate"
        reason = f"repr_sep:O={opp['O_local']:.2f}"

        # park novel / high-mem for later expand
        n_buffered = buffer_novel(
            state, X_work, y_live, wt["sample_weight"],
            priv, mem, mem_thr=mem_thr,
        )

        # protect known: light replay on *raw* known rows (no sep distortion)
        known_mask = ~np.isin(y_live, priv) if priv else np.ones(len(y_live), dtype=bool)
        if known_mask.any() and known_mask.sum() >= 4:
            classes = getattr(clf, "classes_", None)
            if classes is None:
                classes = np.unique(np.concatenate([y_ref, y_live]))
            Xk = X_live[known_mask]
            yk = y_live[known_mask]
            wk = np.ones(len(yk), dtype=float)
            X_fit, y_fit, w_fit = mix_known_replay(
                Xk, yk, wk, X_ref, y_ref,
                replay_frac=min(float(replay_frac), 0.40), seed=seed + 5,
            )
            clf = cold_start_partial_fit(
                clf, X_fit, y_fit, sample_weight=w_fit,
                feature_mask=wt["feature_mask"], classes=classes, n_passes=1,
            )

        # --- Step 3: dynamic head expand ---
        new_priv = [c for c in priv if int(c) not in state.expanded_labels]
        do_exp, exp_why = should_expand_head(
            state, priv, mem, vimp_top,
            expand_min=expand_min,
            n_private_batch=n_priv_batch,
            new_private=new_priv,
        )
        # opportunity must also clear expand bar (mass × O_local)
        if do_exp and not opp["do_expand_ok"] and not force_adapt:
            do_exp, exp_why = False, "O_local_below_exp"
        if do_exp:
            classes = getattr(clf, "classes_", None)
            if classes is None:
                classes = np.unique(np.concatenate([y_ref, y_live]))
            rf = min(float(replay_frac), 0.25) if priv else float(replay_frac)
            clf, expand_info = expand_and_fit(
                clf, state, X_ref, y_ref,
                classes=classes, feature_mask=wt["feature_mask"],
                replay_frac=rf, seed=seed + 7, n_passes=4,
            )
            expand_info["why"] = exp_why
            phase = "expand"
            reason = f"repr_sep+expand:{exp_why}|O={opp['O_local']:.2f}"
        else:
            reason = f"repr_sep+buffer(n={state.buffer_size()})|O={opp['O_local']:.2f}"

        # --- Step 4: optional TTA on S ---
        if do_tta:
            tta_info = tta_entropy_subspace(
                clf, X_work, feat, steps=tta_steps, lr=5e-3,
            )
            if tta_info.get("applied"):
                reason = reason + "+tta_S"
    elif trigger and not opp["do_sep"]:
        reason = f"hold:O_local={opp['O_local']:.2f}<tau_sep"
        phase = "hold"
        # low opportunity: do not carry stale geometry into unknown next batch
        if not opp["keep_offset"]:
            state.sep_offset = None

    return {
        "clf": clf,
        "state": state,
        "phase": phase,
        "gate": gate,
        "boundary": bound,
        "vimp": vimp,
        "membership": mem,
        "membership_raw": mem_raw,
        "feat_idx": feat,
        "weights": wt,
        "opportunity": opp,
        "adapted": phase in ("separate", "expand"),
        "reason": reason,
        "vimp_top_mass": vimp_top,
        "use_private_labels": bool(use_private_labels),
        "weight_mode": wt["weight_mode"],
        "weight_mode_requested": mode_req,
        "subset_frac": wt["subset_frac"],
        "mean_membership": float(np.mean(mem)),
        "membership_conc": _membership_conc(mem),
        "n_replay": int(expand_info.get("n_replay", 0)),
        "n_buffered": int(n_buffered),
        "buffer_size": state.buffer_size(),
        "sep": {
            "n_pushed": sep_info.get("n_pushed"),
            "gap_before": sep_info.get("gap_before"),
            "gap_after": sep_info.get("gap_after"),
            "gap_lift": sep_info.get("gap_lift"),
        } if sep_info else {},
        "expand": expand_info,
        "tta": tta_info,
        "X_separated": X_work,
    }


def predict_with_state(
    clf: SGDClassifier,
    X: np.ndarray,
    state: ColdStartState,
    *,
    membership: Optional[np.ndarray] = None,
    mem_thr: float = 0.55,
) -> np.ndarray:
    """Inference: apply persistent S-offset on high-mem rows, then clf.predict."""
    X = np.asarray(X, float)
    if membership is None or state.sep_offset is None:
        return clf.predict(X)
    Xs = apply_sep_offset(X, membership, state, mem_thr=mem_thr)
    return clf.predict(Xs)
