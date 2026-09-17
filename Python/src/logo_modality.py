"""Leave-one-group-out localization for modalities.

Groups are a partition of X (video / audio / text, or any feature blocks).
LOGO shares are localization proxies, not a unique causal decomposition of
the serving drop, and not a Shapley split of concept vs covariate.

Two ratio layers, both vs the same T=0 reference batch:
  1. serving loss — MSE if Y is continuous, Brier if Y is binary
  2. PO-risk vs MMD²(X_new, X_ref) mix, per group

T is the batch label, not a treatment whose CATE we identify.
PO-risk is a hop of P(Y|X). MMD is P(X) vs D_ref, never pairwise history.

No online-bootstrap.
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
from sklearn.linear_model import LogisticRegression, Ridge

from streaming_po_risk import (
    ACTION_FREEZE,
    ACTION_KEEP,
    ACTION_TRICKY,
    ACTION_WATCH,
    ACTION_XSHIFT,
    _is_binary,
    batch_mse,
    large_deviation,
    mmd_vs_reference,
    po_mse_action,
    rbf_bandwidth,
    streaming_po_and_mse,
)

TOWER_FULL = "full_train"
TOWER_INFER = "infer_only"
TOWER_TOP = "train_top"
TOWER_STEM = "train_stem"
TOWER_FREEZE = "freeze_tower"
FUSION_FULL = "fusion_full"
FUSION_INFER = "fusion_infer"
FUSION_HEAD = "fusion_head"

MIX_LOUD = 0.6


def as_groups(groups: Mapping[str, Sequence[int]] | Mapping[str, slice]) -> dict[str, np.ndarray]:
    """name → column index array."""
    out = {}
    for name, idx in groups.items():
        if isinstance(idx, slice):
            if idx.start is None or idx.stop is None:
                raise ValueError(f"group {name!r} slice needs start and stop")
            out[str(name)] = np.arange(int(idx.start), int(idx.stop))
        else:
            out[str(name)] = np.asarray(idx, dtype=int)
    return out


def drop_group(X, groups: Mapping[str, np.ndarray], leave_out: str) -> np.ndarray:
    """X without one group's columns. Remaining groups stay in name order."""
    X = np.asarray(X, dtype=float)
    keep = []
    for name, idx in groups.items():
        if name == leave_out:
            continue
        keep.append(X[:, np.asarray(idx, dtype=int)])
    if not keep:
        raise ValueError("dropping the only group leaves an empty X")
    return np.hstack(keep)


def serving_mu(model, X, binary: bool) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    if binary:
        proba = model.predict_proba(X)
        if proba.shape[1] == 2:
            return proba[:, 1]
        return proba
    return np.asarray(model.predict(X), dtype=float).ravel()


def fit_serving(X, Y, seed: int = 2026):
    """Cheap serving readout for LOGO. Not the PO-risk RF, not a ViT."""
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    if _is_binary(Y):
        m = LogisticRegression(max_iter=400, random_state=int(seed))
        m.fit(X, Y.astype(int))
        return m, True
        m = Ridge(alpha=1.0)
        m.fit(X, Y)
        return m, False


def brier_or_mse(Y, mu, binary: bool) -> float:
    """Layer-1 loss: Brier when Y is binary, MSE otherwise."""
    Y = np.asarray(Y, dtype=float)
    mu = np.asarray(mu, dtype=float)
    if binary and mu.ndim == 2:
        K = mu.shape[1]
        onehot = np.eye(K)[np.clip(Y.ravel().astype(int), 0, K - 1)]
        return float(np.mean(np.sum((mu - onehot) ** 2, axis=1)))
    return batch_mse(Y, mu)


def _pos_share(deltas: Mapping[str, float]) -> dict[str, float]:
    """Normalize ReLU(delta) to a simplex. Negative deltas stay visible as 0 share."""
    names = list(deltas)
    raw = np.array([max(float(deltas[n]), 0.0) for n in names], dtype=float)
    s = float(raw.sum())
    if s <= 1e-15:
        return {n: 0.0 for n in names}
    return {n: float(raw[i] / s) for i, n in enumerate(names)}


def two_layer_ratios(pi_loss, pi_po, pi_mmd) -> dict[str, dict]:
    """Per group: layer-1 π_loss, layer-2 mix of π_PO vs π_MMD."""
    names = list(pi_loss)
    out = {}
    for g in names:
        po = float(pi_po.get(g, 0.0))
        mmd = float(pi_mmd.get(g, 0.0))
        mix_den = po + mmd
        quiet = mix_den <= 1e-12 and float(pi_loss.get(g, 0.0)) <= 1e-12
        out[g] = {
            "pi_loss": float(pi_loss.get(g, 0.0)),
            "pi_po": po,
            "pi_mmd": mmd,
            "mix_po": 0.0 if mix_den <= 1e-12 else po / mix_den,
            "mix_mmd": 0.0 if mix_den <= 1e-12 else mmd / mix_den,
            "quiet": bool(quiet),
        }
    return out


def dominant_group(ratios: Mapping[str, Mapping]) -> str | None:
    """Argmax of π_loss, else of π_PO+π_MMD. None if everything is quiet."""
    if not ratios:
        return None
    loss = {g: float(r["pi_loss"]) for g, r in ratios.items()}
    if max(loss.values()) > 1e-12:
        return max(loss, key=loss.get)
    shift = {g: float(r["pi_po"]) + float(r["pi_mmd"]) for g, r in ratios.items()}
    if max(shift.values()) > 1e-12:
        return max(shift, key=shift.get)
    return None


def plan_next_batch(global_action: str, ratios: Mapping[str, Mapping]) -> dict:
    """Heuristic for the next incoming batch. Not an autostop.

    keep        → full train every tower + fusion
    watch       → infer only (PO hopped, serving still pays rent)
    x_shift     → stem-adapt towers whose mix is MMD-loud; do not freeze as concept
    freeze      → train_top on PO-loud towers; freeze quiet towers
    tricky      → infer only; no freeze-from-layer
    """
    action = str(global_action)
    towers = {}
    if action == ACTION_KEEP:
        for g in ratios:
            towers[g] = {
                "tower": TOWER_FULL,
                "i_star": None,
                "reason": "global quiet: full backprop",
            }
        fusion = FUSION_FULL
        update = "full_train"
    elif action in (ACTION_WATCH, ACTION_TRICKY):
        why = (
            "PO hopped, MSE holds — watch"
            if action == ACTION_WATCH
            else "MSE broke, PO and MMD quiet — tricky"
        )
        for g in ratios:
            towers[g] = {"tower": TOWER_INFER, "i_star": None, "reason": why}
        fusion = FUSION_INFER
        update = "infer_only"
    elif action == ACTION_XSHIFT:
        for g, r in ratios.items():
            if r["quiet"]:
                towers[g] = {
                    "tower": TOWER_FREEZE,
                    "i_star": None,
                    "reason": "quiet modality under X-shift",
                }
            elif r["mix_mmd"] >= MIX_LOUD:
                towers[g] = {
                    "tower": TOWER_STEM,
                    "i_star": 0,
                    "reason": "MMD-loud: stem-adapt, not concept freeze",
                }
            else:
                towers[g] = {
                    "tower": TOWER_INFER,
                    "i_star": None,
                    "reason": "X-shift but this tower is not MMD-loud",
                }
        fusion = FUSION_INFER
        update = "stem_adapt"
    elif action == ACTION_FREEZE:
        for g, r in ratios.items():
            if r["quiet"]:
                towers[g] = {
                    "tower": TOWER_FREEZE,
                    "i_star": None,
                    "reason": "quiet tower: spend no gradient",
                }
            elif r["mix_po"] >= MIX_LOUD:
                towers[g] = {
                    "tower": TOWER_TOP,
                    "i_star": "top",
                    "reason": "PO-loud: freeze bottom, train top + fusion",
                }
            else:
                towers[g] = {
                    "tower": TOWER_INFER,
                    "i_star": None,
                    "reason": "both-broken globally, this tower is not PO-loud",
                }
        fusion = FUSION_HEAD
        update = "selective_freeze"
    else:
        for g in ratios:
            towers[g] = {"tower": TOWER_INFER, "i_star": None, "reason": f"unknown action {action}"}
        fusion = FUSION_INFER
        update = "infer_only"
    return {
        "global_action": action,
        "update": update,
        "fusion": fusion,
        "towers": towers,
        "dominant": dominant_group(ratios),
    }


def logo_batch(
    X_ref,
    Y_ref,
    X_new,
    Y_new,
    groups,
    *,
    seed: int = 2026,
    po_broken: bool | None = None,
    mse_broken: bool | None = None,
    mmd_broken: bool | None = None,
    po_base: float | None = None,
    mse_base: float | None = None,
    mmd_base: float | None = None,
    ratio: float = 2.0,
) -> dict:
    """One incoming batch: full metrics + LOGO deltas + two-layer ratios + plan."""
    groups = as_groups(groups)
    X_ref = np.asarray(X_ref, dtype=float)
    X_new = np.asarray(X_new, dtype=float)
    Y_ref = np.asarray(Y_ref, dtype=float).ravel()
    Y_new = np.asarray(Y_new, dtype=float).ravel()

    sigma = rbf_bandwidth(X_ref, seed=seed)
    po_full, mse_full = streaming_po_and_mse(X_ref, Y_ref, X_new, Y_new, seed=seed)
    mmd_full = mmd_vs_reference(X_ref, X_new, sigma=sigma, seed=seed)

    serve, binary = fit_serving(X_ref, Y_ref, seed=seed)
    loss_name = "brier" if binary else "mse"
    loss_full = brier_or_mse(Y_new, serving_mu(serve, X_new, binary), binary)

    d_loss, d_po, d_mmd = {}, {}, {}
    dropped = {}
    for g in groups:
        Xr = drop_group(X_ref, groups, g)
        Xn = drop_group(X_new, groups, g)
        po_g, _ = streaming_po_and_mse(Xr, Y_ref, Xn, Y_new, seed=seed)
        sig_g = rbf_bandwidth(Xr, seed=seed)
        mmd_g = mmd_vs_reference(Xr, Xn, sigma=sig_g, seed=seed)
        serve_g, _ = fit_serving(Xr, Y_ref, seed=seed)
        loss_g = brier_or_mse(Y_new, serving_mu(serve_g, Xn, binary), binary)
        d_loss[g] = float(loss_full - loss_g)
        d_po[g] = float(po_full - po_g)
        d_mmd[g] = float(mmd_full - mmd_g)
        dropped[g] = {"po": float(po_g), "mmd": float(mmd_g), "loss": float(loss_g)}

    pi_loss, pi_po, pi_mmd = _pos_share(d_loss), _pos_share(d_po), _pos_share(d_mmd)
    ratios = two_layer_ratios(pi_loss, pi_po, pi_mmd)

    if po_base is not None:
        po_broken = large_deviation(po_full, po_base, ratio=ratio)
    elif po_broken is None:
        po_broken = False
    if mse_base is not None:
        mse_broken = large_deviation(mse_full, mse_base, ratio=ratio)
    elif mse_broken is None:
        mse_broken = False
    if mmd_base is not None:
        mmd_broken = large_deviation(mmd_full, mmd_base, ratio=ratio)
    elif mmd_broken is None:
        mmd_broken = False
    action = po_mse_action(bool(po_broken), bool(mse_broken), bool(mmd_broken))
    plan = plan_next_batch(action, ratios)
    return {
        "loss_name": loss_name,
        "binary": bool(binary),
        "po_full": float(po_full),
        "mmd_full": float(mmd_full),
        "mse_full": float(mse_full),
        "loss_full": float(loss_full),
        "delta_loss": d_loss,
        "delta_po": d_po,
        "delta_mmd": d_mmd,
        "dropped": dropped,
        "pi_loss": pi_loss,
        "pi_po": pi_po,
        "pi_mmd": pi_mmd,
        "ratios": ratios,
        "plan": plan,
        "global_action": action,
    }


def subset_excess(
    X_ref,
    Y_ref,
    X_new,
    Y_new,
    labels_new,
    *,
    labels_ref=None,
    seed: int = 2026,
    min_n: int = 20,
) -> list[dict]:
    """Slice the new batch; rank excess serving loss vs D_ref. Localization, not CATE."""
    X_ref = np.asarray(X_ref, dtype=float)
    X_new = np.asarray(X_new, dtype=float)
    Y_ref = np.asarray(Y_ref, dtype=float).ravel()
    Y_new = np.asarray(Y_new, dtype=float).ravel()
    labels_new = np.asarray(labels_new)
    serve, binary = fit_serving(X_ref, Y_ref, seed=seed)
    loss_ref = brier_or_mse(Y_ref, serving_mu(serve, X_ref, binary), binary)
    sigma = rbf_bandwidth(X_ref, seed=seed)
    labels_ref = None if labels_ref is None else np.asarray(labels_ref)
    rows = []
    for lab in sorted(set(labels_new.tolist()), key=str):
        mask = labels_new == lab
        n = int(mask.sum())
        if n < int(min_n):
            continue
        loss_s = brier_or_mse(Y_new[mask], serving_mu(serve, X_new[mask], binary), binary)
        if labels_ref is not None and int((labels_ref == lab).sum()) >= int(min_n):
            Xr = X_ref[labels_ref == lab]
            sig_s = rbf_bandwidth(Xr, seed=seed)
            mmd_s = mmd_vs_reference(Xr, X_new[mask], sigma=sig_s, seed=seed)
        else:
            mmd_s = mmd_vs_reference(X_ref, X_new[mask], sigma=sigma, seed=seed)
        po_s, mse_s = streaming_po_and_mse(X_ref, Y_ref, X_new[mask], Y_new[mask], seed=seed)
        rows.append(
            {
                "subset": lab,
                "n": n,
                "loss": float(loss_s),
                "excess_loss": float(loss_s - loss_ref),
                "mse": float(mse_s),
                "po": float(po_s),
                "mmd": float(mmd_s),
            }
        )
    rows.sort(key=lambda r: r["excess_loss"], reverse=True)
    return rows


def cumulative_regret(loss_policy, loss_oracle) -> np.ndarray:
    """R_t = Σ_{s≤t} (L_s(policy) − L_s(oracle)). Oracle must be feasible, not a unique decomp."""
    a = np.asarray(loss_policy, dtype=float).ravel()
    b = np.asarray(loss_oracle, dtype=float).ravel()
    if a.shape != b.shape:
        raise ValueError("policy and oracle series must align")
    return np.cumsum(a - b)
