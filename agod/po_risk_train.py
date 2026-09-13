"""PO-risk metric versions → next-step training actuators.

Philosophy
----------
PO / residual-concept is an *external* attributor for P(Y|X) drift.
Several metric logics turn PO-risk (optionally + MMD / proto / ΔPO /
noise-gates) into a simplex ``α``, then map ``α`` into *next-window*
training knobs — towers FWD stay on:

  LR_m · step budget · hard BWD freeze · stacking prior.

Metric versions (sensor → α)
----------------------------
  equal         uniform (no PO)
  po_soft       Softmax(PO / τ)
  po_minus_cov  Softmax((PO − λ·MMD) / τ)
  po_gated      drift-vs-noise gate on PO, then Softmax
  po_proto      Softmax((PO·(1+proto) − λ·MMD) / τ)
  po_delta      Softmax((EMA(PO) + γ·ΔPO) / τ)     anticipatory
  po_budget     floor + (1−floor)·Softmax(PO)      risk budget w/ floor
  po_next       EMA(α) ⊕ Softmax(ΔPO)              next-step forecast
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping, Sequence

import numpy as np

from .lr_controller import alpha_to_lr, softmax_scores, z_norm
from .next_step import predict_next_alpha
from .smooth_router import DriftNoiseGateConfig, drift_vs_noise_gate

METRIC_VERSIONS = (
    "equal",
    "po_soft",
    "po_minus_cov",
    "po_gated",
    "po_proto",
    "po_delta",
    "po_budget",
    "po_next",
)

ACTUATORS = ("lr_mult", "step_alloc", "freeze_mask", "stack_prior")


@dataclass
class PORiskMetricConfig:
    tau: float = 0.30
    lam_cov: float = 1.0
    gamma_delta: float = 0.75
    ema_po: float = 0.55
    budget_floor: float = 0.15
    gate: DriftNoiseGateConfig = field(default_factory=DriftNoiseGateConfig)


@dataclass
class NextStepActuatorConfig:
    beta_lr: float = 0.10
    freeze_theta: float = 0.12
    total_steps: int = 28
    min_steps: int = 2


def _as(
    mods: Sequence[str], d: Mapping[str, float] | None, default: float = 0.0
) -> dict[str, float]:
    d = d or {}
    return {m: float(d.get(m, default)) for m in mods}


def _ema_dict(
    prev: Mapping[str, float] | None,
    cur: Mapping[str, float],
    mods: Sequence[str],
    ema: float,
) -> dict[str, float]:
    if prev is None:
        return {m: float(cur[m]) for m in mods}
    return {
        m: float(ema) * float(prev[m]) + (1.0 - float(ema)) * float(cur[m])
        for m in mods
    }


def metric_to_alpha(
    version: str,
    mods: Sequence[str],
    *,
    po: Mapping[str, float],
    mmd: Mapping[str, float] | None = None,
    proto: Mapping[str, float] | None = None,
    uni_acc: Mapping[str, float] | None = None,
    vimp: Mapping[str, float] | None = None,
    po_prev: Mapping[str, float] | None = None,
    po_ema: Mapping[str, float] | None = None,
    alpha_hist: Sequence[Mapping[str, float]] | None = None,
    cfg: PORiskMetricConfig | None = None,
) -> dict:
    """Map sensors → α under a named PO-risk metric version."""
    cfg = cfg or PORiskMetricConfig()
    mods = list(mods)
    po_a = _as(mods, po)
    mmd_a = _as(mods, mmd)
    proto_a = _as(mods, proto)
    uni = _as(mods, uni_acc, 0.55)
    vim = _as(mods, vimp, 0.0)
    diag: dict = {"version": version}

    if version == "equal":
        score = {m: 0.0 for m in mods}
        alpha = {m: 1.0 / len(mods) for m in mods}
        return {"alpha": alpha, "score": score, "diag": diag}

    if version == "po_soft":
        score = dict(po_a)
        alpha = softmax_scores(score, mods, cfg.tau)
        return {"alpha": alpha, "score": score, "diag": diag}

    if version == "po_minus_cov":
        score = {m: po_a[m] - float(cfg.lam_cov) * mmd_a[m] for m in mods}
        alpha = softmax_scores(score, mods, cfg.tau)
        return {"alpha": alpha, "score": score, "diag": diag}

    if version == "po_gated":
        gate = drift_vs_noise_gate(po_a, uni, vim, mods, cfg=cfg.gate)
        po_g = gate["po_gated"]
        score = {m: float(po_g[m]) for m in mods}
        alpha = softmax_scores(score, mods, cfg.tau)
        diag["gate"] = {
            "is_signal": gate["is_signal"],
            "adapter_only": gate["adapter_only"],
            "frac_signal": gate["frac_signal"],
            "reasons": gate["reasons"],
        }
        return {"alpha": alpha, "score": score, "diag": diag}

    if version == "po_proto":
        score = {
            m: po_a[m] * (1.0 + proto_a[m]) - float(cfg.lam_cov) * mmd_a[m]
            for m in mods
        }
        alpha = softmax_scores(score, mods, cfg.tau)
        return {"alpha": alpha, "score": score, "diag": diag}

    if version == "po_delta":
        base = _ema_dict(po_ema, po_a, mods, cfg.ema_po)
        if po_prev is None:
            delta = {m: 0.0 for m in mods}
        else:
            delta = {m: po_a[m] - float(po_prev[m]) for m in mods}
        score = {m: base[m] + float(cfg.gamma_delta) * delta[m] for m in mods}
        alpha = softmax_scores(score, mods, cfg.tau)
        diag["po_ema"] = base
        diag["delta"] = delta
        return {"alpha": alpha, "score": score, "diag": diag}

    if version == "po_budget":
        soft = softmax_scores(po_a, mods, cfg.tau)
        floor = min(float(cfg.budget_floor), 0.9 / len(mods))
        rem = 1.0 - floor * len(mods)
        alpha = {m: floor + rem * soft[m] for m in mods}
        s = sum(alpha.values())
        alpha = {m: alpha[m] / s for m in mods}
        score = dict(po_a)
        diag["floor"] = floor
        return {"alpha": alpha, "score": score, "diag": diag}

    if version == "po_next":
        hist = list(alpha_hist or [])
        if not hist:
            hist = [softmax_scores(po_a, mods, cfg.tau)]
        sensors_now = {"g": po_a}
        sensors_prev = {"g": _as(mods, po_prev)} if po_prev is not None else None
        alpha = predict_next_alpha(
            hist,
            mods,
            sensors_now=sensors_now,
            sensors_prev=sensors_prev,
            ema=cfg.ema_po,
            tau=cfg.tau,
            sensor_key="g",
        )
        score = dict(po_a)
        return {"alpha": alpha, "score": score, "diag": diag}

    raise ValueError(
        f"unknown PO-risk metric version {version!r}; expected {METRIC_VERSIONS}"
    )


def next_step_actuators(
    alpha: Mapping[str, float],
    mods: Sequence[str],
    *,
    cfg: NextStepActuatorConfig | None = None,
) -> dict:
    """Map α → next-window training adjustments (all four actuators)."""
    cfg = cfg or NextStepActuatorConfig()
    mods = list(mods)
    a = z_norm(alpha, mods)

    lr = alpha_to_lr(a, mods, beta=cfg.beta_lr)

    raw = np.array([max(float(a[m]), 0.0) for m in mods], float)
    if raw.sum() <= 1e-12:
        raw = np.ones(len(mods)) / len(mods)
    else:
        raw = raw / raw.sum()
    steps = np.maximum(np.round(raw * cfg.total_steps), cfg.min_steps).astype(int)
    while int(steps.sum()) > cfg.total_steps and (steps > cfg.min_steps).any():
        j = int(np.argmax(steps))
        if steps[j] > cfg.min_steps:
            steps[j] -= 1
        else:
            break
    while int(steps.sum()) < cfg.total_steps:
        j = int(np.argmax(raw))
        steps[j] += 1
    step_alloc = {m: int(steps[i]) for i, m in enumerate(mods)}

    freeze = {m: bool(a[m] < cfg.freeze_theta) for m in mods}
    if all(freeze.values()):
        freeze[max(mods, key=lambda m: a[m])] = False

    p = np.array([a[m] for m in mods], float)
    p = p[p > 1e-12]
    ent = float(-(p * np.log(p)).sum()) if len(p) else 0.0
    ent_norm = ent / max(np.log(len(mods)), 1e-6)

    return {
        "alpha": a,
        "lr_mult": {m: float(lr[m]) for m in mods},
        "lr_shared": float(lr["shared"]),
        "step_alloc": step_alloc,
        "freeze_mask": freeze,
        "stack_prior": dict(a),
        "active_mods": [m for m in mods if not freeze[m]],
        "flops_rel": float(sum(1.0 for m in mods if not freeze[m]) / len(mods)),
        "alpha_entropy": ent_norm,
        "top_mod": max(mods, key=lambda m: a[m]),
    }


def opportunity_rank(
    cells: Sequence[Mapping],
    *,
    baseline: str = "equal",
) -> list[dict]:
    """Rank metric versions by holdout lift vs equal + side signals.

    Opportunity axes:
      Acc lift / MSE drop vs equal;
      FLOPs saved via freeze (if Acc not hurt);
      mid entropy preferred;
      forecast MAE when present.
    """
    by = {c["version"]: c for c in cells}
    if baseline not in by:
        return []
    b = by[baseline]
    rows = []
    for v, c in by.items():
        d_acc = float(c["mean_acc_lift"] - b["mean_acc_lift"])
        d_mse = float(c["mean_mse_drop"] - b["mean_mse_drop"])
        flops = float(c.get("mean_flops_rel", 1.0))
        ent = float(c.get("mean_alpha_entropy", 1.0))
        mae = float(c.get("mean_prop_mae", float("nan")))
        score = d_acc * 2.0 + d_mse
        if d_acc >= -0.005:
            score += 0.15 * (1.0 - flops)
        score -= 0.05 * abs(ent - 0.65)
        if np.isfinite(mae):
            score -= 0.10 * mae
        rows.append(
            {
                "version": v,
                "opp_score": float(score),
                "d_acc_vs_equal": d_acc,
                "d_mse_vs_equal": d_mse,
                "mean_flops_rel": flops,
                "mean_alpha_entropy": ent,
                "mean_prop_mae": mae,
                "mean_acc_post": float(c["mean_acc_post"]),
                "mean_mse_drop": float(c["mean_mse_drop"]),
                "mean_acc_lift": float(c["mean_acc_lift"]),
            }
        )
    rows.sort(key=lambda r: r["opp_score"], reverse=True)
    return rows
