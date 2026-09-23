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
  po_fuse       long⊗short concept emphasis         (see fuse_long_short)
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
    "po_fuse",
)

ACTUATORS = ("lr_mult", "step_alloc", "freeze_mask", "stack_prior")


@dataclass
class PORiskMetricConfig:
    tau: float = 0.30
    lam_cov: float = 1.0
    gamma_delta: float = 0.75
    ema_po: float = 0.55
    budget_floor: float = 0.15
    # long⊗short fusion (po_fuse): concept-modality emphasis at step t
    omega_long: float = 0.55
    omega_short0: float = 0.45
    spike_gain: float = 1.25
    freeze_from_long: bool = True
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


def _zscore_dict(d: Mapping[str, float], mods: Sequence[str]) -> dict[str, float]:
    vals = np.array([float(d[m]) for m in mods], float)
    mu = float(vals.mean())
    sd = float(vals.std())
    if sd < 1e-12:
        return {m: 0.0 for m in mods}
    return {m: (float(d[m]) - mu) / sd for m in mods}


def structured_epsilon_alpha(
    soft: Mapping[str, float],
    mods: Sequence[str],
    *,
    epsilon: float,
) -> dict[str, float]:
    """Structured ε-mix (ε-greedy *analogue*, not coin-flip).

    ``alpha_m = f + (1 - |M| f) * soft_m`` with ``f = ε/|M|``.
    Matches ``po_budget`` floor semantics: guaranteed warm mass on every
    tower so sensors can still fire, without randomizing α.
    """
    mods = list(mods)
    m = max(len(mods), 1)
    eps = float(np.clip(epsilon, 0.0, 1.0))
    floor = eps / m
    # keep floor feasible
    floor = min(floor, 0.9 / m)
    rem = 1.0 - floor * m
    s = np.array([max(float(soft.get(k, 0.0)), 0.0) for k in mods], float)
    if s.sum() <= 1e-12:
        s = np.ones(m) / m
    else:
        s = s / s.sum()
    out = {k: float(floor + rem * s[i]) for i, k in enumerate(mods)}
    z = sum(out.values())
    return {k: out[k] / z for k in mods}


def fuse_long_short(
    mods: Sequence[str],
    *,
    po: Mapping[str, float],
    po_prev: Mapping[str, float] | None = None,
    po_ema: Mapping[str, float] | None = None,
    mmd: Mapping[str, float] | None = None,
    proto: Mapping[str, float] | None = None,
    alpha_hist: Sequence[Mapping[str, float]] | None = None,
    cfg: PORiskMetricConfig | None = None,
) -> dict:
    """Fuse long-term + short-term sensors for concept-modality emphasis.

    Long track L_m (chronic residual-concept / institutional focus)
      L_m = EMA_ρ(PO_m) · (1 + proto_m)  [+ mild pull from last α]
    Short track S_m (this-step spike / anticipatory tilt)
      S_m = ΔPO_m = PO_m − PO_{m,t−1}   (fallback: PO − EMA)

    Adaptive mix (more short when the spike is large vs long MAD)::

        ω_S = clip( ω_S0 · (1 + spike_gain · ‖S‖_∞ / (MAD(L)+ε)) , 0, 1 )
        ω_L = 1 − ω_S   (or fixed omega_long if spike tiny)

    Score → Softmax::

        score_m = ω_L z(L_m) + ω_S z(S_m) − λ MMD_m

    Hierarchical actuators (returned in diag for next_step_actuators callers)::

        freeze_hint_m = 1{L_m < quantile_θ}   # slow — don't thrash on S
        step_tilt_m   ∝ Softmax(S)_active     # dump steps on short spike
    """
    cfg = cfg or PORiskMetricConfig()
    mods = list(mods)
    po_a = _as(mods, po)
    mmd_a = _as(mods, mmd)
    proto_a = _as(mods, proto)
    L_po = _ema_dict(po_ema, po_a, mods, cfg.ema_po)
    L = {m: L_po[m] * (1.0 + proto_a[m]) for m in mods}
    if alpha_hist:
        last_a = _as(mods, alpha_hist[-1])
        L = {m: 0.85 * L[m] + 0.15 * last_a[m] for m in mods}

    if po_prev is not None:
        S = {m: po_a[m] - float(po_prev[m]) for m in mods}
    else:
        S = {m: po_a[m] - L_po[m] for m in mods}

    L_vals = np.array([L[m] for m in mods], float)
    med = float(np.median(L_vals))
    mad = float(np.median(np.abs(L_vals - med))) + 1e-8
    spike = float(np.max(np.abs([S[m] for m in mods])))
    # cap so a single frame cannot force ω_S=1 from numerical MAD collapse
    spike_ratio = float(np.clip(spike / mad, 0.0, 4.0))
    w_s = float(cfg.omega_short0) * (1.0 + float(cfg.spike_gain) * spike_ratio)
    w_s = float(np.clip(w_s, 0.0, 1.0))
    # if spike is tiny, prefer configured long weight
    if spike_ratio < 0.25:
        w_l = float(cfg.omega_long)
        w_s = 1.0 - w_l
    else:
        w_l = 1.0 - w_s

    zL = _zscore_dict(L, mods)
    zS = _zscore_dict(S, mods)
    score = {
        m: w_l * zL[m] + w_s * zS[m] - float(cfg.lam_cov) * mmd_a[m] for m in mods
    }
    alpha = softmax_scores(score, mods, cfg.tau)

    # slow freeze from long track; keep at least one active
    thr = float(np.quantile(L_vals, 0.30)) if len(mods) > 1 else -1e9
    freeze_hint = {m: bool(cfg.freeze_from_long and L[m] < thr) for m in mods}
    if all(freeze_hint.values()):
        freeze_hint[max(mods, key=lambda m: L[m])] = False

    step_tilt = softmax_scores(S, mods, cfg.tau)
    # zero tilt mass on frozen (hierarchical)
    for m in mods:
        if freeze_hint[m]:
            step_tilt[m] = 0.0
    z = sum(step_tilt.values())
    if z <= 1e-12:
        step_tilt = {m: 1.0 / len(mods) for m in mods}
    else:
        step_tilt = {m: step_tilt[m] / z for m in mods}

    return {
        "alpha": alpha,
        "score": score,
        "long": L,
        "short": S,
        "omega_long": w_l,
        "omega_short": w_s,
        "spike_ratio": spike_ratio,
        "freeze_hint": freeze_hint,
        "step_tilt": step_tilt,
        "top_concept_mod": max(mods, key=lambda m: L[m]),
        "top_spike_mod": max(mods, key=lambda m: S[m]),
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
        # structured-ε mix: f=ε/|M|, ε := budget_floor * |M| (legacy knobs)
        eps = min(float(cfg.budget_floor) * len(mods), 0.9)
        alpha = structured_epsilon_alpha(soft, mods, epsilon=eps)
        score = dict(po_a)
        diag["floor"] = float(eps / max(len(mods), 1))
        diag["structured_epsilon"] = float(eps)
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

    if version == "po_fuse":
        fused = fuse_long_short(
            mods,
            po=po_a,
            po_prev=po_prev,
            po_ema=po_ema,
            mmd=mmd_a,
            proto=proto_a,
            alpha_hist=alpha_hist,
            cfg=cfg,
        )
        diag.update(
            {
                "omega_long": fused["omega_long"],
                "omega_short": fused["omega_short"],
                "spike_ratio": fused["spike_ratio"],
                "long": fused["long"],
                "short": fused["short"],
                "freeze_hint": fused["freeze_hint"],
                "step_tilt": fused["step_tilt"],
                "top_concept_mod": fused["top_concept_mod"],
                "top_spike_mod": fused["top_spike_mod"],
            }
        )
        return {"alpha": fused["alpha"], "score": fused["score"], "diag": diag}

    raise ValueError(
        f"unknown PO-risk metric version {version!r}; expected {METRIC_VERSIONS}"
    )


def next_step_actuators_fused(
    fuse_out: Mapping,
    mods: Sequence[str],
    *,
    cfg: NextStepActuatorConfig | None = None,
) -> dict:
    """Actuators with hierarchical long→freeze, short→step tilt.

    Use after ``fuse_long_short`` / ``metric_to_alpha(..., version='po_fuse')``.
    LR still follows fused α; freeze prefers long-track hint; step_alloc
    follows short-track tilt on the active set (acceleration focus).
    """
    cfg = cfg or NextStepActuatorConfig()
    mods = list(mods)
    alpha = fuse_out.get("alpha") or fuse_out
    # accept either fuse_long_short dict or metric_to_alpha return
    if "diag" in fuse_out and "freeze_hint" in dict(fuse_out.get("diag") or {}):
        diag = dict(fuse_out["diag"])
        alpha = fuse_out["alpha"]
        freeze_hint = diag.get("freeze_hint")
        step_tilt = diag.get("step_tilt")
    else:
        freeze_hint = fuse_out.get("freeze_hint")
        step_tilt = fuse_out.get("step_tilt")

    base = next_step_actuators(alpha, mods, cfg=cfg)
    if freeze_hint:
        freeze = {m: bool(freeze_hint.get(m, False)) for m in mods}
        if all(freeze.values()):
            freeze[max(mods, key=lambda m: float(alpha[m]))] = False
        base["freeze_mask"] = freeze
        base["active_mods"] = [m for m in mods if not freeze[m]]
        base["flops_rel"] = float(len(base["active_mods"]) / max(len(mods), 1))

    if step_tilt:
        active = base["active_mods"]
        raw = np.array(
            [float(step_tilt.get(m, 0.0)) if m in active else 0.0 for m in mods],
            float,
        )
        if raw.sum() <= 1e-12:
            raw = np.array(
                [1.0 if m in active else 0.0 for m in mods], float
            )
        raw = raw / raw.sum()
        steps = np.maximum(np.round(raw * cfg.total_steps), 0).astype(int)
        # ensure each active mod gets ≥ min_steps when possible
        for i, m in enumerate(mods):
            if m in active and steps[i] < cfg.min_steps:
                steps[i] = cfg.min_steps
        # repair sum
        while int(steps.sum()) > cfg.total_steps and (steps > cfg.min_steps).any():
            j = int(np.argmax(steps))
            if steps[j] > cfg.min_steps:
                steps[j] -= 1
            else:
                break
        while int(steps.sum()) < cfg.total_steps:
            j = int(np.argmax(raw))
            steps[j] += 1
        base["step_alloc"] = {m: int(steps[i]) for i, m in enumerate(mods)}
        base["step_tilt"] = {m: float(step_tilt.get(m, 0.0)) for m in mods}

    base["fuse"] = {
        "omega_long": fuse_out.get("omega_long")
        or (fuse_out.get("diag") or {}).get("omega_long"),
        "omega_short": fuse_out.get("omega_short")
        or (fuse_out.get("diag") or {}).get("omega_short"),
        "spike_ratio": fuse_out.get("spike_ratio")
        or (fuse_out.get("diag") or {}).get("spike_ratio"),
        "top_concept_mod": fuse_out.get("top_concept_mod")
        or (fuse_out.get("diag") or {}).get("top_concept_mod"),
        "top_spike_mod": fuse_out.get("top_spike_mod")
        or (fuse_out.get("diag") or {}).get("top_spike_mod"),
    }
    return base


def realize_step_alloc(
    step_alloc: Mapping[str, int],
    freeze_mask: Mapping[str, bool],
    mods: Sequence[str],
    *,
    redistribute: bool = False,
) -> dict[str, int]:
    """Apply freeze to step budgets (R2).

    Frozen mods get 0 steps (FLOPs save). If ``redistribute``, remaining
    budget is re-pumped onto active mods (same total compute, sharper dump).
    Default ``False``: freeze actually cuts update FLOPs.
    """
    mods = list(mods)
    out = {
        m: 0 if freeze_mask.get(m, False) else int(step_alloc.get(m, 0))
        for m in mods
    }
    if not redistribute:
        return out
    active = [m for m in mods if not freeze_mask.get(m, False)]
    if not active:
        return out
    budget = int(sum(int(step_alloc.get(m, 0)) for m in mods))
    raw = np.array([max(float(step_alloc.get(m, 0)), 0.0) for m in active], float)
    if raw.sum() <= 1e-12:
        raw = np.ones(len(active))
    raw = raw / raw.sum()
    steps = np.maximum(np.round(raw * budget), 0).astype(int)
    while int(steps.sum()) > budget and steps.sum() > 0:
        j = int(np.argmax(steps))
        if steps[j] > 0:
            steps[j] -= 1
        else:
            break
    while int(steps.sum()) < budget:
        j = int(np.argmax(raw))
        steps[j] += 1
    return {m: 0 for m in mods} | {m: int(steps[i]) for i, m in enumerate(active)}


def expand_step_schedule(
    step_alloc: Mapping[str, int],
    mods: Sequence[str],
    *,
    mode: str = "block",
) -> list[str]:
    """Expand per-mod budgets into an ordered optimizer-step schedule (R2).

    ``block`` (default): dump highest-budget modality first — matches
    short-track spike dump before Acc collapses.
    ``round_robin``: interleave remaining slots (smoother, less dump-like).
    """
    mods = list(mods)
    rem = {m: max(0, int(step_alloc.get(m, 0))) for m in mods}
    schedule: list[str] = []
    if mode == "round_robin":
        while any(rem[m] > 0 for m in mods):
            for m in sorted(mods, key=lambda x: -rem[x]):
                if rem[m] > 0:
                    schedule.append(m)
                    rem[m] -= 1
        return schedule
    # block: highest remaining budget first, all its steps contiguous
    order = sorted(mods, key=lambda m: (-rem[m], m))
    for m in order:
        schedule.extend([m] * rem[m])
    return schedule


def step_flops_rel(
    realized: Mapping[str, int],
    *,
    total_steps: int,
) -> float:
    """Update-FLOPs proxy under real step_alloc: steps_used / budget."""
    tot = max(int(total_steps), 1)
    used = int(sum(max(0, int(v)) for v in realized.values()))
    return float(min(1.0, used / tot))


ROW_WEIGHT_MODES = ("sqrt", "cbrt", "prop", "uniform")


def row_po_residual(
    y: np.ndarray,
    mu: np.ndarray,
) -> np.ndarray:
    """Observation-level PO proxy: |Y − μ| (residual under control fit)."""
    y = np.asarray(y, dtype=float).reshape(-1)
    mu = np.asarray(mu, dtype=float).reshape(-1)
    if mu.shape[0] != y.shape[0]:
        raise ValueError("y and mu length mismatch")
    return np.abs(y - mu)


def po_iptw_weights(
    po_i: np.ndarray,
    *,
    mode: str = "sqrt",
    rejected: bool = True,
    clip: tuple[float, float] = (0.25, 4.0),
) -> np.ndarray:
    """Reject-gated row weights (Logic B / R3).

    Calm (``rejected=False``) or ``mode='uniform'`` → ones.
    Default soft ``sqrt``; alts ``cbrt`` / ``prop``. Mean-normalized then clipped.
    High PO_i = residual under fit — not intrinsic label difficulty.
    """
    po = np.asarray(po_i, dtype=float).reshape(-1)
    if po.size == 0:
        return po
    if (not rejected) or mode == "uniform":
        return np.ones_like(po, dtype=float)
    if mode == "sqrt":
        raw = np.sqrt(np.maximum(po, 0.0) + 1e-12)
    elif mode == "cbrt":
        raw = np.cbrt(np.maximum(po, 0.0) + 1e-12)
    elif mode == "prop":
        raw = np.maximum(po, 0.0) + 1e-12
    else:
        raise ValueError(f"unknown row weight mode: {mode}")
    w = raw / max(float(raw.mean()), 1e-12)
    lo, hi = float(clip[0]), float(clip[1])
    return np.clip(w, lo, hi).astype(float)


def stream_reject_proxy(
    *,
    po_mods: Mapping[str, float],
    mmd_mods: Mapping[str, float] | None = None,
    po_prev: Mapping[str, float] | None = None,
    mods: Sequence[str] | None = None,
    po_mean_thresh: float = 0.05,
    delta_thresh: float = 0.04,
    mmd_mean_thresh: float = 0.08,
) -> dict:
    """Lightweight reject event when OnlineRFPerm is not in the loop.

    Fires if mean PO, mean |ΔPO|, or mean MMD clears a threshold.
    Replace with real RFPerm/CFPerm reject flag in production; same
    ``po_iptw_weights(..., rejected=)`` hook either way.
    """
    mods = list(mods or po_mods.keys())
    po_vals = np.array([float(po_mods.get(m, 0.0)) for m in mods], float)
    mean_po = float(po_vals.mean()) if len(po_vals) else 0.0
    if po_prev is not None:
        d = np.array(
            [abs(float(po_mods.get(m, 0.0)) - float(po_prev.get(m, 0.0))) for m in mods],
            float,
        )
        mean_d = float(d.mean()) if len(d) else 0.0
    else:
        mean_d = 0.0
    if mmd_mods:
        mmd_vals = np.array([float(mmd_mods.get(m, 0.0)) for m in mods], float)
        mean_mmd = float(mmd_vals.mean()) if len(mmd_vals) else 0.0
    else:
        mean_mmd = 0.0
    reasons = []
    if mean_po >= po_mean_thresh:
        reasons.append("mean_po")
    if mean_d >= delta_thresh:
        reasons.append("mean_delta_po")
    if mean_mmd >= mmd_mean_thresh:
        reasons.append("mean_mmd")
    return {
        "rejected": bool(reasons),
        "reasons": reasons,
        "mean_po": mean_po,
        "mean_delta_po": mean_d,
        "mean_mmd": mean_mmd,
        "proxy": True,
    }


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


def freeze_jaccard(
    freeze_a: Mapping[str, bool],
    freeze_b: Mapping[str, bool],
    mods: Sequence[str],
) -> float:
    """Jaccard similarity of frozen sets (stability of long-track freeze)."""
    a = {m for m in mods if freeze_a.get(m)}
    b = {m for m in mods if freeze_b.get(m)}
    if not a and not b:
        return 1.0
    inter = len(a & b)
    union = len(a | b)
    return float(inter / union) if union else 1.0


def continuous_gain_metrics(
    rows: Sequence[Mapping],
    *,
    mods: Sequence[str] | None = None,
    acc_star: float = 0.55,
) -> dict:
    """R1 continuous-time gains from a compare trajectory.

    Returns
    -------
    mean_flops_rel, cum_flops (sum of flops_rel),
    t_to_acc_star (first t with acc>=acc_star, else None),
    cum_flops_to_acc_star,
    mean_freeze_jaccard (adjacent windows),
    n_windows
    """
    rows = list(rows)
    if not rows:
        return {
            "mean_flops_rel": 1.0,
            "cum_flops": 0.0,
            "t_to_acc_star": None,
            "cum_flops_to_acc_star": None,
            "mean_freeze_jaccard": 1.0,
            "n_windows": 0,
        }
    if mods is None:
        mods = list((rows[0].get("freeze") or rows[0].get("alpha") or {}).keys())
    flops = [float(r.get("flops_rel", 1.0)) for r in rows]
    accs = [float(r.get("acc", 0.0)) for r in rows]
    cum = float(np.cumsum(flops)[-1]) if flops else 0.0
    t_hit = None
    cum_hit = None
    running = 0.0
    for t, (a, f) in enumerate(zip(accs, flops)):
        running += f
        if t_hit is None and a >= float(acc_star):
            t_hit = t
            cum_hit = running
    jacs = []
    for i in range(1, len(rows)):
        fa = rows[i - 1].get("freeze") or {}
        fb = rows[i].get("freeze") or {}
        jacs.append(freeze_jaccard(fa, fb, mods))
    return {
        "mean_flops_rel": float(np.mean(flops)),
        "cum_flops": cum,
        "t_to_acc_star": t_hit,
        "cum_flops_to_acc_star": cum_hit,
        "mean_freeze_jaccard": float(np.mean(jacs)) if jacs else 1.0,
        "n_windows": len(rows),
        "acc_star": float(acc_star),
    }