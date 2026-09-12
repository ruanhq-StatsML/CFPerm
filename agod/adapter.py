"""Layered next-stage adapter + dynamic hard-gate policies.

Design target (not a new backbone):
  L0 sensors → L1 EMA α state → L2 soft LR adapter → L3 hard adapt-gate.

Hard-gate only zeros modality *projection BWD*; FWD inference always runs.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Mapping, Sequence

import numpy as np

from .lr_controller import alpha_to_lr

ADAPTER_LAYERS = (
    "L0_sensor",
    "L1_state_alpha",
    "L2_soft_lr",
    "L3_hard_gate",
)

GATE_POLICIES = (
    "none",       # all modalities active (soft LR only)
    "fixed",      # α_m >= θ
    "quantile",   # keep top mass / drop bottom quantile
    "ema_theta",  # θ tracks EMA of mean(α) · ratio
    "hysteresis", # enter/exit thresholds to reduce chatter
    "random",     # non-attribution sparsity control
)


@dataclass
class AdapterConfig:
    """Knobs for soft LR + hard gate."""

    beta: float = 0.10
    theta: float = 0.28
    theta_on: float = 0.32   # hysteresis: turn ON when α >= theta_on
    theta_off: float = 0.24  # hysteresis: turn OFF when α < theta_off
    q_drop: float = 0.33     # quantile: drop bottom fraction of α
    ema_theta_ratio: float = 0.85  # θ_t = ratio * EMA(mean α)
    p_keep_rand: float = 0.67
    ensure_one: bool = True


@dataclass
class AdapterDecision:
    alpha: dict[str, float]
    lr_mult: dict[str, float]
    active: dict[str, bool]
    theta_used: float
    gate_policy: str
    layer_trace: dict[str, object] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)


class DynamicGate:
    """Stateful L3 hard-gate (needed for ema_theta / hysteresis)."""

    def __init__(self, mods: Sequence[str], cfg: AdapterConfig | None = None):
        self.mods = list(mods)
        self.cfg = cfg or AdapterConfig()
        self.theta_ema = float(self.cfg.theta)
        self.prev_active = {m: True for m in self.mods}

    def reset(self) -> None:
        self.theta_ema = float(self.cfg.theta)
        self.prev_active = {m: True for m in self.mods}

    def _ensure_one(self, active: dict[str, bool], alpha: Mapping[str, float]) -> dict[str, bool]:
        if self.cfg.ensure_one and not any(active.values()):
            mstar = max(self.mods, key=lambda m: float(alpha[m]))
            active = {m: m == mstar for m in self.mods}
        return active

    def decide(
        self,
        alpha: Mapping[str, float],
        *,
        policy: str,
        seed: int = 0,
        t: int = 0,
    ) -> tuple[dict[str, bool], float]:
        if policy not in GATE_POLICIES:
            raise ValueError(f"unknown gate policy {policy!r}; expected {GATE_POLICIES}")

        if policy == "none":
            active = {m: True for m in self.mods}
            return active, 0.0

        if policy == "fixed":
            th = float(self.cfg.theta)
            active = {m: float(alpha[m]) >= th for m in self.mods}
            return self._ensure_one(active, alpha), th

        if policy == "quantile":
            vals = np.array([float(alpha[m]) for m in self.mods], float)
            # keep modalities above the q_drop quantile of α
            th = float(np.quantile(vals, self.cfg.q_drop))
            active = {m: float(alpha[m]) >= th for m in self.mods}
            # if all equal, keep all
            if len(np.unique(np.round(vals, 6))) == 1:
                active = {m: True for m in self.mods}
            return self._ensure_one(active, alpha), th

        if policy == "ema_theta":
            mean_a = float(np.mean([float(alpha[m]) for m in self.mods]))
            self.theta_ema = 0.7 * self.theta_ema + 0.3 * mean_a
            th = float(self.cfg.ema_theta_ratio * self.theta_ema)
            active = {m: float(alpha[m]) >= th for m in self.mods}
            return self._ensure_one(active, alpha), th

        if policy == "hysteresis":
            active = {}
            for m in self.mods:
                a = float(alpha[m])
                if self.prev_active[m]:
                    active[m] = a >= self.cfg.theta_off
                else:
                    active[m] = a >= self.cfg.theta_on
            active = self._ensure_one(active, alpha)
            self.prev_active = dict(active)
            # report mid threshold for telemetry
            th = 0.5 * (self.cfg.theta_on + self.cfg.theta_off)
            return active, float(th)

        # random
        rng = np.random.default_rng(seed + 101 * t + 17)
        active = {m: bool(rng.random() < self.cfg.p_keep_rand) for m in self.mods}
        return self._ensure_one(active, alpha), float(self.cfg.p_keep_rand)


class LayeredAdapter:
    """L2 soft LR + L3 hard gate over EMA α (L1 provided by caller)."""

    def __init__(
        self,
        mods: Sequence[str],
        cfg: AdapterConfig | None = None,
        *,
        gate_policy: str = "fixed",
    ):
        self.mods = list(mods)
        self.cfg = cfg or AdapterConfig()
        self.gate_policy = gate_policy
        self.gate = DynamicGate(self.mods, self.cfg)

    def reset(self) -> None:
        self.gate.reset()

    def step(
        self,
        alpha: Mapping[str, float],
        *,
        seed: int = 0,
        t: int = 0,
        gain: Mapping[str, float] | None = None,
    ) -> AdapterDecision:
        lr_mult = alpha_to_lr(alpha, self.mods, beta=self.cfg.beta, gain=gain)
        active, th = self.gate.decide(
            alpha, policy=self.gate_policy, seed=seed, t=t
        )
        # soft LR still applied on active mods; inactive get LR but grads zeroed upstream
        return AdapterDecision(
            alpha={m: float(alpha[m]) for m in self.mods},
            lr_mult=lr_mult,
            active=active,
            theta_used=th,
            gate_policy=self.gate_policy,
            layer_trace={
                "L2_soft_lr": {m: lr_mult[m] for m in self.mods},
                "L3_hard_gate": {
                    "policy": self.gate_policy,
                    "theta_used": th,
                    "active": active,
                },
            },
        )


def flops_rel_proj(
    active: Mapping[str, bool],
    mods: Sequence[str],
    *,
    c_pf: float = 1.0,
    c_pb: float = 2.0,
    c_sf: float = 0.5,
    c_sb: float = 1.0,
) -> float:
    """Relative projection adapt FLOPs (FWD always on)."""
    fwd = len(mods) * c_pf + c_sf
    bwd = sum(c_pb for m in mods if active[m]) + c_sb
    full = len(mods) * c_pf + c_sf + len(mods) * c_pb + c_sb
    return float((fwd + bwd) / full)


def characterize_gate_traj(traj: Sequence[Mapping], mods: Sequence[str]) -> dict:
    """Summarize dynamic-gate behavior over windows."""
    if not traj:
        return {}
    fracs = {m: float(np.mean([bool(r["active"][m]) for r in traj])) for m in mods}
    switches = 0
    for i in range(1, len(traj)):
        for m in mods:
            if bool(traj[i]["active"][m]) != bool(traj[i - 1]["active"][m]):
                switches += 1
    thetas = [float(r.get("theta_used", np.nan)) for r in traj]
    return {
        "frac_active": fracs,
        "mean_active_mods": float(np.mean([sum(bool(r["active"][m]) for m in mods) for r in traj])),
        "gate_switches": int(switches),
        "switches_per_window": float(switches / max(len(traj) - 1, 1)),
        "mean_theta_used": float(np.nanmean(thetas)),
        "mean_flops_rel": float(np.mean([r["flops_rel"] for r in traj])),
        "mean_acc_lift": float(np.nanmean([r["acc_lift"] for r in traj])),
        "mean_acc_post": float(np.nanmean([r["acc_post"] for r in traj])),
        "mean_cost_utility": float(
            np.nanmean(
                [
                    (r["acc_lift"] / r["flops_rel"]) if r["flops_rel"] > 0 else np.nan
                    for r in traj
                ]
            )
        ),
        "n_windows": len(traj),
    }


def soft_lr_dispersion(lr_mult: Mapping[str, float], mods: Sequence[str]) -> dict:
    """Characterize continuous next-stage step sizes (L2 soft LR)."""
    vals = np.array([float(lr_mult[m]) for m in mods], float)
    vmin = float(vals.min()) if len(vals) else float("nan")
    vmax = float(vals.max()) if len(vals) else float("nan")
    mean = float(vals.mean()) if len(vals) else float("nan")
    std = float(vals.std()) if len(vals) else float("nan")
    # relative range vs equal-LR baseline (=1)
    return {
        "lr_mean": mean,
        "lr_std": std,
        "lr_min": vmin,
        "lr_max": vmax,
        "lr_range": float(vmax - vmin) if len(vals) else float("nan"),
        "lr_ratio": float(vmax / max(vmin, 1e-8)) if len(vals) else float("nan"),
        "lr_cv": float(std / max(abs(mean), 1e-8)) if len(vals) else float("nan"),
    }


def alpha_entropy(alpha: Mapping[str, float], mods: Sequence[str]) -> float:
    """Shannon entropy of routing mass (nats). Equal → log|M|."""
    p = np.array([max(float(alpha[m]), 1e-12) for m in mods], float)
    p = p / p.sum()
    return float(-(p * np.log(p)).sum())


def characterize_soft_lr_traj(traj: Sequence[Mapping], mods: Sequence[str]) -> dict:
    """Summarize continuous soft-LR adapter (no gate zoo)."""
    if not traj:
        return {}
    lifts = np.array([float(r["acc_lift"]) for r in traj], float)
    flops = np.array([float(r.get("flops_rel", 1.0)) for r in traj], float)
    ents = [alpha_entropy(r["alpha"], mods) for r in traj]
    disp = [soft_lr_dispersion(r["lr_mult"], mods) for r in traj]
    n_pos = int(np.sum(lifts > 0))
    return {
        "n_windows": len(traj),
        "mean_acc_lift": float(np.nanmean(lifts)),
        "median_acc_lift": float(np.nanmedian(lifts)),
        "frac_windows_positive_lift": float(n_pos / max(len(traj), 1)),
        "mean_acc_post": float(np.nanmean([float(r["acc_post"]) for r in traj])),
        "mean_flops_rel": float(np.nanmean(flops)),
        "mean_cost_utility": float(
            np.nanmean(
                [
                    (float(r["acc_lift"]) / float(r["flops_rel"]))
                    if float(r.get("flops_rel", 0)) > 0
                    else np.nan
                    for r in traj
                ]
            )
        ),
        "mean_alpha_entropy": float(np.nanmean(ents)),
        "equal_alpha_entropy": float(np.log(len(mods))),
        "mean_lr_std": float(np.nanmean([d["lr_std"] for d in disp])),
        "mean_lr_range": float(np.nanmean([d["lr_range"] for d in disp])),
        "mean_lr_ratio": float(np.nanmean([d["lr_ratio"] for d in disp])),
        "mean_lr_cv": float(np.nanmean([d["lr_cv"] for d in disp])),
    }


def n_star_for_lift(
    curve: Sequence[Mapping],
    *,
    n_key: str = "n_adapt",
    lift_key: str = "mean_acc_lift",
    eps: float = 0.0,
) -> dict:
    """Smallest sample size on an ascending N-curve with mean Acc lift > eps.

    ``curve`` rows must already be sorted by ``n_key`` ascending.
    Returns n_star / first_positive + full curve diagnostics.
    """
    ordered = sorted(curve, key=lambda r: float(r[n_key]))
    n_star = None
    first = None
    for row in ordered:
        lift = float(row[lift_key])
        if first is None and lift > eps:
            first = float(row[n_key])
        if n_star is None and lift > eps:
            # require non-decreasing positivity: keep first crossing
            n_star = float(row[n_key])
            break
    return {
        "eps": float(eps),
        "n_star": n_star,
        "first_positive_n": first,
        "max_lift": float(max((float(r[lift_key]) for r in ordered), default=float("nan"))),
        "best_n": (
            float(max(ordered, key=lambda r: float(r[lift_key]))[n_key])
            if ordered
            else None
        ),
    }
