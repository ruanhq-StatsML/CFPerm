"""Smooth drift-vs-noise router for AGOD modality weights / LR.

Problem: high PO-risk can be true concept drift *or* noisy / low-info
modalities. Blindly up-weighting then memorizes noise.

Gate (per modality, per window):
  boost iff  PO high  AND  unimodal task metric not bad  AND  Fisher/VIMP not low
  else       damp PO contribution (down-weight / adapter-only hint)

Smooth control law (estimator has variance — do not hard-step):

  g_m     = normalize(w1·PO_gated_m + w2·MMD_m + w3·VIMP_m)
  α_raw   = Softmax(g / τ)
  α       = EMA(α_raw)
  w_m     = w0 · (β + (1-β) · α_m · |M|)     # bounded
  loss    = Σ_m w_m · loss_m   (or LR_m ∝ w_m on a shared CE)

FWD stays on for all modalities; the gate only reshapes adapt spend.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping, MutableMapping, Sequence

import numpy as np


def _as_dict(mods: Sequence[str], values: Mapping[str, float] | None, default: float = 0.0) -> dict[str, float]:
    values = values or {}
    return {m: float(values.get(m, default)) for m in mods}


def normalize_nonneg(d: Mapping[str, float], mods: Sequence[str]) -> dict[str, float]:
    v = np.array([max(float(d[m]), 0.0) for m in mods], float)
    s = float(v.sum())
    if s <= 1e-12:
        return {m: 1.0 / len(mods) for m in mods}
    v = v / s
    return {m: float(v[i]) for i, m in enumerate(mods)}


def softmax_tau(d: Mapping[str, float], mods: Sequence[str], tau: float) -> dict[str, float]:
    z = np.array([float(d[m]) for m in mods], float) / max(float(tau), 1e-6)
    z = z - z.max()
    e = np.exp(z)
    e = e / e.sum()
    return {m: float(e[i]) for i, m in enumerate(mods)}


@dataclass
class DriftNoiseGateConfig:
    """Thresholds for 'true drift' vs 'noise drifting'."""

    po_high: float = 0.02          # absolute PO-risk floor to even consider boost
    po_quantile: float = 0.55      # also require PO ≥ this cross-mod quantile
    uni_acc_floor: float = 0.52    # unimodal hold Acc must clear this
    uni_drop_tol: float = 0.05     # or Acc drop vs ref ≤ this
    vimp_floor: float = 0.02       # Fisher/importance proxy floor
    vimp_quantile: float = 0.35    # and not in the bottom quantile
    damp: float = 0.15             # multiply PO by this when gated as noise
    adapter_only_damp: float = 0.05  # even stronger damp → near adapter-only


@dataclass
class SmoothRouterConfig:
    w_po: float = 1.0
    w_mmd: float = 0.75
    w_vimp: float = 0.50
    tau: float = 0.30
    ema: float = 0.40
    beta: float = 0.10          # LR / loss-weight floor mix
    w0: float = 1.0
    gate: DriftNoiseGateConfig = field(default_factory=DriftNoiseGateConfig)


def fisher_proxy_from_vimp(vimp: Mapping[str, float], mods: Sequence[str]) -> dict[str, float]:
    """Use mean RF VIMP as a cheap Fisher/importance proxy (scale-free later)."""
    return _as_dict(mods, vimp, 0.0)


def drift_vs_noise_gate(
    po: Mapping[str, float],
    uni_acc: Mapping[str, float],
    vimp: Mapping[str, float],
    mods: Sequence[str],
    *,
    uni_acc_ref: Mapping[str, float] | None = None,
    cfg: DriftNoiseGateConfig | None = None,
) -> dict:
    """Decide which modalities look like *signal drift* vs *noise drift*.

    Returns gated PO, per-mod flags, and adapter_only hints.
    """
    cfg = cfg or DriftNoiseGateConfig()
    mods = list(mods)
    po_a = _as_dict(mods, po)
    uni = _as_dict(mods, uni_acc, 0.5)
    vim = fisher_proxy_from_vimp(vimp, mods)
    ref = _as_dict(mods, uni_acc_ref, 0.5) if uni_acc_ref is not None else None

    po_vals = np.array([po_a[m] for m in mods], float)
    v_vals = np.array([vim[m] for m in mods], float)
    po_q = float(np.quantile(po_vals, cfg.po_quantile)) if len(po_vals) else 0.0
    v_q = float(np.quantile(v_vals, cfg.vimp_quantile)) if len(v_vals) else 0.0

    gated = {}
    is_signal = {}
    adapter_only = {}
    reasons = {}
    for m in mods:
        high_po = (po_a[m] >= cfg.po_high) and (po_a[m] >= po_q)
        uni_ok = uni[m] >= cfg.uni_acc_floor
        if ref is not None:
            uni_ok = uni_ok or ((ref[m] - uni[m]) <= cfg.uni_drop_tol)
        imp_ok = (vim[m] >= cfg.vimp_floor) and (vim[m] >= v_q)
        sig = bool(high_po and uni_ok and imp_ok)
        is_signal[m] = sig
        if sig:
            gated[m] = po_a[m]
            adapter_only[m] = False
            reasons[m] = "signal_drift"
        elif high_po and not uni_ok:
            # PO high but unimodal task already bad → noise / weak predictor
            gated[m] = cfg.adapter_only_damp * po_a[m]
            adapter_only[m] = True
            reasons[m] = "noise_po_weak_uni"
        elif high_po and not imp_ok:
            gated[m] = cfg.damp * po_a[m]
            adapter_only[m] = True
            reasons[m] = "noise_po_low_fisher"
        else:
            gated[m] = po_a[m]  # low PO: leave as-is (no boost path)
            adapter_only[m] = False
            reasons[m] = "low_po_pass"

    return {
        "po_gated": gated,
        "is_signal": is_signal,
        "adapter_only": adapter_only,
        "reasons": reasons,
        "po_quantile": po_q,
        "vimp_quantile": v_q,
        "frac_signal": float(np.mean([1.0 if is_signal[m] else 0.0 for m in mods])),
    }


def compose_g(
    po_gated: Mapping[str, float],
    mmd: Mapping[str, float],
    vimp: Mapping[str, float],
    mods: Sequence[str],
    *,
    w_po: float = 1.0,
    w_mmd: float = 0.75,
    w_vimp: float = 0.50,
) -> dict[str, float]:
    """g_m = normalize(w1·PO_gated + w2·MMD + w3·VIMP)."""
    mods = list(mods)
    raw = {
        m: float(w_po) * max(float(po_gated.get(m, 0.0)), 0.0)
        + float(w_mmd) * max(float(mmd.get(m, 0.0)), 0.0)
        + float(w_vimp) * max(float(vimp.get(m, 0.0)), 0.0)
        for m in mods
    }
    return normalize_nonneg(raw, mods)


def weights_from_alpha(
    alpha: Mapping[str, float],
    mods: Sequence[str],
    *,
    beta: float = 0.10,
    w0: float = 1.0,
) -> dict[str, float]:
    """w_m = w0 · (β + (1-β) · α_m · |M|) with shared mean."""
    mods = list(mods)
    inv = float(len(mods))
    out = {
        m: float(w0) * float(beta + (1.0 - beta) * float(alpha[m]) * inv) for m in mods
    }
    out["shared"] = float(np.mean([out[m] for m in mods]))
    return out


class SmoothDriftNoiseRouter:
    """Stateful EMA router with drift-vs-noise PO gating."""

    def __init__(self, mods: Sequence[str], cfg: SmoothRouterConfig | None = None):
        self.mods = list(mods)
        self.cfg = cfg or SmoothRouterConfig()
        n = len(self.mods)
        self.alpha = {m: 1.0 / n for m in self.mods}
        self.last_gate: dict = {}
        self.last_g: dict[str, float] = dict(self.alpha)

    def update(
        self,
        *,
        po: Mapping[str, float],
        mmd: Mapping[str, float],
        vimp: Mapping[str, float],
        uni_acc: Mapping[str, float],
        uni_acc_ref: Mapping[str, float] | None = None,
    ) -> dict:
        gate = drift_vs_noise_gate(
            po,
            uni_acc,
            vimp,
            self.mods,
            uni_acc_ref=uni_acc_ref,
            cfg=self.cfg.gate,
        )
        g = compose_g(
            gate["po_gated"],
            mmd,
            vimp,
            self.mods,
            w_po=self.cfg.w_po,
            w_mmd=self.cfg.w_mmd,
            w_vimp=self.cfg.w_vimp,
        )
        raw = softmax_tau(g, self.mods, self.cfg.tau)
        ema = float(self.cfg.ema)
        for m in self.mods:
            self.alpha[m] = ema * self.alpha[m] + (1.0 - ema) * raw[m]
        s = sum(self.alpha.values())
        self.alpha = {m: self.alpha[m] / s for m in self.mods}
        w = weights_from_alpha(
            self.alpha, self.mods, beta=self.cfg.beta, w0=self.cfg.w0
        )
        # near-adapter-only: shrink LR further when gated as noise
        for m in self.mods:
            if gate["adapter_only"].get(m):
                w[m] = float(self.cfg.beta * self.cfg.w0)
        w["shared"] = float(np.mean([w[m] for m in self.mods]))

        self.last_gate = gate
        self.last_g = g
        return {
            "alpha": dict(self.alpha),
            "alpha_raw": raw,
            "g": g,
            "weights": w,
            "lr_mult": w,  # same bounded map; caller may use as LR multipliers
            "gate": gate,
        }


def weighted_modality_loss(
    losses: Mapping[str, float],
    weights: Mapping[str, float],
    mods: Sequence[str],
) -> float:
    """loss = Σ_m w_m · loss_m  (for logging / numpy side)."""
    return float(sum(float(weights[m]) * float(losses[m]) for m in mods))
