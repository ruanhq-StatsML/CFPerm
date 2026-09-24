"""Soft IPTW burn policy + FLOPs ledger (√ / ∛).

Question
--------
Should we *burn* soft PO-weights (w∝PO^α) on reject batches?

Key separation (算力账)
-----------------------
- **Adaptation FLOPs** (big): OnlineRFPerm reject → RF / μ0 re-fit.
  Accounted by ``expected_adapt_flops`` / duty·n_control.
- **Weighting FLOPs** (tiny): w_i = PO_i^α is O(n) — **does not move the bill**.

So √ vs ∛ is **not** a compute choice.  It is a bias–variance / MSE-risk
choice *after* you already paid for the reject.

Burn logic (该不该烧)
--------------------
1. Default **KEEP_UNIFORM** — evidence: uniform best on 5/6 packs.
2. **BURN_α** only if gated_α beats uniform on sig-MSE (rel < 1).
3. If forced to gate (ops already rejecting) but gated hurts:
   **SOFTEN_ONLY** → prefer ∛ over √ (less damage; same FLOPs).
4. Never burn soft weights to \"save FLOPs\" — they don't.

Decision codes
--------------
- ``KEEP_UNIFORM``
- ``BURN_SQRT`` / ``BURN_CBRT``
- ``SOFTEN_ONLY``  (gate already on; use ∛)
- ``INDETERMINATE``
"""
from __future__ import annotations

from typing import Any, Dict, Mapping, Optional, Sequence

import numpy as np


def weighting_flops_proxy(n: int) -> float:
    """O(n) reweight — negligible vs RF fit; kept for ledger completeness."""
    return float(max(int(n), 0))


def soft_flops_ledger(
    *,
    adapt_flops: float,
    n_rows: int = 0,
) -> Dict[str, float]:
    """Split adaptation vs weighting so soft-α is not mistaken for a FLOPs knob."""
    w = weighting_flops_proxy(n_rows)
    adapt = float(adapt_flops) if np.isfinite(adapt_flops) else float("nan")
    total = (
        float(adapt + w)
        if np.isfinite(adapt)
        else float("nan")
    )
    share = float(w / total) if np.isfinite(total) and total > 0 else float("nan")
    return {
        "adapt_flops": adapt,
        "weighting_flops": w,
        "total_flops": total,
        "weighting_share": share,  # typically ≪ 1%
    }


def burn_soft_weights(
    *,
    rel_gated_sqrt: float,
    rel_gated_cbrt: float,
    soft_win_cbrt_le_sqrt: Optional[bool] = None,
    gate_already_on: bool = True,
    improve_eps: float = 1e-6,
) -> Dict[str, Any]:
    """Decide whether to burn soft IPTW given rel-MSE vs uniform.

    ``rel_*`` = mse_gated / mse_uniform ( <1 ⇒ improves ).
    """
    rs = float(rel_gated_sqrt) if rel_gated_sqrt is not None else float("nan")
    rc = float(rel_gated_cbrt) if rel_gated_cbrt is not None else float("nan")
    sqrt_helps = np.isfinite(rs) and rs < 1.0 - improve_eps
    cbrt_helps = np.isfinite(rc) and rc < 1.0 - improve_eps

    if sqrt_helps or cbrt_helps:
        # Pick the better improver
        if sqrt_helps and cbrt_helps:
            if rs <= rc:
                code, alpha, why = (
                    "BURN_SQRT",
                    0.5,
                    "both beat uniform; √ has lower rel MSE",
                )
            else:
                code, alpha, why = (
                    "BURN_CBRT",
                    1.0 / 3.0,
                    "both beat uniform; ∛ has lower rel MSE",
                )
        elif sqrt_helps:
            code, alpha, why = "BURN_SQRT", 0.5, "only √ beats uniform on sig-MSE"
        else:
            code, alpha, why = "BURN_CBRT", 1.0 / 3.0, "only ∛ beats uniform on sig-MSE"
        return {
            "decision": code,
            "burn": True,
            "alpha": alpha,
            "reason": why,
            "flops_note": "weighting FLOPs≈0; bill already paid by gate/refit",
        }

    # Neither improves vs uniform
    if gate_already_on:
        prefer_cbrt = soft_win_cbrt_le_sqrt is not False  # default soften
        return {
            "decision": "SOFTEN_ONLY",
            "burn": False,  # do not claim IPTW win
            "alpha": (1.0 / 3.0) if prefer_cbrt else 0.5,
            "reason": (
                "gated IPTW does not beat uniform; √ hurts less than ∛ here"
                if not prefer_cbrt
                else (
                    "gated IPTW does not beat uniform; if reject path is mandatory, "
                    "use ∛ to limit damage (same FLOPs as √)"
                )
            ),
            "flops_note": "α choice ≠ FLOPs; keep adaptation budget on refit/duty",
        }

    return {
        "decision": "KEEP_UNIFORM",
        "burn": False,
        "alpha": None,
        "reason": "no evidence gated soft IPTW beats uniform — do not burn",
        "flops_note": "saving weighting FLOPs is irrelevant; skip the risk",
    }


def attach_burn_to_power_card(card: Mapping[str, Any]) -> Dict[str, Any]:
    """Enrich one ``power_card_from_dataset`` row with burn decision + ledger."""
    out = dict(card)
    by = {m["mode"]: m for m in (card.get("modes") or [])}
    gs = by.get("gated_sqrt") or {}
    gc = by.get("gated_cbrt") or {}
    decision = burn_soft_weights(
        rel_gated_sqrt=float(gs.get("rel_vs_uniform", float("nan"))),
        rel_gated_cbrt=float(gc.get("rel_vs_uniform", float("nan"))),
        soft_win_cbrt_le_sqrt=card.get("soft_win_cbrt_le_sqrt"),
        gate_already_on=True,
    )
    # adapt flops from gated mse_eff denominator path: use expected refit scale if present
    adapt = float("nan")
    for m in (gs, gc):
        # recover adapt flops from mse_eff definition if possible — optional
        break
    duty = card.get("gate_duty")
    # ledger without needing n_rows: weighting share → 0 message
    ledger = soft_flops_ledger(adapt_flops=1.0, n_rows=0)  # placeholder scale
    ledger = {
        "adapt_flops": "duty·n·fit (see po_eff)",
        "weighting_flops": "~O(n) negligible",
        "weighting_share": "~0",
        "duty": duty,
    }
    out["burn"] = decision
    out["flops_ledger"] = ledger
    return out


def summarize_burn_decisions(
    cards: Sequence[Mapping[str, Any]],
) -> Dict[str, Any]:
    from collections import Counter

    enriched = [attach_burn_to_power_card(c) for c in cards if c.get("ok")]
    counts = Counter((c.get("burn") or {}).get("decision") for c in enriched)
    n_burn = sum(1 for c in enriched if (c.get("burn") or {}).get("burn"))
    return {
        "n_packs": len(enriched),
        "n_burn": n_burn,
        "decision_counts": dict(counts),
        "cards": enriched,
        "headline": (
            f"burn soft IPTW on {n_burn}/{len(enriched)} packs; "
            f"decisions={dict(counts)}. "
            "Default KEEP_UNIFORM / SOFTEN_ONLY — α is not a FLOPs knob."
        ),
        "policy": (
            "1) Do not burn soft weights to save FLOPs (they don't). "
            "2) Burn only if gated_α beats uniform on sig-MSE. "
            "3) If gate is mandatory and hurts, SOFTEN_ONLY with ∛. "
            "4) Keep the real bill on duty × refit (算力账)."
        ),
    }
