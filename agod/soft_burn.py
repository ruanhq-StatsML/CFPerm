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
        "alpha_changes_flops": 0.0,  # √ vs ∛ identical O(n)
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
        # Soften is a *policy*: always ∛ to shrink weight tails (same FLOPs).
        # Do not flip to √ just because √ hurts slightly less among losers.
        return {
            "decision": "SOFTEN_ONLY",
            "burn": False,  # do not claim IPTW win
            "alpha": 1.0 / 3.0,
            "reason": (
                "gated IPTW does not beat uniform; if reject path is mandatory, "
                "use ∛ to limit damage (same FLOPs as √)"
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


def _resolve_adapt_flops(
    card: Mapping[str, Any],
    *,
    n_batches: Optional[int] = None,
    batch_size: Optional[int] = None,
    n_control: Optional[int] = None,
) -> float:
    """Prefer card.expected_adapt_flops; else recompute from duty · refit."""
    raw = card.get("expected_adapt_flops")
    try:
        v = float(raw) if raw is not None else float("nan")
    except Exception:
        v = float("nan")
    if np.isfinite(v):
        return v
    duty = card.get("gate_duty")
    try:
        duty_f = float(duty) if duty is not None else float("nan")
    except Exception:
        duty_f = float("nan")
    if not np.isfinite(duty_f):
        return float("nan")
    from agod.po_eff import expected_adapt_flops

    nb = int(n_batches if n_batches is not None else card.get("n_batches") or 40)
    bs = int(batch_size if batch_size is not None else card.get("batch_size") or 100)
    nc = int(n_control if n_control is not None else card.get("n_control") or 1)
    return float(
        expected_adapt_flops(
            "refit",
            gate_duty=duty_f,
            n_batches=nb,
            batch_size=bs,
            n_control=nc,
        )
    )


def attach_burn_to_power_card(
    card: Mapping[str, Any],
    *,
    n_batches: Optional[int] = None,
    batch_size: Optional[int] = None,
    n_control: Optional[int] = None,
    n_rows: Optional[int] = None,
) -> Dict[str, Any]:
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
    adapt = _resolve_adapt_flops(
        card, n_batches=n_batches, batch_size=batch_size, n_control=n_control
    )
    # Default n_rows ≈ stream length if caller omits: n_batches · batch_size
    if n_rows is None:
        nb = int(n_batches if n_batches is not None else card.get("n_batches") or 0)
        bs = int(batch_size if batch_size is not None else card.get("batch_size") or 0)
        n_rows = int(nb * bs) if nb > 0 and bs > 0 else 0
    ledger = soft_flops_ledger(adapt_flops=adapt, n_rows=int(n_rows))
    ledger["duty"] = card.get("gate_duty")
    out["burn"] = decision
    out["flops_ledger"] = ledger
    return out


def summarize_burn_decisions(
    cards: Sequence[Mapping[str, Any]],
    *,
    n_batches: Optional[int] = None,
    batch_size: Optional[int] = None,
    n_control: Optional[int] = None,
    n_rows: Optional[int] = None,
) -> Dict[str, Any]:
    from collections import Counter

    enriched = [
        attach_burn_to_power_card(
            c,
            n_batches=n_batches,
            batch_size=batch_size,
            n_control=n_control,
            n_rows=n_rows,
        )
        for c in cards
        if c.get("ok")
    ]
    counts = Counter((c.get("burn") or {}).get("decision") for c in enriched)
    n_burn = sum(1 for c in enriched if (c.get("burn") or {}).get("burn"))
    shares = [
        float((c.get("flops_ledger") or {}).get("weighting_share", float("nan")))
        for c in enriched
    ]
    finite_shares = [s for s in shares if np.isfinite(s)]
    mean_w_share = float(np.mean(finite_shares)) if finite_shares else float("nan")
    share_s = f"{mean_w_share:.2e}" if np.isfinite(mean_w_share) else "n/a"
    return {
        "n_packs": len(enriched),
        "n_burn": n_burn,
        "decision_counts": dict(counts),
        "mean_weighting_share": mean_w_share,
        "cards": enriched,
        "headline": (
            f"burn soft IPTW on {n_burn}/{len(enriched)} packs; "
            f"decisions={dict(counts)}; mean weighting_share={share_s}. "
            "Default KEEP_UNIFORM / SOFTEN_ONLY — α is not a FLOPs knob."
        ),
        "policy": (
            "1) Do not burn soft weights to save FLOPs (they don't). "
            "2) Burn only if gated_α beats uniform on sig-MSE. "
            "3) If gate is mandatory and hurts, SOFTEN_ONLY with ∛. "
            "4) Keep the real bill on duty × refit (算力账)."
        ),
    }
