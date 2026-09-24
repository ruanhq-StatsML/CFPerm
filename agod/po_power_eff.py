"""IPTW power softness: sqrt vs cbrt under the same adaptation FLOPs.

When OnlineRFPerm rejects, reweight with w∝PO^α.  α=1/2 (sqrt) vs α=1/3
(cbrt) have **essentially identical adaptation FLOPs** (weighting is O(n);
the RF re-fit dominates).  So mse_eff collapses to comparing next-MSE —
the soft knob is pure **bias–variance of IPTW**, not compute.

Justify
-------
High PO variance ⇒ heavy sqrt tails over-weight a few hard rows ⇒ next-batch
MSE can *worsen* vs uniform.  cbrt shrinks weight CV → closer to uniform,
often ``gated_cbrt ≤ gated_sqrt`` on sig MSE when the PO ranker is imperfect.

Report: relative MSE vs uniform, pairwise soft win rate, and (for gated
modes) mse_eff using reject-only FLOPs so it sits on the same ruler as
``po_eff.scorecard_from_summary``.
"""
from __future__ import annotations

from collections import Counter
from typing import Any, Dict, List, Mapping, Optional

import numpy as np

from agod.po_eff import expected_adapt_flops, mse_efficiency


def _mse(block: Mapping[str, Any], mode: str) -> float:
    r = (block.get("results") or {}).get(mode) or {}
    v = r.get("mse_mean_sig", r.get("mse_mean"))
    try:
        return float(v)
    except Exception:
        return float("nan")


def _duty(block: Mapping[str, Any]) -> float:
    r = (block.get("results") or {}).get("uniform") or {}
    try:
        return float(r.get("gate_duty", float("nan")))
    except Exception:
        return float("nan")


def power_card_from_dataset(
    name: str,
    block: Mapping[str, Any],
    *,
    n_batches: int,
    batch_size: int,
    n_control: int = 1,
) -> Dict[str, Any]:
    modes = ["uniform", "sqrt", "cbrt", "gated_sqrt", "gated_cbrt", "dre"]
    unif = _mse(block, "uniform")
    duty = _duty(block)
    # gated IPTW pays reject-only refit-like cost (same order as po_eff refit)
    e_flops = expected_adapt_flops(
        "refit",
        gate_duty=duty if np.isfinite(duty) else 0.0,
        n_batches=n_batches,
        batch_size=batch_size,
        n_control=n_control,
    )
    rows = []
    for m in modes:
        mse = _mse(block, m)
        rel = float(mse / unif) if np.isfinite(mse) and np.isfinite(unif) and unif != 0 else float("nan")
        rows.append(
            {
                "mode": m,
                "mse_mean_sig": mse,
                "rel_vs_uniform": rel,
                "mse_eff": mse_efficiency(mse, unif, e_flops)
                if m.startswith("gated_")
                else float("nan"),
            }
        )
    by = {r["mode"]: r for r in rows}
    gs, gc = by.get("gated_sqrt", {}), by.get("gated_cbrt", {})
    soft_win = None
    if np.isfinite(gs.get("mse_mean_sig", np.nan)) and np.isfinite(
        gc.get("mse_mean_sig", np.nan)
    ):
        soft_win = bool(gc["mse_mean_sig"] <= gs["mse_mean_sig"])
    # best among gated + uniform on sig MSE
    cand = {
        m: by[m]["mse_mean_sig"]
        for m in ("uniform", "gated_sqrt", "gated_cbrt")
        if np.isfinite(by.get(m, {}).get("mse_mean_sig", np.nan))
    }
    best = min(cand, key=cand.get) if cand else None
    reading = _power_reading(soft_win, best, duty)
    return {
        "dataset": name,
        "ok": True,
        "gate_duty": duty,
        "modes": rows,
        "soft_win_cbrt_le_sqrt": soft_win,
        "best_sig_mode": best,
        "reading": reading,
    }


def _power_reading(
    soft_win: Optional[bool], best: Optional[str], duty: float
) -> str:
    duty_s = f"duty={duty:.2f}" if np.isfinite(duty) else "duty=?"
    if best == "uniform":
        base = f"{duty_s}: keep uniform — gated IPTW does not buy sig MSE"
    elif best == "gated_cbrt":
        base = f"{duty_s}: gated∛ best sig MSE (softer weights)"
    elif best == "gated_sqrt":
        base = f"{duty_s}: gated√ best sig MSE"
    else:
        base = f"{duty_s}: best={best}"
    if soft_win is True:
        return base + "; ∛≤√ on this pack"
    if soft_win is False:
        return base + "; √ beats ∛ here"
    return base


def power_scorecard_from_summary(summary: Mapping[str, Any]) -> Dict[str, Any]:
    n_batches = int(summary.get("n_batches") or 40)
    batch_size = int(summary.get("batch_size") or 100)
    n_control = int(summary.get("n_control") or 1)
    cards: List[Dict[str, Any]] = []
    for name, block in (summary.get("datasets") or {}).items():
        if isinstance(block, dict):
            cards.append(
                power_card_from_dataset(
                    name,
                    block,
                    n_batches=n_batches,
                    batch_size=batch_size,
                    n_control=n_control,
                )
            )
    ok = [c for c in cards if c.get("ok")]
    soft = [c["soft_win_cbrt_le_sqrt"] for c in ok if c.get("soft_win_cbrt_le_sqrt") is not None]
    soft_rate = float(np.mean(soft)) if soft else float("nan")
    best_counts = Counter(c["best_sig_mode"] for c in ok if c.get("best_sig_mode"))
    return {
        "n_datasets": len(ok),
        "n_batches": n_batches,
        "batch_size": batch_size,
        "soft_win_rate_cbrt_le_sqrt": soft_rate,
        "best_sig_counts": dict(best_counts),
        "cards": cards,
        "headline": (
            f"gated∛≤gated√ on {soft_rate:.0%} of packs; "
            f"best_sig counts {dict(best_counts)}. "
            "Same FLOPs — softness is the only knob."
        ),
        "note": (
            "sqrt vs cbrt IPTW share adaptation FLOPs; compare rel MSE / soft wins, "
            "not compute. Prefer ∛ when PO ranks well but √ IPTW hurts."
        ),
    }
