"""Grad-OnlineRFPerm freeze closed-loop: MSE vs FLOPs efficiency.

After the first Grad reject, policies freeze some MLP layers while adapting.
Relative to ``always_adapt``:

  mse_rel   = MSE_post(policy) / MSE_post(always)
  flops_rel = FLOPs_proxy(policy) / FLOPs_proxy(always)

``freeze_eff = (1 - mse_rel) / flops_rel``
  >0  ⇒ better MSE *and* we normalize by remaining compute
  <0  ⇒ MSE worse than always_adapt (FLOPs save alone is not a win)

Pareto: mse_rel < 1 and flops_rel < 1 ⇒ dominates always_adapt on both axes.
"""
from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional

import numpy as np

# Documented aggregates from Grad_OnlineRFPerm_extras (multi-seed means)
CANONICAL_FREEZE: Dict[str, Dict[str, Dict[str, float]]] = {
    "electricity": {
        "always_adapt": {"mse_rel": 1.0, "flops_rel": 1.0},
        "freeze_early": {"mse_rel": 0.86, "flops_rel": 0.70},
        "freeze_low_share": {"mse_rel": 0.86, "flops_rel": 0.70},
        "no_adapt": {"mse_rel": float("nan"), "flops_rel": 0.0},  # collapses
    },
    "synthetic": {
        "always_adapt": {"mse_rel": 1.0, "flops_rel": 1.0},
        "freeze_early": {"mse_rel": 1.25, "flops_rel": 0.70},  # worse than low_share
        "freeze_low_share": {"mse_rel": 1.15, "flops_rel": 0.57},
        "no_adapt": {"mse_rel": float("nan"), "flops_rel": 0.0},
    },
}


def freeze_efficiency(mse_rel: float, flops_rel: float) -> float:
    if not np.isfinite(mse_rel) or not np.isfinite(flops_rel) or flops_rel <= 0:
        return float("nan")
    return float((1.0 - mse_rel) / flops_rel)


def dominates_always(mse_rel: float, flops_rel: float) -> bool:
    return bool(
        np.isfinite(mse_rel)
        and np.isfinite(flops_rel)
        and mse_rel < 1.0
        and flops_rel < 1.0
    )


def rows_from_rel_table(
    dataset: str, table: Mapping[str, Mapping[str, float]]
) -> Dict[str, Any]:
    rows: List[Dict[str, Any]] = []
    for pol, vals in table.items():
        mr = float(vals.get("mse_rel", float("nan")))
        fr = float(vals.get("flops_rel", float("nan")))
        rows.append(
            {
                "policy": pol,
                "mse_rel": mr,
                "flops_rel": fr,
                "freeze_eff": freeze_efficiency(mr, fr),
                "dominates_always": dominates_always(mr, fr),
            }
        )
    spend = [
        r
        for r in rows
        if r["policy"] not in ("always_adapt", "no_adapt")
        and np.isfinite(r["freeze_eff"])
    ]
    best = max(spend, key=lambda r: r["freeze_eff"])["policy"] if spend else None
    dom = [r["policy"] for r in rows if r["dominates_always"]]
    return {
        "dataset": dataset,
        "ok": True,
        "policies": rows,
        "best_freeze_eff": best,
        "dominates_always": dom,
        "reading": _read(dataset, best, dom),
    }


def _read(dataset: str, best: Optional[str], dom: List[str]) -> str:
    if dom:
        return f"{dataset}: {', '.join(dom)} Pareto-beat always_adapt (MSE↓ & FLOPs↓)"
    if best:
        return f"{dataset}: best freeze_eff={best} (may trade MSE for FLOPs)"
    return f"{dataset}: no usable freeze policy"


def scorecard_from_freeze_bundle(bundle: Mapping[str, Any]) -> Dict[str, Any]:
    """Parse ``exp_freeze_loop`` output: {ds: {_agg: {pol: {mse_post_mean, flops_mean}}}}."""
    cards = []
    for ds, block in bundle.items():
        if not isinstance(block, dict):
            continue
        agg = block.get("_agg") or block
        if "always_adapt" not in agg:
            continue
        base_m = float(agg["always_adapt"].get("mse_post_mean", float("nan")))
        base_f = float(agg["always_adapt"].get("flops_mean", float("nan")))
        table = {}
        for pol, v in agg.items():
            if not isinstance(v, dict):
                continue
            m = float(v.get("mse_post_mean", float("nan")))
            f = float(v.get("flops_mean", float("nan")))
            table[pol] = {
                "mse_rel": float(m / base_m) if np.isfinite(m) and base_m else float("nan"),
                "flops_rel": float(f / base_f) if np.isfinite(f) and base_f else float("nan"),
            }
        cards.append(rows_from_rel_table(ds, table))
    return _finalize(cards, source="bundle")


def scorecard_canonical() -> Dict[str, Any]:
    cards = [rows_from_rel_table(ds, tab) for ds, tab in CANONICAL_FREEZE.items()]
    return _finalize(cards, source="canonical_extras_doc")


def _finalize(cards: List[Dict[str, Any]], *, source: str) -> Dict[str, Any]:
    from collections import Counter

    ok = [c for c in cards if c.get("ok")]
    best = Counter(c["best_freeze_eff"] for c in ok if c.get("best_freeze_eff"))
    n_dom = sum(len(c.get("dominates_always") or []) for c in ok)
    return {
        "source": source,
        "n_datasets": len(ok),
        "best_freeze_eff_counts": dict(best),
        "n_pareto_dominances": n_dom,
        "cards": cards,
        "headline": (
            f"freeze_eff best {dict(best)}; "
            f"{n_dom} Pareto wins vs always_adapt (MSE↓&FLOPs↓). "
            "no_adapt collapses — freeze ≠ stop learning."
        ),
    }
