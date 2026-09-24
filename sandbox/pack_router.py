"""FW-router: pack → classical model dispatch for the forecast flywheel.

Opportunity
-----------
Bakeoff winners differ by pack (metro→HGB, waymo→ridge, pm25→naive).
A single global HGB wastes compute on RW packs and loses on linear packs.

Router agent
------------
1. Read bakeoff card ``best_rmse`` (or lift thresholds).
2. Dispatch that model into ``run_flywheel``.
3. Score vs a forced global-HGB baseline (MAE / surprise_rate).

Stays in sandbox — not AGOD.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence

import numpy as np

from sandbox.forecast_bakeoff import Pack, bakeoff_all, load_packs, run_bakeoff
from sandbox.forecast_flywheel import run_flywheel


def route_model_from_card(card: Dict[str, Any]) -> str:
    """Map one bakeoff card → specialist model name."""
    if not card.get("ok"):
        return "ridge"  # cheap default
    best = str(card.get("best_rmse") or "ridge")
    if best in ("hgb", "ridge", "naive_last"):
        return best
    lift = float(card.get("rmse_lift_vs_naive", {}).get("hgb", float("nan")))
    if np.isfinite(lift) and lift > 0.15:
        return "hgb"
    if np.isfinite(lift) and lift < 0.02:
        return "naive_last"
    return "ridge"


def build_routing_table(bakeoff: Dict[str, Any]) -> Dict[str, str]:
    """pack name → model."""
    table = {}
    for c in bakeoff.get("cards") or []:
        if not c.get("ok"):
            continue
        table[str(c["dataset"])] = route_model_from_card(c)
    return table


def run_routed_vs_global(
    packs: Sequence[Pack],
    *,
    bakeoff: Optional[Dict[str, Any]] = None,
    seed: int = 0,
    warm: int = 160,
    max_steps: int = 180,
) -> Dict[str, Any]:
    """Compare FW-router specialists vs forcing HGB on every pack."""
    if bakeoff is None:
        # light per-pack bakeoff if caller didn't pass one
        cards = [run_bakeoff(p, seed=seed) for p in packs]
        bakeoff = {"cards": cards, "n_packs": sum(1 for c in cards if c.get("ok"))}
    table = build_routing_table(bakeoff)
    rows = []
    for p in packs:
        routed = table.get(p.name, "ridge")
        r_spec = run_flywheel(
            p, model_name=routed, warm=warm, max_steps=max_steps, seed=seed
        )
        r_glob = run_flywheel(
            p, model_name="hgb", warm=warm, max_steps=max_steps, seed=seed
        )
        if not (r_spec.get("ok") and r_glob.get("ok")):
            continue
        mae_gain = float(r_glob["mae"] - r_spec["mae"])  # >0 ⇒ router better
        rows.append(
            {
                "pack": p.name,
                "routed_model": routed,
                "mae_routed": r_spec["mae"],
                "mae_global_hgb": r_glob["mae"],
                "mae_gain_vs_global_hgb": mae_gain,
                "surprise_routed": r_spec["surprise_rate"],
                "surprise_global_hgb": r_glob["surprise_rate"],
                "fallback_routed": r_spec.get("fell_back_to_naive"),
                "fallback_global_hgb": r_glob.get("fell_back_to_naive"),
            }
        )
    gains = [r["mae_gain_vs_global_hgb"] for r in rows]
    n_win = sum(1 for g in gains if g > 1e-9)
    return {
        "routing_table": table,
        "rows": rows,
        "n_packs": len(rows),
        "n_router_mae_win": n_win,
        "mean_mae_gain_vs_global_hgb": float(np.mean(gains)) if gains else float("nan"),
        "headline": _headline(table, n_win, len(rows), gains),
        "opportunity": "FW-router",
        "note": "specialist dispatch from bakeoff vs forced global HGB",
    }


def _headline(
    table: Dict[str, str], n_win: int, n: int, gains: List[float]
) -> str:
    mean_g = float(np.mean(gains)) if gains else float("nan")
    return (
        f"FW-router table={table}; "
        f"MAE better than global-HGB on {n_win}/{n} packs "
        f"(mean gain={mean_g:.4g})"
    )


def router_suite(*, max_n: int = 5000, seed: int = 0) -> Dict[str, Any]:
    packs = load_packs(max_n=max_n)
    bakeoff = bakeoff_all(max_n=max_n, seed=seed)
    return run_routed_vs_global(packs, bakeoff=bakeoff, seed=seed)
