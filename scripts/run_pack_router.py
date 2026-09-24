#!/usr/bin/env python3
"""FW-router scorecard: pack→model vs global HGB.

  PYTHONPATH=. python3 scripts/run_pack_router.py --out results/sandbox_router
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from sandbox.pack_router import router_suite

ROOT = Path(__file__).resolve().parents[1]


def _clean(o):
    if isinstance(o, float) and (o != o or o in (float("inf"), float("-inf"))):
        return None
    if isinstance(o, (np.floating,)):
        v = float(o)
        return None if v != v else v
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, list):
        return [_clean(v) for v in o]
    return o


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=ROOT / "results/sandbox_router")
    ap.add_argument("--max-n", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rep = router_suite(max_n=args.max_n, seed=args.seed)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "router_scorecard.json").write_text(
        json.dumps(_clean(rep), indent=2) + "\n"
    )

    def f(x, nd=4):
        try:
            v = float(x)
            return "—" if v != v else f"{v:.{nd}g}"
        except Exception:
            return "—"

    lines = [
        "# FW-router scorecard",
        "",
        f"**Headline:** {rep.get('headline')}",
        "",
        f"Routing table: `{rep.get('routing_table')}`",
        "",
        "| pack | routed | MAE routed | MAE global HGB | gain |",
        "|---|---|---:|---:|---:|",
    ]
    for r in rep.get("rows") or []:
        lines.append(
            f"| `{r['pack']}` | `{r['routed_model']}` | {f(r['mae_routed'])} | "
            f"{f(r['mae_global_hgb'])} | {f(r['mae_gain_vs_global_hgb'])} |"
        )
    lines += [
        "",
        "## Reading",
        "",
        "- gain>0 ⇒ specialist beats forced HGB on that pack.",
        "- This is the P0 flywheel opportunity: route first, deepen later.",
        "",
    ]
    md = "\n".join(lines) + "\n"
    (args.out / "ROUTER_SCORECARD.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
