#!/usr/bin/env python3
"""FW-pm25-naive regime watcher scorecard.

  PYTHONPATH=. python3 scripts/run_naive_regime_watch.py --out results/sandbox_naive_watch
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from sandbox.naive_regime_watch import suite_from_loaders

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
    ap.add_argument("--out", type=Path, default=ROOT / "results/sandbox_naive_watch")
    ap.add_argument("--max-n", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rep = suite_from_loaders(max_n=args.max_n, seed=args.seed)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "naive_watch.json").write_text(json.dumps(_clean(rep), indent=2) + "\n")

    def f(x, nd=4):
        try:
            v = float(x)
            return "—" if v != v else f"{v:.{nd}g}"
        except Exception:
            return "—"

    lines = [
        "# FW-pm25-naive regime watch",
        "",
        f"**Headline:** {rep.get('headline')}",
        "",
        "| pack | router | MAE | surprise_rate | regime_alerts |",
        "|---|---|---:|---:|---:|",
    ]
    for r in rep.get("rows") or []:
        lines.append(
            f"| `{r['dataset']}` | `{r['router_model']}` | {f(r['mae'])} | "
            f"{f(r['surprise_rate'])} | {r['n_regime_alerts']} |"
        )
    lines += [
        "",
        "## Reading",
        "",
        "- On naive-routed packs: keep last-value; `regime_alert` ⇒ re-run bakeoff/router.",
        "- Do not promote HGB just because surprises exist — check streak first.",
        "",
    ]
    md = "\n".join(lines) + "\n"
    (args.out / "NAIVE_WATCH.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
