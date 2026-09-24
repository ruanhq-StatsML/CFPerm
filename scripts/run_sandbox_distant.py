#!/usr/bin/env python3
"""Far-from-AGOD sandbox: forecast bakeoff + DiffusionDB themes.

  PYTHONPATH=. python3 scripts/run_sandbox_distant.py \\
    --out results/sandbox_distant
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from sandbox.forecast_bakeoff import bakeoff_all
from sandbox.prompt_themes import run_theme_job

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


def render_md(forecast: dict, themes: dict) -> str:
    lines = [
        "# Sandbox distant bakeoff (away from AGOD methods)",
        "",
        "Two unrelated classical ML jobs — **no** excess / PO / soft-burn / claim router.",
        "",
        "## A. Next-step stream forecast",
        "",
        f"- packs ok: {forecast.get('n_packs')}",
        f"- best_rmse counts: `{forecast.get('best_counts')}`",
        f"- mean HGB RMSE lift vs naive: **{forecast.get('mean_hgb_rmse_lift_vs_naive')}**",
        "",
        "| dataset | best | naive RMSE | ridge RMSE | HGB RMSE | HGB lift |",
        "|---|---|---:|---:|---:|---:|",
    ]

    def f(x, nd=4):
        try:
            v = float(x)
            if v != v:
                return "—"
            return f"{v:.{nd}g}"
        except Exception:
            return "—"

    for c in forecast.get("cards") or []:
        if not c.get("ok"):
            continue
        m = c["models"]
        lines.append(
            f"| `{c['dataset']}` | `{c['best_rmse']}` | "
            f"{f(m['naive_last']['rmse'])} | {f(m['ridge']['rmse'])} | "
            f"{f(m['hgb']['rmse'])} | {f(c['rmse_lift_vs_naive']['hgb'])} |"
        )
    lines += [
        "",
        "## B. DiffusionDB prompt themes",
        "",
    ]
    if not themes.get("ok"):
        lines.append(f"- skipped: `{themes.get('reason')}`")
    else:
        lines.append(f"- n_prompts: {themes.get('n_prompts')}")
        lines.append(f"- best k: **{themes.get('best_k')}** (silhouette={themes.get('best_silhouette')})")
        lines.append("- top terms (first 3 clusters):")
        for i, terms in enumerate(themes.get("best_top_terms_head") or []):
            lines.append(f"  - C{i}: {', '.join(terms[:6])}")
    lines += [
        "",
        "## Effect reading (plain)",
        "",
        "- Forecast: if HGB lift ≫ 0 on a pack, lag+X features beat last-value; if ≈0, series is near random-walk.",
        "- Themes: silhouette picks a usable k; terms are descriptive clusters only.",
        "- Explicitly **not** an AGOD efficiency / causal / transfer claim.",
        "",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=ROOT / "results/sandbox_distant")
    ap.add_argument("--max-n", type=int, default=8000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    forecast = bakeoff_all(max_n=args.max_n, seed=args.seed)
    themes = run_theme_job(max_n=min(4000, args.max_n), seed=args.seed)
    blob = {"forecast": forecast, "themes": themes}
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "sandbox_distant.json").write_text(json.dumps(_clean(blob), indent=2) + "\n")
    md = render_md(forecast, themes)
    (args.out / "SANDBOX_DISTANT.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
