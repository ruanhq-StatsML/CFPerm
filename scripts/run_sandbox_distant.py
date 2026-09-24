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

from sandbox.forecast_bakeoff import bakeoff_all, load_packs
from sandbox.forecast_flywheel import flywheel_suite, opportunity_map_from_bakeoff
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


def render_md(forecast: dict, themes: dict, flywheel: dict | None = None) -> str:
    lines = [
        "# Sandbox distant bakeoff (away from AGOD methods)",
        "",
        "Classical ML jobs — **no** excess / PO / soft-burn / claim router.",
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
        lines.append(
            f"- best k: **{themes.get('best_k')}** "
            f"(silhouette={themes.get('best_silhouette')})"
        )
        lines.append("- top terms (first 3 clusters):")
        for i, terms in enumerate(themes.get("best_top_terms_head") or []):
            lines.append(f"  - C{i}: {', '.join(terms[:6])}")

    fw = flywheel or {}
    lines += ["", "## C. Classical forecast → agent flywheel", ""]
    if fw:
        lines.append(f"**Headline:** {fw.get('headline')}")
        lines.append("")
        lines.append("| id | pack | opportunity | priority |")
        lines.append("|---|---|---|---|")
        for o in fw.get("opportunities") or []:
            lines.append(
                f"| `{o.get('id')}` | `{o.get('pack')}` | {o.get('opportunity')} | "
                f"{o.get('priority')} |"
            )
        lines.append("")
        lines.append("| pack | model | MAE | surprise_rate | retrains | fallback |")
        lines.append("|---|---|---:|---:|---:|:---:|")
        for r in fw.get("runs") or []:
            if not r.get("ok"):
                continue
            lines.append(
                f"| `{r.get('dataset')}` | `{r.get('model')}` | {f(r.get('mae'))} | "
                f"{f(r.get('surprise_rate'))} | {r.get('n_retrain')} | "
                f"{'Y' if r.get('fell_back_to_naive') else 'N'} |"
            )
    lines += [
        "",
        "## Effect reading (plain)",
        "",
        "- Forecast: HGB lift ≫ 0 ⇒ learnable pack; ≈0 ⇒ random-walk (naive agent).",
        "- Flywheel: opportunity is **pack→model routing + surprise/retrain/fallback**, not deeper nets.",
        "- Themes: silhouette/terms only — context sticker for other agents.",
        "- Not an AGOD efficiency / causal claim.",
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
    packs = load_packs(max_n=args.max_n)
    flywheel = flywheel_suite(forecast, packs, seed=args.seed)
    blob = {"forecast": forecast, "themes": themes, "flywheel": flywheel}
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "sandbox_distant.json").write_text(
        json.dumps(_clean(blob), indent=2) + "\n"
    )
    md = render_md(forecast, themes, flywheel)
    (args.out / "SANDBOX_DISTANT.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
