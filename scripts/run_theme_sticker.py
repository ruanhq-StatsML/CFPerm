#!/usr/bin/env python3
"""FW-theme sticker scorecard.

  PYTHONPATH=. python3 scripts/run_theme_sticker.py --out results/sandbox_theme
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from sandbox.theme_sticker import theme_sticker_scorecard

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
    ap.add_argument("--out", type=Path, default=ROOT / "results/sandbox_theme")
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--max-n", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rep = theme_sticker_scorecard(max_n=args.max_n, k=args.k, seed=args.seed)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "theme_sticker.json").write_text(json.dumps(_clean(rep), indent=2) + "\n")
    if not rep.get("ok"):
        print(rep)
        return
    lines = [
        "# FW-theme sticker scorecard",
        "",
        f"**Headline:** {rep.get('headline')}",
        "",
        f"- k={rep.get('k')} · fit prompts={rep.get('n_prompts_fit')} · "
        f"stuck steps={rep.get('n_steps_stuck')}",
        f"- paths with theme tag: {rep.get('n_paths_with_theme')}",
        f"- cluster_counts: `{rep.get('cluster_counts')}`",
        "",
        "## Top terms",
        "",
    ]
    for i, terms in enumerate(rep.get("top_terms") or []):
        lines.append(f"- C{i}: {', '.join(terms)}")
    ex = rep.get("example") or {}
    lines += [
        "",
        "## Example annotated step",
        "",
        f"- theme_cluster: {ex.get('theme_cluster')}",
        f"- decision_path: `{ex.get('decision_path')}`",
        f"- text head: {ex.get('theme_text_head')}",
        "",
        "Reading: sticker is context for other agents — not a quality claim.",
        "",
    ]
    md = "\n".join(lines) + "\n"
    (args.out / "THEME_STICKER.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
