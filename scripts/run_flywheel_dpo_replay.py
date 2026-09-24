#!/usr/bin/env python3
"""Pollute labels → shuffle decision_path replay → DPO feed flywheel.

  PYTHONPATH=. python3 scripts/run_flywheel_dpo_replay.py \\
    --out results/sandbox_dpo_replay
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from sandbox.forecast_bakeoff import load_packs
from sandbox.flywheel_dpo_replay import run_pollute_replay_dpo

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
    ap.add_argument("--out", type=Path, default=ROOT / "results/sandbox_dpo_replay")
    ap.add_argument("--max-n", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--pack", type=str, default="metro_interstate")
    args = ap.parse_args()
    packs = {p.name: p for p in load_packs(max_n=args.max_n)}
    if args.pack not in packs:
        # fallback first pack
        pack = next(iter(packs.values()))
    else:
        pack = packs[args.pack]
    # model hint from name
    model = "hgb" if "metro" in pack.name else ("ridge" if "waymo" in pack.name else "naive_last")
    if model == "naive_last":
        model = "hgb"  # still run a learner to see pollution bite
    out_pack = args.out / pack.name
    rep = run_pollute_replay_dpo(
        pack, model_name=model, seed=args.seed, out_dir=out_pack
    )
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(_clean(rep), indent=2) + "\n")
    md = [
        "# Flywheel DPO + label pollution + path replay",
        "",
        f"**Pack:** `{rep['dataset']}` · model=`{rep['model']}`",
        "",
        f"**Reading:** {rep.get('reading')}",
        "",
        "## Pollution",
        "",
        f"- method: `{rep['pollution'].get('method')}`",
        f"- n_polluted indices: {rep['pollution'].get('n_polluted')}",
        f"- sampled periods: {rep['sampled_periods']}",
        "",
        "## Rewards",
        "",
        f"- clean: {rep['rewards']['clean']}",
        f"- polluted: {rep['rewards']['polluted']}",
        f"- replay shuffled path: {rep['rewards']['replay_shuffled_path']}",
        "",
        "## DPO",
        "",
        f"- n_pairs: {rep['dpo'].get('n_pairs')}",
        f"- mean_loss: {rep['dpo'].get('mean_loss')}",
        f"- mean_margin: {rep['dpo'].get('mean_margin')}",
        f"- frac_chosen_better: {rep['dpo'].get('frac_chosen_better')}",
        "",
        "## Decision path replay example",
        "",
        f"- orig: `{rep['replay_path_example'].get('orig')}`",
        f"- shuffled: `{rep['replay_path_example'].get('shuffled')}`",
        "",
        f"Step JSONL under `{out_pack}/`.",
        "",
    ]
    text = "\n".join(md) + "\n"
    (args.out / "DPO_REPLAY.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
