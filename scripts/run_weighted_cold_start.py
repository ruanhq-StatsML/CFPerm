#!/usr/bin/env python3
"""Weighted cold-start under category shift (CIFAR C-GCD + Office-Home PDA proxies).

  python3 scripts/run_weighted_cold_start.py --quick
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from category_shift_boards import run_suite  # noqa: E402

OUT = ROOT / "results" / "weighted_cold_start"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")


def _readme(payload: dict) -> str:
    cg_w, cg_u, oh = payload["cgcd_weighted"], payload["cgcd_uniform"], payload["officehome"]
    lines = [
        "# Weighted cold-start under category shift",
        "",
        payload["takeaway"],
        "",
        "## CIFAR-100 C-GCD proxy (known → staged novel influx)",
        "",
        "| method | mean known acc | mean novel acc | mean forgetting | feat hit |",
        "| --- | --- | --- | --- | --- |",
        "| weighted | %.3f | %.3f | %.3f | %.3f |" % (
            cg_w["mean_known_acc"], cg_w["mean_novel_acc"], cg_w["mean_forgetting"], cg_w["mean_feat_hit"],
        ),
        "| uniform | %.3f | %.3f | %.3f | — |" % (
            cg_u["mean_known_acc"], cg_u["mean_novel_acc"], cg_u["mean_forgetting"],
        ),
        "",
        "Δ novel acc (weighted − uniform): **%.3f**" % payload["delta_novel_acc"],
        "",
        "Δ forgetting (uniform − weighted; higher = weighted forgets less): **%.3f**" % payload["delta_forgetting"],
        "",
        "## Office-Home PDA proxy (partial label overlap)",
        "",
        "| overlap | shared acc | incomparable rate | adapt rate |",
        "| --- | --- | --- | --- |",
    ]
    for ov, v in oh["by_overlap"].items():
        lines.append("| %s | %.3f | %.3f | %.3f |" % (ov, v["mean_shared_acc"], v["incomparable_rate"], v["adapt_rate"]))
    lines.extend([
        "",
        "Weights: sample `w=1+λ·membership` in top-K VIMP subspace; soft-freeze non-selected features in `partial_fit`.",
        "",
        "Boundary: Jaccard(label_ref, label_live) < 0.2 → **incomparable** (skip directed adapt).",
        "",
        "```bash",
        "python3 scripts/run_weighted_cold_start.py --quick",
        "```",
        "",
    ])
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    payload = run_suite(seed=int(args.seed), quick=bool(args.quick))
    (OUT / "weighted_cold_start.json").write_text(json.dumps(payload, indent=2, default=str))
    readme = _readme(payload)
    (OUT / "README.md").write_text(readme)
    DOCS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "README.md", DOCS / "Weighted_cold_start_category_shift.md")
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "README.md", ARTIFACTS / "weighted_cold_start_README.md")
    shutil.copy(OUT / "weighted_cold_start.json", ARTIFACTS / "weighted_cold_start.json")
    print("[wcs] Δnovel=%.3f  Δforget=%.3f" % (payload["delta_novel_acc"], payload["delta_forgetting"]), flush=True)
    for ov, v in payload["officehome"]["by_overlap"].items():
        print("[wcs] overlap=%s  shared=%.3f  incomparable=%.2f" % (ov, v["mean_shared_acc"], v["incomparable_rate"]), flush=True)
    print("[wcs] wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
