#!/usr/bin/env python3
"""Benchmark FSDS/CFPerm sample + observation-subset weights on real/synth boards.

  python3 scripts/run_fsds_cfperm_weight_benchmark.py --quick
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from fsds_cfperm_weight_benchmark import run_benchmark  # noqa: E402

OUT = ROOT / "results" / "fsds_cfperm_weight_bench"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")


def _fmt(x):
    if x is None:
        return "—"
    try:
        return "%.3f" % float(x)
    except Exception:
        return "—"


def _readme(payload: dict) -> str:
    lines = [
        "# FSDS/CFPerm weight benchmark",
        "",
        payload["takeaway"],
        "",
        "## Aggregate (vs uniform)",
        "",
        "| channel | mean Δ novel | mean Δ acc |",
        "| --- | --- | --- |",
        "| sample_w (membership) | %s | — |" % _fmt(payload.get("mean_delta_novel_sample")),
        "| subset_w (hard obs) | %s | — |" % _fmt(payload.get("mean_delta_novel_subset")),
        "| sample × subset | %s | %s |" % (
            _fmt(payload.get("mean_delta_novel_both")),
            _fmt(payload.get("mean_delta_acc_both")),
        ),
        "| sample_w fused | — | %s |" % _fmt(payload.get("mean_delta_acc_sample_fused")),
        "| subset_w fused | — | %s |" % _fmt(payload.get("mean_delta_acc_subset_fused")),
        "| both fused | %s | %s |" % (
            _fmt(payload.get("mean_delta_novel_both_fused")),
            _fmt(payload.get("mean_delta_acc_both_fused")),
        ),
        "",
        "## Per dataset",
        "",
        "| dataset | n_ref/n_live | subset_frac | vimp_top | Δacc both | Δacc both_fused | Δnovel both_fused |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for r in payload["boards"]:
        lines.append(
            "| %s | %d/%d | %s | %s | %s | %s | %s |" % (
                r["dataset"], r["n_ref"], r["n_live"],
                _fmt(r["subset_frac"]), _fmt(r["vimp_top_mass"]),
                _fmt(r.get("delta_acc_both")), _fmt(r.get("delta_acc_both_fused")),
                _fmt(r.get("delta_novel_both_fused")),
            )
        )
    lines.extend([
        "",
        "### Scheme detail (acc / novel / forgetting)",
        "",
    ])
    for r in payload["boards"]:
        lines.append("**%s** — %s" % (r["dataset"], r.get("note") or ""))
        lines.append("")
        lines.append("| scheme | acc | novel | known | forgetting | weight_gini |")
        lines.append("| --- | --- | --- | --- | --- | --- |")
        for name, s in r["schemes"].items():
            lines.append("| %s | %s | %s | %s | %s | %s |" % (
                name, _fmt(s.get("acc")), _fmt(s.get("novel_acc")),
                _fmt(s.get("known_acc")), _fmt(s.get("forgetting")), _fmt(s.get("weight_gini")),
            ))
        lines.append("")
    lines.extend([
        "Weights: `sample_w = 1+λ·membership`; `subset_w` boosts observations with membership ≥ q.",
        "",
        "```bash",
        "python3 scripts/run_fsds_cfperm_weight_benchmark.py --quick",
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
    payload = run_benchmark(seed=int(args.seed), quick=bool(args.quick))
    (OUT / "fsds_cfperm_weight_bench.json").write_text(json.dumps(payload, indent=2, default=str))
    readme = _readme(payload)
    (OUT / "README.md").write_text(readme)
    DOCS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "README.md", DOCS / "FSDS_CFPerm_weight_benchmark.md")
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "README.md", ARTIFACTS / "fsds_cfperm_weight_bench_README.md")
    shutil.copy(OUT / "fsds_cfperm_weight_bench.json", ARTIFACTS / "fsds_cfperm_weight_bench.json")
    print("[wbench] boards=%d  Δacc both_fused=%s  Δnovel both_fused=%s" % (
        payload["n_boards"], _fmt(payload.get("mean_delta_acc_both_fused")),
        _fmt(payload.get("mean_delta_novel_both_fused")),
    ), flush=True)
    for r in payload["boards"]:
        print("[wbench] %-10s  Δacc_fused=%s  Δnovel_fused=%s  subset_frac=%s" % (
            r["dataset"], _fmt(r.get("delta_acc_both_fused")),
            _fmt(r.get("delta_novel_both_fused")), _fmt(r["subset_frac"]),
        ), flush=True)
    print("[wbench] wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
