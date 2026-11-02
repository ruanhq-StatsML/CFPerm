#!/usr/bin/env python3
"""Downsample vs FSDS weight ablation.

  python3 scripts/run_downsample_vs_weight_ablation.py --quick
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from downsample_vs_weight_ablation import run_downsample_ablation  # noqa: E402

OUT = ROOT / "results" / "downsample_vs_weight"
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
        "# Downsample vs FSDS sample-weight ablation",
        "",
        payload["takeaway"],
        "",
        "## Aggregate (novel acc deltas)",
        "",
        "| contrast | mean Δ novel |",
        "| --- | --- |",
        "| ds_known − uniform | %s |" % _fmt(payload.get("mean_delta_ds_known_vs_uniform")),
        "| ds_random − uniform | %s |" % _fmt(payload.get("mean_delta_ds_random_vs_uniform")),
        "| ds_mem − uniform | %s |" % _fmt(payload.get("mean_delta_ds_mem_vs_uniform")),
        "| best_weight − uniform | %s |" % _fmt(payload.get("mean_delta_best_weight_vs_uniform")),
        "| best_weight − ds_known (**incremental**) | %s |" % _fmt(payload.get("mean_delta_best_weight_vs_ds_known")),
        "| ds_known+sample_f − ds_known | %s |" % _fmt(payload.get("mean_delta_combo_vs_ds_known")),
        "",
        "Weights beat ds_known on **%d / %d** boards." % (
            payload.get("n_weights_beat_ds", 0), payload.get("n_boards", 0),
        ),
        "",
        "## Per board",
        "",
        "| dataset | uniform | ds_known | ds_mem | best_w | best−ds | knownΔ@best | verdict |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for r in payload["boards"]:
        u = r["schemes"]["uniform"]
        dsk = r["schemes"]["ds_known"]
        dsm = r["schemes"]["ds_mem"]
        bw = r["schemes"][r["best_weight_scheme"]]
        verd = "weights>ds" if r["weights_beat_downsample"] else (
            "ds enough" if r["delta_novel_ds_known_vs_uniform"] > 0.02 else "neither"
        )
        lines.append("| %s | %s | %s | %s | %s (%s) | %s | %s | %s |" % (
            r["dataset"], _fmt(u["novel_acc"]), _fmt(dsk["novel_acc"]),
            _fmt(dsm["novel_acc"]), _fmt(bw["novel_acc"]), r["best_weight_scheme"],
            _fmt(r["delta_novel_best_weight_vs_ds_known"]),
            _fmt(bw["known_acc"] - u["known_acc"]), verd,
        ))
    lines.extend([
        "",
        "### Scheme detail",
        "",
    ])
    for r in payload["boards"]:
        lines.append("**%s** — %s" % (r["dataset"], r.get("note") or ""))
        lines.append("")
        lines.append("| scheme | novel | known | forgetting | n_train | novel_frac |")
        lines.append("| --- | --- | --- | --- | --- | --- |")
        for name, s in r["schemes"].items():
            lines.append("| %s | %s | %s | %s | %d | %s |" % (
                name, _fmt(s.get("novel_acc")), _fmt(s.get("known_acc")),
                _fmt(s.get("forgetting")), s.get("n_train", 0),
                _fmt(s.get("novel_frac_train")),
            ))
        lines.append("")
    lines.extend([
        "```bash",
        "python3 scripts/run_downsample_vs_weight_ablation.py --quick",
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
    payload = run_downsample_ablation(seed=int(args.seed), quick=bool(args.quick))
    (OUT / "downsample_vs_weight.json").write_text(json.dumps(payload, indent=2, default=str))
    readme = _readme(payload)
    (OUT / "README.md").write_text(readme)
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "README.md", ARTIFACTS / "downsample_vs_weight_README.md")
    shutil.copy(OUT / "downsample_vs_weight.json", ARTIFACTS / "downsample_vs_weight.json")
    print("[ds-ablate] boards=%d  mean(best_w−ds_known)=%s  beat_ds=%d/%d" % (
        payload["n_boards"], _fmt(payload.get("mean_delta_best_weight_vs_ds_known")),
        payload.get("n_weights_beat_ds", 0), payload.get("n_boards", 0),
    ), flush=True)
    for r in payload["boards"]:
        print("[ds-ablate] %-8s  ds−u=%s  best_w−u=%s  best_w−ds=%s  best=%s" % (
            r["dataset"],
            _fmt(r["delta_novel_ds_known_vs_uniform"]),
            _fmt(r["delta_novel_best_weight_vs_uniform"]),
            _fmt(r["delta_novel_best_weight_vs_ds_known"]),
            r["best_weight_scheme"],
        ), flush=True)
    print("[ds-ablate] wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
