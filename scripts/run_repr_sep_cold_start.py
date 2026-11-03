#!/usr/bin/env python3
"""Representation-first cold-start board (4-step pipeline) + streaming performance.

  python3 scripts/run_repr_sep_cold_start.py --quick
  python3 scripts/run_repr_sep_cold_start.py --multiseed   # streaming performance
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from category_shift_boards import cold_start_efficiency, run_cgcd  # noqa: E402

OUT = ROOT / "results" / "repr_sep_cold_start"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")


def _fmt(x):
    try:
        return "%.3f" % float(x)
    except Exception:
        return "—"


def _run_streaming_performance(seeds, *, quick: bool, expand_min: int) -> dict:
    """Full streaming C-GCD performance table (primary deliverable)."""
    modes = ("repr_sep", "weighted", "membership", "uniform")
    agg = {m: {"known": [], "novel": [], "forget": [], "stages": []} for m in modes}
    eff_ns = {m: [] for m in modes}
    for seed in seeds:
        for mode in modes:
            r = run_cgcd(
                seed=int(seed), mode=mode, quick=bool(quick), expand_min=int(expand_min),
            )
            agg[mode]["known"].append(r["mean_known_acc"])
            agg[mode]["novel"].append(r["mean_novel_acc"])
            agg[mode]["forget"].append(r["mean_forgetting"])
            agg[mode]["stages"].append(r["stages"])
        eff = cold_start_efficiency(seed=int(seed), quick=bool(quick), target=0.60)
        for m, cur in eff["curves"].items():
            eff_ns[m].append(cur["n_star"])

    def mean_std(xs):
        xs = [x for x in xs if x == x]
        return float(np.mean(xs)), float(np.std(xs))

    summary = []
    for m in modes:
        kn, ks = mean_std(agg[m]["known"])
        nn, ns = mean_std(agg[m]["novel"])
        fn, fs = mean_std(agg[m]["forget"])
        n_stages = len(agg[m]["stages"][0])
        stage_novel = []
        for s in range(n_stages):
            vals = [agg[m]["stages"][i][s].get("novel_acc") for i in range(len(seeds))]
            vals = [v for v in vals if v is not None and v == v]
            stage_novel.append(float(np.mean(vals)) if vals else float("nan"))
        nstars = [x for x in eff_ns[m] if x is not None]
        summary.append({
            "mode": m,
            "known_mean": kn, "known_std": ks,
            "novel_mean": nn, "novel_std": ns,
            "forget_mean": fn, "forget_std": fs,
            "stage_novel_mean": stage_novel,
            "n_star_mean": float(np.mean(nstars)) if nstars else None,
            "n_star_list": eff_ns[m],
        })
    u = next(s for s in summary if s["mode"] == "uniform")
    for s in summary:
        s["delta_novel_vs_uniform"] = s["novel_mean"] - u["novel_mean"]
        s["delta_forget_vs_uniform"] = u["forget_mean"] - s["forget_mean"]
    return {
        "protocol": "C-GCD streaming proxy",
        "quick": bool(quick),
        "seeds": list(seeds),
        "kpi": "primary=novel_acc; secondary=known/forgetting; efficiency=n* to novel>=0.60",
        "summary": summary,
    }


def _streaming_md(perf: dict) -> str:
    lines = [
        "# Streaming C-GCD proxy — performance",
        "",
        "Protocol: %s (quick=%s). Seeds %s." % (
            perf["protocol"], perf["quick"], perf["seeds"],
        ),
        "KPI: **novel_acc** primary; known retention / forgetting; efficiency n*.",
        "",
        "| mode | known | novel | forgetting | Δnovel vs uniform | n* (mean) |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for s in perf["summary"]:
        nst = "—" if s["n_star_mean"] is None else "%.1f" % s["n_star_mean"]
        lines.append(
            "| %s | %.3f±%.3f | %.3f±%.3f | %.3f±%.3f | %+.3f | %s |" % (
                s["mode"], s["known_mean"], s["known_std"],
                s["novel_mean"], s["novel_std"],
                s["forget_mean"], s["forget_std"],
                s["delta_novel_vs_uniform"], nst,
            )
        )
    lines.extend([
        "", "## Per-stage novel_acc (mean over seeds)", "",
        "| mode | " + " | ".join("s%d" % i for i in range(len(perf["summary"][0]["stage_novel_mean"]))) + " |",
        "| --- | " + " | ".join("---" for _ in perf["summary"][0]["stage_novel_mean"]) + " |",
    ])
    for s in perf["summary"]:
        lines.append(
            "| %s | %s |" % (s["mode"], " | ".join("%.3f" % v for v in s["stage_novel_mean"]))
        )
    lines.append("")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--expand-min", type=int, default=10)
    ap.add_argument("--tta", action="store_true")
    ap.add_argument("--multiseed", action="store_true", help="streaming performance over 3 seeds")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    modes = {}
    for mode in ("repr_sep", "weighted", "uniform"):
        modes[mode] = run_cgcd(
            seed=args.seed, mode=mode, quick=bool(args.quick),
            expand_min=args.expand_min, do_tta=bool(args.tta) and mode == "repr_sep",
        )

    r, w, u = modes["repr_sep"], modes["weighted"], modes["uniform"]
    payload = {
        "seed": args.seed,
        "quick": bool(args.quick),
        "expand_min": args.expand_min,
        "tta": bool(args.tta),
        "modes": modes,
        "delta_novel_repr_vs_uniform": r["mean_novel_acc"] - u["mean_novel_acc"],
        "delta_novel_repr_vs_weighted": r["mean_novel_acc"] - w["mean_novel_acc"],
        "delta_forget_repr_vs_uniform": u["mean_forgetting"] - r["mean_forgetting"],
        "takeaway": (
            "Streaming KPI = novel_acc (+ known/forgetting). "
            "repr_sep = opportunity-gated attribute → sep on S → expand."
        ),
    }

    # always emit single-seed per-stage streaming performance
    stage_rows = []
    for row in r.get("stages") or []:
        stage_rows.append({
            "stage": row.get("stage"),
            "novel_acc": row.get("novel_acc"),
            "known_acc": row.get("known_acc"),
            "forgetting": row.get("forgetting"),
            "phase": row.get("phase"),
            "O_local": row.get("O_local"),
            "O_persist": row.get("O_persist"),
            "gap_lift": row.get("gap_lift"),
            "expand_fitted": row.get("expand_fitted"),
            "n_novel_train": row.get("n_novel_train"),
        })
    payload["streaming_performance"] = {
        "seed": args.seed,
        "summary": {
            "repr_sep": {
                "known": r["mean_known_acc"], "novel": r["mean_novel_acc"],
                "forgetting": r["mean_forgetting"],
            },
            "weighted": {
                "known": w["mean_known_acc"], "novel": w["mean_novel_acc"],
                "forgetting": w["mean_forgetting"],
            },
            "uniform": {
                "known": u["mean_known_acc"], "novel": u["mean_novel_acc"],
                "forgetting": u["mean_forgetting"],
            },
        },
        "repr_sep_stages": stage_rows,
    }
    (OUT / "repr_sep_cold_start.json").write_text(json.dumps(payload, indent=2))
    (OUT / "streaming_performance.json").write_text(
        json.dumps(payload["streaming_performance"], indent=2)
    )

    rows_diag = []
    for row in stage_rows:
        rows_diag.append(
            "| %s | %s | %s | %s | %s | %s | %s |" % (
                row.get("stage"), _fmt(row.get("novel_acc")), _fmt(row.get("known_acc")),
                _fmt(row.get("forgetting")), row.get("phase"),
                _fmt(row.get("O_local")), _fmt(row.get("gap_lift")),
            )
        )
    md = "\n".join([
        "# Streaming cold-start performance",
        "",
        payload["takeaway"],
        "",
        "## Summary (this seed)",
        "",
        "| mode | known | novel | forgetting | feat_hit |",
        "| --- | --- | --- | --- | --- |",
        "| repr_sep | %s | %s | %s | %s |" % (
            _fmt(r["mean_known_acc"]), _fmt(r["mean_novel_acc"]),
            _fmt(r["mean_forgetting"]), _fmt(r["mean_feat_hit"]),
        ),
        "| weighted | %s | %s | %s | %s |" % (
            _fmt(w["mean_known_acc"]), _fmt(w["mean_novel_acc"]),
            _fmt(w["mean_forgetting"]), _fmt(w["mean_feat_hit"]),
        ),
        "| uniform | %s | %s | %s | — |" % (
            _fmt(u["mean_known_acc"]), _fmt(u["mean_novel_acc"]),
            _fmt(u["mean_forgetting"]),
        ),
        "",
        "Δ novel (repr − uniform): **%s**" % _fmt(payload["delta_novel_repr_vs_uniform"]),
        "Δ novel (repr − weighted): **%s**" % _fmt(payload["delta_novel_repr_vs_weighted"]),
        "",
        "## Streaming per-stage (repr_sep)",
        "",
        "| stage | novel | known | forget | phase | O_local | gap_lift |",
        "| --- | --- | --- | --- | --- | --- | --- |",
        *rows_diag,
        "",
    ])

    if args.multiseed:
        seeds = [args.seed, args.seed + 1, args.seed + 2]
        # full stream when not --quick
        perf = _run_streaming_performance(seeds, quick=bool(args.quick), expand_min=args.expand_min)
        (OUT / "streaming_performance_multiseed.json").write_text(json.dumps(perf, indent=2))
        smd = _streaming_md(perf)
        (OUT / "STREAMING_PERFORMANCE.md").write_text(smd)
        (DOCS / "Streaming_cold_start_performance.md").write_text(smd)
        md = smd + "\n" + md
        payload["streaming_performance_multiseed"] = perf

    (OUT / "README.md").write_text(md)
    DOCS.mkdir(parents=True, exist_ok=True)
    (DOCS / "Repr_sep_cold_start.md").write_text(md)
    try:
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        shutil.copy2(OUT / "repr_sep_cold_start.json", ARTIFACTS / "repr_sep_cold_start.json")
        shutil.copy2(OUT / "streaming_performance.json", ARTIFACTS / "streaming_performance.json")
        shutil.copy2(OUT / "README.md", ARTIFACTS / "Repr_sep_cold_start.md")
        if (OUT / "STREAMING_PERFORMANCE.md").exists():
            shutil.copy2(OUT / "STREAMING_PERFORMANCE.md", ARTIFACTS / "STREAMING_PERFORMANCE.md")
            shutil.copy2(
                OUT / "streaming_performance_multiseed.json",
                ARTIFACTS / "streaming_performance_multiseed.json",
            )
    except Exception:
        pass
    print(md)
    print("wrote", OUT / "streaming_performance.json")


if __name__ == "__main__":
    main()
