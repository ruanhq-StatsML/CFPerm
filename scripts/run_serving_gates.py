#!/usr/bin/env python3
"""Run the two serving gates on tables that are already on disk.

Same object, same gates, different Y. Live facets today:

  审核 / judge     llm_audit/xy_hh_helpful_{consistent,hop}.csv
                   llm_audit/xy_wildguard.csv
                   llm_audit/xy_hh_multistep_{consistent,hop}.csv
  Graph-RAG        graph_rag_batches/xy_graph_query{,_hop}.csv
  混合检索         hybrid_retrieval/xy_hotpot_hybrid{,_hop}.csv

Graph-RAG gates run on the query rows (Y in {0,1}); the one-row-per-window
files are the serving-table formulation, not the probe sample.

Hotpot file order is not wall-clock time.

Usage::

    PYTHONPATH=. python3 scripts/run_serving_gates.py
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.llm_audit_online_bootstrap_prototype import (  # noqa: E402
    CUT_BATCH,
    GATE,
    N_REF_BATCHES,
    run_stream,
)

AUDIT = ROOT / "results" / "manuscript" / "llm_audit"
GRAPH = ROOT / "results" / "manuscript" / "graph_rag_batches"
HYBRID = ROOT / "results" / "manuscript" / "hybrid_retrieval"
OUT = ROOT / "results" / "manuscript" / "serving_gates"
SEED = 2026
LEAK = ("chosen", "rejected")

STREAMS = [
    {
        "facet": "审核 / judge",
        "y_meaning": "过 / 不过",
        "fire_means": "政策包或 judge 换代",
        "not_means": "拒绝率；HH chosen",
        "name": "judge_helpful_quiet",
        "title": "HH helpful, consistent",
        "path": AUDIT / "xy_hh_helpful_consistent.csv",
        "regime": "quiet",
        "labeled_cut": False,
    },
    {
        "facet": "审核 / judge",
        "y_meaning": "过 / 不过",
        "fire_means": "政策包或 judge 换代",
        "not_means": "拒绝率；HH chosen",
        "name": "judge_helpful_hop",
        "title": "HH helpful, policy hop",
        "path": AUDIT / "xy_hh_helpful_hop.csv",
        "regime": "hop",
        "labeled_cut": True,
    },
    {
        "facet": "审核 / judge",
        "y_meaning": "过 / 不过",
        "fire_means": "政策包或 judge 换代",
        "not_means": "拒绝率；HH chosen",
        "name": "judge_wildguard_native",
        "title": "WildGuard native is_unharmful",
        "path": AUDIT / "xy_wildguard.csv",
        "regime": "quiet",
        "labeled_cut": False,
    },
    {
        "facet": "审核 / judge",
        "y_meaning": "过 / 不过（这一跳）",
        "fire_means": "政策包或 judge 换代",
        "not_means": "拒绝率；HH chosen；整段对话成功",
        "name": "judge_multistep_quiet",
        "title": "HH multi-step, consistent",
        "path": AUDIT / "xy_hh_multistep_consistent.csv",
        "regime": "quiet",
        "labeled_cut": False,
    },
    {
        "facet": "审核 / judge",
        "y_meaning": "过 / 不过（这一跳）",
        "fire_means": "政策包或 judge 换代",
        "not_means": "拒绝率；HH chosen；整段对话成功",
        "name": "judge_multistep_hop",
        "title": "HH multi-step, policy hop",
        "path": AUDIT / "xy_hh_multistep_hop.csv",
        "regime": "hop",
        "labeled_cut": True,
    },
    {
        "facet": "Graph-RAG 子图",
        "y_meaning": "支撑节点在图包里",
        "fire_means": "图包或 community 换代",
        "not_means": "单点相关性",
        "name": "graphrag_local_quiet",
        "title": "local pack (seed ∪ 1-hop)",
        "path": GRAPH / "xy_graph_query.csv",
        "regime": "quiet",
        "labeled_cut": False,
    },
    {
        "facet": "Graph-RAG 子图",
        "y_meaning": "支撑节点在图包里",
        "fire_means": "图包或 community 换代",
        "not_means": "单点相关性",
        "name": "graphrag_community_hop",
        "title": "rewire + largest-CC pack",
        "path": GRAPH / "xy_graph_query_hop.csv",
        "regime": "hop",
        "labeled_cut": True,
    },
    {
        "facet": "混合检索",
        "y_meaning": "融合后的答案成不成立",
        "fire_means": "某一路索引或融合坏了",
        "not_means": "单路 Recall",
        "name": "hybrid_fused_quiet",
        "title": "RRF fused gold-in-topk",
        "path": HYBRID / "xy_hotpot_hybrid.csv",
        "regime": "quiet",
        "labeled_cut": False,
    },
    {
        "facet": "混合检索",
        "y_meaning": "融合后的答案成不成立",
        "fire_means": "某一路索引或融合坏了",
        "not_means": "单路 Recall",
        "name": "hybrid_dense_hop",
        "title": "dense channel flipped after cut",
        "path": HYBRID / "xy_hotpot_hybrid_hop.csv",
        "regime": "hop",
        "labeled_cut": True,
    },
]


def x_columns(fieldnames) -> list[str]:
    return [c for c in fieldnames if c.startswith("x_")]


def load_xy_table(path: Path):
    """Any (y, batch, x_*) table. chosen/rejected is not Y."""
    with path.open() as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"empty table {path}")
    keys = {k.lower() for k in rows[0]}
    for bad in LEAK:
        if bad in keys:
            raise RuntimeError(f"{path} still has {bad} — that is not Y")
    cols = x_columns(rows[0].keys())
    if not cols:
        raise RuntimeError(f"{path} has no x_* columns")
    y = np.asarray([int(float(r["y"])) for r in rows], dtype=int)
    batch = np.asarray([int(r["batch"]) for r in rows], dtype=int)
    X = np.asarray([[float(r[c]) for c in cols] for r in rows], dtype=float)
    return X, y, batch, cols


def compact(rec: dict, spec: dict) -> dict:
    boot = rec["bootstrap"]
    return {
        "facet": spec["facet"],
        "name": rec["name"],
        "title": rec["title"],
        "regime": rec["regime"],
        "y_meaning": spec["y_meaning"],
        "fire_means": spec["fire_means"],
        "not_means": spec["not_means"],
        "n": rec["n"],
        "n_batches": rec["n_batches"],
        "y_pass_rate": rec["y_pass_rate"],
        "mu_ref": rec["mu_ref"],
        "mean_delta": boot["mean_delta"],
        "n_boot_fires": boot["n_fires"],
        "first_boot_fire": boot["first_sig_abs_batch"],
        "delay_batches": boot["delay_batches"],
        "hop_fire_at_cut": rec["hop_fire_at_cut"],
        "n_hop_fires": rec["n_hop_fires"],
        "first_hop_abs_batch": rec["first_hop_abs_batch"],
        "path": str(spec["path"].relative_to(ROOT)),
    }


def plot_deltas(runs: list[tuple[dict, dict]], path: Path) -> None:
    pairs = [
        ("judge", "judge_helpful_quiet", "judge_helpful_hop"),
        ("judge-multistep", "judge_multistep_quiet", "judge_multistep_hop"),
        ("Graph-RAG", "graphrag_local_quiet", "graphrag_community_hop"),
        ("hybrid", "hybrid_fused_quiet", "hybrid_dense_hop"),
    ]
    by_name = {rec["name"]: rec for rec, _spec in runs}
    fig, axes = plt.subplots(len(pairs), 2, figsize=(9.4, 9.4), sharex=True)
    for i, (facet, quiet, hop) in enumerate(pairs):
        for j, name in enumerate((quiet, hop)):
            ax = axes[i, j]
            rec = by_name.get(name)
            if rec is None:
                ax.set_axis_off()
                continue
            xs = [r["abs_batch"] for r in rec["rows"]]
            ys = [r["s"] for r in rec["rows"]]
            ax.axhline(0.0, color="#888", lw=0.8)
            if rec["regime"] == "hop":
                ax.axvline(CUT_BATCH - 0.5, color="#9b2c2c", ls="--", lw=1.0)
            ax.plot(xs, ys, marker="o", color="#1f4e79")
            ax.set_title(f"{facet} · {rec['regime']}", fontsize=10)
            if i == len(pairs) - 1:
                ax.set_xlabel("batch")
            if j == 0:
                ax.set_ylabel(r"$\Delta_t$")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def _fmt(x) -> str:
    if x is None:
        return "—"
    if isinstance(x, bool):
        return "yes" if x else "no"
    if isinstance(x, float):
        if not np.isfinite(x):
            return "—"
        return f"{x:.3f}"
    return str(x)


def render_report(rows: list[dict]) -> str:
    lines = [
        "# Serving gates — same map, tables already on disk",
        "",
        "One command. Two gates. Different `Y`. Fire is a hop of `P(Y|X)`, not a quality score.",
        "",
        "```bash",
        "PYTHONPATH=. python3 scripts/run_serving_gates.py",
        "```",
        "",
        "| Facet | Stream | Y | regime | pass | mean Δ | boot fires | hop@cut | n hop |",
        "|---|---|---|---|---:|---:|---:|---|---:|",
    ]
    for r in rows:
        lines.append(
            "| {facet} | {title} | {y} | {regime} | {pass_} | {delta} | {nf} | {hop} | {nh} |".format(
                facet=r["facet"],
                title=r["title"],
                y=r["y_meaning"],
                regime=r["regime"],
                pass_=_fmt(r["y_pass_rate"]),
                delta=_fmt(r["mean_delta"]),
                nf=r["n_boot_fires"],
                hop=_fmt(r["hop_fire_at_cut"]),
                nh=r["n_hop_fires"],
            )
        )
    lines += [
        "",
        "Quiet streams should stay near Δ = 0. Labeled hops should lift Δ and, when the adjacent-window ratio clears γ, fire last-two at the cut.",
        "Multi-step audit is the same map: one row per assistant hop, Y = pass/fail for that hop, not the whole thread.",
        "Graph-RAG and hybrid Hotpot tables are pool snapshots: `batch` is file/window index, not wall-clock time.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    runs = []
    skipped = []
    for spec in STREAMS:
        if not spec["path"].exists():
            skipped.append(str(spec["path"].relative_to(ROOT)))
            continue
        X, y, batch, _cols = load_xy_table(spec["path"])
        rec = run_stream(
            spec["name"],
            spec["title"],
            X,
            y,
            batch,
            regime=spec["regime"],
            n_ref_batches=N_REF_BATCHES,
            cut_batch=CUT_BATCH,
            gate=GATE,
            seed=SEED,
            labeled_cut=spec["labeled_cut"],
        )
        runs.append((rec, spec))
    if not runs:
        raise RuntimeError("no serving tables on disk")
    summary = [compact(rec, spec) for rec, spec in runs]
    plot_deltas(runs, OUT / "delta_by_facet.png")
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    report = render_report(summary)
    if skipped:
        report += "Skipped missing files:\n\n" + "\n".join(f"- `{p}`" for p in skipped) + "\n"
    (OUT / "REPORT.md").write_text(report, encoding="utf-8")
    print(report)
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
