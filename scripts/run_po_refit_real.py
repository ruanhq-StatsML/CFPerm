#!/usr/bin/env python3
"""Gated PO-refit on real consecutive-batch streams (no K-fold).

Train on recent batches, score the next batch. Residual gate compares
two consecutive hops:

  ρ = err(fit B_{t-1} → B_t) / err(fit B_{t-2} → B_{t-1})

  PYTHONPATH=. python3 scripts/run_po_refit_real.py --learners rf,xgb
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from agod.real_data import iter_real_streams
from agod.po_refit import (
    run_dre_hop,
    run_resid_stream,
    run_refit_stream,
    run_switch_stream,
    run_uniform_on_same_rows,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "po_refit_real"
DOCS = ROOT / "docs" / "agod"

METHODS = (
    "uniform_pair",
    "gated_pair",
    "switch",
    "resid",
    "resid_po",
    "uniform_hop",
    "dre_hop",
)


def run_one(stream, learner, gate, seed):
    kw = dict(learner=learner, seed=seed)
    return {
        "uniform_pair": run_uniform_on_same_rows(stream, assign="pair", **kw),
        "gated_pair": run_refit_stream(
            stream, assign="pair", gate=gate, always=False, **kw
        ),
        "switch": run_switch_stream(stream, gate=gate, **kw),
        "resid": run_resid_stream(stream, gate=2.0, po_on_fire=False, **kw),
        "resid_po": run_resid_stream(stream, gate=2.0, po_on_fire=True, **kw),
        "uniform_hop": run_uniform_on_same_rows(stream, assign="hop", **kw),
        "dre_hop": run_dre_hop(stream, **kw),
    }


def write_md(rows, path):
    lines = [
        "# Real-data consecutive-batch PO-refit (no K-fold)",
        "",
        "Each hop trains on the latest consecutive batches and scores the",
        "**next** batch. Residual gate:",
        "",
        "`ρ = err(B_{t-1} → B_t) / err(B_{t-2} → B_{t-1})`",
        "",
        "MSE ↓ for continuous, Acc ↑ for discrete. Uniform is the default.",
        "",
    ]
    learners = sorted({r["learner"] for r in rows})
    datasets = []
    for r in rows:
        if r["dataset"] not in datasets:
            datasets.append(r["dataset"])
    for learner in learners:
        lines += [f"## `{learner}`", ""]
        lines += [
            "| dataset | task | n_batches | uniform_pair | gated_pair | switch | resid | resid_po | uniform_hop | dre_hop | fire_resid |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for ds in datasets:
            sub = [r for r in rows if r["learner"] == learner and r["dataset"] == ds]
            if not sub:
                continue
            by_m = {r["method"]: r for r in sub}
            task = sub[0]["task"]
            nb = sub[0]["n_batches"]

            def sc(m):
                return by_m[m]["online_score"]

            lines.append(
                "| `%s` | %s | %d | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.2f |"
                % (
                    ds,
                    task,
                    nb,
                    sc("uniform_pair"),
                    sc("gated_pair"),
                    sc("switch"),
                    sc("resid"),
                    sc("resid_po"),
                    sc("uniform_hop"),
                    sc("dre_hop"),
                    by_m["resid"]["fire_rate"],
                )
            )
        lines.append("")
    lines += [
        "- **resid** drops the old batch only when the consecutive hop error jumps.",
        "- **resid_po** = resid train-set + √PO on the new batch.",
        "- **gated_pair / switch** use PO-ratio, which misses a global label map flip.",
        "- No K-fold; no shuffle. Row order is the stream clock.",
        "",
    ]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--learners", default="rf,xgb")
    ap.add_argument("--n-batches", type=int, default=10)
    ap.add_argument("--n-per", type=int, default=100)
    ap.add_argument("--gate", type=float, default=1.25)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--pca-d", type=int, default=32)
    args = ap.parse_args()
    learners = [s.strip() for s in args.learners.split(",") if s.strip()]
    rows = []
    for stream in iter_real_streams(
        n_per=args.n_per,
        n_batches=args.n_batches,
        pca_d=args.pca_d,
        seed=args.seed,
    ):
        task = stream.task
        nb = int(stream.batch.max()) + 1
        print(
            f"=== {stream.name} task={task} batches={nb} "
            f"n={len(stream.y)} d={stream.X.shape[1]} ===",
            flush=True,
        )
        for learner in learners:
            print(f"  learner={learner}", flush=True)
            recs = run_one(stream, learner, args.gate, args.seed)
            for method, rec in recs.items():
                rows.append(
                    {
                        "dataset": stream.name,
                        "task": task,
                        "learner": learner,
                        "method": method,
                        "n_batches": nb,
                        "online_score": rec["online_mse"],
                        "fire_rate": rec["fire_rate"],
                        "path": rec["path"],
                        "history": [
                            {
                                "t": h["t"],
                                "score": h["next_mse"],
                                "fired": h["fired"],
                                "ratio": h.get("ratio"),
                            }
                            for h in rec["history"]
                        ],
                    }
                )
                print(
                    f"    {method}: {rec['metric']}={rec['online_mse']:.4f} "
                    f"fire={rec['fire_rate']:.2f}",
                    flush=True,
                )
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    payload = {
        "gate": args.gate,
        "learners": learners,
        "n_per": args.n_per,
        "n_batches": args.n_batches,
        "rows": rows,
    }
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_md(rows, OUT / "PO_REFIT_REAL_REPORT.md")
    write_md(rows, DOCS / "AGOD_po_refit_real.md")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
