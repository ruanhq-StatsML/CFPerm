#!/usr/bin/env python3
"""Gated online-RF √po_risk0 reweighting on consecutive-batch real streams.

Uniform on the last two batches is the default. √PO IPTW only
fires when consecutive OOS probe error jumps:

  e_now  = err(μ0 fit B_{t-1} → B_t)
  e_prev = err(μ0 fit B_{t-2} → B_{t-1})
  fire iff e_now / e_prev ≥ γ   (default 1.5)

  PYTHONPATH=. python3 scripts/run_po_refit_real.py --learners rf --n-per 200
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from agod.online_rfperm import run_online_rfperm
from agod.performance_tex import write_all_tex, write_real_tex
from agod.po_refit import run_resid_stream, run_uniform_last_two
from agod.real_data import iter_real_streams

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "po_refit_real"
DOCS = ROOT / "docs" / "agod"

METHODS = (
    "uniform_pair",
    "rfperm",
    "resid",
)


def run_one(stream, learner, gate, seed):
    kw = dict(learner=learner, seed=seed)
    return {
        "uniform_pair": run_uniform_last_two(stream, **kw),
        "rfperm": run_online_rfperm(stream, gate=gate, **kw),
        "resid": run_resid_stream(stream, gate=2.0, po_on_fire=False, **kw),
    }


def write_md(rows, path, *, n_per, gate):
    lines = [
        "# Real consecutive-batch PO-risk adaptation",
        "",
        "Overview: `docs/agod/AGOD_overview.md`.",
        "LaTeX: `docs/agod/AGOD_po_refit_real_tables.tex`,",
        "`docs/agod/AGOD_performance_tables.tex`.",
        "",
        f"Batch size **{n_per}**, row order is the stream clock (no shuffle,",
        f"no K-fold). Gate γ={gate}: fire only when consecutive OOS probe",
        "error jumps. Always-on DRE / always-on √PO are off this board.",
        "",
        "- **uniform_pair**: last two batches, w=1.",
        "- **rfperm**: same rows; on fire, T=1 gets `w=√po_risk0`.",
        "- **resid**: residual hop-gate (drop old batch) as a reference.",
        "",
        "Continuous tasks report **RMSE** (↓); discrete report **Acc** (↑).",
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
            "| dataset | clock | task | n_batches | uniform_pair | rfperm | resid | fire_rfperm | fire_resid |",
            "|---|---|---|---:|---:|---:|---:|---:|---:|",
        ]
        for ds in datasets:
            sub = [r for r in rows if r["learner"] == learner and r["dataset"] == ds]
            if not sub:
                continue
            by_m = {r["method"]: r for r in sub}
            task = sub[0]["task"]
            nb = sub[0]["n_batches"]
            clock = sub[0].get("clock") or ""

            def sc(m):
                v = by_m[m]["online_score"]
                if task == "mse":
                    return float(v) ** 0.5
                return float(v)

            lines.append(
                "| `%s` | %s | %s | %d | %.4f | %.4f | %.4f | %.2f | %.2f |"
                % (
                    ds,
                    clock,
                    "RMSE" if task == "mse" else "acc",
                    nb,
                    sc("uniform_pair"),
                    sc("rfperm"),
                    sc("resid"),
                    by_m["rfperm"]["fire_rate"],
                    by_m["resid"]["fire_rate"],
                )
            )
        lines.append("")
    lines += [
        "Quiet streams should match uniform. A real P(Y|X) hop should fire",
        "rfperm; RMSE/Acc then shows whether reweighting helped the",
        "**next** batch. Subset localization is off this board — too thin",
        "at batch size 200.",
        "",
    ]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--learners", default="rf,xgb")
    ap.add_argument("--n-batches", type=int, default=24)
    ap.add_argument("--n-per", type=int, default=200)
    ap.add_argument("--gate", type=float, default=1.5)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--pca-d", type=int, default=32)
    ap.add_argument(
        "--clocks",
        default="time,shift,spatial",
        help="comma list: time, shift, spatial, local",
    )
    ap.add_argument(
        "--download",
        action="store_true",
        default=True,
        help="fetch remote OpenML/UCI time tables if no cache",
    )
    ap.add_argument("--no-download", action="store_true")
    args = ap.parse_args()
    learners = [s.strip() for s in args.learners.split(",") if s.strip()]
    clocks = [s.strip() for s in args.clocks.split(",") if s.strip()]
    download = bool(args.download) and not bool(args.no_download)
    rows = []
    for stream in iter_real_streams(
        n_per=args.n_per,
        n_batches=args.n_batches,
        pca_d=args.pca_d,
        seed=args.seed,
        max_n=max(8000, args.n_per * args.n_batches),
        download=download,
        clocks=clocks,
    ):
        task = stream.task
        nb = int(stream.batch.max()) + 1
        clock = (stream.meta or {}).get("clock", "")
        print(
            f"=== {stream.name} clock={clock} task={task} batches={nb} "
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
                        "clock": clock,
                        "task": task,
                        "learner": learner,
                        "method": method,
                        "n_batches": nb,
                        "n_per": int(args.n_per),
                        "online_score": rec["online_mse"],
                        "fire_rate": rec["fire_rate"],
                        "path": rec["path"],
                        "history": [
                            {
                                "t": h["t"],
                                "score": h["next_mse"],
                                "fired": h["fired"],
                                "ratio": h.get("ratio"),
                                "n_train": h.get("n_train"),
                                "po_t1": h.get("po_t1"),
                                "w_t1": h.get("w_t1"),
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
        "clocks": clocks,
        "rows": rows,
    }
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_md(
        rows,
        OUT / "PO_REFIT_REAL_REPORT.md",
        n_per=args.n_per,
        gate=args.gate,
    )
    write_md(
        rows,
        DOCS / "AGOD_po_refit_real.md",
        n_per=args.n_per,
        gate=args.gate,
    )
    write_real_tex(
        rows,
        args.gate,
        args.n_per,
        OUT / "PO_refit_real.tex",
    )
    write_real_tex(
        rows,
        args.gate,
        args.n_per,
        DOCS / "AGOD_po_refit_real_tables.tex",
    )
    gated_json = ROOT / "results" / "po_refit_gated" / "summary.json"
    if gated_json.is_file():
        gated = json.loads(gated_json.read_text(encoding="utf-8"))
        write_all_tex(gated, payload, DOCS / "AGOD_performance_tables.tex")
        write_all_tex(gated, payload, OUT / "AGOD_performance_tables.tex")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
