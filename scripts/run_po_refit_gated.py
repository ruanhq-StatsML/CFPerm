#!/usr/bin/env python3
"""Uniform vs gated rolling PO-learner refit (next-batch MSE).

Claim: uniform is slightly better when batches look alike; re-adjust
(√PO weights from a *refit* PO-learner) only when the new batch is
clearly different.

  T=0  most recent control (batch t-2)
  T=1  上一批 ∪ 这一批
  train on T=1 → score batch t+1

  PYTHONPATH=. python3 scripts/run_po_refit_gated.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from agod.po_refit import make_batch_stream, run_refit_stream, run_uniform_on_same_rows

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "po_refit_gated"
DOCS = ROOT / "docs" / "agod"

SCENES = (
    ("similar", dict(cov=0.0, concept_at=None, concept=0.0)),
    ("covariate", dict(cov=0.55, concept_at=None, concept=0.0)),
    ("concept", dict(cov=0.0, concept_at=4, concept=1.0)),
    ("mixed", dict(cov=0.40, concept_at=4, concept=1.0)),
)


def run_scene(name, spec, seeds, n_batches, n_per, p, gate):
    rows = []
    for seed in seeds:
        stream = make_batch_stream(
            n_batches=n_batches,
            n_per=n_per,
            p=p,
            seed=seed,
            **spec,
        )
        methods = {
            "uniform_pair": run_uniform_on_same_rows(stream, assign="pair"),
            "gated_pair": run_refit_stream(stream, assign="pair", gate=gate, always=False),
            "always_pair": run_refit_stream(stream, assign="pair", gate=gate, always=True),
            "uniform_hop": run_uniform_on_same_rows(stream, assign="hop"),
            "gated_hop": run_refit_stream(stream, assign="hop", gate=gate, always=False),
            "always_hop": run_refit_stream(stream, assign="hop", gate=gate, always=True),
        }
        for m, rec in methods.items():
            rows.append(
                {
                    "scene": name,
                    "seed": int(seed),
                    "method": m,
                    "online_mse": rec["online_mse"],
                    "fire_rate": rec["fire_rate"],
                    "path": rec["path"],
                }
            )
    return rows


def summarize(rows):
    out = {}
    for scene, spec in SCENES:
        sub = [r for r in rows if r["scene"] == scene]
        methods = sorted({r["method"] for r in sub})
        out[scene] = {"spec": {k: spec[k] for k in spec}, "methods": {}}
        for m in methods:
            ms = np.array([r["online_mse"] for r in sub if r["method"] == m], float)
            fr = np.array([r["fire_rate"] for r in sub if r["method"] == m], float)
            out[scene]["methods"][m] = {
                "mse_mean": float(ms.mean()),
                "mse_std": float(ms.std()),
                "fire_mean": float(fr.mean()),
            }
    return out


def write_md(summary, gate, path):
    lines = [
        "# Rolling PO-learner refit, gated re-adjustment",
        "",
        "Uniform stays the default. When a new batch arrives, **refit** the",
        "PO-learner: `T=0` = most recent control, `T=1` = 上一批 ∪ 这一批.",
        "√PO weights apply **only if** mean PO-risk on T=1 / T=0 ≥ "
        f"`gate={gate}`.",
        "",
        "| scene | uniform_pair | gated_pair | always_pair | fire_gated | uniform_hop | gated_hop |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for scene, cell in summary.items():
        m = cell["methods"]
        def mse(k):
            return m[k]["mse_mean"]
        lines.append(
            "| `%s` | %.3f | %.3f | %.3f | %.2f | %.3f | %.3f |"
            % (
                scene,
                mse("uniform_pair"),
                mse("gated_pair"),
                mse("always_pair"),
                m["gated_pair"]["fire_mean"],
                mse("uniform_hop"),
                mse("gated_hop"),
            )
        )
    lines += [
        "",
        "- **similar**: batches share P(Y|X) → uniform should be slightly better;",
        "  gated fire-rate should stay low.",
        "- **concept / mixed**: a clear hop → gated re-adjusts; always-on PO",
        "  can overfit the current contrast.",
        "- pair = train on last two batches; hop = train on the new batch only.",
        "",
    ]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def write_tex(summary, gate, path):
    lines = [
        r"% Rolling PO-learner refit. Auto-generated.",
        r"\begin{table}[ht]\centering",
        r"\caption{Next-batch Ridge MSE. Uniform is the default. The PO-learner is",
        r"refit at every hop with $T{=}0$ the most recent control batch and",
        r"$T{=}1$ the previous$\cup$current batch. $\sqrt{\mathrm{PO}}$ weights fire",
        r"only when the T=1 / T=0 PO-risk ratio exceeds $%.2f$.}" % gate,
        r"\label{tab:po-refit-gated}",
        r"\small",
        r"\begin{tabular}{@{}lcccc@{}}\toprule",
        r"Scene & uniform (pair) & gated PO & always PO & gate fire \\",
        r"\midrule",
    ]
    for scene, cell in summary.items():
        m = cell["methods"]
        lines.append(
            r"%s & $%.3f$ & $%.3f$ & $%.3f$ & $%.2f$ \\"
            % (
                scene,
                m["uniform_pair"]["mse_mean"],
                m["gated_pair"]["mse_mean"],
                m["always_pair"]["mse_mean"],
                m["gated_pair"]["fire_mean"],
            )
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--n-batches", type=int, default=8)
    ap.add_argument("--n-per", type=int, default=100)
    ap.add_argument("--p", type=int, default=12)
    ap.add_argument("--gate", type=float, default=1.25)
    args = ap.parse_args()
    seeds = list(range(2026, 2026 + int(args.seeds)))
    rows = []
    for name, spec in SCENES:
        print("scene", name, flush=True)
        rows.extend(
            run_scene(
                name,
                spec,
                seeds,
                args.n_batches,
                args.n_per,
                args.p,
                args.gate,
            )
        )
    summary = summarize(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    payload = {"gate": args.gate, "seeds": seeds, "summary": summary, "rows": rows}
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_md(summary, args.gate, OUT / "PO_REFIT_GATED_REPORT.md")
    write_md(summary, args.gate, DOCS / "AGOD_po_refit_gated.md")
    write_tex(summary, args.gate, OUT / "PO_refit_gated.tex")
    write_tex(summary, args.gate, DOCS / "AGOD_po_refit_gated_tables_only.tex")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
