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

from agod.po_refit import (
    make_batch_stream,
    run_adaptive_stream,
    run_dre_hop,
    run_oracle_switch,
    run_resid_stream,
    run_refit_stream,
    run_switch_stream,
    run_uniform_on_same_rows,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "po_refit_gated"
DOCS = ROOT / "docs" / "agod"

SCENES = (
    ("similar", dict(cov=0.0, concept_at=None, concept=0.0)),
    ("covariate", dict(cov=0.55, concept_at=None, concept=0.0)),
    ("concept", dict(cov=0.0, concept_at=4, concept=1.0)),
    ("mixed", dict(cov=0.40, concept_at=4, concept=1.0)),
)

MAIN_COLS = (
    "uniform_pair",
    "gated_pair",
    "soft_pair",
    "always_pair",
    "switch",
    "adaptive",
    "oracle",
    "uniform_hop",
    "gated_hop",
    "soft_hop",
    "dre_hop",
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
            "soft_pair": run_refit_stream(
                stream, assign="pair", gate=gate, always=False, soft=True
            ),
            "always_pair": run_refit_stream(stream, assign="pair", gate=gate, always=True),
            "uniform_hop": run_uniform_on_same_rows(stream, assign="hop"),
            "gated_hop": run_refit_stream(stream, assign="hop", gate=gate, always=False),
            "soft_hop": run_refit_stream(
                stream, assign="hop", gate=gate, always=False, soft=True
            ),
            "always_hop": run_refit_stream(stream, assign="hop", gate=gate, always=True),
            "adaptive": run_adaptive_stream(stream, gate=gate),
            "switch": run_switch_stream(stream, gate=gate),
            "resid": run_resid_stream(stream, gate=2.0, po_on_fire=False),
            "resid_po": run_resid_stream(stream, gate=2.0, po_on_fire=True),
            "dre_hop": run_dre_hop(stream),
            "oracle": run_oracle_switch(stream),
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
                    "history": [
                        {
                            "t": h["t"],
                            "next_mse": h["next_mse"],
                            "fired": h["fired"],
                            "ratio": h.get("ratio"),
                        }
                        for h in rec["history"]
                    ],
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
            cell = {
                "mse_mean": float(ms.mean()),
                "mse_std": float(ms.std()),
                "fire_mean": float(fr.mean()),
            }
            hops = {}
            for r in sub:
                if r["method"] != m:
                    continue
                for h in r.get("history") or []:
                    hops.setdefault(int(h["t"]), []).append(float(h["next_mse"]))
            if hops:
                cell["by_t"] = {
                    str(t): float(np.mean(vs)) for t, vs in sorted(hops.items())
                }
            out[scene]["methods"][m] = cell
    return out


def _mse(m, k):
    return m[k]["mse_mean"]


def write_md(summary, gate, path):
    lines = [
        "# Rolling PO-learner refit, gated re-adjustment",
        "",
        "Uniform stays the default. When a new batch arrives, **refit** the",
        "PO-learner: `T=0` = most recent control (`B_{t-2}`), `T=1` = 上一批 ∪ 这一批.",
        "√PO weights apply **only if** mean PO-risk on T=1 / T=0 ≥ "
        f"`gate={gate}`.",
        "",
        "## Board (next-batch Ridge MSE)",
        "",
        "| scene | uniform_pair | gated_pair | switch | adaptive | resid | resid_po | oracle | uniform_hop | gated_hop | dre_hop |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for scene, cell in summary.items():
        m = cell["methods"]

        def mse(k):
            return m[k]["mse_mean"]

        lines.append(
            "| `%s` | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |"
            % (
                scene,
                mse("uniform_pair"),
                mse("gated_pair"),
                mse("switch"),
                mse("adaptive"),
                mse("resid"),
                mse("resid_po"),
                mse("oracle"),
                mse("uniform_hop"),
                mse("gated_hop"),
                mse("dre_hop"),
            )
        )
    lines += [
        "",
        "Fire rates (PO-ratio γ=1.25 vs residual γ=2.0):",
        "",
        "| scene | gated_pair | gated_hop | switch | resid | resid_po | soft_pair |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for scene, cell in summary.items():
        m = cell["methods"]
        lines.append(
            "| `%s` | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f |"
            % (
                scene,
                m["gated_pair"]["fire_mean"],
                m["gated_hop"]["fire_mean"],
                m["switch"]["fire_mean"],
                m["resid"]["fire_mean"],
                m["resid_po"]["fire_mean"],
                m["soft_pair"]["fire_mean"],
            )
        )
    if "concept" in summary:
        m = summary["concept"]["methods"]
        hops = sorted(int(t) for t in m["uniform_pair"].get("by_t", {}))
        if hops:
            lines += [
                "",
                "## Concept stream, MSE by hop (predict `B_{t+1}`)",
                "",
                "Cut is `concept_at=4`. Predicting `B_4` uses hop `t=3` (both",
                "pre-cut) — that first post-shift batch is structurally",
                "unpredictable from labels. Gate can fire from `t=4`.",
                "",
                "| t (train) | test | uniform_pair | switch | resid | resid_po | oracle | gated_hop |",
                "|---|---|---:|---:|---:|---:|---:|---:|",
            ]
            for t in hops:
                test = f"B_{t + 1}"
                mark = " ← cut" if t + 1 == 4 else ""
                def at(method):
                    return m[method]["by_t"][str(t)]
                lines.append(
                    "| %d | `%s`%s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |"
                    % (
                        t,
                        test,
                        mark,
                        at("uniform_pair"),
                        at("switch"),
                        at("resid"),
                        at("resid_po"),
                        at("oracle"),
                        at("gated_hop"),
                    )
                )
    lines += [
        "",
        "- **similar / covariate**: `P(Y|X)` stable → uniform slightly better;",
        "  always-on √PO hurts. PO-ratio gate fire-rate should stay low.",
        "- **PO-ratio vs residual gate**: a *global* concept flip raises φ² on",
        "  both T=0 and T=1, so ρ=r̄1/r̄0 stays near 1 and switch/adaptive",
        "  keep pooling. Residual gate = MSE(fit 上一批 → 这一批) / train MSE;",
        "  that is the 'batches clearly different' detector.",
        "- **switch vs resid**: same train-set idea; differ only in the gate.",
        "- **resid vs resid_po**: drop old batch, then optional √PO on the new one.",
        "- **oracle**: knows `concept_at`; train-set upper bound, not deployable.",
        "- **dre_hop**: X-only density ratio. Misses label shift; overreacts to X-hop.",
        "- pair = last two batches; hop = new batch only.",
        "",
        f"{len(MAIN_COLS)} methods, Ridge next-batch MSE, gate γ={gate}.",
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
        r"only when the T=1 / T=0 PO-risk ratio exceeds $%.2f$. switch drops the old" % gate,
        r"batch without reweighting; oracle knows the concept cut.}",
        r"\label{tab:po-refit-gated}",
        r"\small",
        r"\setlength{\tabcolsep}{3.5pt}",
        r"\begin{tabular}{@{}lcccccc@{}}\toprule",
        r"Scene & unif-pair & gated-PO & resid & resid+PO & oracle & hop \\",
        r"\midrule",
    ]
    for scene, cell in summary.items():
        m = cell["methods"]
        lines.append(
            r"%s & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\"
            % (
                scene,
                m["uniform_pair"]["mse_mean"],
                m["gated_pair"]["mse_mean"],
                m["resid"]["mse_mean"],
                m["resid_po"]["mse_mean"],
                m["oracle"]["mse_mean"],
                m["uniform_hop"]["mse_mean"],
            )
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=8)
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
