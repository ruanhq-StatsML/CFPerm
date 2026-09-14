#!/usr/bin/env python3
"""Uniform vs gated online-RF √po_risk0 (next-batch MSE).

Uniform on the last two batches is the default. √PO / localization
only fire when consecutive OOS probe error jumps (distribution shift),
not every hop.

  PYTHONPATH=. python3 scripts/run_po_refit_gated.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from agod.online_rfperm import run_rfperm_stream
from agod.po_refit import (
    make_batch_stream,
    run_oracle_switch,
    run_resid_stream,
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

BOARD_METHODS = (
    "uniform_pair",
    "rfperm",
    "local",
    "resid",
    "oracle",
)


def run_scene(name, spec, seeds, n_batches, n_per, p, gate, learner):
    rows = []
    for seed in seeds:
        stream = make_batch_stream(
            n_batches=n_batches,
            n_per=n_per,
            p=p,
            seed=seed,
            **spec,
        )
        kw = dict(learner=learner, seed=seed)
        methods = {
            "uniform_pair": run_uniform_on_same_rows(stream, assign="pair", **kw),
            "rfperm": run_rfperm_stream(stream, gate=gate, localize=False, **kw),
            "local": run_rfperm_stream(
                stream, gate=gate, localize=True, q=0.30, **kw
            ),
            "resid": run_resid_stream(stream, gate=2.0, po_on_fire=False, **kw),
            "oracle": run_oracle_switch(stream, **kw),
        }
        for m, rec in methods.items():
            rows.append(
                {
                    "scene": name,
                    "learner": learner,
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
                            "n_train": h.get("n_train"),
                        }
                        for h in rec["history"]
                    ],
                }
            )
    return rows


def summarize(rows, learners):
    out = {}
    for learner in learners:
        out[learner] = {}
        for scene, spec in SCENES:
            sub = [r for r in rows if r["scene"] == scene and r["learner"] == learner]
            methods = sorted({r["method"] for r in sub})
            out[learner][scene] = {"spec": {k: spec[k] for k in spec}, "methods": {}}
            for m in methods:
                ms = np.array([r["online_mse"] for r in sub if r["method"] == m], float)
                fr = np.array([r["fire_rate"] for r in sub if r["method"] == m], float)
                cell = {
                    "mse_mean": float(ms.mean()) if ms.size else float("nan"),
                    "mse_std": float(ms.std()) if ms.size else float("nan"),
                    "fire_mean": float(fr.mean()) if fr.size else float("nan"),
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
                out[learner][scene]["methods"][m] = cell
    return out


def _fmt_board(learner, scenes):
    lines = [
        f"### `{learner}`",
        "",
        "| scene | uniform_pair | rfperm | local | resid | oracle |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for scene, cell in scenes.items():
        m = cell["methods"]

        def mse(k):
            return m[k]["mse_mean"]

        lines.append(
            "| `%s` | %.3f | %.3f | %.3f | %.3f | %.3f |"
            % (
                scene,
                mse("uniform_pair"),
                mse("rfperm"),
                mse("local"),
                mse("resid"),
                mse("oracle"),
            )
        )
    lines += [
        "",
        f"Fire rates (`{learner}`):",
        "",
        "| scene | rfperm | local | resid |",
        "|---|---:|---:|---:|",
    ]
    for scene, cell in scenes.items():
        m = cell["methods"]
        lines.append(
            "| `%s` | %.2f | %.2f | %.2f |"
            % (
                scene,
                m["rfperm"]["fire_mean"],
                m["local"]["fire_mean"],
                m["resid"]["fire_mean"],
            )
        )
    if "concept" in scenes:
        m = scenes["concept"]["methods"]
        hops = sorted(int(t) for t in m["uniform_pair"].get("by_t", {}))
        if hops:
            lines += [
                "",
                f"Concept hop path (`{learner}`, cut at `B_4`):",
                "",
                "| t (train) | test | uniform_pair | rfperm | local | resid | oracle |",
                "|---|---|---:|---:|---:|---:|---:|",
            ]
            for t in hops:
                test = f"B_{t + 1}"
                mark = " ← cut" if t + 1 == 4 else ""

                def at(method):
                    return m[method]["by_t"][str(t)]

                lines.append(
                    "| %d | `%s`%s | %.3f | %.3f | %.3f | %.3f | %.3f |"
                    % (
                        t,
                        test,
                        mark,
                        at("uniform_pair"),
                        at("rfperm"),
                        at("local"),
                        at("resid"),
                        at("oracle"),
                    )
                )
    lines.append("")
    return lines


def write_md(summary, gate, learners, path):
    lines = [
        "# Gated online-RF √po_risk0 vs PO-tail localization",
        "",
        "Probe is the shallow IPTW RF (`n_estimators=20`, `max_depth=4`)",
        "on 上一批 as T=0. Instance `po_risk0` is `|Y−μ0(X)|` mixed with the",
        "batch gap (`instance_po_risk`, mix=0.5). The next-batch model is",
        f"`rf` / `xgb` / `mlp` (not Ridge). Gate γ={gate}:",
        "",
        "`e_now = err(μ0 fit B_{t-1} → B_t)`, "
        "`e_prev = err(μ0 fit B_{t-2} → B_{t-1})`",
        "",
        "Fire only if `e_now / e_prev ≥ γ` (skip first hop). In-sample",
        "`e1>e0` would fire every hop on trees — do not use it.",
        "",
        "- **uniform_pair**: last two batches, w=1 (default).",
        "- **rfperm**: on fire, train current batch with `w=√po_risk0`.",
        "- **local**: same gate; train the high-`po_risk0` tail (`q=0.3`).",
        "- **resid**: consecutive residual hop-gate (full learner MSE).",
        "- **oracle**: knows the concept cut (train-set upper bound).",
        "",
        "Always-on DRE / always-on √PO are off this board.",
        "",
        "## Board (next-batch MSE)",
        "",
    ]
    for learner in learners:
        lines.extend(_fmt_board(learner, summary[learner]))
    lines += [
        "- **similar / covariate**: `P(Y|X)` stable → uniform; rfperm/local",
        "  should stay quiet.",
        "- **concept**: fire at the cut hop, then drop the old batch.",
        "- **local** vs **rfperm**: same gate; localization keeps the high-PO",
        "  tail instead of IPTW on the whole new batch.",
        "",
    ]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def write_tex(summary, gate, learners, path):
    lines = [
        r"% Gated online-RF sqrt(po_risk0) vs PO-tail localization.",
    ]
    for learner in learners:
        lines += [
            r"\begin{table}[ht]\centering",
            r"\caption{Next-batch MSE with \texttt{%s}. Uniform default." % learner,
            r"RFPerm fires on consecutive OOS probe jump "
            r"$(\gamma=%.2f)$.}" % gate,
            r"\label{tab:po-refit-%s}" % learner,
            r"\small",
            r"\setlength{\tabcolsep}{3.5pt}",
            r"\begin{tabular}{@{}lccccc@{}}\toprule",
            r"Scene & unif-pair & rfperm & local & resid & oracle \\",
            r"\midrule",
        ]
        for scene, cell in summary[learner].items():
            m = cell["methods"]
            lines.append(
                r"%s & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ \\"
                % (
                    scene,
                    m["uniform_pair"]["mse_mean"],
                    m["rfperm"]["mse_mean"],
                    m["local"]["mse_mean"],
                    m["resid"]["mse_mean"],
                    m["oracle"]["mse_mean"],
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
    ap.add_argument("--gate", type=float, default=1.5)
    ap.add_argument(
        "--learners",
        default="rf,xgb",
        help="comma list: rf, xgb, mlp",
    )
    args = ap.parse_args()
    seeds = list(range(2026, 2026 + int(args.seeds)))
    learners = [s.strip() for s in str(args.learners).split(",") if s.strip()]
    rows = []
    for learner in learners:
        for name, spec in SCENES:
            print("learner", learner, "scene", name, flush=True)
            rows.extend(
                run_scene(
                    name,
                    spec,
                    seeds,
                    args.n_batches,
                    args.n_per,
                    args.p,
                    args.gate,
                    learner,
                )
            )
    summary = summarize(rows, learners)
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    payload = {
        "gate": args.gate,
        "seeds": seeds,
        "learners": learners,
        "summary": summary,
        "rows": rows,
    }
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_md(summary, args.gate, learners, OUT / "PO_REFIT_GATED_REPORT.md")
    write_md(summary, args.gate, learners, DOCS / "AGOD_po_refit_gated.md")
    write_tex(summary, args.gate, learners, OUT / "PO_refit_gated.tex")
    write_tex(summary, args.gate, learners, DOCS / "AGOD_po_refit_gated_tables_only.tex")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
