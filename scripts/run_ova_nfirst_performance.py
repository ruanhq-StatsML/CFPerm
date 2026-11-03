#!/usr/bin/env python3
"""OVA AUC/acc on first n_first=100 of each new category (CE head).

  python3 scripts/run_ova_nfirst_performance.py --n-first 100
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

from first_category_performance import (  # noqa: E402
    LOSS_NAME,
    evaluate_board_ova_n_first,
    load_boards_n_first,
    run_ova_n_first_benchmark,
)
from weighted_cold_start import SEED  # noqa: E402

OUT = ROOT / "results" / "ova_nfirst"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")


def _fmt(x):
    try:
        if x is None or x != x:
            return "—"
        return "%.3f" % float(x)
    except Exception:
        return "—"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-first", type=int, default=100)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--multiseed", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    seeds = [args.seed, args.seed + 1, args.seed + 2] if args.multiseed else [args.seed]
    per_seed = []
    for s in seeds:
        print("=== seed", s, "n_first", args.n_first, "loss", LOSS_NAME, flush=True)
        per_seed.append(run_ova_n_first_benchmark(seed=s, n_first=args.n_first))

    # aggregate by dataset across seeds
    from collections import defaultdict
    agg = defaultdict(lambda: {
        "uniform_auc": [], "weighted_auc": [], "fused_auc": [],
        "uniform_acc": [], "weighted_acc": [], "fused_acc": [],
    })
    for run in per_seed:
        for b in run["boards"]:
            name = b["dataset"]
            sc = b["schemes"]
            agg[name]["uniform_auc"].append(sc["uniform"]["mean_ova_auc"])
            agg[name]["weighted_auc"].append(sc["weighted"]["mean_ova_auc"])
            agg[name]["fused_auc"].append(sc["both_fused"]["mean_ova_auc"])
            agg[name]["uniform_acc"].append(sc["uniform"]["mean_acc_on_first"])
            agg[name]["weighted_acc"].append(sc["weighted"]["mean_acc_on_first"])
            agg[name]["fused_acc"].append(sc["both_fused"]["mean_acc_on_first"])

    rows = []
    for name, d in agg.items():
        ua, wa, fa = map(np.asarray, (d["uniform_auc"], d["weighted_auc"], d["fused_auc"]))
        uc, wc, fc = map(np.asarray, (d["uniform_acc"], d["weighted_acc"], d["fused_acc"]))
        rows.append({
            "dataset": name,
            "ova_auc_uniform": float(np.nanmean(ua)),
            "ova_auc_weighted": float(np.nanmean(wa)),
            "ova_auc_fused": float(np.nanmean(fa)),
            "delta_auc_weighted": float(np.nanmean(wa - ua)),
            "delta_auc_fused": float(np.nanmean(fa - ua)),
            "acc_first_uniform": float(np.nanmean(uc)),
            "acc_first_weighted": float(np.nanmean(wc)),
            "acc_first_fused": float(np.nanmean(fc)),
            "delta_acc_weighted": float(np.nanmean(wc - uc)),
            "delta_acc_fused": float(np.nanmean(fc - uc)),
            "P_fused_auc_gt_uniform": float(np.nanmean(fa > ua)),
        })

    payload = {
        "loss": LOSS_NAME,
        "n_first": args.n_first,
        "seeds": seeds,
        "runs": per_seed,
        "by_dataset": rows,
        "mean_delta_auc_fused": float(np.nanmean([r["delta_auc_fused"] for r in rows])) if rows else None,
        "mean_delta_acc_fused": float(np.nanmean([r["delta_acc_fused"] for r in rows])) if rows else None,
        "takeaway": (
            "Method still uses CE (log_loss). Empowerment locus = new category "
            "first-n=%d window; KPI = OVA AUC + acc_on_first; weights vs uniform."
            % args.n_first
        ),
    }
    (OUT / "ova_nfirst.json").write_text(json.dumps(payload, indent=2))

    lines = [
        "# OVA first-n performance (CE head)",
        "",
        payload["takeaway"],
        "",
        f"n_first=**{args.n_first}**, seeds={seeds}, loss=`{LOSS_NAME}`.",
        "",
        "| dataset | OVA-AUC uniform | weighted | fused | ΔAUC fused | acc@first U | fused | Δacc |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for r in rows:
        lines.append(
            "| {dataset} | {ua} | {wa} | {fa} | {da} | {uc} | {fc} | {dc} |".format(
                dataset=r["dataset"],
                ua=_fmt(r["ova_auc_uniform"]),
                wa=_fmt(r["ova_auc_weighted"]),
                fa=_fmt(r["ova_auc_fused"]),
                da=_fmt(r["delta_auc_fused"]),
                uc=_fmt(r["acc_first_uniform"]),
                fc=_fmt(r["acc_first_fused"]),
                dc=_fmt(r["delta_acc_fused"]),
            )
        )
    lines.extend([
        "",
        f"mean Δ OVA-AUC (fused−uniform): **{_fmt(payload['mean_delta_auc_fused'])}**",
        f"mean Δ acc@first (fused−uniform): **{_fmt(payload['mean_delta_acc_fused'])}**",
        "",
        "## Per-category detail (last seed, fused)",
        "",
    ])
    last = per_seed[-1]
    for b in last["boards"]:
        lines.append(f"### {b['dataset']}")
        lines.append("")
        lines.append("| category | n_pos | OVA-AUC fused | OVA-AUC uniform | acc@first fused | acc@first U |")
        lines.append("| --- | --- | --- | --- | --- | --- |")
        fu = b["schemes"]["both_fused"]["per_category"]
        un = b["schemes"]["uniform"]["per_category"]
        for c in sorted(fu.keys(), key=int):
            lines.append(
                "| %s | %s | %s | %s | %s | %s |" % (
                    c, fu[c]["n_pos"],
                    _fmt(fu[c]["ova_auc"]), _fmt(un[c]["ova_auc"]),
                    _fmt(fu[c]["acc_on_first"]), _fmt(un[c]["acc_on_first"]),
                )
            )
        lines.append("")

    lines.extend([
        "```bash",
        "python3 scripts/run_ova_nfirst_performance.py --n-first 100 --multiseed",
        "```",
        "",
    ])
    md = "\n".join(lines)
    (OUT / "README.md").write_text(md)
    DOCS.mkdir(parents=True, exist_ok=True)
    (DOCS / "OVA_nfirst_performance.md").write_text(md)
    try:
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        shutil.copy2(OUT / "ova_nfirst.json", ARTIFACTS / "ova_nfirst.json")
        shutil.copy2(OUT / "README.md", ARTIFACTS / "OVA_nfirst_performance.md")
    except Exception:
        pass
    print(md)
    print("wrote", OUT / "ova_nfirst.json")


if __name__ == "__main__":
    main()
