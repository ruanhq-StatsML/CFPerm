#!/usr/bin/env python3
"""Text multiclass OVA first-n (no tabular).

  python3 scripts/run_text_ova_nfirst.py --n-first 100 --multiseed
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from text_ova_nfirst import LOSS_NAME, run_text_ova_nfirst  # noqa: E402
from weighted_cold_start import SEED  # noqa: E402

OUT = ROOT / "results" / "text_ova_nfirst"
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
    runs = []
    for s in seeds:
        print("=== text OVA seed", s, "n_first", args.n_first, flush=True)
        runs.append(run_text_ova_nfirst(seed=s, n_first=args.n_first))

    agg = defaultdict(lambda: {
        "u_auc": [], "w_auc": [], "f_auc": [],
        "u_acc": [], "w_acc": [], "f_acc": [],
    })
    for run in runs:
        for b in run["boards"]:
            name = b["dataset"]
            sc = b["schemes"]
            agg[name]["u_auc"].append(sc["uniform"]["mean_ova_auc"])
            agg[name]["w_auc"].append(sc["weighted"]["mean_ova_auc"])
            agg[name]["f_auc"].append(sc["both_fused"]["mean_ova_auc"])
            agg[name]["u_acc"].append(sc["uniform"]["mean_acc_on_first"])
            agg[name]["w_acc"].append(sc["weighted"]["mean_acc_on_first"])
            agg[name]["f_acc"].append(sc["both_fused"]["mean_acc_on_first"])

    rows = []
    for name, d in agg.items():
        ua, wa, fa = map(np.asarray, (d["u_auc"], d["w_auc"], d["f_auc"]))
        uc, wc, fc = map(np.asarray, (d["u_acc"], d["w_acc"], d["f_acc"]))
        rows.append({
            "dataset": name,
            "ova_auc_uniform": float(np.nanmean(ua)),
            "ova_auc_weighted": float(np.nanmean(wa)),
            "ova_auc_fused": float(np.nanmean(fa)),
            "delta_auc_fused": float(np.nanmean(fa - ua)),
            "delta_auc_weighted": float(np.nanmean(wa - ua)),
            "acc_first_uniform": float(np.nanmean(uc)),
            "acc_first_fused": float(np.nanmean(fc)),
            "delta_acc_fused": float(np.nanmean(fc - uc)),
            "P_fused_auc_gt_U": float(np.nanmean(fa > ua)),
        })

    payload = {
        "domain": "text_multiclass_only",
        "tabular_excluded": True,
        "loss": LOSS_NAME,
        "weight_channel": "sample-level (FSDS localizes; w_i in partial_fit)",
        "n_first": args.n_first,
        "seeds": seeds,
        "runs": runs,
        "by_dataset": rows,
        "mean_delta_auc_fused": float(np.nanmean([r["delta_auc_fused"] for r in rows])) if rows else None,
        "mean_delta_acc_fused": float(np.nanmean([r["delta_acc_fused"] for r in rows])) if rows else None,
        "takeaway": (
            "Text multiclass only — tabular skipped. CE head + sample-level weights; "
            "KPI = OVA AUC / acc on first n_first of each new category."
        ),
    }
    (OUT / "text_ova_nfirst.json").write_text(json.dumps(payload, indent=2))

    lines = [
        "# Text multiclass OVA first-n (no tabular)",
        "",
        payload["takeaway"],
        "",
        f"n_first=**{args.n_first}**, seeds={seeds}, loss=`{LOSS_NAME}`, weights=**sample-level**.",
        "",
        "| dataset | OVA-AUC U | weighted | fused | ΔAUC fused | acc@first U | fused | Δacc |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for r in rows:
        lines.append(
            "| {d} | {ua} | {wa} | {fa} | {da} | {uc} | {fc} | {dc} |".format(
                d=r["dataset"],
                ua=_fmt(r["ova_auc_uniform"]), wa=_fmt(r["ova_auc_weighted"]),
                fa=_fmt(r["ova_auc_fused"]), da=_fmt(r["delta_auc_fused"]),
                uc=_fmt(r["acc_first_uniform"]), fc=_fmt(r["acc_first_fused"]),
                dc=_fmt(r["delta_acc_fused"]),
            )
        )
    lines.extend([
        "",
        f"mean Δ OVA-AUC fused: **{_fmt(payload['mean_delta_auc_fused'])}**",
        f"mean Δ acc@first fused: **{_fmt(payload['mean_delta_acc_fused'])}**",
        "",
        "## Per-category (last seed, fused vs uniform)",
        "",
    ])
    for b in runs[-1]["boards"]:
        lines.append(f"### {b['dataset']}")
        lines.append("")
        lines.append("| cat | n_pos | AUC fused | AUC U | acc fused | acc U |")
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
        "python3 scripts/run_text_ova_nfirst.py --n-first 100 --multiseed",
        "```",
        "",
    ])
    md = "\n".join(lines)
    (OUT / "README.md").write_text(md)
    DOCS.mkdir(parents=True, exist_ok=True)
    (DOCS / "Text_OVA_nfirst.md").write_text(md)

    # plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        names = [r["dataset"] for r in rows]
        x = np.arange(len(names))
        w = 0.35
        fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.0), constrained_layout=True)
        ax = axes[0]
        ax.bar(x - w/2, [r["ova_auc_uniform"] for r in rows], w, color="#6B7280", label="uniform")
        ax.bar(x + w/2, [r["ova_auc_fused"] for r in rows], w, color="#DC2626", label="sample-level fused")
        ax.set_xticks(x); ax.set_xticklabels(names, rotation=25, ha="right")
        ax.set_ylabel("OVA AUC"); ax.set_ylim(0.4, 1.05)
        ax.set_title(f"Text multiclass — OVA AUC on first n={args.n_first}")
        ax.legend(frameon=False, fontsize=8)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.grid(True, axis="y", alpha=0.25, ls="--")
        ax = axes[1]
        deltas = [r["delta_auc_fused"] for r in rows]
        ax.bar(x, deltas, color=["#059669" if v >= 0 else "#DC2626" for v in deltas])
        ax.axhline(0, color="#374151", lw=1)
        ax.set_xticks(x); ax.set_xticklabels(names, rotation=25, ha="right")
        ax.set_ylabel("Δ OVA-AUC (fused−uniform)")
        ax.set_title("Sample-level weights lift on new-category first-n")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.grid(True, axis="y", alpha=0.25, ls="--")
        fig.suptitle(f"{LOSS_NAME} · text only · no tabular", fontsize=12, fontweight="bold")
        figdir = OUT / "figures"; figdir.mkdir(parents=True, exist_ok=True)
        p = figdir / "text_ova_nfirst100.png"
        fig.savefig(p, dpi=170, bbox_inches="tight"); plt.close(fig)
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, ARTIFACTS / p.name)
        shutil.copy2(p, DOCS / p.name)
    except Exception as e:
        print("plot warn:", e)

    try:
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        shutil.copy2(OUT / "text_ova_nfirst.json", ARTIFACTS / "text_ova_nfirst.json")
        shutil.copy2(OUT / "README.md", ARTIFACTS / "Text_OVA_nfirst.md")
    except Exception:
        pass
    print(md)
    print("wrote", OUT / "text_ova_nfirst.json")


if __name__ == "__main__":
    main()
