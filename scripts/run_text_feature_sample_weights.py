#!/usr/bin/env python3
"""Feature-level vs sample-level weights — FULL held-out LIVE accuracy.

  python3 scripts/run_text_feature_sample_weights.py --multiseed
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

from text_feature_sample_weights import (  # noqa: E402
    LOSS_NAME,
    _FULL_N_PER_NOVEL,
    run_text_full_acc_benchmark,
)
from weighted_cold_start import SEED  # noqa: E402

OUT = ROOT / "results" / "text_feature_sample_weights"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")
SCHEMES = ("uniform", "sample_w", "feature_w", "both")


def _fmt(x):
    try:
        if x is None or x != x:
            return "—"
        return "%.3f" % float(x)
    except Exception:
        return "—"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--n-per-novel", type=int, default=_FULL_N_PER_NOVEL)
    ap.add_argument("--multiseed", action="store_true")
    ap.add_argument("--alpha-feat", type=float, default=3.0)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    seeds = [args.seed, args.seed + 1, args.seed + 2] if args.multiseed else [args.seed]
    runs = []
    for s in seeds:
        print("=== text FULL-acc seed", s, "n_per_novel", args.n_per_novel, flush=True)
        runs.append(
            run_text_full_acc_benchmark(
                seed=s,
                n_per_novel=args.n_per_novel,
                schemes=SCHEMES,
                alpha_feat=args.alpha_feat,
            )
        )

    # aggregate absolute + delta per dataset across seeds
    agg = defaultdict(lambda: {sch: {"acc": [], "novel": [], "known": []} for sch in SCHEMES})
    for run in runs:
        for b in run["boards"]:
            name = b["dataset"]
            for sch in SCHEMES:
                row = b["schemes"].get(sch, {})
                agg[name][sch]["acc"].append(row.get("acc", float("nan")))
                agg[name][sch]["novel"].append(row.get("novel_acc", float("nan")))
                agg[name][sch]["known"].append(row.get("known_acc", float("nan")))

    rows = []
    for name, d in agg.items():
        u_acc = np.asarray(d["uniform"]["acc"], float)
        row = {
            "dataset": name,
            "n_seeds": len(u_acc),
        }
        for sch in SCHEMES:
            a = np.asarray(d[sch]["acc"], float)
            n = np.asarray(d[sch]["novel"], float)
            k = np.asarray(d[sch]["known"], float)
            row[f"acc_{sch}"] = float(np.nanmean(a))
            row[f"novel_acc_{sch}"] = float(np.nanmean(n))
            row[f"known_acc_{sch}"] = float(np.nanmean(k))
            if sch != "uniform":
                row[f"delta_acc_{sch}"] = float(np.nanmean(a - u_acc))
                row[f"delta_novel_{sch}"] = float(
                    np.nanmean(n - np.asarray(d["uniform"]["novel"], float))
                )
                row[f"P_acc_gt_U_{sch}"] = float(np.nanmean(a > u_acc))
        rows.append(row)

    def _mean_key(key):
        vals = [r[key] for r in rows if key in r and r[key] == r[key]]
        return float(np.mean(vals)) if vals else None

    payload = {
        "domain": "text_multiclass_only",
        "tabular_excluded": True,
        "loss": LOSS_NAME,
        "kpi": "FULL held-out LIVE overall / known / novel accuracy (NOT first-n subset)",
        "weight_channels": {
            "sample_w": "FSDS membership → w_i in partial_fit",
            "feature_w": "RF-domain VIMP → x_j ← x_j·(1+αṽ_j)",
            "both": "feature_w + sample_w",
        },
        "n_per_novel": args.n_per_novel,
        "alpha_feat": args.alpha_feat,
        "seeds": seeds,
        "runs": runs,
        "by_dataset": rows,
        "mean_acc": {sch: _mean_key(f"acc_{sch}") for sch in SCHEMES},
        "mean_novel_acc": {sch: _mean_key(f"novel_acc_{sch}") for sch in SCHEMES},
        "mean_delta_acc": {
            "sample_w": _mean_key("delta_acc_sample_w"),
            "feature_w": _mean_key("delta_acc_feature_w"),
            "both": _mean_key("delta_acc_both"),
        },
        "mean_delta_novel_acc": {
            "sample_w": _mean_key("delta_novel_sample_w"),
            "feature_w": _mean_key("delta_novel_feature_w"),
            "both": _mean_key("delta_novel_both"),
        },
        "takeaway": (
            "Compare feature-level VIMP scales vs sample-level CE weights on FULL "
            "held-out LIVE accuracy for text multiclass category shift."
        ),
    }
    (OUT / "text_feature_sample_weights.json").write_text(json.dumps(payload, indent=2))

    lines = [
        "# Feature-level vs sample-level weights — FULL accuracy",
        "",
        payload["takeaway"],
        "",
        (
            f"KPI = **overall / known / novel accuracy on full held-out LIVE** "
            f"(not first-n subset). n_per_novel≈**{args.n_per_novel}**, "
            f"seeds={seeds}, loss=`{LOSS_NAME}`, α_feat={args.alpha_feat}."
        ),
        "",
        "## Overall accuracy (full eval)",
        "",
        "| dataset | uniform | sample_w | feature_w | both | Δfeat | Δsample | Δboth |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for r in rows:
        lines.append(
            "| {d} | {u} | {s} | {f} | {b} | {df} | {ds} | {db} |".format(
                d=r["dataset"],
                u=_fmt(r["acc_uniform"]),
                s=_fmt(r["acc_sample_w"]),
                f=_fmt(r["acc_feature_w"]),
                b=_fmt(r["acc_both"]),
                df=_fmt(r.get("delta_acc_feature_w")),
                ds=_fmt(r.get("delta_acc_sample_w")),
                db=_fmt(r.get("delta_acc_both")),
            )
        )
    lines.extend([
        "",
        f"mean overall acc: U={_fmt(payload['mean_acc']['uniform'])}  "
        f"sample={_fmt(payload['mean_acc']['sample_w'])}  "
        f"feature={_fmt(payload['mean_acc']['feature_w'])}  "
        f"both={_fmt(payload['mean_acc']['both'])}",
        f"mean Δacc vs U: feature={_fmt(payload['mean_delta_acc']['feature_w'])}  "
        f"sample={_fmt(payload['mean_delta_acc']['sample_w'])}  "
        f"both={_fmt(payload['mean_delta_acc']['both'])}",
        "",
        "## Novel-class accuracy (full eval)",
        "",
        "| dataset | uniform | sample_w | feature_w | both | Δfeat | Δsample | Δboth |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ])
    for r in rows:
        lines.append(
            "| {d} | {u} | {s} | {f} | {b} | {df} | {ds} | {db} |".format(
                d=r["dataset"],
                u=_fmt(r["novel_acc_uniform"]),
                s=_fmt(r["novel_acc_sample_w"]),
                f=_fmt(r["novel_acc_feature_w"]),
                b=_fmt(r["novel_acc_both"]),
                df=_fmt(r.get("delta_novel_feature_w")),
                ds=_fmt(r.get("delta_novel_sample_w")),
                db=_fmt(r.get("delta_novel_both")),
            )
        )
    lines.extend([
        "",
        f"mean novel_acc: U={_fmt(payload['mean_novel_acc']['uniform'])}  "
        f"sample={_fmt(payload['mean_novel_acc']['sample_w'])}  "
        f"feature={_fmt(payload['mean_novel_acc']['feature_w'])}  "
        f"both={_fmt(payload['mean_novel_acc']['both'])}",
        f"mean Δnovel vs U: feature={_fmt(payload['mean_delta_novel_acc']['feature_w'])}  "
        f"sample={_fmt(payload['mean_delta_novel_acc']['sample_w'])}  "
        f"both={_fmt(payload['mean_delta_novel_acc']['both'])}",
        "",
        "## Known retention (full eval)",
        "",
        "| dataset | uniform | sample_w | feature_w | both |",
        "| --- | --- | --- | --- | --- |",
    ])
    for r in rows:
        lines.append(
            "| {d} | {u} | {s} | {f} | {b} |".format(
                d=r["dataset"],
                u=_fmt(r["known_acc_uniform"]),
                s=_fmt(r["known_acc_sample_w"]),
                f=_fmt(r["known_acc_feature_w"]),
                b=_fmt(r["known_acc_both"]),
            )
        )
    lines.extend([
        "",
        "```bash",
        "python3 scripts/run_text_feature_sample_weights.py --multiseed",
        "```",
        "",
    ])
    md = "\n".join(lines)
    (OUT / "README.md").write_text(md)
    DOCS.mkdir(parents=True, exist_ok=True)
    (DOCS / "Text_feature_sample_weights.md").write_text(md)

    # plot overall + novel Δ
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        names = [r["dataset"] for r in rows]
        x = np.arange(len(names))
        w = 0.2
        fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.2), constrained_layout=True)

        ax = axes[0]
        for i, (sch, color) in enumerate([
            ("uniform", "#6B7280"),
            ("sample_w", "#2563EB"),
            ("feature_w", "#DC2626"),
            ("both", "#059669"),
        ]):
            ax.bar(
                x + (i - 1.5) * w,
                [r[f"acc_{sch}"] for r in rows],
                w, color=color, label=sch,
            )
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=25, ha="right")
        ax.set_ylabel("Overall accuracy (full held-out LIVE)")
        ax.set_ylim(0.0, 1.05)
        ax.set_title("Text multiclass — FULL overall accuracy")
        ax.legend(frameon=False, fontsize=8, ncol=2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, axis="y", alpha=0.25, ls="--")

        ax = axes[1]
        for i, (sch, color, label) in enumerate([
            ("sample_w", "#2563EB", "Δ sample_w"),
            ("feature_w", "#DC2626", "Δ feature_w"),
            ("both", "#059669", "Δ both"),
        ]):
            ax.bar(
                x + (i - 1) * w,
                [r.get(f"delta_acc_{sch}", 0.0) for r in rows],
                w, color=color, label=label,
            )
        ax.axhline(0, color="#374151", lw=1)
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=25, ha="right")
        ax.set_ylabel("Δ overall acc vs uniform")
        ax.set_title("Lift on FULL accuracy (not first-n)")
        ax.legend(frameon=False, fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, axis="y", alpha=0.25, ls="--")

        fig.suptitle(
            f"{LOSS_NAME} · feature-level VIMP vs sample-level · text only",
            fontsize=12, fontweight="bold",
        )
        figdir = OUT / "figures"
        figdir.mkdir(parents=True, exist_ok=True)
        p = figdir / "text_feature_sample_full_acc.png"
        fig.savefig(p, dpi=170, bbox_inches="tight")
        plt.close(fig)

        # novel-acc companion figure
        fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.2), constrained_layout=True)
        ax = axes[0]
        for i, (sch, color) in enumerate([
            ("uniform", "#6B7280"),
            ("sample_w", "#2563EB"),
            ("feature_w", "#DC2626"),
            ("both", "#059669"),
        ]):
            ax.bar(
                x + (i - 1.5) * w,
                [r[f"novel_acc_{sch}"] for r in rows],
                w, color=color, label=sch,
            )
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=25, ha="right")
        ax.set_ylabel("Novel-class accuracy (full held-out)")
        ax.set_ylim(0.0, 1.05)
        ax.set_title("Novel accuracy on FULL LIVE eval")
        ax.legend(frameon=False, fontsize=8, ncol=2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, axis="y", alpha=0.25, ls="--")

        ax = axes[1]
        for i, (sch, color, label) in enumerate([
            ("sample_w", "#2563EB", "Δ sample_w"),
            ("feature_w", "#DC2626", "Δ feature_w"),
            ("both", "#059669", "Δ both"),
        ]):
            ax.bar(
                x + (i - 1) * w,
                [r.get(f"delta_novel_{sch}", 0.0) for r in rows],
                w, color=color, label=label,
            )
        ax.axhline(0, color="#374151", lw=1)
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=25, ha="right")
        ax.set_ylabel("Δ novel_acc vs uniform")
        ax.set_title("Novel lift on FULL eval")
        ax.legend(frameon=False, fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, axis="y", alpha=0.25, ls="--")
        fig.suptitle(
            f"{LOSS_NAME} · novel-class FULL accuracy · text only",
            fontsize=12, fontweight="bold",
        )
        p2 = figdir / "text_feature_sample_novel_acc.png"
        fig.savefig(p2, dpi=170, bbox_inches="tight")
        plt.close(fig)

        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        for src in (p, p2):
            shutil.copy2(src, ARTIFACTS / src.name)
            shutil.copy2(src, DOCS / src.name)
    except Exception as e:
        print("plot warn:", e)

    try:
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        shutil.copy2(OUT / "text_feature_sample_weights.json", ARTIFACTS / "text_feature_sample_weights.json")
        shutil.copy2(OUT / "README.md", ARTIFACTS / "Text_feature_sample_weights.md")
    except Exception:
        pass

    print(md)
    print("wrote", OUT / "text_feature_sample_weights.json")
    print("mean_acc", payload["mean_acc"])
    print("mean_delta_acc", payload["mean_delta_acc"])


if __name__ == "__main__":
    main()
