#!/usr/bin/env python3
"""Feature-level / sample-level / EWMA decay weights — FULL held-out LIVE accuracy.

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
    _DEFAULT_EWMA_LAM,
    _FULL_N_PER_NOVEL,
    run_text_full_acc_benchmark,
)
from weighted_cold_start import SEED  # noqa: E402

OUT = ROOT / "results" / "text_feature_sample_weights"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")
SCHEMES = ("uniform", "sample_w", "feature_w", "both", "ewma", "ewma_feat")


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
    ap.add_argument("--ewma-lam", type=float, default=_DEFAULT_EWMA_LAM)
    ap.add_argument("--no-incremental", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    seeds = [args.seed, args.seed + 1, args.seed + 2] if args.multiseed else [args.seed]
    runs = []
    for s in seeds:
        print(
            "=== text FULL-acc seed", s,
            "n_per_novel", args.n_per_novel,
            "ewma_lam", args.ewma_lam,
            flush=True,
        )
        runs.append(
            run_text_full_acc_benchmark(
                seed=s,
                n_per_novel=args.n_per_novel,
                schemes=SCHEMES,
                alpha_feat=args.alpha_feat,
                ewma_lam=args.ewma_lam,
                run_incremental=not args.no_incremental,
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

    # aggregate incremental curves: mean Δacc(k) across boards×seeds
    incr_agg = defaultdict(lambda: defaultdict(list))  # sch -> k -> [delta]
    incr_abs = defaultdict(lambda: defaultdict(list))  # sch -> k -> [acc]
    for run in runs:
        for inc in run.get("incremental") or []:
            for sch, pts in (inc.get("deltas_vs_uniform") or {}).items():
                for p in pts:
                    incr_agg[sch][p["k"]].append(p["delta_acc"])
            for sch, pts in (inc.get("curve") or {}).items():
                for p in pts:
                    incr_abs[sch][p["k"]].append(p["acc"])

    incr_summary = {}
    for sch, by_k in incr_agg.items():
        incr_summary[sch] = [
            {
                "k": int(k),
                "mean_delta_acc": float(np.mean(vs)),
                "mean_acc": float(np.mean(incr_abs[sch][k])) if incr_abs[sch][k] else None,
                "mean_acc_uniform": (
                    float(np.mean(incr_abs["uniform"][k])) if incr_abs["uniform"][k] else None
                ),
            }
            for k, vs in sorted(by_k.items())
        ]

    def _mean_key(key):
        vals = [r[key] for r in rows if key in r and r[key] == r[key]]
        return float(np.mean(vals)) if vals else None

    payload = {
        "domain": "text_multiclass_only",
        "tabular_excluded": True,
        "loss": LOSS_NAME,
        "kpi": "FULL held-out LIVE overall / known / novel accuracy (NOT first-n subset)",
        "ewma_formula": "w_i = exp(-lam * n_idx_novel_category)",
        "ewma_lam": args.ewma_lam,
        "weight_channels": {
            "sample_w": "FSDS membership → w_i in partial_fit",
            "feature_w": "RF-domain VIMP → x_j ← x_j·(1+αṽ_j)",
            "both": "feature_w + sample_w",
            "ewma": "exp(-λ·n_idx) on novel-category arrival index",
            "ewma_feat": "ewma + feature_w",
        },
        "n_per_novel": args.n_per_novel,
        "alpha_feat": args.alpha_feat,
        "seeds": seeds,
        "runs": runs,
        "by_dataset": rows,
        "incremental_mean_delta": incr_summary,
        "mean_acc": {sch: _mean_key(f"acc_{sch}") for sch in SCHEMES},
        "mean_novel_acc": {sch: _mean_key(f"novel_acc_{sch}") for sch in SCHEMES},
        "mean_delta_acc": {
            sch: _mean_key(f"delta_acc_{sch}") for sch in SCHEMES if sch != "uniform"
        },
        "mean_delta_novel_acc": {
            sch: _mean_key(f"delta_novel_{sch}") for sch in SCHEMES if sch != "uniform"
        },
        "takeaway": (
            "EWMA decay w=exp(-λ·n_idx) puts CE mass on earliest novel-category hits; "
            "compare FULL overall/novel accuracy and incremental Δacc(k) vs uniform."
        ),
    }
    (OUT / "text_feature_sample_weights.json").write_text(json.dumps(payload, indent=2))

    lines = [
        "# Feature / sample / EWMA weights — FULL accuracy",
        "",
        payload["takeaway"],
        "",
        (
            f"KPI = **overall / known / novel accuracy on full held-out LIVE** "
            f"(not first-n). n_per_novel≈**{args.n_per_novel}**, seeds={seeds}, "
            f"loss=`{LOSS_NAME}`, α_feat={args.alpha_feat}, "
            f"EWMA `w=exp(-λ·n_idx)` with λ=**{args.ewma_lam}**."
        ),
        "",
        "## Overall accuracy (full eval)",
        "",
        "| dataset | U | sample | feature | both | **ewma** | ewma+feat | Δewma | Δfeat |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for r in rows:
        lines.append(
            "| {d} | {u} | {s} | {f} | {b} | {e} | {ef} | {de} | {df} |".format(
                d=r["dataset"],
                u=_fmt(r["acc_uniform"]),
                s=_fmt(r["acc_sample_w"]),
                f=_fmt(r["acc_feature_w"]),
                b=_fmt(r["acc_both"]),
                e=_fmt(r["acc_ewma"]),
                ef=_fmt(r["acc_ewma_feat"]),
                de=_fmt(r.get("delta_acc_ewma")),
                df=_fmt(r.get("delta_acc_feature_w")),
            )
        )
    lines.extend([
        "",
        "mean overall acc: "
        + "  ".join(f"{sch}={_fmt(payload['mean_acc'][sch])}" for sch in SCHEMES),
        "mean Δacc vs U: "
        + "  ".join(
            f"{sch}={_fmt(payload['mean_delta_acc'][sch])}"
            for sch in SCHEMES if sch != "uniform"
        ),
        "",
        "## Novel-class accuracy (full eval)",
        "",
        "| dataset | U | sample | feature | both | **ewma** | ewma+feat | Δewma | Δfeat |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ])
    for r in rows:
        lines.append(
            "| {d} | {u} | {s} | {f} | {b} | {e} | {ef} | {de} | {df} |".format(
                d=r["dataset"],
                u=_fmt(r["novel_acc_uniform"]),
                s=_fmt(r["novel_acc_sample_w"]),
                f=_fmt(r["novel_acc_feature_w"]),
                b=_fmt(r["novel_acc_both"]),
                e=_fmt(r["novel_acc_ewma"]),
                ef=_fmt(r["novel_acc_ewma_feat"]),
                de=_fmt(r.get("delta_novel_ewma")),
                df=_fmt(r.get("delta_novel_feature_w")),
            )
        )
    lines.extend([
        "",
        "mean novel_acc: "
        + "  ".join(f"{sch}={_fmt(payload['mean_novel_acc'][sch])}" for sch in SCHEMES),
        "mean Δnovel vs U: "
        + "  ".join(
            f"{sch}={_fmt(payload['mean_delta_novel_acc'][sch])}"
            for sch in SCHEMES if sch != "uniform"
        ),
        "",
        "## Incremental Δacc(k) — EWMA early-stage lift",
        "",
        "Train on known + first-k novel hits/category; score fixed FULL held-out eval.",
        "",
        "| k | Δewma | Δewma+feat | Δfeature | acc U | acc ewma |",
        "| --- | --- | --- | --- | --- | --- |",
    ])
    # align k rows
    ks = sorted({p["k"] for pts in incr_summary.values() for p in pts})
    by_sch_k = {
        sch: {p["k"]: p for p in pts} for sch, pts in incr_summary.items()
    }
    for k in ks:
        pe = by_sch_k.get("ewma", {}).get(k, {})
        pef = by_sch_k.get("ewma_feat", {}).get(k, {})
        pf = by_sch_k.get("feature_w", {}).get(k, {})
        lines.append(
            "| {k} | {de} | {def_} | {df} | {u} | {e} |".format(
                k=k,
                de=_fmt(pe.get("mean_delta_acc")),
                def_=_fmt(pef.get("mean_delta_acc")),
                df=_fmt(pf.get("mean_delta_acc")),
                u=_fmt(pe.get("mean_acc_uniform")),
                e=_fmt(pe.get("mean_acc")),
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

    # plots
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        names = [r["dataset"] for r in rows]
        x = np.arange(len(names))
        w = 0.13
        fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.4), constrained_layout=True)

        ax = axes[0]
        scheme_style = [
            ("uniform", "#6B7280"),
            ("sample_w", "#2563EB"),
            ("feature_w", "#DC2626"),
            ("both", "#059669"),
            ("ewma", "#D97706"),
            ("ewma_feat", "#7C3AED"),
        ]
        for i, (sch, color) in enumerate(scheme_style):
            ax.bar(
                x + (i - 2.5) * w,
                [r[f"acc_{sch}"] for r in rows],
                w, color=color, label=sch,
            )
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=25, ha="right")
        ax.set_ylabel("Overall accuracy (full held-out LIVE)")
        ax.set_ylim(0.0, 1.05)
        ax.set_title("FULL overall accuracy — +EWMA")
        ax.legend(frameon=False, fontsize=7, ncol=3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, axis="y", alpha=0.25, ls="--")

        ax = axes[1]
        delta_style = [
            ("sample_w", "#2563EB", "Δ sample"),
            ("feature_w", "#DC2626", "Δ feature"),
            ("ewma", "#D97706", "Δ ewma"),
            ("ewma_feat", "#7C3AED", "Δ ewma+feat"),
        ]
        wd = 0.18
        for i, (sch, color, label) in enumerate(delta_style):
            ax.bar(
                x + (i - 1.5) * wd,
                [r.get(f"delta_acc_{sch}", 0.0) for r in rows],
                wd, color=color, label=label,
            )
        ax.axhline(0, color="#374151", lw=1)
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=25, ha="right")
        ax.set_ylabel("Δ overall acc vs uniform")
        ax.set_title("Incremental lift on FULL accuracy")
        ax.legend(frameon=False, fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, axis="y", alpha=0.25, ls="--")

        fig.suptitle(
            f"{LOSS_NAME} · EWMA w=exp(-λ·n_idx) λ={args.ewma_lam} · text only",
            fontsize=12, fontweight="bold",
        )
        figdir = OUT / "figures"
        figdir.mkdir(parents=True, exist_ok=True)
        p = figdir / "text_feature_sample_full_acc.png"
        fig.savefig(p, dpi=170, bbox_inches="tight")
        plt.close(fig)

        # incremental curve figure
        fig, ax = plt.subplots(figsize=(7.5, 4.2), constrained_layout=True)
        for sch, color, label in [
            ("ewma", "#D97706", "Δ ewma"),
            ("ewma_feat", "#7C3AED", "Δ ewma+feat"),
            ("feature_w", "#DC2626", "Δ feature"),
        ]:
            pts = incr_summary.get(sch) or []
            if not pts:
                continue
            ks_ = [p["k"] for p in pts]
            ys = [p["mean_delta_acc"] for p in pts]
            ax.plot(ks_, ys, "-o", color=color, lw=2, ms=5, label=label)
        ax.axhline(0, color="#374151", lw=1)
        ax.set_xlabel("k = first novel hits / category (train)")
        ax.set_ylabel("mean Δ overall acc vs uniform")
        ax.set_title(f"Incremental EWMA lift  (λ={args.ewma_lam})")
        ax.legend(frameon=False, fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, alpha=0.25, ls="--")
        p_inc = figdir / "text_ewma_incremental_delta.png"
        fig.savefig(p_inc, dpi=170, bbox_inches="tight")
        plt.close(fig)

        # novel-acc companion
        fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.2), constrained_layout=True)
        ax = axes[0]
        for i, (sch, color) in enumerate(scheme_style):
            ax.bar(
                x + (i - 2.5) * w,
                [r[f"novel_acc_{sch}"] for r in rows],
                w, color=color, label=sch,
            )
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=25, ha="right")
        ax.set_ylabel("Novel-class accuracy (full held-out)")
        ax.set_ylim(0.0, 1.05)
        ax.set_title("Novel accuracy on FULL LIVE eval")
        ax.legend(frameon=False, fontsize=7, ncol=3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, axis="y", alpha=0.25, ls="--")

        ax = axes[1]
        for i, (sch, color, label) in enumerate(delta_style):
            ax.bar(
                x + (i - 1.5) * wd,
                [r.get(f"delta_novel_{sch}", 0.0) for r in rows],
                wd, color=color, label=label,
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
            f"{LOSS_NAME} · novel-class FULL accuracy · +EWMA",
            fontsize=12, fontweight="bold",
        )
        p2 = figdir / "text_feature_sample_novel_acc.png"
        fig.savefig(p2, dpi=170, bbox_inches="tight")
        plt.close(fig)

        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        for src in (p, p2, p_inc):
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
    print("incremental_mean_delta", incr_summary)


if __name__ == "__main__":
    main()
