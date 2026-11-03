#!/usr/bin/env python3
"""Feature vs sample vs both weights — FULL holdout acc / OVA-AUC curves.

No EWMA. Primary KPI = overall accuracy + mean novel OVA-AUC on full LIVE holdout.

  python3 scripts/run_text_full_weight_curves.py --multiseed
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

from text_full_weight_curves import (  # noqa: E402
    LOSS_NAME,
    SCHEMES,
    run_text_full_weight_curves,
)
from weighted_cold_start import SEED  # noqa: E402

OUT = ROOT / "results" / "text_full_weight_curves"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")

SCHEME_STYLE = [
    ("uniform", "#6B7280", "uniform"),
    ("feature_w", "#DC2626", "feature_w"),
    ("sample_w", "#2563EB", "sample_w"),
    ("both", "#059669", "both"),
]


def _fmt(x):
    try:
        if x is None or x != x:
            return "—"
        return "%.3f" % float(x)
    except Exception:
        return "—"


def _aggregate(runs):
    """Mean curves / endpoint metrics across seeds × boards."""
    # board -> scheme -> budget -> list of metric dicts
    by = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    endpoints = defaultdict(lambda: defaultdict(list))  # board -> scheme -> [endpoint]
    for run in runs:
        for board in run["boards"]:
            name = board["dataset"]
            for sch, pts in board["curves"].items():
                for p in pts:
                    by[name][sch][int(p["budget"])].append(p)
            for sch, ep in board["endpoint"].items():
                endpoints[name][sch].append(ep)

    rows = []
    for name, sch_map in by.items():
        budgets = sorted({b for sch in sch_map.values() for b in sch})
        curves = {}
        for sch in SCHEMES:
            curves[sch] = []
            for b in budgets:
                pts = sch_map[sch].get(b) or []
                if not pts:
                    continue
                curves[sch].append({
                    "budget": b,
                    "acc": float(np.nanmean([p["acc"] for p in pts])),
                    "ova_auc_novel_mean": float(
                        np.nanmean([p["ova_auc_novel_mean"] for p in pts])
                    ),
                    "novel_acc": float(np.nanmean([p.get("novel_acc", np.nan) for p in pts])),
                    "known_acc": float(np.nanmean([p.get("known_acc", np.nan) for p in pts])),
                })
        ep = {}
        for sch in SCHEMES:
            eps = endpoints[name][sch]
            if not eps:
                continue
            ep[sch] = {
                "acc": float(np.nanmean([e.get("acc", np.nan) for e in eps])),
                "ova_auc_novel_mean": float(
                    np.nanmean([e.get("ova_auc_novel_mean", np.nan) for e in eps])
                ),
                "novel_acc": float(np.nanmean([e.get("novel_acc", np.nan) for e in eps])),
                "known_acc": float(np.nanmean([e.get("known_acc", np.nan) for e in eps])),
            }
        u = ep.get("uniform", {})
        deltas = {}
        for sch in SCHEMES:
            if sch == "uniform" or sch not in ep:
                continue
            deltas[sch] = {
                "delta_acc": float(ep[sch]["acc"] - u["acc"]) if u else None,
                "delta_ova_auc": (
                    float(ep[sch]["ova_auc_novel_mean"] - u["ova_auc_novel_mean"])
                    if u and ep[sch]["ova_auc_novel_mean"] == ep[sch]["ova_auc_novel_mean"]
                    and u.get("ova_auc_novel_mean") == u.get("ova_auc_novel_mean")
                    else None
                ),
            }
        rows.append({
            "dataset": name,
            "curves": curves,
            "endpoint": ep,
            "deltas_vs_uniform": deltas,
        })
    return rows


def _mean_curve(rows, scheme, metric):
    """Average metric over datasets at each budget."""
    by_b = defaultdict(list)
    for r in rows:
        for p in r["curves"].get(scheme) or []:
            by_b[int(p["budget"])].append(p[metric])
    return sorted(
        ({"budget": b, metric: float(np.nanmean(vs))} for b, vs in by_b.items()),
        key=lambda d: d["budget"],
    )


def _plot(rows, payload, out_dir: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    names = [r["dataset"] for r in rows]

    # --- 1) Acc + AUC curves (mean over datasets) ---
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), constrained_layout=True)
    for ax, metric, ylab, title in (
        (axes[0], "acc", "Overall accuracy (full LIVE holdout)", "Accuracy vs novel budget"),
        (axes[1], "ova_auc_novel_mean", "Mean novel OVA-AUC (full LIVE holdout)", "OVA-AUC vs novel budget"),
    ):
        for sch, color, label in SCHEME_STYLE:
            pts = _mean_curve(rows, sch, metric)
            if not pts:
                continue
            ax.plot(
                [p["budget"] for p in pts],
                [p[metric] for p in pts],
                "-o", color=color, lw=2.2, ms=5, label=label,
            )
        ax.set_xlabel("Novel adapt samples (budget)")
        ax.set_ylabel(ylab)
        ax.set_title(title)
        ax.legend(frameon=False, fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, alpha=0.25, ls="--")
    fig.suptitle(
        f"{LOSS_NAME} · feature vs sample vs both · NO EWMA · text only",
        fontsize=12, fontweight="bold",
    )
    p_curves = out_dir / "text_full_weight_acc_auc_curves.png"
    fig.savefig(p_curves, dpi=170, bbox_inches="tight")
    plt.close(fig)

    # --- 2) Per-dataset endpoint bars ---
    x = np.arange(len(names))
    w = 0.18
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4), constrained_layout=True)
    for ax, metric, ylab, title in (
        (axes[0], "acc", "Overall accuracy", "Endpoint overall accuracy"),
        (axes[1], "ova_auc_novel_mean", "Mean novel OVA-AUC", "Endpoint novel OVA-AUC"),
    ):
        for i, (sch, color, label) in enumerate(SCHEME_STYLE):
            vals = [r["endpoint"].get(sch, {}).get(metric, np.nan) for r in rows]
            ax.bar(x + (i - 1.5) * w, vals, w, color=color, label=label)
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=25, ha="right")
        ax.set_ylabel(ylab)
        ax.set_ylim(0.0, 1.05)
        ax.set_title(title)
        ax.legend(frameon=False, fontsize=8, ncol=2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, axis="y", alpha=0.25, ls="--")
    fig.suptitle(
        f"{LOSS_NAME} · FULL holdout endpoint · NO EWMA",
        fontsize=12, fontweight="bold",
    )
    p_bars = out_dir / "text_full_weight_endpoint_bars.png"
    fig.savefig(p_bars, dpi=170, bbox_inches="tight")
    plt.close(fig)

    # --- 3) Δacc / ΔAUC vs uniform ---
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), constrained_layout=True)
    delta_style = [
        ("feature_w", "#DC2626", "Δ feature_w"),
        ("sample_w", "#2563EB", "Δ sample_w"),
        ("both", "#059669", "Δ both"),
    ]
    wd = 0.22
    for ax, key, ylab, title in (
        (axes[0], "delta_acc", "Δ overall acc vs uniform", "Δ accuracy (endpoint)"),
        (axes[1], "delta_ova_auc", "Δ mean novel OVA-AUC vs uniform", "Δ OVA-AUC (endpoint)"),
    ):
        for i, (sch, color, label) in enumerate(delta_style):
            vals = [
                (r["deltas_vs_uniform"].get(sch) or {}).get(key, 0.0) for r in rows
            ]
            ax.bar(x + (i - 1) * wd, vals, wd, color=color, label=label)
        ax.axhline(0, color="#374151", lw=1)
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=25, ha="right")
        ax.set_ylabel(ylab)
        ax.set_title(title)
        ax.legend(frameon=False, fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, axis="y", alpha=0.25, ls="--")
    fig.suptitle("Lift vs uniform on FULL LIVE holdout (no EWMA)", fontsize=12, fontweight="bold")
    p_delta = out_dir / "text_full_weight_delta.png"
    fig.savefig(p_delta, dpi=170, bbox_inches="tight")
    plt.close(fig)

    # also per-dataset small multiples for curves
    n = len(rows)
    fig, axes = plt.subplots(n, 2, figsize=(10.5, 2.6 * n), constrained_layout=True)
    if n == 1:
        axes = np.asarray([axes])
    for i, r in enumerate(rows):
        for j, (metric, ylab) in enumerate((
            ("acc", "overall acc"),
            ("ova_auc_novel_mean", "OVA-AUC"),
        )):
            ax = axes[i, j]
            for sch, color, label in SCHEME_STYLE:
                pts = r["curves"].get(sch) or []
                if not pts:
                    continue
                ax.plot(
                    [p["budget"] for p in pts],
                    [p[metric] for p in pts],
                    "-o", color=color, lw=1.8, ms=3.5,
                    label=label if i == 0 else None,
                )
            ax.set_title(f"{r['dataset']} · {ylab}", fontsize=10)
            if i == n - 1:
                ax.set_xlabel("novel budget")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.grid(True, alpha=0.25, ls="--")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False, fontsize=9)
    fig.suptitle("Per-dataset FULL holdout curves · NO EWMA", fontsize=12, fontweight="bold", y=1.01)
    p_per = out_dir / "text_full_weight_per_dataset_curves.png"
    fig.savefig(p_per, dpi=160, bbox_inches="tight")
    plt.close(fig)

    return [p_curves, p_bars, p_delta, p_per]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--multiseed", action="store_true")
    ap.add_argument("--n-first-board", type=int, default=400,
                    help="max novel samples / class when building FULL LIVE boards")
    ap.add_argument("--feat-alpha", type=float, default=3.0)
    ap.add_argument(
        "--budgets", type=int, nargs="+",
        default=[20, 40, 60, 80, 100, 150, 200],
    )
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    seeds = [args.seed, args.seed + 1, args.seed + 2] if args.multiseed else [args.seed]
    runs = []
    for s in seeds:
        print(
            "=== full-weight-curves seed", s,
            "budgets", args.budgets,
            "feat_alpha", args.feat_alpha,
            flush=True,
        )
        runs.append(
            run_text_full_weight_curves(
                seed=s,
                n_first_board=args.n_first_board,
                budgets=args.budgets,
                feat_alpha=args.feat_alpha,
            )
        )

    rows = _aggregate(runs)

    def _mean_ep(metric, sch):
        vals = [
            r["endpoint"].get(sch, {}).get(metric)
            for r in rows
            if r["endpoint"].get(sch, {}).get(metric) is not None
        ]
        vals = [v for v in vals if v == v]
        return float(np.mean(vals)) if vals else None

    def _mean_delta(key, sch):
        vals = [
            (r["deltas_vs_uniform"].get(sch) or {}).get(key)
            for r in rows
        ]
        vals = [v for v in vals if v is not None and v == v]
        return float(np.mean(vals)) if vals else None

    payload = {
        "domain": "text_multiclass_only",
        "tabular_excluded": True,
        "loss": LOSS_NAME,
        "ewma": False,
        "kpi": "FULL held-out LIVE overall accuracy + mean novel OVA-AUC (NOT first-n subset)",
        "schemes": {
            "uniform": "CE, flat sample weights, raw features",
            "feature_w": "X *= (1 + α · vimp) from RF-domain VIMP",
            "sample_w": "sample_weight from FSDS membership",
            "both": "feature_w + sample_w",
        },
        "feat_alpha": args.feat_alpha,
        "budgets": list(args.budgets),
        "n_first_board": args.n_first_board,
        "seeds": seeds,
        "runs": runs,
        "by_dataset": rows,
        "mean_acc": {sch: _mean_ep("acc", sch) for sch in SCHEMES},
        "mean_ova_auc": {sch: _mean_ep("ova_auc_novel_mean", sch) for sch in SCHEMES},
        "mean_delta_acc": {
            sch: _mean_delta("delta_acc", sch) for sch in SCHEMES if sch != "uniform"
        },
        "mean_delta_ova_auc": {
            sch: _mean_delta("delta_ova_auc", sch) for sch in SCHEMES if sch != "uniform"
        },
        "mean_curves": {
            sch: {
                "acc": _mean_curve(rows, sch, "acc"),
                "ova_auc": _mean_curve(rows, sch, "ova_auc_novel_mean"),
            }
            for sch in SCHEMES
        },
        "takeaway": (
            "Feature-level = X·(1+α·vimp); sample-level = membership w_i; both = product. "
            "Curves on FULL LIVE holdout overall acc + novel OVA-AUC. No EWMA."
        ),
    }
    (OUT / "text_full_weight_curves.json").write_text(json.dumps(payload, indent=2))

    lines = [
        "# Feature vs sample weights — FULL holdout acc / AUC curves",
        "",
        payload["takeaway"],
        "",
        (
            f"KPI = **overall accuracy + mean novel OVA-AUC on full held-out LIVE** "
            f"(not first-n). budgets={list(args.budgets)}, seeds={seeds}, "
            f"loss=`{LOSS_NAME}`, α_feat={args.feat_alpha}, **EWMA=off**."
        ),
        "",
        "## Endpoint overall accuracy (full eval)",
        "",
        "| dataset | uniform | feature_w | sample_w | both | Δfeat | Δsample | Δboth |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for r in rows:
        ep = r["endpoint"]
        d = r["deltas_vs_uniform"]
        lines.append(
            "| {ds} | {u} | {f} | {s} | {b} | {df} | {ds_} | {db} |".format(
                ds=r["dataset"],
                u=_fmt(ep.get("uniform", {}).get("acc")),
                f=_fmt(ep.get("feature_w", {}).get("acc")),
                s=_fmt(ep.get("sample_w", {}).get("acc")),
                b=_fmt(ep.get("both", {}).get("acc")),
                df=_fmt((d.get("feature_w") or {}).get("delta_acc")),
                ds_=_fmt((d.get("sample_w") or {}).get("delta_acc")),
                db=_fmt((d.get("both") or {}).get("delta_acc")),
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
        "## Endpoint mean novel OVA-AUC (full eval)",
        "",
        "| dataset | uniform | feature_w | sample_w | both | Δfeat | Δsample | Δboth |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ])
    for r in rows:
        ep = r["endpoint"]
        d = r["deltas_vs_uniform"]
        lines.append(
            "| {ds} | {u} | {f} | {s} | {b} | {df} | {ds_} | {db} |".format(
                ds=r["dataset"],
                u=_fmt(ep.get("uniform", {}).get("ova_auc_novel_mean")),
                f=_fmt(ep.get("feature_w", {}).get("ova_auc_novel_mean")),
                s=_fmt(ep.get("sample_w", {}).get("ova_auc_novel_mean")),
                b=_fmt(ep.get("both", {}).get("ova_auc_novel_mean")),
                df=_fmt((d.get("feature_w") or {}).get("delta_ova_auc")),
                ds_=_fmt((d.get("sample_w") or {}).get("delta_ova_auc")),
                db=_fmt((d.get("both") or {}).get("delta_ova_auc")),
            )
        )
    lines.extend([
        "",
        "mean OVA-AUC: "
        + "  ".join(f"{sch}={_fmt(payload['mean_ova_auc'][sch])}" for sch in SCHEMES),
        "mean ΔOVA-AUC vs U: "
        + "  ".join(
            f"{sch}={_fmt(payload['mean_delta_ova_auc'][sch])}"
            for sch in SCHEMES if sch != "uniform"
        ),
        "",
        "## Mean curves (acc / OVA-AUC vs novel budget)",
        "",
        "| budget | acc U | acc feat | acc sample | acc both | AUC U | AUC feat | AUC sample | AUC both |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ])
    mc = payload["mean_curves"]
    budgets = sorted({p["budget"] for p in mc["uniform"]["acc"]})
    acc_maps = {
        sch: {p["budget"]: p["acc"] for p in mc[sch]["acc"]} for sch in SCHEMES
    }
    auc_maps = {
        sch: {p["budget"]: p["ova_auc_novel_mean"] for p in mc[sch]["ova_auc"]}
        for sch in SCHEMES
    }
    for b in budgets:
        lines.append(
            "| {b} | {au} | {af} | {as_} | {ab} | {uu} | {uf} | {us} | {ub} |".format(
                b=b,
                au=_fmt(acc_maps["uniform"].get(b)),
                af=_fmt(acc_maps["feature_w"].get(b)),
                as_=_fmt(acc_maps["sample_w"].get(b)),
                ab=_fmt(acc_maps["both"].get(b)),
                uu=_fmt(auc_maps["uniform"].get(b)),
                uf=_fmt(auc_maps["feature_w"].get(b)),
                us=_fmt(auc_maps["sample_w"].get(b)),
                ub=_fmt(auc_maps["both"].get(b)),
            )
        )
    lines.extend([
        "",
        "```bash",
        "python3 scripts/run_text_full_weight_curves.py --multiseed",
        "```",
        "",
    ])
    md = "\n".join(lines)
    (OUT / "README.md").write_text(md)
    DOCS.mkdir(parents=True, exist_ok=True)
    (DOCS / "Text_full_weight_curves.md").write_text(md)

    figdir = OUT / "figures"
    try:
        paths = _plot(rows, payload, figdir)
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        for src in paths:
            shutil.copy2(src, ARTIFACTS / src.name)
            shutil.copy2(src, DOCS / src.name)
    except Exception as e:
        print("plot warn:", e)
        paths = []

    try:
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        shutil.copy2(OUT / "text_full_weight_curves.json", ARTIFACTS / "text_full_weight_curves.json")
        shutil.copy2(OUT / "README.md", ARTIFACTS / "Text_full_weight_curves.md")
    except Exception:
        pass

    print(md)
    print("wrote", OUT / "text_full_weight_curves.json")
    print("mean_acc", payload["mean_acc"])
    print("mean_ova_auc", payload["mean_ova_auc"])
    print("mean_delta_acc", payload["mean_delta_acc"])
    print("mean_delta_ova_auc", payload["mean_delta_ova_auc"])
    print("figures", [str(p) for p in paths])


if __name__ == "__main__":
    main()
