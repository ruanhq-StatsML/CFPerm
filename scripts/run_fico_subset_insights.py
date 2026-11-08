#!/usr/bin/env python3
"""Auto-run post-hoc subset localization package on FICO designs + insights.

Uses the notebook package:
  feature_localize → subset_on_driver → subset_localize (PSI/TVD/Hellinger)

  python3 scripts/run_fico_subset_insights.py
  python3 scripts/run_fico_subset_insights.py --quick
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from fico_heloc_cfperm import (  # noqa: E402
    clean_heloc,
    load_heloc,
    median_split_batches,
    top_lr_features,
)
from fico_heloc_weak_cd import make_weak_cd_batches  # noqa: E402
from posthoc_subset_localization import (  # noqa: E402
    feature_localize,
    run_posthoc_on_batches,
    subset_localize,
    subset_on_driver,
)

OUT = ROOT / "results" / "fico_heloc_cfperm" / "subset_insights"
CACHE = ROOT / "data" / "heloc_dataset_v1.csv"


def _insights_from_tables(
    *,
    design: str,
    summary: pd.DataFrame,
    subset: pd.DataFrame,
    zcut: pd.DataFrame,
    feat: pd.DataFrame,
) -> list[str]:
    lines = [f"### {design}", ""]
    if summary.empty:
        lines.append("- No summary rows.")
        return lines

    # concentration
    conc = summary["concentrated"].dropna().astype(str)
    if len(conc):
        mode = conc.value_counts().idxmax()
        lines.append(
            f"- Driver-cut concentration mostly on **`{mode}`** "
            f"({int((conc == mode).sum())}/{len(conc)} cells)."
        )
    if "contrast" in summary and summary["contrast"].notna().any():
        c = summary["contrast"].dropna()
        lines.append(
            f"- |contrast| mean=`{c.abs().mean():.4f}` "
            f"(max=`{c.abs().max():.4f}`) — larger ⇒ stronger low/high asymmetry."
        )
    if "dmse_high" in summary.columns and "dmse_low" in summary.columns:
        hi = summary["dmse_high"].dropna()
        lo = summary["dmse_low"].dropna()
        if len(hi) and len(lo):
            lines.append(
                f"- Mean ΔMSE: low=`{lo.mean():.4f}`, high=`{hi.mean():.4f}` "
                f"(high/low ratio≈`{(hi.mean() / (lo.mean() + 1e-12)):.2f}`)."
            )
    if "hit_at_k" in summary and summary["hit_at_k"].notna().any():
        lines.append(
            f"- PFI hit@k mean=`{summary['hit_at_k'].mean():.2f}` "
            f"(AUC mean=`{summary['auc_abs'].mean():.2f}`)."
        )

    # z-cut: where PSI / Hellinger peak
    if zcut is not None and len(zcut):
        z = zcut.copy()
        if "Hellinger" in z:
            top = z.sort_values("Hellinger", ascending=False).head(5)
            lines.append("- Top Hellinger localizations (z-cut × feature):")
            for _, r in top.iterrows():
                lines.append(
                    f"  - `{r.get('z_cut', '?')}` grp={r.get('grp', '?')} / "
                    f"`{r.get('y_from_x', '?')}`: H={r['Hellinger']:.3f}, "
                    f"PSI={r.get('PSI', float('nan')):.3f}, SMD={r.get('SMD', float('nan')):.3f} "
                    f"(n0={r.get('n0')}, n1={r.get('n1')})"
                )
        # groups with largest |SMD|
        if "SMD" in z:
            top_s = z.assign(abs_smd=z["SMD"].abs()).sort_values("abs_smd", ascending=False).head(3)
            lines.append("- Largest |SMD| pockets:")
            for _, r in top_s.iterrows():
                lines.append(
                    f"  - `{r.get('z_cut')}` grp={r.get('grp')} → `{r.get('y_from_x')}` "
                    f"SMD={r['SMD']:.3f}"
                )

    # feature rank: are true drivers rising?
    if feat is not None and len(feat) and "is_true" in feat.columns:
        true_top = feat[(feat["rank"] <= 5) & (feat["is_true"] == 1)]
        lines.append(
            f"- True drivers appearing in PFI top-5: "
            f"`{len(true_top)}` row-hits across draws "
            f"(unique features: {sorted(true_top['feature'].unique().tolist())})."
        )
    lines.append("")
    return lines


def run_median_design(df, *, max_n, seed, B, n_perm, top_k):
    splits = top_lr_features(df, top_k=3, seed=seed)
    summaries, subs, zcuts, feats = [], [], [], []
    for i, feat in enumerate(splits):
        X, Y, W, names = median_split_batches(
            df, feat, max_n_per_arm=max_n, seed=seed + i
        )
        X0, Y0, X1, Y1 = X[W == 0], Y[W == 0], X[W == 1], Y[W == 1]
        out = run_posthoc_on_batches(
            X0,
            Y0,
            X1,
            Y1,
            names,
            split_feature=feat,
            seed=seed + 100 * i,
            B=B,
            top_k=top_k,
            n_perm=n_perm,
            true_feature=feat,
            z_cols=[feat] + [c for c in names if c != feat][:2],
        )
        summaries.append(out["summary"])
        feats.append(out["feature_rank"])
        subs.append(out["subset_driver"])
        if out["zcut"] is not None and len(out["zcut"]):
            zcuts.append(out["zcut"])
    return (
        pd.DataFrame(summaries),
        pd.concat(feats, ignore_index=True),
        pd.concat(subs, ignore_index=True),
        pd.concat(zcuts, ignore_index=True) if zcuts else pd.DataFrame(),
    )


def run_weak_cd_design(df, *, n_draws, max_n, seed, intensity, B, n_perm, top_k):
    summaries, subs, zcuts, feats = [], [], [], []
    for d in range(n_draws):
        draw_seed = seed + 101 * d
        X, Y, W, names, meta = make_weak_cd_batches(
            df,
            k_drivers=5,
            intensity=intensity,
            max_n_per_arm=max_n,
            seed=draw_seed,
        )
        X0, Y0, X1, Y1 = X[W == 0], Y[W == 0], X[W == 1], Y[W == 1]
        print(
            f"[subset] weak-cd draw={d} drivers={meta.driver_features} ΔȲ={meta.delta_y:.3f}",
            flush=True,
        )
        feat_df, fmeta, _ = feature_localize(
            X0,
            Y0,
            X1,
            Y1,
            seed=draw_seed,
            names=names,
            true_idx=meta.driver_idx,
            top_k=top_k,
            n_perm=n_perm,
        )
        sub_df, contrast, flag = subset_on_driver(
            X0,
            Y0,
            X1,
            Y1,
            B=B,
            j_star=int(fmeta["top_j"]),
            driver_idx=list(fmeta["driver_idx"]),
            seed=draw_seed + 3,
            n_perm=n_perm,
        )
        # z-cuts on PFI top + true drivers
        selected = [fmeta["top_name"]] + [
            r.feature
            for _, r in feat_df.head(top_k).iterrows()
            if r.feature != fmeta["top_name"]
        ]
        selected = list(dict.fromkeys(selected + meta.driver_features))[: top_k + 2]
        z_cols = list(
            dict.fromkeys([fmeta["top_name"]] + meta.driver_features[:2])
        )
        df0 = pd.DataFrame(X0, columns=names)
        df1 = pd.DataFrame(X1, columns=names)
        df0["_Y"] = Y0
        df1["_Y"] = Y1
        loc = subset_localize(
            df0,
            df1,
            x_cols=names,
            y_col="_Y",
            z_cols=z_cols,
            selected=selected,
            top_k=top_k,
        )
        dmse_low = (
            float(sub_df.loc[sub_df.subset == "low", "dmse"].iloc[0])
            if (sub_df.subset == "low").any()
            else np.nan
        )
        dmse_high = (
            float(sub_df.loc[sub_df.subset == "high", "dmse"].iloc[0])
            if (sub_df.subset == "high").any()
            else np.nan
        )
        summaries.append(
            dict(
                draw=d,
                design="weak5_nonlinear_cd",
                drivers=",".join(meta.driver_features),
                pfi_top=fmeta["top_name"],
                hit_at_k=fmeta["hit_at_k"],
                auc_abs=fmeta["auc_abs"],
                contrast=contrast,
                concentrated=flag,
                delta_y=meta.delta_y,
                dmse_low=dmse_low,
                dmse_high=dmse_high,
            )
        )
        feats.append(feat_df.assign(draw=d, design="weak5_nonlinear_cd"))
        subs.append(
            sub_df.assign(
                draw=d,
                design="weak5_nonlinear_cd",
                pfi_top=fmeta["top_name"],
                drivers=",".join(meta.driver_features),
            )
        )
        if loc["table"] is not None and len(loc["table"]):
            zcuts.append(
                loc["table"].assign(
                    draw=d,
                    design="weak5_nonlinear_cd",
                    pfi_top=fmeta["top_name"],
                )
            )
    return (
        pd.DataFrame(summaries),
        pd.concat(feats, ignore_index=True),
        pd.concat(subs, ignore_index=True),
        pd.concat(zcuts, ignore_index=True) if zcuts else pd.DataFrame(),
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--n-draws", type=int, default=5)
    ap.add_argument("--intensity", type=float, default=2.5)
    ap.add_argument("--max-n-per-arm", type=int, default=1500)
    ap.add_argument("--B", type=int, default=30)
    ap.add_argument("--n-perm", type=int, default=25)
    ap.add_argument("--top-k", type=int, default=5)
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    if args.quick:
        n_draws, max_n, B, n_perm = 2, 250, 8, 8
    else:
        n_draws, max_n, B, n_perm = (
            int(args.n_draws),
            int(args.max_n_per_arm),
            int(args.B),
            int(args.n_perm),
        )

    df = clean_heloc(load_heloc(cache_path=CACHE))
    OUT.mkdir(parents=True, exist_ok=True)

    print("subset package: median-split design…", flush=True)
    med_sum, med_feat, med_sub, med_z = run_median_design(
        df,
        max_n=max_n,
        seed=int(args.seed),
        B=B,
        n_perm=n_perm,
        top_k=int(args.top_k),
    )
    print("subset package: weak-5 nonlinear CD…", flush=True)
    cd_sum, cd_feat, cd_sub, cd_z = run_weak_cd_design(
        df,
        n_draws=n_draws,
        max_n=max_n,
        seed=int(args.seed),
        intensity=float(args.intensity),
        B=B,
        n_perm=n_perm,
        top_k=int(args.top_k),
    )

    # write tables
    med_sum.to_csv(OUT / "median_summary.csv", index=False)
    med_feat.to_csv(OUT / "median_feature_rank.csv", index=False)
    med_sub.to_csv(OUT / "median_subset_driver.csv", index=False)
    med_z.to_csv(OUT / "median_zcut_dist.csv", index=False)
    cd_sum.to_csv(OUT / "weakcd_summary.csv", index=False)
    cd_feat.to_csv(OUT / "weakcd_feature_rank.csv", index=False)
    cd_sub.to_csv(OUT / "weakcd_subset_driver.csv", index=False)
    cd_z.to_csv(OUT / "weakcd_zcut_dist.csv", index=False)

    # side-by-side subset_on_driver
    both_sub = pd.concat(
        [
            med_sub.assign(design="median_split_covshift")
            if "design" not in med_sub.columns
            else med_sub,
            cd_sub,
        ],
        ignore_index=True,
        sort=False,
    )
    both_sub.to_csv(OUT / "subset_driver_comparison.csv", index=False)

    lines = [
        "# FICO subset localization insights (auto)",
        "",
        "Package: `feature_localize` → `subset_on_driver` → `subset_localize` "
        "(SMD + PSI / TVD / Hellinger).",
        "",
        f"Settings: max_n_per_arm=`{max_n}`, B_boot=`{B}`, n_perm_pfi=`{n_perm}`, "
        f"weak-cd draws=`{n_draws}`, intensity=`{args.intensity}`.",
        "",
        "## Comparison snapshot",
        "",
        "| design | mean |contrast| | mean ΔMSE low | mean ΔMSE high | mean hit@k |",
        "|---|---:|---:|---:|---:|",
    ]

    def _snap(sum_df):
        if sum_df.empty:
            return (np.nan, np.nan, np.nan, np.nan)
        c = sum_df["contrast"].dropna().abs().mean() if "contrast" in sum_df else np.nan
        lo = sum_df["dmse_low"].dropna().mean() if "dmse_low" in sum_df else np.nan
        hi = sum_df["dmse_high"].dropna().mean() if "dmse_high" in sum_df else np.nan
        hit = sum_df["hit_at_k"].dropna().mean() if "hit_at_k" in sum_df else np.nan
        return c, lo, hi, hit

    for name, sdf in (
        ("median_split_covshift", med_sum),
        ("weak5_nonlinear_cd", cd_sum),
    ):
        c, lo, hi, hit = _snap(sdf)
        lines.append(
            f"| {name} | {c:.4f} | {lo:.4f} | {hi:.4f} | {hit:.2f} |"
        )

    lines += ["", "## Insights", ""]
    lines += _insights_from_tables(
        design="median_split_covshift",
        summary=med_sum,
        subset=med_sub,
        zcut=med_z,
        feat=med_feat,
    )
    lines += _insights_from_tables(
        design="weak5_nonlinear_cd",
        summary=cd_sum,
        subset=cd_sub,
        zcut=cd_z,
        feat=cd_feat,
    )
    lines += [
        "## Cross-design justification",
        "",
        "- **Median-split** localizes a *single strong axis* of covariate shift; "
        "ΔVIMP/ΔMSE often sit on one arm of that axis, and Hellinger peaks along "
        "the split feature's own bins.",
        "- **Weak-5 nonlinear CD** keeps `P(X)` shared, so localization is about "
        "*where in X-space the concept map bites*: `subset_on_driver` finds the "
        "PFI-top cut with larger |ΔVIMP| / ΔMSE, and z-cuts with high Hellinger/"
        "PSI mark the pockets where selected features diverge across batches "
        "under the new labeling.",
        "- Operational takeaway: after CFPerm flags a shift, run this package to "
        "**name the subset (low/high of j★)** that concentrates the damage — "
        "that is the actionable post-hoc localization, not a CS/CD decomposition.",
        "",
        "Tables: `*_summary.csv`, `*_subset_driver.csv`, `*_zcut_dist.csv`, "
        "`subset_driver_comparison.csv`.",
        "",
    ]
    (OUT / "INSIGHTS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (OUT / "INSIGHTS.json").write_text(
        json.dumps(
            {
                "median_summary": med_sum.to_dict(orient="records"),
                "weakcd_summary": cd_sum.to_dict(orient="records"),
            },
            indent=2,
            default=str,
        ),
        encoding="utf-8",
    )
    print("wrote", OUT / "INSIGHTS.md", flush=True)
    print((OUT / "INSIGHTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
