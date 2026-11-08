#!/usr/bin/env python3
"""Post-hoc subset localization tables on FICO HELOC median-split batches.

  python3 scripts/run_fico_posthoc_localization.py
  python3 scripts/run_fico_posthoc_localization.py --quick
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from fico_heloc_cfperm import (  # noqa: E402
    clean_heloc,
    load_heloc,
    median_split_batches,
    top_lr_features,
)
from posthoc_subset_localization import run_posthoc_on_batches  # noqa: E402

OUT = ROOT / "results" / "fico_heloc_cfperm"
CACHE = ROOT / "data" / "heloc_dataset_v1.csv"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--max-n-per-arm", type=int, default=1500)
    ap.add_argument("--B", type=int, default=30, help="bootstrap reps for subset_on_driver")
    ap.add_argument("--n-perm", type=int, default=25, help="PFI inner permutations")
    ap.add_argument("--top-k", type=int, default=5)
    args = ap.parse_args()

    max_n = 250 if args.quick else int(args.max_n_per_arm)
    B = 8 if args.quick else int(args.B)
    n_perm = 8 if args.quick else int(args.n_perm)

    df = clean_heloc(load_heloc(cache_path=CACHE))
    split_feats = top_lr_features(df, top_k=3, seed=int(args.seed))
    print("posthoc splits:", split_feats, "max_n=", max_n, "B=", B, flush=True)

    feat_rows = []
    sub_rows = []
    zcut_rows = []
    summaries = []

    for i, feat in enumerate(split_feats):
        X, Y, W, names = median_split_batches(
            df, feat, max_n_per_arm=max_n, seed=int(args.seed) + i
        )
        X0, Y0 = X[W == 0], Y[W == 0]
        X1, Y1 = X[W == 1], Y[W == 1]
        print(f"[posthoc] split={feat} n0={len(Y0)} n1={len(Y1)}", flush=True)
        out = run_posthoc_on_batches(
            X0,
            Y0,
            X1,
            Y1,
            names,
            split_feature=feat,
            seed=int(args.seed) + 100 * i,
            B=B,
            top_k=int(args.top_k),
            n_perm=n_perm,
            true_feature=feat,
        )
        feat_rows.append(out["feature_rank"])
        sub_rows.append(out["subset_driver"])
        if out["zcut"] is not None and len(out["zcut"]):
            zcut_rows.append(out["zcut"])
        summaries.append(out["summary"])
        print(
            f"  driver={out['summary']['driver']} "
            f"concentrated={out['summary']['concentrated']} "
            f"contrast={out['summary']['contrast']:.4f}",
            flush=True,
        )

    OUT.mkdir(parents=True, exist_ok=True)
    feat_df = pd.concat(feat_rows, ignore_index=True)
    sub_df = pd.concat(sub_rows, ignore_index=True)
    zcut_df = pd.concat(zcut_rows, ignore_index=True) if zcut_rows else pd.DataFrame()
    sum_df = pd.DataFrame(summaries)

    feat_df.to_csv(OUT / "posthoc_feature_rank.csv", index=False)
    sub_df.to_csv(OUT / "posthoc_subset_driver.csv", index=False)
    zcut_df.to_csv(OUT / "posthoc_zcut_dist.csv", index=False)
    sum_df.to_csv(OUT / "posthoc_summary.csv", index=False)
    with open(OUT / "posthoc_summary.json", "w", encoding="utf-8") as f:
        json.dump(summaries, f, indent=2)

    # Markdown appendix
    lines = [
        "# Post-hoc subset localization (FICO HELOC)",
        "",
        "After CFPerm / median-split batches: PFI degree → driver `j★`, then "
        "`subset_on_driver` (low/high median cut) with bootstrap ΔVIMP CI + ΔMSE; "
        "z-cuts report SMD + PSI / TVD / Hellinger.",
        "",
        "## Summary",
        "",
        "| split | driver | concentrated | contrast | ΔMSE low | ΔMSE high | hit@k |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for s in summaries:
        lines.append(
            f"| `{s['split_feature']}` | `{s['driver']}` | {s['concentrated']} | "
            f"{s['contrast']:.4f} | {s['dmse_low']:.4f} | {s['dmse_high']:.4f} | "
            f"{s['hit_at_k'] if s['hit_at_k'] == s['hit_at_k'] else float('nan'):.2f} |"
        )
    lines += [
        "",
        "## Files",
        "",
        "- `posthoc_feature_rank.csv` — |ΔPFI| ranking per split",
        "- `posthoc_subset_driver.csv` — low/high ΔVIMP (±1.96 SE) and ΔMSE",
        "- `posthoc_zcut_dist.csv` — group-wise SMD + PSI/TVD/Hellinger",
        "- `posthoc_summary.csv` / `.json`",
        "",
    ]
    (OUT / "POSTHOC.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
