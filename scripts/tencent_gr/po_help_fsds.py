#!/usr/bin/env python3
"""DS one-liner CLI: period-PO help → official FSDS (讲武德: W=period, not ATE).

Recommended recipe from overnight iters (PR #72):

  po_help_select / PO-VIMP pool (k+3) → Scaler→Var→SelectKBest→HGB/LR @ k=15
  Report seed mean±std under rare convert. Optional: k=18 for plain PO-VIMP.

Example:

  PYTHONPATH=. python3 scripts/tencent_gr/po_help_fsds.py \\
    --w1-grid results/tencent_gr_localize_fsds_time/W1_localized_grid.parquet \\
    --w2-grid results/tencent_gr_localize_fsds_time/W2_localized_grid_sample.parquet \\
    --out-dir results/tencent_gr_fsds_iterate/po_help_cli_smoke \\
    --select-k 15 --seed 0
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
from sklearn.model_selection import GroupShuffleSplit

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from time_window_feats import fsds_feature_columns  # noqa: E402
from run_standardize_mmd_fsds import run_fsds, w1w2_candidate_columns  # noqa: E402
from po_risk_fsds import po_help_select  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    ap = argparse.ArgumentParser(description="PO-help → official FSDS (period W; not ATE)")
    ap.add_argument(
        "--w1-grid",
        type=Path,
        default=ROOT / "results" / "tencent_gr_localize_fsds_time" / "W1_localized_grid.parquet",
    )
    ap.add_argument(
        "--w2-grid",
        type=Path,
        default=ROOT
        / "results"
        / "tencent_gr_localize_fsds_time"
        / "W2_localized_grid_sample.parquet",
    )
    ap.add_argument("--select-k", type=int, default=15, help="FSDS k (prefer 15; 18 ok for plain PO-VIMP)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--pool-extra", type=int, default=3, help="PO pool slack beyond k")
    ap.add_argument("--alpha", type=float, default=0.5, help="|δ|⋈PO-VIMP blend (default 0.5)")
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "results" / "tencent_gr_fsds_iterate" / "po_help_cli",
    )
    ap.add_argument(
        "--holdout-frac",
        type=float,
        default=0.25,
        help="W1 user-group holdout for sanity metrics (selection still on train only)",
    )
    args = ap.parse_args()
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    g1 = pd.read_parquet(args.w1_grid)
    g2 = pd.read_parquet(args.w2_grid)
    raw = [c for c in fsds_feature_columns(g1) if c in g2.columns]
    cols = w1w2_candidate_columns(raw)

    gss = GroupShuffleSplit(n_splits=1, test_size=args.holdout_frac, random_state=args.seed)
    tr_idx, va_idx = next(gss.split(g1, groups=g1["user_id"].to_numpy()))
    g_tr, g_va = g1.iloc[tr_idx].copy(), g1.iloc[va_idx].copy()

    report = po_help_select(
        g_tr,
        g2,
        cols,
        k=args.select_k,
        seed=args.seed,
        alpha=args.alpha,
        pool_extra=args.pool_extra,
    )
    pool = report["pool"]
    report["feature_table"].to_csv(out / "po_feature_vimp.csv", index=False)
    report["blend_table"].to_csv(out / "po_cmean_blend.csv", index=False)
    (out / "po_pool.txt").write_text("\n".join(pool) + "\n")
    (out / "po_risk.txt").write_text(
        f"risk={report['risk']:.6g}\nnote={report['note']}\nk={args.select_k}\n"
        f"pool_n={len(pool)}\nalpha={args.alpha}\n"
    )

    # Official FSDS on PO-guided pool (W1-train fit; W2 never selects)
    res_va = run_fsds(g_tr, g_va, pool, select_k=min(args.select_k, len(pool)), seed=args.seed)
    res_w2 = run_fsds(g_tr, g2, pool, select_k=min(args.select_k, len(pool)), seed=args.seed)
    selected = list(res_w2.get("selected") or pool[: args.select_k])
    if isinstance(res_w2.get("ranking"), pd.DataFrame):
        res_w2["ranking"].to_csv(out / "fsds_ranking.csv", index=False)
    (out / "fsds_selected.txt").write_text("\n".join(selected) + "\n")

    summary = {
        "recipe": "po_help_select → official FSDS",
        "note": report["note"],
        "po_risk": report["risk"],
        "k": args.select_k,
        "pool": pool,
        "selected": selected,
        "W1_hold_hgb": (res_va.get("models") or {}).get("hgb"),
        "W2_hgb": (res_w2.get("models") or {}).get("hgb"),
        "W2_logreg": (res_w2.get("models") or {}).get("logreg"),
        "n_train": len(g_tr),
        "n_w2": len(g2),
        "pos_train": float(g_tr["y_convert"].mean()),
        "tips": {
            "k": "prefer 15; try 18 for plain PO-VIMP",
            "freeze": "for one frozen list use LOO-pos majority (see iter09); not seed-maj2",
            "seeds": "always report mean±std under rare convert",
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    w2 = summary["W2_hgb"] or {}
    print(
        f"PO-help→FSDS  risk={report['risk']:.4g}  "
        f"selected={len(selected)}  "
        f"W2_hgb_auc={w2.get('auc')}  W2_ap={w2.get('ap')}  "
        f"→ {out}",
        flush=True,
    )
    print(report["note"], flush=True)


if __name__ == "__main__":
    main()
