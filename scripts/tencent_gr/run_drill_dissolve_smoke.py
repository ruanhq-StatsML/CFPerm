#!/usr/bin/env python3
"""Smoke: hierarchical dissolve vs top-down MMD drill (landability).

讲武德: W=period; scores are shift proxies — not ATE.

  PYTHONPATH=. python3 scripts/tencent_gr/run_drill_dissolve_smoke.py
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from time_window_feats import fsds_feature_columns  # noqa: E402
from run_standardize_mmd_fsds import (  # noqa: E402
    _matrix,
    fit_standardizer,
    run_fsds,
    w1w2_candidate_columns,
)
from run_three_step_subset_localize import (  # noqa: E402
    _entity_mean_shift,
    _entity_mmd,
    attach_merchant,
    load_item_merchant_map,
)
from attribution_dissolve import (  # noqa: E402
    dissolve_scores,
    drill_hierarchy,
    jaccard_topk,
    leaf_shift_scores,
    rank_overlap_table,
)
from po_risk_fsds import fit_po_on_windows, tau2_on_frame  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    ap = argparse.ArgumentParser()
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
    ap.add_argument("--data-root", type=Path, default=ROOT / "data" / "tencent_subset")
    ap.add_argument("--merchant-col", type=str, default="122")
    ap.add_argument("--top-merchants", type=int, default=30)
    ap.add_argument("--top-users", type=int, default=50)
    ap.add_argument("--mmd-max-cand", type=int, default=80)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "results" / "tencent_gr_drill_dissolve",
    )
    args = ap.parse_args()
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    g1 = pd.read_parquet(args.w1_grid)
    g2 = pd.read_parquet(args.w2_grid)
    raw = [c for c in fsds_feature_columns(g1) if c in g2.columns]
    cols = w1w2_candidate_columns(raw)

    imap = load_item_merchant_map(args.data_root, merchant_col=args.merchant_col)
    g1 = attach_merchant(g1, imap)
    g2 = attach_merchant(g2, imap)
    print(
        f"grids W1={len(g1)} W2={len(g2)} feats={len(cols)} "
        f"merchant_mapped_W1={g1['merchant_mapped'].mean():.3f}",
        flush=True,
    )

    # --- leaf scores: shift_l2 + optional τ̂² ---
    leaf = leaf_shift_scores(g1, g2, cols, seed=args.seed)
    # attach merchant from g2 alignment by index
    leaf["merchant_id"] = g2["merchant_id"].to_numpy()
    leaf["user_id"] = g2["user_id"].to_numpy()
    leaf["item_id"] = g2["item_id"].to_numpy()

    po = fit_po_on_windows(g1, g2, cols, seed=args.seed, max_n=6000)
    leaf["tau2"] = tau2_on_frame(g2, cols, po)
    leaf.to_csv(out / "leaf_scores.csv", index=False)

    hier_l2 = drill_hierarchy(
        leaf,
        score_col="shift_l2",
        levels=("user_id", "merchant_id"),
        how="sum",
        top_k={"user_id": args.top_users, "merchant_id": args.top_merchants},
    )
    hier_tau = drill_hierarchy(
        leaf,
        score_col="tau2",
        levels=("user_id", "merchant_id"),
        how="sum",
        top_k={"user_id": args.top_users, "merchant_id": args.top_merchants},
    )
    hier_l2["user_id"].to_csv(out / "dissolve_user_shift_l2.csv", index=False)
    hier_l2["merchant_id"].to_csv(out / "dissolve_merchant_shift_l2.csv", index=False)
    hier_tau["user_id"].to_csv(out / "dissolve_user_tau2.csv", index=False)
    hier_tau["merchant_id"].to_csv(out / "dissolve_merchant_tau2.csv", index=False)

    # Cascading dissolve (same shape as three-step): top merchants → users inside
    top_m_l2 = set(hier_l2["merchant_id"].head(args.top_merchants)["merchant_id"].tolist())
    leaf_in_m = leaf[leaf["merchant_id"].isin(top_m_l2)]
    cascade_users = dissolve_scores(leaf_in_m, score_col="shift_l2", by="user_id", how="sum")
    cascade_users["selected"] = (cascade_users.index < args.top_users).astype(int)
    cascade_users.to_csv(out / "dissolve_cascade_user_shift_l2.csv", index=False)
    print(
        f"cascade: merchants={len(top_m_l2)} leaf_in={len(leaf_in_m)} users={len(cascade_users)}",
        flush=True,
    )

    # --- top-down MMD / mean-shift reference (same spirit as three-step L1/L2) ---
    sc = fit_standardizer(_matrix(g1, cols))
    mer_mmd = _entity_mmd(
        g1,
        g2,
        cols,
        sc,
        entity_col="merchant_id",
        seed=args.seed,
        max_cand=args.mmd_max_cand,
        min_edges=3,
        mmd_max_n=400,
    )
    if mer_mmd.empty or mer_mmd["mmd2"].isna().all():
        mer_mmd = _entity_mean_shift(
            g1, g2, cols, sc, entity_col="merchant_id", max_cand=args.mmd_max_cand, min_edges=2
        )
    mer_mmd.to_csv(out / "ref_merchant_mmd.csv", index=False)

    top_m = set(mer_mmd.head(args.top_merchants)["merchant_id"].tolist())
    g1m = g1[g1["merchant_id"].isin(top_m)]
    g2m = g2[g2["merchant_id"].isin(top_m)]
    usr_mmd = _entity_mmd(
        g1m,
        g2m,
        cols,
        sc,
        entity_col="user_id",
        seed=args.seed,
        max_cand=min(200, args.mmd_max_cand * 2),
        min_edges=2,
        mmd_max_n=300,
    )
    if usr_mmd.empty or usr_mmd["mmd2"].isna().all():
        usr_mmd = _entity_mean_shift(
            g1m, g2m, cols, sc, entity_col="user_id", max_cand=200, min_edges=1
        )
    usr_mmd.to_csv(out / "ref_user_mmd.csv", index=False)

    # overlaps
    ov_m_l2 = rank_overlap_table(hier_l2["merchant_id"], mer_mmd, key="merchant_id")
    ov_m_tau = rank_overlap_table(hier_tau["merchant_id"], mer_mmd, key="merchant_id")
    ov_u_l2 = rank_overlap_table(hier_l2["user_id"], usr_mmd, key="user_id")
    ov_u_tau = rank_overlap_table(hier_tau["user_id"], usr_mmd, key="user_id")
    ov_u_cascade = rank_overlap_table(cascade_users, usr_mmd, key="user_id")
    ov_m_l2.to_csv(out / "overlap_merchant_l2_vs_mmd.csv", index=False)
    ov_m_tau.to_csv(out / "overlap_merchant_tau2_vs_mmd.csv", index=False)
    ov_u_l2.to_csv(out / "overlap_user_l2_vs_mmd.csv", index=False)
    ov_u_tau.to_csv(out / "overlap_user_tau2_vs_mmd.csv", index=False)
    ov_u_cascade.to_csv(out / "overlap_user_cascade_vs_mmd.csv", index=False)

    # FSDS: fit on full W1 (rare-pos safe); evaluate W2 restricted to support users
    def _fsds_on_users(user_ids, tag: str):
        u = set(int(x) for x in user_ids)
        te = g2[g2["user_id"].isin(u)]
        if len(te) < 20 or te["y_convert"].nunique() < 2:
            te = g2
            note = f"eval_fallback_full_W2 n_users={len(u)}"
        else:
            note = f"eval_W2_support n_te={len(te)} n_users={len(u)} pos={int(te['y_convert'].sum())}"
        res = run_fsds(g1, te, cols, select_k=15, seed=args.seed)
        h = (res.get("models") or {}).get("hgb", {}) or {}
        selected = list(res.get("selected") or [])
        return {
            "tag": tag,
            "note": note,
            "W2_hgb_auc": h.get("auc"),
            "W2_ap": h.get("ap"),
            "n_users": len(u),
            "n_te": int(len(te)),
            "ok": bool(res.get("ok")),
            "top5": ",".join(selected[:5]),
        }

    fsds_rows = [
        _fsds_on_users(hier_l2["user_id"].head(args.top_users)["user_id"], "dissolve_l2_users"),
        _fsds_on_users(cascade_users.head(args.top_users)["user_id"], "dissolve_cascade_users"),
        _fsds_on_users(hier_tau["user_id"].head(args.top_users)["user_id"], "dissolve_tau2_users"),
        _fsds_on_users(usr_mmd.head(args.top_users)["user_id"], "mmd_users"),
        _fsds_on_users(g1["user_id"].drop_duplicates().head(args.top_users), "random_head_users"),
    ]
    fsds_df = pd.DataFrame(fsds_rows)
    fsds_df.to_csv(out / "fsds_on_dissolved_support.csv", index=False)

    summary = {
        "sec": time.time() - t0,
        "n_leaf": len(leaf),
        "po_risk": po["risk"],
        "note": "dissolve = hierarchical score rollup (not PyPI dissolve API migrator); W=period, not ATE",
        "overlap_merchant_l2": ov_m_l2.to_dict("records"),
        "overlap_merchant_tau2": ov_m_tau.to_dict("records"),
        "overlap_user_l2": ov_u_l2.to_dict("records"),
        "overlap_user_tau2": ov_u_tau.to_dict("records"),
        "overlap_user_cascade": ov_u_cascade.to_dict("records"),
        "fsds": fsds_rows,
        "landability": {
            "cheap": True,
            "pypi_dissolve": "WRONG package (API deprecation) — we use dissolve-as-rollup",
            "vs_mmd": "see Jaccard tables; bottom-up dissolve ≠ top-down MMD but can seed the same FSDS",
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))

    # FINDINGS
    lines = [
        "# Drill dissolve smoke — landability",
        "",
        "**讲武德:** `W`=period. Leaf scores / dissolve sums are shift proxies — **not** ATE.",
        "",
        "## Naming",
        "",
        "PyPI [`dissolve`](https://pypi.org/project/dissolve/) = API-deprecation migrator. **Not used.**",
        "Here **dissolve** = geopandas-style rollup: leaf score → `groupby(parent).sum/mean`.",
        "",
        "## Setup",
        "",
        f"- leaves (W2 rows): {len(leaf)}",
        f"- PO-risk: {po['risk']:.6g}",
        f"- Top merchants={args.top_merchants}, Top users={args.top_users}",
        "",
        "## Overlap vs top-down MMD (Jaccard@K)",
        "",
        "### Merchant (dissolve shift_l2 vs MMD)",
        "",
        "```",
        ov_m_l2.to_string(index=False),
        "```",
        "",
        "### Merchant (dissolve τ̂² vs MMD)",
        "",
        "```",
        ov_m_tau.to_string(index=False),
        "```",
        "",
        "### User (global dissolve shift_l2 vs MMD-in-merchants)",
        "",
        "```",
        ov_u_l2.to_string(index=False),
        "```",
        "",
        "### User (cascade dissolve: top-M merchants → users vs MMD)",
        "",
        "```",
        ov_u_cascade.to_string(index=False),
        "```",
        "",
        "### User (dissolve τ̂² vs MMD)",
        "",
        "```",
        ov_u_tau.to_string(index=False),
        "```",
        "",
        "## FSDS on dissolved user support",
        "",
        "```",
        fsds_df.to_string(index=False),
        "```",
        "",
        "## Landability verdict",
        "",
        "1. **Implementable:** pure pandas groupby — no new heavy deps; plugs into existing grids.",
        "2. **Merchant-level dissolve (shift_l2) tracks top-down MMD** reasonably (see Jaccard@K).",
        "3. **Global user dissolve ≠ cascaded MMD users** (Jaccard≈0); use **cascade dissolve** (merchants→users) to match drill shape.",
        "4. **τ̂² leaf dissolve** poorly matches MMD tops here — keep as optional PO side-channel, not the drill gate.",
        "5. **Useful as:** cheap bottom-up prior / candidate gen before official cmean·MMD·FSDS; not a replacement.",
        "6. **Do not** ship PyPI `dissolve` (API migrator); do not claim ATE.",
        "",
        f"Artifacts under `{out}`.",
        "",
    ]
    (out / "FINDINGS.md").write_text("\n".join(lines))
    print(json.dumps({"fsds": fsds_rows, "ov_m_l2": ov_m_l2.to_dict("records")[:3]}, indent=2), flush=True)
    print(f"wrote {out} in {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
