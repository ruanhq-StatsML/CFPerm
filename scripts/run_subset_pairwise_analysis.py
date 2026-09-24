#!/usr/bin/env python3
"""Subset pairwise analysis on intensity baseline + credit (Tencent edges).

Given: edge grain + pairwise chunks + ops intensity as buy-volume baseline.
Subset-analysis becomes: cut the *same* sorted edge stream into cohorts, then
re-run adjacent pairwise transfer inside each cohort.

Cohorts (judgment):
  - intensity_high / intensity_low : tertiles of **item** buy-volume
        ``i_log1p_n_exp`` (edge ``e_n_exp`` is mostly 1 — too flat to slice)
  - credit_high / credit_low       : top/bottom by ``i_item_credit_rank``
        (1=highest share_linear; more separable than raw credit_last zeros)
  - all                            : full stream (reference)

Also dumps a **pack pairwise scorecard** (from board summary) and a
**calendar-vs-count** note (equal-N spans vs DiffusionDB calendar-width board).

  PYTHONPATH=. python3 scripts/run_subset_pairwise_analysis.py \\
    --chunk-size 1000 --out results/subset_pairwise_analysis
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from scripts.run_sample_chunk_adjacent_board import (
    TENCENT_CONVERT_LEAK,
    TENCENT_VOLUME,
    chunk_by_n,
    fs_adjacent,
    interpret_pack,
    mean_jaccard,
    mean_key,
)

ROOT = Path(__file__).resolve().parents[1]
GRID = ROOT / "results/tencent_gr_tabular_ui/feature_grid.parquet"
BOARD_SUMMARY = ROOT / "results/sample_chunk_adjacent_board/summary.json"
CAL_SUMMARY = ROOT / "results/diffusiondb_adjacent_board/summary.json"


def _tencent_xy(
    df: pd.DataFrame, *, panel: str
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    y_col = "y_convert"
    drop = {
        "user_id",
        "item_id",
        y_col,
        "e_last_ts",
        "last_ts",
        *TENCENT_CONVERT_LEAK,
    }
    if panel == "content":
        drop |= set(TENCENT_VOLUME)
    feats = [
        c
        for c in df.columns
        if c not in drop and pd.api.types.is_numeric_dtype(df[c])
    ]
    X = np.nan_to_num(df[feats].to_numpy(float), nan=0.0, posinf=0.0, neginf=0.0)
    y = df[y_col].to_numpy(float)
    return X, y, feats


def _tertile_masks(s: pd.Series) -> Dict[str, np.ndarray]:
    q1, q2 = s.quantile([1 / 3, 2 / 3])
    v = s.to_numpy(dtype=float)
    # if degenerate (almost constant), caller should pick another column
    return {
        "low": v <= float(q1),
        "mid": (v > float(q1)) & (v <= float(q2)),
        "high": v > float(q2),
    }


def _rank_tail_masks(rank: pd.Series, *, frac: float = 0.25) -> Dict[str, np.ndarray]:
    """Low rank number = higher credit share in our tables."""
    v = rank.to_numpy(dtype=float)
    lo, hi = np.nanquantile(v, [frac, 1.0 - frac])
    return {
        "high": v <= float(lo),  # best ranks
        "low": v >= float(hi),  # worst ranks
    }


def run_subset(
    df_sorted: pd.DataFrame,
    mask: np.ndarray,
    *,
    name: str,
    panel: str,
    chunk_size: int,
    select_k: int,
    seed: int,
    max_pairs: int,
    stride: Optional[int] = None,
) -> Dict[str, Any]:
    sub = df_sorted.loc[mask].reset_index(drop=True)
    n = len(sub)
    stride_eff = int(stride) if stride and stride > 0 else int(chunk_size)
    need = chunk_size + (stride_eff if stride_eff < chunk_size else chunk_size)
    if n < need:
        return {
            "ok": False,
            "name": name,
            "panel": panel,
            "n": n,
            "reason": f"too_small(need>={need})",
        }
    X, y, names = _tencent_xy(sub, panel=panel)
    if stride_eff >= chunk_size:
        # non-overlap: classic chunk ids
        T = chunk_by_n(n, chunk_size)
        spans = []
        for cid in sorted(np.unique(T)):
            sl = sub.iloc[T == cid]
            spans.append(
                float((sl["e_last_ts"].max() - sl["e_last_ts"].min()) / 3600.0)
            )
        rows = fs_adjacent(
            X,
            y,
            T,
            names,
            select_k=select_k,
            y_quantile=0.7,
            seed=seed,
            max_pairs=max_pairs,
            min_n=min(80, chunk_size // 2),
        )
        mode = "nonoverlap"
    else:
        # overlapping windows of length N, step=stride; pairwise consecutive starts
        from itertools import pairwise as _pw

        starts = list(range(0, n - chunk_size + 1, stride_eff))
        if len(starts) < 2:
            return {
                "ok": False,
                "name": name,
                "panel": panel,
                "n": n,
                "reason": "stride_too_few_windows",
            }
        pairs = list(_pw(starts))
        if max_pairs > 0 and len(pairs) > max_pairs:
            idx = np.linspace(0, len(pairs) - 1, num=max_pairs, dtype=int)
            pairs = [pairs[i] for i in idx]
        # map to pseudo chunk ids 0..W-1 on a synthetic T for reporting spans only
        # run FS manually per pair
        rows = []
        spans = []
        for s0, s1 in pairs:
            m0 = np.zeros(n, dtype=bool)
            m1 = np.zeros(n, dtype=bool)
            m0[s0 : s0 + chunk_size] = True
            m1[s1 : s1 + chunk_size] = True
            T = np.full(n, -1, dtype=int)
            T[m0] = 0
            T[m1] = 1
            # only keep the two windows' rows for a tiny adjacent call
            keep = m0 | m1
            # Remap to contiguous for fs_adjacent with ids 0,1
            Tk = np.where(m0[keep], 0, 1)
            pair_rows = fs_adjacent(
                X[keep],
                y[keep],
                Tk,
                names,
                select_k=select_k,
                y_quantile=0.7,
                seed=seed,
                max_pairs=1,
                min_n=min(80, chunk_size // 2),
            )
            rows.extend(pair_rows)
            sl0 = sub.iloc[s0 : s0 + chunk_size]
            spans.append(
                float((sl0["e_last_ts"].max() - sl0["e_last_ts"].min()) / 3600.0)
            )
        mode = f"stride_{stride_eff}"

    m_auc = mean_key(rows, "hgb_auc")
    m_lr = mean_key(rows, "logreg_auc")
    m_dy = mean_key(rows, "delta_Y")
    m_j = mean_jaccard(rows)
    return {
        "ok": True,
        "name": name,
        "panel": panel,
        "n": n,
        "window_mode": mode,
        "n_pairs": len(rows),
        "mean_auc": m_auc,
        "mean_logreg_auc": m_lr,
        "mean_delta_Y": m_dy,
        "mean_jaccard": m_j,
        "mean_span_hours": float(np.mean(spans)) if spans else None,
        "span_hours_by_chunk": spans[:8],
        "reading": interpret_pack(m_auc, m_j),
        "top0": (rows[0]["top_fsds"][:4] if rows else []),
        "y_rate": float(y.mean()),
    }


def pack_pairwise_scorecard(board: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Flatten cross-pack sample-count pairwise metrics."""
    out = []
    for row in board.get("cross_pack") or []:
        out.append(
            {
                "pack": row.get("dataset"),
                "grain": "equal_count_edges_or_rows",
                "N": row.get("chunk"),
                "mean_auc": row.get("mean_auc"),
                "mean_logreg_auc": row.get("mean_logreg_auc"),
                "mean_jaccard": row.get("mean_jaccard"),
                "mean_delta_Y": row.get("mean_delta_Y"),
                "mean_fsds_cmean_jaccard": row.get("mean_fsds_cmean_jaccard"),
                "reading": row.get("reading"),
            }
        )
    return out


def calendar_width_scorecard(cal: Optional[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """DiffusionDB calendar-width adjacent board (different grain)."""
    if not cal:
        return []
    out = []
    for b in cal.get("boards") or []:
        rows = b.get("rows") or []
        aucs = [r["hgb_auc"] for r in rows]
        dys = [r["delta_Y"] for r in rows]
        out.append(
            {
                "pack": "diffusiondb",
                "grain": "calendar_width_min",
                "width_min": b.get("width_min"),
                "n_bins": b.get("n_bins"),
                "n_reported": b.get("n_reported"),
                "mean_auc": float(np.mean(aucs)) if aucs else None,
                "mean_delta_Y": float(np.mean(dys)) if dys else None,
                "note": (
                    "calendar bins; n per bin unbalanced — contrast with equal-count "
                    "sample-chunk board on same pack"
                ),
            }
        )
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--chunk-size", type=int, default=1000)
    ap.add_argument("--n-cap", type=int, default=12000, help="cap edges after sort")
    ap.add_argument("--select-k", type=int, default=15)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-pairs", type=int, default=8)
    ap.add_argument(
        "--cross",
        action="store_true",
        help="also run intensity×credit 2×2 cross cohorts",
    )
    ap.add_argument(
        "--stride",
        type=int,
        default=0,
        help="if 0 < stride < chunk-size, use overlapping windows (robustness)",
    )
    ap.add_argument(
        "--by-advertiser",
        type=int,
        default=0,
        help="if >0, also slice top-K merchants (item_feat.122) as advertiser packs",
    )
    ap.add_argument(
        "--out", type=Path, default=ROOT / "results/subset_pairwise_analysis"
    )
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    df = pd.read_parquet(GRID).sort_values("e_last_ts").reset_index(drop=True)
    if args.n_cap < len(df):
        idx = np.linspace(0, len(df) - 1, num=args.n_cap, dtype=int)
        df = df.iloc[idx].reset_index(drop=True)

    # Buy-volume: prefer item-level exp (edge e_n_exp ≈ 1 for most rows)
    inten_col = (
        "i_log1p_n_exp"
        if "i_log1p_n_exp" in df.columns
        else ("i_n_exp" if "i_n_exp" in df.columns else "e_log1p_exp")
    )
    inten = _tertile_masks(df[inten_col].astype(float))
    # Last-touch: item_credit_rank (1=best); fallback credit_last tertile among >0
    if "i_item_credit_rank" in df.columns:
        credit = _rank_tail_masks(df["i_item_credit_rank"].astype(float), frac=0.25)
        cred_col = "i_item_credit_rank"
    else:
        cred_col = "i_credit_last" if "i_credit_last" in df.columns else "i_share_last"
        credit = _tertile_masks(df[cred_col].astype(float))

    cohorts = [
        ("all", np.ones(len(df), dtype=bool)),
        ("intensity_high", inten["high"]),
        ("intensity_low", inten["low"]),
        ("credit_high", credit["high"]),
        ("credit_low", credit["low"]),
    ]
    if args.cross:
        cohorts += [
            ("I_high_C_high", inten["high"] & credit["high"]),
            ("I_high_C_low", inten["high"] & credit["low"]),
            ("I_low_C_high", inten["low"] & credit["high"]),
            ("I_low_C_low", inten["low"] & credit["low"]),
        ]

    # Optional advertiser / account-like slice
    advertiser_note = None
    if args.by_advertiser > 0:
        try:
            import sys

            sys.path.insert(0, str(ROOT / "scripts/tencent_gr"))
            from run_three_step_subset_localize import (  # type: ignore
                attach_merchant,
                load_item_merchant_map,
            )

            imap = load_item_merchant_map(ROOT / "data/tencent_subset", merchant_col="122")
            df = attach_merchant(df, imap)
            vc = df["merchant_id"].value_counts()
            top = vc.head(int(args.by_advertiser)).index.tolist()
            advertiser_note = {
                "key": "merchant_id",
                "top": [str(x) for x in top],
            }
            for mid in top[: min(3, len(top))]:
                cohorts.append((f"adv_{mid}", (df["merchant_id"] == mid).to_numpy()))
        except Exception as e:
            # Fallback: top items as creative/account proxy when merchant parquet missing
            vc = df["item_id"].value_counts()
            top = vc.head(int(args.by_advertiser)).index.tolist()
            advertiser_note = {
                "ok": False,
                "merchant_error": str(e),
                "fallback_key": "item_id",
                "note": (
                    "item_feat merchant map unavailable; using top item_id by edge "
                    "count as creative/account proxy (not true advertiser)."
                ),
                "top": [str(x) for x in top],
            }
            for iid in top[: min(3, len(top))]:
                cohorts.append((f"item_{iid}", (df["item_id"] == iid).to_numpy()))

    stride = args.stride if args.stride > 0 else None
    subset_rows: List[Dict[str, Any]] = []
    for name, mask in cohorts:
        for panel in ("ops", "content"):
            print(
                f"subset {name} panel={panel} n={int(mask.sum())} stride={stride}",
                flush=True,
            )
            subset_rows.append(
                run_subset(
                    df,
                    mask,
                    name=name,
                    panel=panel,
                    chunk_size=args.chunk_size,
                    select_k=args.select_k,
                    seed=args.seed,
                    max_pairs=args.max_pairs,
                    stride=stride,
                )
            )

    board = (
        json.loads(BOARD_SUMMARY.read_text()) if BOARD_SUMMARY.is_file() else {}
    )
    cal = json.loads(CAL_SUMMARY.read_text()) if CAL_SUMMARY.is_file() else None
    pack_card = pack_pairwise_scorecard(board)
    cal_card = calendar_width_scorecard(cal)

    # equal-count → calendar side-effect on full tencent stream
    T = chunk_by_n(len(df), args.chunk_size)
    eq_spans = []
    for cid in sorted(np.unique(T))[:8]:
        sl = df.iloc[T == cid]
        eq_spans.append(
            {
                "chunk": int(cid),
                "n": int((T == cid).sum()),
                "span_hours": float(
                    (sl["e_last_ts"].max() - sl["e_last_ts"].min()) / 3600.0
                ),
                "y_rate": float(sl["y_convert"].mean()),
            }
        )

    report = {
        "note": (
            "Subset pairwise = same edge sort + pairwise(N), but only edges in cohort. "
            "Intensity tertiles = buy-volume baseline slices; credit tertiles = "
            "last-touch slices. Pack scorecard = sample-count pairwise; calendar "
            "scorecard = DiffDB width-min board (different grain)."
        ),
        "chunk_size": args.chunk_size,
        "n_edges_used": len(df),
        "intensity_col": inten_col,
        "credit_col": cred_col,
        "stride": stride,
        "cross": bool(args.cross),
        "advertiser": advertiser_note,
        "subsets": subset_rows,
        "pack_pairwise_scorecard": pack_card,
        "calendar_width_scorecard": cal_card,
        "equal_count_calendar_side_effect": eq_spans,
        "sec": float(time.time() - t0),
    }
    (args.out / "summary.json").write_text(
        json.dumps(report, indent=2, default=str) + "\n"
    )

    md = [
        "# Subset pairwise analysis",
        "",
        report["note"],
        "",
        "## Tencent subsets (edge pairwise)",
        "",
        "| subset | panel | n | pairs | HGB | LogReg | ΔȲ | Jac | mean_span_h | reading | top0 |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    for r in subset_rows:
        if not r.get("ok"):
            md.append(
                f"| {r['name']} | {r['panel']} | {r.get('n')} | — | — | — | — | — | — | {r.get('reason')} | — |"
            )
            continue
        md.append(
            f"| {r['name']} | {r['panel']} | {r['n']} | {r['n_pairs']} | "
            f"{(r['mean_auc'] or float('nan')):.3f} | "
            f"{(r['mean_logreg_auc'] or float('nan')):.3f} | "
            f"{(r['mean_delta_Y'] or 0):.4g} | "
            f"{(r['mean_jaccard'] if r['mean_jaccard'] is not None else float('nan')):.2f} | "
            f"{(r['mean_span_hours'] or float('nan')):.1f} | "
            f"{r['reading']} | {', '.join(r['top0'][:3])} |"
        )

    md += [
        "",
        "## Pack pairwise scorecard (sample-count grain)",
        "",
        "| pack | N | HGB | LogReg | Jac | ΔȲ | reading |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for r in pack_card:
        md.append(
            f"| {r['pack']} | {r['N']} | "
            f"{(r['mean_auc'] or float('nan')):.3f} | "
            f"{(r['mean_logreg_auc'] or float('nan')):.3f} | "
            f"{(r['mean_jaccard'] if r['mean_jaccard'] is not None else float('nan')):.2f} | "
            f"{abs(r['mean_delta_Y'] or 0):.4g} | {r['reading']} |"
        )

    md += [
        "",
        "## Calendar-width scorecard (DiffDB — other grain)",
        "",
        "| width_min | n_bins | n_reported | mean_auc | mean_ΔȲ |",
        "|---:|---:|---:|---:|---:|",
    ]
    for r in cal_card:
        md.append(
            f"| {r['width_min']} | {r['n_bins']} | {r['n_reported']} | "
            f"{(r['mean_auc'] if r['mean_auc'] is not None else float('nan'))} | "
            f"{(r['mean_delta_Y'] if r['mean_delta_Y'] is not None else float('nan'))} |"
        )

    md += [
        "",
        "## Equal-count → calendar side effect (Tencent edges)",
        "",
        "Same N edges ≠ same hours — burstiness:",
        "",
        "| chunk | n | span_h | y_rate |",
        "|---:|---:|---:|---:|",
    ]
    for s in eq_spans:
        md.append(
            f"| {s['chunk']} | {s['n']} | {s['span_hours']:.1f} | {s['y_rate']:.4f} |"
        )

    md += [
        "",
        "## How to read",
        "",
        "1. **intensity_high vs low (ops)**: if both stay ~1.0 AUC, intensity baseline "
        "is cohort-robust; if only high stays high, baseline is buy-volume concentrated.",
        "2. **content on intensity_low**: last-touch signal without heavy spend — "
        "stronger creative/path story.",
        "3. **credit_high content**: should surface `i_credit_*` / `i_share_*` tops.",
        "4. **Pack scorecard**: compare DiffDB weak vs Tencent ops strong vs Metro shifting.",
        "5. **Calendar scorecard**: width-min grain; do not mix cells with equal-count N.",
        "6. **credit_low y_rate→0**: worst credit-rank edges often have no converts — "
        "pairwise FS has no pairs (expected); credit_high is the actionable slice.",
        "",
    ]
    (args.out / "SUBSET_PAIRWISE.md").write_text("\n".join(md) + "\n")
    # refresh md table was already written — rewrite full file including insight
    print(
        json.dumps(
            {
                "n_subsets": len(subset_rows),
                "ok": sum(1 for r in subset_rows if r.get("ok")),
                "intensity_col": inten_col,
                "credit_col": cred_col,
                "highlights": {
                    "intensity_ops_auc": {
                        r["name"]: r.get("mean_auc")
                        for r in subset_rows
                        if r.get("panel") == "ops"
                        and r["name"].startswith("intensity")
                        and r.get("ok")
                    },
                    "credit_high_content_top": next(
                        (
                            r.get("top0")
                            for r in subset_rows
                            if r.get("name") == "credit_high"
                            and r.get("panel") == "content"
                            and r.get("ok")
                        ),
                        None,
                    ),
                    "credit_low_y_rate": next(
                        (
                            r.get("y_rate")
                            for r in subset_rows
                            if r.get("name") == "credit_low" and r.get("panel") == "ops"
                        ),
                        None,
                    ),
                },
                "out": str(args.out),
                "sec": report["sec"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
