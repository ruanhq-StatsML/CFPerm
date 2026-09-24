#!/usr/bin/env python3
"""Multi-method feature panel (cmean / MMD-LOCO / PO-VIMP / FSDS).

Not PO-only: surfaces agreement & tension across shift descriptors and
supervised FSDS tips. Display / JSON only — does not change Drill gates.

  PYTHONPATH=. python3 scripts/tencent_gr/feature_methods_panel.py \\
    --feat-diag results/tencent_gr_w1w2_mmd_po_fsds/feature_shift_diagnostics.csv \\
    --ranking results/tencent_gr_w1w2_mmd_po_fsds/fsds_feature_ranking.csv
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

METHOD_COLS = {
    "cmean": "cmean_abs",
    "mmd_loco": "mmd_loco",
    "po_vimp": "po_vimp",
    "fsds": "f_score",
}

# Cross-domain aliases → canonical columns (same statistical roles).
COL_ALIASES: Dict[str, List[str]] = {
    "cmean_abs": ["cmean_abs", "abs_delta", "abs_mean_delta"],
    "mmd_loco": ["mmd_loco", "mmd_delta"],
    "po_vimp": ["po_vimp", "vimp_cov", "vimp", "cov_vimp"],
    "f_score": ["f_score", "fsds_f", "score"],
    "feature": ["feature", "token", "col", "name"],
    "mean_W1": ["mean_W1", "mean_early", "mean_ref"],
    "mean_W2": ["mean_W2", "mean_late", "mean_cur"],
}

METHOD_ROLE = {
    "cmean": "shift_mean",  # |μ_cur − μ_ref|
    "mmd_loco": "shift_mmd",  # ΔMMD when dropping j
    "po_vimp": "shift_po_or_cov",  # PO residual VIMP or X→T cov VIMP
    "fsds": "supervised_y",  # SelectKBest F on y|support
}


def _rename_aliases(df: pd.DataFrame) -> pd.DataFrame:
    """Map domain-specific column names onto the TencentGR canonical schema."""
    if df is None or df.empty:
        return pd.DataFrame() if df is None else df.copy()
    out = df.copy()
    lower = {c.lower(): c for c in out.columns}
    renames = {}
    for canon, aliases in COL_ALIASES.items():
        if canon in out.columns:
            continue
        for a in aliases:
            if a in out.columns:
                renames[a] = canon
                break
            if a.lower() in lower:
                renames[lower[a.lower()]] = canon
                break
    if renames:
        out = out.rename(columns=renames)
    return out


def normalize_feature_tables(
    feat_diag: Optional[pd.DataFrame] = None,
    ranking: Optional[pd.DataFrame] = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Normalize any domain's method CSVs to canonical feature / score columns."""
    diag = _rename_aliases(feat_diag) if feat_diag is not None else pd.DataFrame()
    rank = _rename_aliases(ranking) if ranking is not None else pd.DataFrame()
    if not diag.empty and "feature" in diag.columns:
        diag["feature"] = diag["feature"].astype(str)
        if "cmean_abs" not in diag.columns and "delta" in diag.columns:
            diag["cmean_abs"] = pd.to_numeric(diag["delta"], errors="coerce").abs()
    if not rank.empty and "feature" in rank.columns:
        rank["feature"] = rank["feature"].astype(str)
    return diag, rank


def _top_features(df: pd.DataFrame, col: str, k: int) -> List[str]:
    if col not in df.columns or df.empty:
        return []
    sub = df[["feature", col]].dropna()
    if sub.empty:
        return []
    # higher magnitude = more important; for mmd_loco use abs (can be neg)
    vals = pd.to_numeric(sub[col], errors="coerce")
    if col == "mmd_loco":
        sub = sub.assign(_score=vals.abs())
    else:
        sub = sub.assign(_score=vals)
    sub = sub.dropna(subset=["_score"]).sort_values("_score", ascending=False)
    return [str(x) for x in sub["feature"].head(k).tolist()]


def _rank_map(df: pd.DataFrame, col: str) -> Dict[str, float]:
    if col not in df.columns or df.empty:
        return {}
    s = pd.to_numeric(df[col], errors="coerce")
    if col == "mmd_loco":
        s = s.abs()
    tmp = df.assign(_s=s).dropna(subset=["_s"])
    if tmp.empty:
        return {}
    # rank 1 = largest
    r = tmp["_s"].rank(ascending=False, method="average")
    return {str(f): float(rk) for f, rk in zip(tmp["feature"], r)}


def build_feature_methods_panel(
    feat_diag: Optional[pd.DataFrame],
    ranking: Optional[pd.DataFrame] = None,
    *,
    tip_features: Optional[Sequence[str]] = None,
    tip_signs: Optional[Mapping[str, str]] = None,
    top_k: int = 8,
    domain: Optional[str] = None,
) -> Dict[str, Any]:
    """Build multi-method feature view for summary / review card."""
    diag, rank = normalize_feature_tables(feat_diag, ranking)

    # unify on feature
    if not diag.empty and "feature" in diag.columns and not rank.empty and "feature" in rank.columns:
        merged = diag.merge(
            rank[["feature"] + [c for c in ("f_score", "rank", "selected") if c in rank.columns]],
            on="feature",
            how="outer",
        )
    elif not diag.empty:
        merged = diag.copy()
    elif not rank.empty:
        merged = rank.copy()
    else:
        return {
            "domain": domain,
            "methods": {},
            "tops": {},
            "consensus_top": [],
            "agree_all": [],
            "fsds_only": [],
            "shift_only": [],
            "po_only": [],
            "tip_method_table": [],
            "read": "无 feature diagnostics / FSDS ranking",
            "note": "multi-method panel; Drill untouched",
        }

    tops: Dict[str, List[str]] = {}
    ranks: Dict[str, Dict[str, float]] = {}
    for name, col in METHOD_COLS.items():
        tops[name] = _top_features(merged, col, top_k)
        ranks[name] = _rank_map(merged, col)

    # rank-average consensus over methods that have scores
    all_feats = set()
    for rm in ranks.values():
        all_feats.update(rm.keys())
    consensus_rows = []
    for f in all_feats:
        present = [ranks[m][f] for m in ranks if f in ranks[m]]
        if not present:
            continue
        consensus_rows.append(
            {
                "feature": f,
                "rank_avg": float(np.mean(present)),
                "n_methods": len(present),
            }
        )
    consensus_rows.sort(key=lambda r: (r["rank_avg"], -r["n_methods"]))
    consensus_top = [r["feature"] for r in consensus_rows[:top_k]]

    set_cmean = set(tops.get("cmean") or [])
    set_mmd = set(tops.get("mmd_loco") or [])
    set_po = set(tops.get("po_vimp") or [])
    set_fsds = set(tops.get("fsds") or [])
    shift_union = set_cmean | set_mmd | set_po
    agree_all = sorted(set_cmean & set_mmd & set_po & set_fsds) if set_fsds else sorted(
        set_cmean & set_mmd & set_po
    )
    # looser: in ≥3 of 4
    vote = {}
    for s in (set_cmean, set_mmd, set_po, set_fsds):
        for f in s:
            vote[f] = vote.get(f, 0) + 1
    agree_ge3 = sorted([f for f, v in vote.items() if v >= 3])

    fsds_only = sorted(set_fsds - shift_union)
    shift_only = sorted(shift_union - set_fsds)
    po_only = sorted(set_po - (set_cmean | set_mmd | set_fsds))

    tips = list(tip_features) if tip_features else list(consensus_top)
    signs = dict(tip_signs or {})
    tip_table: List[Dict[str, Any]] = []
    for j in tips[: max(top_k, 12)]:
        row: Dict[str, Any] = {
            "feature": j,
            "sign": signs.get(j, "0"),
            "in_cmean_top": j in set_cmean,
            "in_mmd_top": j in set_mmd,
            "in_po_top": j in set_po,
            "in_fsds_top": j in set_fsds,
            "votes": int(vote.get(j, 0)),
        }
        # raw scores if available
        hit = merged[merged["feature"].astype(str) == str(j)]
        if len(hit):
            h = hit.iloc[0]
            for name, col in METHOD_COLS.items():
                if col in hit.columns and pd.notna(h.get(col)):
                    row[name] = float(h[col])
        tip_table.append(row)

    methods_meta = {
        name: {
            "col": col,
            "role": METHOD_ROLE[name],
            "top": tops.get(name) or [],
            "available": col in merged.columns
            and pd.to_numeric(merged[col], errors="coerce").notna().any(),
        }
        for name, col in METHOD_COLS.items()
    }

    read_bits = []
    if agree_ge3:
        read_bits.append(f"多法共识(≥3): {', '.join(agree_ge3[:5])}")
    if fsds_only:
        read_bits.append(f"仅FSDS(y判别、未必shift): {', '.join(fsds_only[:4])}")
    if shift_only:
        read_bits.append(f"仅shift描述(cmean/MMD/PO): {', '.join(shift_only[:4])}")
    if po_only:
        read_bits.append(f"仅PO/cov-VIMP: {', '.join(po_only[:3])}")
    if not read_bits:
        read_bits.append("方法面板已挂；暂无交叉共识")

    return {
        "domain": domain,
        "methods": methods_meta,
        "tops": tops,
        "consensus_top": consensus_top,
        "agree_all": agree_all,
        "agree_ge3": agree_ge3,
        "fsds_only": fsds_only,
        "shift_only": shift_only,
        "po_only": po_only,
        "tip_method_table": tip_table,
        "read": "；".join(read_bits),
        "note": (
            "cmean/MMD/PO-or-covVIMP = shift 描述；FSDS = y|support 监督 tip；"
            "列名可跨域别名归一；并列不互相改门"
        ),
        "top_k": top_k,
    }


def panel_from_dir(
    out_dir: Path,
    *,
    tip_features: Optional[Sequence[str]] = None,
    tip_signs: Optional[Mapping[str, str]] = None,
    top_k: int = 8,
) -> Dict[str, Any]:
    diag_p = out_dir / "feature_shift_diagnostics.csv"
    rank_p = out_dir / "fsds_feature_ranking.csv"
    diag = pd.read_csv(diag_p) if diag_p.exists() else None
    rank = pd.read_csv(rank_p) if rank_p.exists() else None
    return build_feature_methods_panel(
        diag, rank, tip_features=tip_features, tip_signs=tip_signs, top_k=top_k
    )


def main() -> None:
    ap = argparse.ArgumentParser(description="Multi-method feature panel")
    ap.add_argument("--feat-diag", type=Path, default=None)
    ap.add_argument("--ranking", type=Path, default=None)
    ap.add_argument("--out-dir", type=Path, default=None, help="dir with both CSVs")
    ap.add_argument("--top-k", type=int, default=8)
    ap.add_argument("--out-json", type=Path, default=None)
    args = ap.parse_args()
    if args.out_dir:
        panel = panel_from_dir(args.out_dir, top_k=args.top_k)
    else:
        diag = pd.read_csv(args.feat_diag) if args.feat_diag and args.feat_diag.exists() else None
        rank = pd.read_csv(args.ranking) if args.ranking and args.ranking.exists() else None
        panel = build_feature_methods_panel(diag, rank, top_k=args.top_k)
    text = json.dumps(panel, indent=2, ensure_ascii=False)
    print(text)
    if args.out_json:
        args.out_json.write_text(text + "\n")


if __name__ == "__main__":
    main()
