#!/usr/bin/env python3
"""Sample-chunk adjacent boards — transfer probe, not a production model.

Judgment (not just "cut every N"):
  The HGB number is a *next-chunk transfer probe* after SelectKBest on chunk t.
  It answers: "features that scored for Y on this batch — do they still rank Y
  on the next batch?"  It does **not** claim causal tip, OOD detection, or
  that HGB is the right production scorer.

  Read three numbers together:
    ΔȲ          — did the outcome level move?
    hgb_auc     — did the selected association transfer?
    jaccard_top — did the *which features* stay the same?

  High AUC + high Jaccard  → persistent association (often operational intensity)
  High AUC + low Jaccard   → transferable predictivity but shifting drivers
  Low AUC                  → association does not travel (DiffusionDB tokens ≈ this)

Datasets / panels:
  diffusiondb, tencent_gr (ops), tencent_gr_content (no volume), waymo_proxy,
  metro_interstate, beijing_pm25

  PYTHONPATH=. python3 scripts/run_sample_chunk_adjacent_board.py \\
    --chunk-sizes 1000,2000 --out results/sample_chunk_adjacent_board
"""
from __future__ import annotations

import argparse
import json
import time
from itertools import pairwise
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from agod.stream_packs import load_beijing_pm25, load_metro_interstate
from agod.transfer_null import (
    enrich_row_with_null,
    pack_excess_ci,
    select_pair_indices,
    summarize_null_pack,
)
from scripts.generate_board_reason_codes import generate_reason_codes, render_md
from scripts.run_diffusiondb_temporal_fsds import (
    build_light_token_X,
    load_subset as load_diffusion_subset,
)

ROOT = Path(__file__).resolve().parents[1]

# Tencent: convert leakage vs volume intensity (judgment: volume drives ~1.0 AUC)
TENCENT_CONVERT_LEAK = {
    "e_n_cnv",
    "u_n_cnv",
    "u_has_convert",
    "u_log1p_n_cnv",
    "u_cvr",
    "u_ctcvr",
    "u_user_ctcvr_rank",
    "i_n_cnv",
    "i_n_as_convert_terminal",
    "i_log1p_n_cnv",
    "i_item_cnv_rank",
    "i_cvr",
    "i_ctcvr",
}
TENCENT_VOLUME = {
    "e_n_exp",
    "e_n_clk",
    "e_log1p_exp",
    "e_log1p_clk",
    "e_ctr",
    "u_n_events",
    "u_n_exp",
    "u_n_clk",
    "u_log1p_n_events",
    "u_log1p_n_exp",
    "u_log1p_n_clk",
    "u_ctr",
    "u_user_activity_rank",
    "i_n_exp",
    "i_n_clk",
    "i_n_users",
    "i_ctr",
    "i_log1p_n_exp",
    "i_log1p_n_clk",
    "i_log1p_n_users",
    "i_item_pop_rank",
}


def chunk_by_n(n_rows: int, chunk_size: int) -> np.ndarray:
    """Integer chunk id for each row after sort (every ``chunk_size`` samples)."""
    if chunk_size <= 0:
        raise ValueError("chunk_size must be > 0")
    return (np.arange(n_rows) // int(chunk_size)).astype(int)


def adjacent_chunk_pairs(T: np.ndarray) -> List[Tuple[int, int]]:
    """Occupied chunk ids as adjacent pairs — same idea as ``itertools.pairwise``.

    Mentally (after ``df = df.sort_values(e_last_ts)``)::

        cuts = [0, N, 2N, ...]
        for i0, i1 in pairwise(range(n_chunks)):
            df1 = df.iloc[i0*N : (i0+1)*N]   # train chunk
            df2 = df.iloc[i1*N : (i1+1)*N]   # test chunk
    """
    occ = sorted(int(t) for t in np.unique(T))
    return list(pairwise(occ))


def jaccard(a: Sequence[str], b: Sequence[str]) -> float:
    sa, sb = set(a), set(b)
    if not sa and not sb:
        return 1.0
    if not sa or not sb:
        return 0.0
    return float(len(sa & sb) / len(sa | sb))


def annotate_stability(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Attach jaccard_top vs previous adjacent pair (None on first)."""
    out = []
    prev_top: Optional[List[str]] = None
    for r in rows:
        rr = dict(r)
        rr["jaccard_top"] = (
            None if prev_top is None else jaccard(prev_top, r.get("top_fsds") or [])
        )
        out.append(rr)
        prev_top = list(r.get("top_fsds") or [])
    return out


def mean_jaccard(rows: List[Dict[str, Any]]) -> Optional[float]:
    vals = [r["jaccard_top"] for r in rows if r.get("jaccard_top") is not None]
    return float(np.mean(vals)) if vals else None


def mean_key(rows: List[Dict[str, Any]], key: str) -> Optional[float]:
    vals = [r[key] for r in rows if r.get(key) is not None and np.isfinite(r[key])]
    return float(np.mean(vals)) if vals else None


def ops_content_gap(report: Dict[str, Any], chunk: str) -> Optional[Dict[str, Any]]:
    """Judgment metric: how much of ops transfer is pure volume."""
    ops = report.get("datasets", {}).get("tencent_gr", {})
    content = report.get("datasets", {}).get("tencent_gr_content", {})
    if not ops.get("ok") or not content.get("ok"):
        return None
    ob = ops["by_chunk"].get(chunk)
    cb = content["by_chunk"].get(chunk)
    if not ob or not cb or ob.get("mean_auc") is None or cb.get("mean_auc") is None:
        return None
    return {
        "chunk": int(chunk),
        "ops_auc": ob["mean_auc"],
        "content_auc": cb["mean_auc"],
        "auc_gap_ops_minus_content": float(ob["mean_auc"] - cb["mean_auc"]),
        "ops_top0": (ob["rows"][0]["top_fsds"][:3] if ob.get("rows") else []),
        "content_top0": (cb["rows"][0]["top_fsds"][:3] if cb.get("rows") else []),
        "reading": (
            "volume explains most transfer"
            if ob["mean_auc"] - cb["mean_auc"] >= 0.1
            else "non-volume feats still carry substantial transfer"
        ),
    }


def ship_gate_from_board(cross_rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Explicit: this board alone never greenlights HGB to production."""
    return {
        "promote_HGB_to_production": False,
        "reason": (
            "Adjacent-chunk transfer probe ≠ online scorer. Missing: label delay, "
            "calibration under serve skew, Acc/latency budget, abstain/rollback, "
            "and causal/PO checks for 'true drivers'."
        ),
        "what_board_can_greenlight": [
            "investigate a feature family (ops intensity vs content/credit)",
            "flag weak-transfer packs (e.g. DiffusionDB tokens for NSFW)",
            "compare grains N=1000 vs 2000 as business windows",
        ],
        "what_needs_other_tools": {
            "true_driver": "ablation / PO-risk / randomized or quasi-exp — not AUC",
            "ship_model": "holdout Acc + calibration + cost/latency + gray rollback",
        },
        "cross_pack_snapshot": cross_rows,
    }


def fs_adjacent(
    X: np.ndarray,
    y: np.ndarray,
    T: np.ndarray,
    feat_names: List[str],
    *,
    select_k: int,
    y_quantile: float,
    seed: int,
    max_pairs: int,
    min_n: int = 80,
    n_null_perm: int = 5,
    pair_sample: str = "reservoir",
) -> List[Dict[str, Any]]:
    pairs = adjacent_chunk_pairs(T)
    n_all = len(pairs)
    if max_pairs > 0 and n_all > max_pairs:
        idx = select_pair_indices(
            n_all, max_pairs, mode=pair_sample, seed=seed
        )
        pairs = [pairs[i] for i in idx]
    rows = []
    for t0, t1 in pairs:
        m0, m1 = T == t0, T == t1
        if int(m0.sum()) < min_n or int(m1.sum()) < min_n:
            continue
        y0, y1 = y[m0], y[m1]
        # binary labels stay binary; continuous → early-quantile binary
        uniq0 = set(np.unique(y0[~np.isnan(y0)]).tolist()) if len(y0) else set()
        if uniq0 and uniq0 <= {0.0, 1.0}:
            yb0, yb1 = y0.astype(int), y1.astype(int)
            thr = 0.5
        else:
            thr = float(np.quantile(y0, y_quantile))
            yb0 = (y0 >= thr).astype(int)
            yb1 = (y1 >= thr).astype(int)
        if len(np.unique(yb0)) < 2 or len(np.unique(yb1)) < 2:
            continue
        k = min(select_k, X.shape[1], max(1, int(m0.sum()) - 1))
        try:
            pre = Pipeline(
                [("sc", StandardScaler()), ("var", VarianceThreshold(1e-10))]
            )
            Xv = pre.fit_transform(X[m0], yb0)
            k_eff = min(k, Xv.shape[1])
            sel = SelectKBest(f_classif, k=k_eff)
            Xt = sel.fit_transform(Xv, yb0)
            Xte = sel.transform(pre.transform(X[m1]))
        except Exception:
            continue
        hgb = HistGradientBoostingClassifier(
            max_depth=3, learning_rate=0.1, max_iter=60, random_state=seed
        )
        hgb.fit(Xt, yb0)
        proba_h = hgb.predict_proba(Xte)[:, 1]
        auc_h = float(roc_auc_score(yb1, proba_h))
        # Second probe: same selected X, linear — if both transfer, not an HGB quirk
        lr = LogisticRegression(max_iter=200, random_state=seed)
        proba_l = None
        try:
            lr.fit(Xt, yb0)
            proba_l = lr.predict_proba(Xte)[:, 1]
            auc_l = float(roc_auc_score(yb1, proba_l))
            brier_l = float(brier_score_loss(yb1, proba_l))
        except Exception:
            auc_l, brier_l = float("nan"), float("nan")
            proba_l = None
        try:
            brier_h = float(brier_score_loss(yb1, proba_h))
        except Exception:
            brier_h = float("nan")
        var_mask = pre.named_steps["var"].get_support()
        cols_var = [c for c, m in zip(feat_names, var_mask) if m]
        ranking = (
            pd.DataFrame({"feature": cols_var, "f_score": sel.scores_})
            .sort_values("f_score", ascending=False)
            .head(5)
        )
        dmu = X[m1].mean(0) - X[m0].mean(0)
        top_cmean = [feat_names[i] for i in np.argsort(-np.abs(dmu))[:5]]
        top_fsds = ranking["feature"].tolist()
        row = {
            "t0": int(t0),
            "t1": int(t1),
            "n0": int(m0.sum()),
            "n1": int(m1.sum()),
            "y_mean_0": float(np.mean(y0)),
            "y_mean_1": float(np.mean(y1)),
            "delta_Y": float(np.mean(y1) - np.mean(y0)),
            "y_threshold": thr,
            "hgb_auc": auc_h,
            "logreg_auc": auc_l,
            "hgb_brier": brier_h,
            "logreg_brier": brier_l,
            "top_fsds": top_fsds,
            "top_cmean": top_cmean,
            "fsds_cmean_jaccard": jaccard(top_fsds, top_cmean),
        }
        # Null on HGB scores + ECE on both probes (ranking vs calibration)
        if n_null_perm > 0:
            row = enrich_row_with_null(
                row,
                yb1,
                proba_h,
                auc_key="hgb_auc",
                n_perm=n_null_perm,
                seed=seed + int(t0) * 17,
                selected_k=k_eff,
                proba_lr=proba_l,
            )
        rows.append(row)
    return annotate_stability(rows)


def _tencent_frame(n_sample: int) -> Tuple[pd.DataFrame, str, Optional[str]]:
    path = ROOT / "results/tencent_gr_tabular_ui/feature_grid.parquet"
    df = pd.read_parquet(path)
    y_col = "y_convert" if "y_convert" in df.columns else "e_ctr"
    ts_col = "e_last_ts" if "e_last_ts" in df.columns else None
    if ts_col:
        df = df.sort_values(ts_col).reset_index(drop=True)
    if n_sample < len(df):
        idx = np.linspace(0, len(df) - 1, num=n_sample, dtype=int)
        df = df.iloc[idx].reset_index(drop=True)
    return df, y_col, ts_col


def _tencent_xy(
    df: pd.DataFrame,
    y_col: str,
    ts_col: Optional[str],
    extra_drop: set,
    panel: str,
) -> Tuple[np.ndarray, np.ndarray, List[str], Dict]:
    drop = {
        "user_id",
        "item_id",
        y_col,
        "e_last_ts",
        "last_ts",
        *TENCENT_CONVERT_LEAK,
        *extra_drop,
    }
    feat_cols = [
        c
        for c in df.columns
        if c not in drop and pd.api.types.is_numeric_dtype(df[c])
    ]
    X = np.nan_to_num(df[feat_cols].to_numpy(dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    y = df[y_col].to_numpy(float)
    return X, y, feat_cols, {
        "dataset": panel,
        "y": y_col,
        "sort": ts_col or "row",
        "n": len(df),
        "d": len(feat_cols),
        "dropped_volume": sorted(extra_drop & set(df.columns)),
    }


def load_diffusion(n_sample: int, seed: int) -> Tuple[np.ndarray, np.ndarray, List[str], Dict]:
    meta = ROOT / "data/diffusiondb/metadata.parquet"
    df = load_diffusion_subset(meta, n_sample=n_sample, seed=seed)
    df = df.sort_values("timestamp").reset_index(drop=True)
    X, names, _ = build_light_token_X(
        df["prompt_clean"].tolist(), method="tfidf", max_features=256, seed=seed
    )
    y = df["image_nsfw"].to_numpy(float)
    return X, y, names, {
        "dataset": "diffusiondb",
        "y": "image_nsfw",
        "sort": "timestamp",
        "n": len(df),
        "read_as": "low AUC expected: prompt tokens weakly transfer for NSFW",
    }


def load_tencent(n_sample: int, seed: int) -> Tuple[np.ndarray, np.ndarray, List[str], Dict]:
    del seed
    df, y_col, ts_col = _tencent_frame(n_sample)
    X, y, names, meta = _tencent_xy(df, y_col, ts_col, set(), "tencent_gr_ops")
    meta["read_as"] = (
        "high AUC usually = exposure/click intensity transfers; not a content tip"
    )
    return X, y, names, meta


def load_tencent_content(
    n_sample: int, seed: int
) -> Tuple[np.ndarray, np.ndarray, List[str], Dict]:
    """Judgment panel: drop volume so board is not just 'more exp → convert'."""
    del seed
    df, y_col, ts_col = _tencent_frame(n_sample)
    X, y, names, meta = _tencent_xy(
        df, y_col, ts_col, set(TENCENT_VOLUME), "tencent_gr_content"
    )
    meta["read_as"] = (
        "volume dropped; AUC should fall if intensity was the only transferable signal"
    )
    return X, y, names, meta


def load_waymo(n_sample: int, seed: int) -> Tuple[np.ndarray, np.ndarray, List[str], Dict]:
    del seed
    blob = np.load(ROOT / "data/stream_packs/waymo_proxy/waymo_proxy_xy.npz")
    X = np.asarray(blob["X"], dtype=float)
    y = np.asarray(blob["y"], dtype=float).ravel()
    if n_sample < len(X):
        X, y = X[:n_sample], y[:n_sample]
    names = [f"x{j}" for j in range(X.shape[1])]
    return X, y, names, {
        "dataset": "waymo_proxy",
        "y": "proxy_y",
        "sort": "row_index",
        "n": len(y),
        "d": X.shape[1],
        "read_as": "synthetic gradual drift; high AUC+Jaccard = planted kinematics persist",
    }


def load_metro(n_sample: int, seed: int) -> Tuple[np.ndarray, np.ndarray, List[str], Dict]:
    del seed
    X, y, meta = load_metro_interstate(ROOT, max_n=n_sample)
    names = [f"m{j}" for j in range(X.shape[1])]
    return X.astype(float), y.astype(float), names, {
        "dataset": "metro_interstate",
        "y": meta.get("target", "traffic_volume"),
        "sort": "date_time",
        "n": int(meta["n"]),
        "d": int(meta["d"]),
        "read_as": "hourly traffic; calendar/weather feats should transfer across adjacent hours",
    }


def load_beijing(n_sample: int, seed: int) -> Tuple[np.ndarray, np.ndarray, List[str], Dict]:
    del seed
    X, y, meta = load_beijing_pm25(ROOT, max_n=n_sample)
    names = [f"p{j}" for j in range(X.shape[1])]
    return X.astype(float), y.astype(float), names, {
        "dataset": "beijing_pm25",
        "y": meta.get("target", "pm2.5"),
        "sort": "stamp",
        "n": int(meta["n"]),
        "d": int(meta["d"]),
        "read_as": "meteo+lag; high transfer expected under smooth pollution regimes",
    }


LOADERS = {
    "diffusiondb": load_diffusion,
    "tencent_gr": load_tencent,
    "tencent_gr_content": load_tencent_content,
    "waymo_proxy": load_waymo,
    "metro_interstate": load_metro,
    "beijing_pm25": load_beijing,
}


def interpret_pack(
    mean_auc: Optional[float],
    mean_j: Optional[float],
    mean_excess: Optional[float] = None,
) -> str:
    from agod.claim_router import route_pack_interpretation

    routed = route_pack_interpretation(
        mean_auc=mean_auc,
        mean_jaccard=mean_j,
        mean_excess=mean_excess,
    )
    return str(routed["text"])


def plot_dataset_board(name: str, by_chunk: Dict[int, List[Dict]], out_png: Path) -> None:
    sizes = sorted(by_chunk.keys())
    fig, axes = plt.subplots(1, len(sizes), figsize=(4.2 * len(sizes), 3.6), squeeze=False)
    for ax, cs in zip(axes[0], sizes):
        rows = by_chunk[cs]
        if not rows:
            ax.set_title(f"chunk={cs} (empty)")
            ax.axis("off")
            continue
        xs = np.arange(len(rows))
        dys = [r["delta_Y"] for r in rows]
        aucs = [r["hgb_auc"] for r in rows]
        jacs = [r["jaccard_top"] if r["jaccard_top"] is not None else np.nan for r in rows]
        ax.bar(xs - 0.2, dys, width=0.25, color="#4C78A8", label="ΔȲ")
        ax2 = ax.twinx()
        ax2.plot(xs, aucs, "o-", color="#F58518", ms=4, label="HGB AUC")
        lr = [r.get("logreg_auc", np.nan) for r in rows]
        ax2.plot(xs, lr, "^:", color="#B279A2", ms=4, label="LogReg AUC")
        ax2.plot(xs, jacs, "s--", color="#54A24B", ms=4, label="Jaccard top5")
        ax.axhline(0, color="gray", lw=0.7)
        ax.set_title(f"every {cs} samples")
        ax.set_xlabel("adjacent chunk idx")
        ax.set_ylabel("ΔȲ")
        ax2.set_ylabel("AUC / Jaccard")
        ax2.set_ylim(0, 1.05)
    fig.suptitle(f"{name}: dual probe (HGB/LogReg) + Jaccard", y=1.02)
    fig.tight_layout()
    fig.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)


def plot_cross_pack(summary_rows: List[Dict[str, Any]], out_png: Path) -> None:
    if not summary_rows:
        return
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for r in summary_rows:
        ax.scatter(
            r["mean_auc"],
            abs(r["mean_delta_Y"]) if r["mean_delta_Y"] is not None else 0.0,
            s=40 + 80 * (r["mean_jaccard"] or 0.0),
            label=f"{r['dataset']}@N={r['chunk']}",
        )
        ax.annotate(
            f"{r['dataset'][:6]}@{r['chunk']}",
            (r["mean_auc"], abs(r["mean_delta_Y"] or 0.0)),
            fontsize=7,
            alpha=0.85,
        )
    ax.axvline(0.65, color="gray", ls=":", lw=0.8)
    ax.set_xlabel("mean next-chunk HGB AUC (transfer probe)")
    ax.set_ylabel("|mean ΔȲ|")
    ax.set_title("Cross-pack: transfer vs level-shift (marker size ∝ Jaccard)")
    ax.set_xlim(0.45, 1.02)
    fig.tight_layout()
    fig.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--datasets",
        type=str,
        default=(
            "diffusiondb,tencent_gr,tencent_gr_content,"
            "waymo_proxy,metro_interstate,beijing_pm25"
        ),
    )
    ap.add_argument("--chunk-sizes", type=str, default="1000,2000")
    ap.add_argument("--n-sample", type=int, default=8000, help="cap rows per dataset")
    ap.add_argument("--select-k", type=int, default=20)
    ap.add_argument("--y-quantile", type=float, default=0.7)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-pairs", type=int, default=12)
    ap.add_argument(
        "--pair-sample",
        type=str,
        default="reservoir",
        choices=("reservoir", "linspace", "head"),
        help="how to thin adjacent pairs under --max-pairs (default: unbiased reservoir)",
    )
    ap.add_argument(
        "--out", type=Path, default=ROOT / "results/sample_chunk_adjacent_board"
    )
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    chunk_sizes = [int(x) for x in args.chunk_sizes.split(",") if x.strip()]
    datasets = [x.strip() for x in args.datasets.split(",") if x.strip()]
    t0 = time.time()
    report: Dict[str, Any] = {
        "note": (
            "HGB AUC = next-chunk transfer probe after SelectKBest on chunk t; "
            "not a causal tip / not a production model. Read with ΔȲ, Jaccard, "
            "excess/ECE/CI. Pair thinning: reservoir (unbiased) by default."
        ),
        "pair_sample": args.pair_sample,
        "max_pairs": args.max_pairs,
        "chunk_sizes": chunk_sizes,
        "datasets": {},
    }

    body_md: List[str] = []
    cross_rows: List[Dict[str, Any]] = []

    for ds in datasets:
        if ds not in LOADERS:
            print("skip unknown", ds)
            continue
        print("load", ds, flush=True)
        try:
            X, y, names, meta = LOADERS[ds](args.n_sample, args.seed)
        except Exception as e:
            report["datasets"][ds] = {"ok": False, "error": str(e)}
            continue
        by_chunk: Dict[int, List[Dict]] = {}
        for cs in chunk_sizes:
            if cs >= len(y):
                by_chunk[cs] = []
                continue
            T = chunk_by_n(len(y), cs)
            rows = fs_adjacent(
                X,
                y,
                T,
                names,
                select_k=args.select_k,
                y_quantile=args.y_quantile,
                seed=args.seed,
                max_pairs=args.max_pairs,
                pair_sample=args.pair_sample,
            )
            by_chunk[cs] = rows
            print(f"  chunk={cs} pairs_reported={len(rows)}", flush=True)

        png = args.out / f"{ds}_board.png"
        plot_dataset_board(ds, by_chunk, png)
        pack_blk: Dict[str, Any] = {}
        for cs, rows in by_chunk.items():
            m_auc = mean_key(rows, "hgb_auc")
            m_lr = mean_key(rows, "logreg_auc")
            m_dy = mean_key(rows, "delta_Y")
            m_j = mean_jaccard(rows)
            m_fc = mean_key(rows, "fsds_cmean_jaccard")
            m_bh = mean_key(rows, "hgb_brier")
            null_sum = summarize_null_pack(rows)
            m_ex = null_sum.get("mean_excess_auc")
            excess_ci = pack_excess_ci(rows, key="excess_auc", seed=args.seed)
            ece_ci = pack_excess_ci(rows, key="hgb_ece", seed=args.seed + 1)
            reading = interpret_pack(m_auc, m_j, m_ex)
            pack_blk[str(cs)] = {
                "n_pairs": len(rows),
                "mean_delta_Y": m_dy,
                "mean_auc": m_auc,
                "mean_logreg_auc": m_lr,
                "mean_hgb_brier": m_bh,
                "mean_jaccard_top": m_j,
                "mean_fsds_cmean_jaccard": m_fc,
                "mean_null_auc": null_sum.get("mean_null_auc"),
                "mean_excess_auc": m_ex,
                "excess_auc_ci90": {
                    "lo": excess_ci.get("lo"),
                    "hi": excess_ci.get("hi"),
                    "block_size": excess_ci.get("block_size"),
                },
                "mean_probe_eff": null_sum.get("mean_probe_eff"),
                "mean_hgb_ece": null_sum.get("mean_hgb_ece"),
                "hgb_ece_ci90": {
                    "lo": ece_ci.get("lo"),
                    "hi": ece_ci.get("hi"),
                },
                "mean_logreg_ece": null_sum.get("mean_logreg_ece"),
                "reading": reading,
                "rows": rows,
            }
            if rows:
                cross_rows.append(
                    {
                        "dataset": ds,
                        "chunk": cs,
                        "mean_auc": m_auc,
                        "mean_logreg_auc": m_lr,
                        "mean_delta_Y": m_dy,
                        "mean_jaccard": m_j,
                        "mean_fsds_cmean_jaccard": m_fc,
                        "mean_excess_auc": m_ex,
                        "excess_auc_ci90_lo": excess_ci.get("lo"),
                        "excess_auc_ci90_hi": excess_ci.get("hi"),
                        "mean_probe_eff": null_sum.get("mean_probe_eff"),
                        "mean_hgb_ece": null_sum.get("mean_hgb_ece"),
                        "mean_logreg_ece": null_sum.get("mean_logreg_ece"),
                        "reading": reading,
                    }
                )
        report["datasets"][ds] = {
            "ok": True,
            "meta": meta,
            "by_chunk": pack_blk,
            "board_png": str(png.name),
        }
        body_md += [
            f"## {ds}",
            f"- meta: `{meta}`",
            f"- ![{ds}]({png.name})",
            "",
        ]
        for cs, rows in by_chunk.items():
            blk = pack_blk[str(cs)]
            body_md += [
                f"### every {cs} samples",
                (
                    f"- pairs={len(rows)} · HGB={blk['mean_auc']} · "
                    f"LogReg={blk['mean_logreg_auc']} · "
                    f"excess={blk.get('mean_excess_auc')} "
                    f"CI90=[{(blk.get('excess_auc_ci90') or {}).get('lo')},"
                    f"{(blk.get('excess_auc_ci90') or {}).get('hi')}] · "
                    f"probe_eff={blk.get('mean_probe_eff')} · "
                    f"ECE(H/L)={blk.get('mean_hgb_ece')}/"
                    f"{blk.get('mean_logreg_ece')} · "
                    f"Jaccard={blk['mean_jaccard_top']} · "
                    f"**{blk['reading']}**"
                ),
                "",
                "| t0→t1 | ΔȲ | HGB | excess | ECE_h | ECE_lr | Jac | top |",
                "|---|---:|---:|---:|---:|---:|---:|---|",
            ]
            for r in rows[:8]:
                jac = (
                    f"{r['jaccard_top']:.2f}"
                    if r.get("jaccard_top") is not None
                    else "—"
                )
                ex = r.get("excess_auc", float("nan"))
                eh = r.get("hgb_ece", float("nan"))
                el = r.get("logreg_ece", float("nan"))
                ex_s = f"{ex:.3f}" if np.isfinite(ex) else "—"
                eh_s = f"{eh:.3f}" if np.isfinite(eh) else "—"
                el_s = f"{el:.3f}" if np.isfinite(el) else "—"
                body_md.append(
                    f"| {r['t0']}→{r['t1']} | {r['delta_Y']:.4f} | "
                    f"{r['hgb_auc']:.3f} | {ex_s} | {eh_s} | {el_s} | {jac} | "
                    f"{', '.join(r['top_fsds'][:3])} |"
                )
            body_md.append("")

    cross_png = args.out / "cross_pack_transfer.png"
    plot_cross_pack(cross_rows, cross_png)
    report["cross_pack"] = cross_rows
    report["cross_pack_png"] = cross_png.name
    report["ops_content_gap"] = [
        g
        for g in (ops_content_gap(report, str(cs)) for cs in chunk_sizes)
        if g is not None
    ]
    report["ship_gate"] = ship_gate_from_board(cross_rows)
    report["sec"] = float(time.time() - t0)

    md = [
        "# Sample-chunk adjacent boards (multi-dataset)",
        "",
        report["note"],
        "",
        f"![cross]({cross_png.name})",
        "",
        "**ship_gate:** `promote_HGB_to_production=false` — board alone never ships.",
        "",
        "| pack@N | HGB | LogReg | |ΔȲ| | Jac | F∩c | reading |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for r in cross_rows:
        def _fmt(x: Optional[float]) -> str:
            return f"{x:.3f}" if x is not None and np.isfinite(x) else "—"

        def _fmt2(x: Optional[float]) -> str:
            return f"{x:.2f}" if x is not None and np.isfinite(x) else "—"

        md.append(
            f"| {r['dataset']}@{r['chunk']} | {_fmt(r['mean_auc'])} | "
            f"{_fmt(r['mean_logreg_auc'])} | "
            f"{abs(r['mean_delta_Y'] or 0):.4g} | "
            f"{_fmt2(r['mean_jaccard'])} | "
            f"{_fmt2(r['mean_fsds_cmean_jaccard'])} | "
            f"{r['reading']} |"
        )
    md += ["", "### ops − content gap (Tencent)", ""]
    if report["ops_content_gap"]:
        md += [
            "| N | ops AUC | content AUC | gap | reading |",
            "|---:|---:|---:|---:|---|",
        ]
        for g in report["ops_content_gap"]:
            md.append(
                f"| {g['chunk']} | {g['ops_auc']:.3f} | {g['content_auc']:.3f} | "
                f"{g['auc_gap_ops_minus_content']:.3f} | {g['reading']} |"
            )
            md.append(
                f"- tops ops `{g['ops_top0']}` vs content `{g['content_top0']}`"
            )
    else:
        md.append("_run both `tencent_gr` and `tencent_gr_content` to populate_")
    md += ["", "---", ""] + body_md

    (args.out / "summary.json").write_text(
        json.dumps(report, indent=2, default=str) + "\n"
    )
    (args.out / "SAMPLE_CHUNK_BOARD.md").write_text("\n".join(md) + "\n")

    # Auto reason-code generation from the same summary
    rc = generate_reason_codes(report)
    report["reason_codes"] = {
        "n_codes": rc["n_codes"],
        "codes": [c["code"] for c in rc["codes"]],
        "catalog_version": rc["catalog_version"],
    }
    rc_dir = args.out / "reason_codes"
    rc_dir.mkdir(parents=True, exist_ok=True)
    (rc_dir / "reason_codes.json").write_text(
        json.dumps(rc, indent=2, ensure_ascii=False, default=str) + "\n"
    )
    (rc_dir / "REASON_CODES.md").write_text(render_md(rc))
    (rc_dir / "paste_for_agent.txt").write_text(rc["paste_for_agent"])
    ad = rc.get("ad_scenario") or {}
    if ad.get("ok"):
        (rc_dir / "ad_scenario.json").write_text(
            json.dumps(ad, indent=2, ensure_ascii=False, default=str) + "\n"
        )
        (rc_dir / "ad_scenario_paste.txt").write_text(ad["paste_for_agent"])
        report["reason_codes"]["ad_scenario"] = {
            "family_code": ad.get("family_code"),
            "sign_Dy_board": ad.get("sign_Dy_board"),
            "read": ad.get("read"),
        }
    # rewrite summary with reason_codes index
    (args.out / "summary.json").write_text(
        json.dumps(report, indent=2, default=str) + "\n"
    )

    print(
        json.dumps(
            {
                "datasets": list(report["datasets"].keys()),
                "chunk_sizes": chunk_sizes,
                "ops_content_gap": report["ops_content_gap"],
                "ship_gate": {
                    "promote_HGB_to_production": report["ship_gate"][
                        "promote_HGB_to_production"
                    ],
                    "reason": report["ship_gate"]["reason"],
                },
                "reason_codes": report["reason_codes"],
                "cross_pack": [
                    {
                        "dataset": r["dataset"],
                        "chunk": r["chunk"],
                        "mean_auc": r["mean_auc"],
                        "mean_logreg_auc": r["mean_logreg_auc"],
                        "mean_jaccard": r["mean_jaccard"],
                        "reading": r["reading"],
                    }
                    for r in cross_rows
                ],
                "out": str(args.out),
                "sec": report["sec"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
