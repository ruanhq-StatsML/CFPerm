#!/usr/bin/env python3
"""TencentGR merchant multi-step attribution board.

Funnel steps: expose (0) → click (1) → convert (2).
Visualizes:
  1) Aggregate funnel counts / rates
  2) Multi-step lag attribution (exp→clk, clk→cnv, exp→cnv)
  3) Time-window attribution rates (clk→cnv within 5m..7d)
  4) Optional: mm_emb coverage (item OID hit-rate for emb_84 / emb_82)

  PYTHONPATH=. python3 scripts/tencent_gr/viz_multistep_attribution.py \
    --root data/tencent_subset --max-users 4000
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[2]
A_EXP, A_CLK, A_CNV = 0, 1, 2
Event = Tuple[int, int, int]
BUCKETS_MIN = [5, 30, 60, 360, 1440, 10080]


def _pq(d: Path) -> List[Path]:
    return sorted(p for p in d.glob("*.parquet") if p.is_file())


def parse_events(seq: list) -> List[Event]:
    evs = []
    for e in seq:
        if not isinstance(e, dict):
            e = dict(e)
        evs.append((int(e["item_id"]), int(e["action_type"]), int(e["timestamp"])))
    evs.sort(key=lambda x: x[2])
    return evs


def iter_users(seq_dir: Path, max_users: int):
    n = 0
    for path in _pq(seq_dir):
        pf = pq.ParquetFile(path)
        for rg in range(pf.num_row_groups):
            df = pf.read_row_group(rg).to_pandas()
            for _, row in df.iterrows():
                yield int(row["user_id"]), parse_events(list(row["seq"]))
                n += 1
                if n >= max_users:
                    return


def accumulate(evs: Sequence[Event], agg: dict) -> None:
    last_exp: Dict[int, int] = {}
    first_exp: Dict[int, int] = {}
    last_clk: Dict[int, int] = {}
    first_clk: Dict[int, int] = {}
    last_any_clk: Optional[int] = None

    for iid, act, ts in evs:
        if act == A_EXP:
            agg["n_exp"] += 1
            last_exp[iid] = ts
            first_exp.setdefault(iid, ts)
        elif act == A_CLK:
            agg["n_clk"] += 1
            if iid in last_exp:
                agg["exp2clk"].append((ts - last_exp[iid]) / 60.0)
            last_clk[iid] = ts
            first_clk.setdefault(iid, ts)
            last_any_clk = ts
        elif act == A_CNV:
            agg["n_cnv"] += 1
            if iid in last_clk:
                dt = (ts - last_clk[iid]) / 60.0
                agg["clk2cnv"].append(dt)
                for b in BUCKETS_MIN:
                    if dt <= b:
                        agg["clk2cnv_bucket"][b] += 1.0
            else:
                agg["cnv_wo_clk"] += 1
            if iid in first_clk:
                agg["firstclk2cnv"].append((ts - first_clk[iid]) / 60.0)
            if iid in first_exp:
                agg["exp2cnv"].append((ts - first_exp[iid]) / 60.0)
            if last_any_clk is not None:
                agg["anyclk2cnv"].append((ts - last_any_clk) / 60.0)


def mm_coverage(root: Path, item_rids: set, indexer: dict, max_rows: int = 400000) -> dict:
    """Fraction of seen item RIDs that have mm_emb (via OID map)."""
    i_map = indexer.get("i", {})  # OID -> RID typically; confirm below
    # indexer['i'] is usually OID->RID; invert to RID->OID
    rid2oid: Dict[int, int] = {}
    sample_items = list(i_map.items())[:5]
    # Heuristic: if values look like dense small ints and keys larger → OID->RID
    keys = [k for k, _ in sample_items]
    vals = [v for _, v in sample_items]
    oid2rid = {int(k): int(v) for k, v in i_map.items()}
    for oid, rid in oid2rid.items():
        rid2oid[rid] = oid

    need_oids = {rid2oid[r] for r in item_rids if r in rid2oid}
    out = {"n_items_seen": len(item_rids), "n_with_oid": len(need_oids)}
    for name, dname in [
        ("emb84_d32", "emb_84_32_parquet"),
        ("emb82_d1024", "emb_82_1024_parquet"),
    ]:
        emb_dir = root / "mm_emb" / dname
        hit = 0
        scanned = 0
        if emb_dir.is_dir():
            want = set(need_oids)
            for p in _pq(emb_dir):
                df = pq.read_table(p, columns=["anonymous_cid"]).to_pandas()
                for cid in df["anonymous_cid"].to_numpy():
                    try:
                        oid = int(cid)
                    except Exception:
                        continue
                    scanned += 1
                    if oid in want:
                        hit += 1
                        want.discard(oid)
                    if max_rows and scanned >= max_rows:
                        break
                if max_rows and scanned >= max_rows:
                    break
                if not want:
                    break
        out[name] = {
            "hit": hit,
            "need": len(need_oids),
            "rate": float(hit) / float(len(need_oids) or 1),
            "scanned": scanned,
        }
    return out


def _pct(xs: List[float], q: float) -> float:
    if not xs:
        return float("nan")
    return float(np.percentile(np.asarray(xs, float), q))


def summarize(agg: dict, n_users: int) -> dict:
    n_exp, n_clk, n_cnv = agg["n_exp"], agg["n_clk"], agg["n_cnv"]
    buckets = {
        str(b): {
            "cnt": float(agg["clk2cnv_bucket"][b]),
            "rate_of_cnv": float(agg["clk2cnv_bucket"][b]) / float(n_cnv or 1),
        }
        for b in BUCKETS_MIN
    }
    return {
        "n_users": n_users,
        "funnel": {
            "n_exp": int(n_exp),
            "n_clk": int(n_clk),
            "n_cnv": int(n_cnv),
            "ctr": float(n_clk) / float(n_exp or 1),
            "cvr": float(n_cnv) / float(n_clk or 1),
            "ctcvr": float(n_cnv) / float(n_exp or 1),
            "cnv_wo_prior_clk": int(agg["cnv_wo_clk"]),
            "cnv_wo_prior_clk_rate": float(agg["cnv_wo_clk"]) / float(n_cnv or 1),
        },
        "lags_minutes": {
            k: {
                "n": len(agg[k]),
                "p50": _pct(agg[k], 50),
                "p90": _pct(agg[k], 90),
                "mean": float(np.mean(agg[k])) if agg[k] else float("nan"),
            }
            for k in ("exp2clk", "clk2cnv", "firstclk2cnv", "exp2cnv", "anyclk2cnv")
        },
        "clk2cnv_within_min": buckets,
    }


def plot_board(summary: dict, mm: Optional[dict], out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(12.5, 8.2))
    gs = fig.add_gridspec(2, 2, hspace=0.38, wspace=0.28)

    # --- (1) funnel ---
    ax0 = fig.add_subplot(gs[0, 0])
    f = summary["funnel"]
    stages = ["expose", "click", "convert"]
    counts = [f["n_exp"], f["n_clk"], f["n_cnv"]]
    colors = ["#4C78A8", "#F58518", "#54A24B"]
    bars = ax0.bar(stages, counts, color=colors, width=0.62)
    for b, c in zip(bars, counts):
        ax0.text(
            b.get_x() + b.get_width() / 2,
            c,
            f"{c:,}",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    ax0.set_title("Multi-step funnel (merchant seq)", fontsize=11)
    ax0.set_ylabel("event count")
    ax0.text(
        0.02,
        0.98,
        f"CTR={f['ctr']:.3f}  CVR={f['cvr']:.3f}  CTCVR={f['ctcvr']:.4f}\n"
        f"cnv w/o prior clk = {f['cnv_wo_prior_clk_rate']:.1%}",
        transform=ax0.transAxes,
        va="top",
        fontsize=8,
        family="monospace",
    )

    # --- (2) lag distributions ---
    ax1 = fig.add_subplot(gs[0, 1])
    lag_keys = ["exp2clk", "clk2cnv", "exp2cnv"]
    labels = ["exp→clk", "clk→cnv", "exp→cnv"]
    data = []
    for k in lag_keys:
        xs = summary["lags_minutes"][k]
        # reconstruct rough samples via log-uniform around p50/p90 is wrong;
        # we only have summary here — plot p50/p90 bars instead
        data.append((xs["p50"], xs["p90"], xs["n"]))
    x = np.arange(len(labels))
    p50 = [d[0] for d in data]
    p90 = [d[1] for d in data]
    ax1.bar(x - 0.18, p50, width=0.36, label="p50 (min)", color="#4C78A8")
    ax1.bar(x + 0.18, p90, width=0.36, label="p90 (min)", color="#E45756")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels)
    ax1.set_yscale("symlog", linthresh=1.0)
    ax1.set_ylabel("minutes (symlog)")
    ax1.set_title("Step-lag attribution (minutes)", fontsize=11)
    ax1.legend(fontsize=8, loc="upper left")
    for i, d in enumerate(data):
        ax1.text(i, max(d[0], d[1]) * 1.05 + 0.1, f"n={d[2]}", ha="center", fontsize=7)

    # --- (3) window rates ---
    ax2 = fig.add_subplot(gs[1, 0])
    buckets = summary["clk2cnv_within_min"]
    order = [str(b) for b in BUCKETS_MIN]
    cnts = [buckets[b]["cnt"] for b in order]
    labels_b = ["5m", "30m", "1h", "6h", "1d", "7d"]
    # same-item clk→cnv is rare in this slice; also show any-click→cnv p50
    ax2.bar(labels_b, cnts, color="#54A24B", width=0.65)
    ax2.set_ylabel("# conversions credited")
    ax2.set_title("same-item clk→cnv within window (counts)", fontsize=11)
    for lab, c in zip(labels_b, cnts):
        ax2.text(lab, c + max(cnts) * 0.03 + 0.05, f"{int(c)}", ha="center", fontsize=8)
    any_lag = summary["lags_minutes"]["anyclk2cnv"]
    ax2.text(
        0.98,
        0.95,
        f"any-click→cnv n={any_lag['n']}\np50={any_lag['p50']:.0f}m  p90={any_lag['p90']:.0f}m\n"
        f"(same-item clk→cnv rare: {summary['funnel']['cnv_wo_prior_clk_rate']:.0%} cnv w/o prior clk)",
        transform=ax2.transAxes,
        ha="right",
        va="top",
        fontsize=7.5,
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#F7F7F7", edgecolor="#CCC"),
    )

    # --- (4) mm emb coverage / step diagram ---
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.axis("off")
    ax3.set_title("Multi-step attribution schema + mm_emb", fontsize=11, pad=8)
    # flow boxes
    boxes = [
        (0.08, 0.62, "① expose\nitem shown"),
        (0.38, 0.62, "② click\nuser engages"),
        (0.68, 0.62, "③ convert\npurchase"),
    ]
    for x0, y0, txt in boxes:
        ax3.add_patch(
            plt.Rectangle(
                (x0, y0),
                0.22,
                0.28,
                fill=True,
                facecolor="#EEF3F8",
                edgecolor="#4C78A8",
                lw=1.5,
                transform=ax3.transAxes,
                clip_on=False,
            )
        )
        ax3.text(
            x0 + 0.11,
            y0 + 0.14,
            txt,
            ha="center",
            va="center",
            fontsize=9,
            transform=ax3.transAxes,
        )
    for x0 in (0.30, 0.60):
        ax3.annotate(
            "",
            xy=(x0 + 0.08, 0.76),
            xytext=(x0, 0.76),
            xycoords=ax3.transAxes,
            textcoords=ax3.transAxes,
            arrowprops=dict(arrowstyle="->", color="#333", lw=1.4),
        )
    ax3.text(
        0.5,
        0.48,
        "credit windows on clk→cnv / exp→cnv lags\n"
        "(first-touch / last-touch / any-click proxies)",
        ha="center",
        va="center",
        fontsize=8,
        transform=ax3.transAxes,
        style="italic",
    )
    if mm:
        lines = [
            f"items seen: {mm['n_items_seen']:,}  (with OID map: {mm['n_with_oid']:,})",
            f"mm_emb84 d=32   hit={mm['emb84_d32']['hit']:,}  "
            f"rate={mm['emb84_d32']['rate']:.1%}",
            f"mm_emb82 d=1024 hit={mm['emb82_d1024']['hit']:,}  "
            f"rate={mm['emb82_d1024']['rate']:.1%}",
        ]
        ax3.text(
            0.05,
            0.05,
            "\n".join(lines),
            ha="left",
            va="bottom",
            fontsize=8,
            family="monospace",
            transform=ax3.transAxes,
            bbox=dict(boxstyle="round,pad=0.35", facecolor="#F7F7F7", edgecolor="#CCC"),
        )
    else:
        ax3.text(0.5, 0.15, "(mm_emb skipped)", ha="center", fontsize=8, transform=ax3.transAxes)

    fig.suptitle(
        f"TencentGR merchant multi-step attribution  (n_users={summary['n_users']})",
        fontsize=13,
        y=0.98,
    )
    out = out_dir / "tencent_merchant_multistep_attribution.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT / "data" / "tencent_subset")
    ap.add_argument("--max-users", type=int, default=4000)
    ap.add_argument("--skip-mm", action="store_true")
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "results" / "tencent_gr_attr_viz",
    )
    args = ap.parse_args()

    agg = {
        "n_exp": 0,
        "n_clk": 0,
        "n_cnv": 0,
        "cnv_wo_clk": 0,
        "exp2clk": [],
        "clk2cnv": [],
        "firstclk2cnv": [],
        "exp2cnv": [],
        "anyclk2cnv": [],
        "clk2cnv_bucket": defaultdict(float),
    }
    item_rids: set = set()
    n_users = 0
    print(f"scanning ≤{args.max_users} users from {args.root / 'seq'} ...", flush=True)
    for _uid, evs in iter_users(args.root / "seq", args.max_users):
        accumulate(evs, agg)
        for iid, _, _ in evs:
            item_rids.add(iid)
        n_users += 1
        if n_users % 500 == 0:
            print(f"  ... {n_users} users", flush=True)

    summary = summarize(agg, n_users)
    mm = None
    if not args.skip_mm:
        print("mm_emb coverage ...", flush=True)
        import pickle

        with open(args.root / "indexer.pkl", "rb") as f:
            indexer = pickle.load(f)
        mm = mm_coverage(args.root, item_rids, indexer)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "multistep_attribution_summary.json").write_text(
        json.dumps({"summary": summary, "mm_emb": mm}, indent=2)
    )
    png = plot_board(summary, mm, args.out_dir)
    md = args.out_dir / "README.md"
    f = summary["funnel"]
    md.write_text(
        "\n".join(
            [
                "# TencentGR merchant multi-step attribution",
                "",
                f"- users: **{summary['n_users']}**",
                f"- funnel: exp={f['n_exp']:,} → clk={f['n_clk']:,} → cnv={f['n_cnv']:,}",
                f"- CTR={f['ctr']:.4f}, CVR={f['cvr']:.4f}, CTCVR={f['ctcvr']:.4f}",
                f"- figure: `{png.name}`",
                "",
            ]
        )
    )
    print(f"wrote {png}", flush=True)
    print(json.dumps(summary["funnel"], indent=2))


if __name__ == "__main__":
    main()
