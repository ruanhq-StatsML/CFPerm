#!/usr/bin/env python3
"""Run gated PO-risk on the TencentGR 150-feat FS matrix.

Clock = last event timestamp in the prefix (``_seq_t_end``). Label =
``future_cnv`` (classification). Same last-two-batch methods as the
real-clock board: uniform / DRE / rfperm / resid.

  PYTHONPATH=. python3 scripts/tencent_gr/run_po_on_fs.py \\
    --feats results/tencent_gr_fs150/user_feats_selected.parquet
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from agod.online_rfperm import run_online_rfperm
from agod.po_refit import (
    run_dre_last_two,
    run_resid_stream,
    run_uniform_last_two,
    stream_from_xy,
)


def load_stream(path: Path, n_per: int, n_batches: int):
    df = pd.read_parquet(path)
    y_col = "_label_future_cnv" if "_label_future_cnv" in df.columns else "future_cnv"
    clock = "_seq_t_end" if "_seq_t_end" in df.columns else None
    if clock:
        df = df.sort_values(clock, kind="mergesort").reset_index(drop=True)
    drop = {"user_id", y_col, "_label_pay_user", "_seq_t_end"}
    num = [
        c
        for c in df.columns
        if c not in drop and np.issubdtype(df[c].dtype, np.number)
    ]
    X = df[num].fillna(0).to_numpy(np.float64)
    y = df[y_col].fillna(0).to_numpy(np.int64)
    need = int(n_per) * int(n_batches)
    if len(X) < need:
        n_batches = max(6, len(X) // int(n_per))
        need = int(n_per) * n_batches
    return stream_from_xy(
        X,
        y,
        n_per=n_per,
        n_batches=n_batches,
        name="tencent_fs150",
        task="acc",
        meta={"clock": "seq_t_end", "n_feat": len(num), "path": str(path)},
    ), {
        "n_rows": int(len(df)),
        "n_used": int(need),
        "n_feat": int(len(num)),
        "pos_rate": float(y[:need].mean()) if need else float("nan"),
        "clock": "seq_t_end" if clock else "file_order",
        "n_per": int(n_per),
        "n_batches": int(n_batches),
    }


def _score(rec):
    task = rec.get("metric") or rec.get("task") or "acc"
    v = rec["online_mse"] if "online_mse" in rec else rec["online_score"]
    return float(v), str(task)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--feats",
        type=Path,
        default=Path("results/tencent_gr_fs150/user_feats_selected.parquet"),
    )
    ap.add_argument("--n-per", type=int, default=250)
    ap.add_argument("--n-batches", type=int, default=12)
    ap.add_argument("--gate", type=float, default=1.5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("results/tencent_gr_fs150"))
    args = ap.parse_args()

    stream, meta = load_stream(args.feats, args.n_per, args.n_batches)
    kw = dict(learner="rf", seed=args.seed)
    recs = {
        "uniform_pair": run_uniform_last_two(stream, **kw),
        "dre": run_dre_last_two(stream, **kw),
        "rfperm": run_online_rfperm(stream, gate=args.gate, **kw),
        "resid": run_resid_stream(stream, gate=2.0, po_on_fire=False, **kw),
    }
    board = {}
    for name, rec in recs.items():
        score, task = _score(rec)
        board[name] = {
            "online_score": score,
            "task": task,
            "fire_rate": float(rec.get("fire_rate", 0.0)),
            "n_hops": int(len(rec.get("history") or [])),
        }
    summary = {"meta": meta, "gate": args.gate, "board": board}
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "PO_STREAM_SUMMARY.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    uni = board["uniform_pair"]["online_score"]
    dre = board["dre"]["online_score"]
    rfp = board["rfperm"]["online_score"]
    res = board["resid"]["online_score"]
    fire = board["rfperm"]["fire_rate"]
    md = f"""# TencentGR 150-feat PO-risk stream

Clock = last prefix event time (`_seq_t_end`). Label = `future_cnv` (Acc ↑).
Batches = {meta['n_batches']} × {meta['n_per']} (used {meta['n_used']} / {meta['n_rows']} users, {meta['n_feat']} feats).
Gate γ={args.gate}.

| method | Acc | fire |
|---|---:|---:|
| uniform last-two | {uni:.4f} | — |
| DRE last-two | {dre:.4f} | — |
| rfperm √PO | {rfp:.4f} | {fire:.2f} |
| resid drop-old | {res:.4f} | {board['resid']['fire_rate']:.2f} |

Quiet hops should match uniform. Fire only when consecutive OOS probe
error jumps with a non-vacuous denominator.
"""
    (args.out / "PO_STREAM_REPORT.md").write_text(md, encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print("wrote", args.out)


if __name__ == "__main__":
    main()
