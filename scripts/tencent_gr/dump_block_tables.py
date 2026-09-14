#!/usr/bin/env python3
"""把 block 中间表落盘。Spark 侧同一套粒。

  python3 scripts/tencent_gr/dump_block_tables.py
  python3 scripts/tencent_gr/dump_block_tables.py --max-users 6000
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from block_tables import build_user_table, events_frame  # noqa: E402


def add_empty_lags(post):
    """空路径档 + 空路径滞后（此前各单 empty_any 的均值）。"""
    import numpy as np

    p = post.sort_values(["user_id", "cnv_ts"]).copy()
    p["empty_any"] = p["dt_any_min"].isna().astype(float)
    p["empty_same"] = p["wo_prior_clk"].fillna(True).astype(float)
    shifted = p.groupby("user_id")["empty_any"].shift(1)
    ok = shifted.notna()
    csum = shifted.fillna(0.0).groupby(p["user_id"]).cumsum()
    ccnt = ok.groupby(p["user_id"]).cumsum()
    p["lag_empty_any"] = np.where(ccnt > 0, csum / ccnt, np.nan)
    return p


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=ROOT / "results/tencent_gr_fs150/tables")
    ap.add_argument("--root", type=Path, default=ROOT / "data/tencent_subset")
    ap.add_argument("--max-users", type=int, default=0, help="0 = 只写 demo 三张表")
    ap.add_argument("--prefix-frac", type=float, default=0.75)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    if args.max_users <= 0:
        from block_tables import CLK, CNV, EXP

        t0 = 1_700_000_000
        evs = [
            (11, EXP, t0, 99.0),
            (11, CLK, t0 + 20, 99.0),
            (11, CNV, t0 + 3600, 99.0),
            (11, CLK, t0 + 3600 + 600, 99.0),
            (22, EXP, t0 + 86400 + 10, 50.0),
            (22, CLK, t0 + 86400 + 30, 50.0),
            (22, CNV, t0 + 86400 + 40, 50.0),
            (99, EXP, t0 + 20 * 86400, None),
        ]
        ev = events_frame(1, evs)
    else:
        from run_post_cnv import load_prefix_ev

        ev = load_prefix_ev(args.root, args.max_users, args.prefix_frac)

    user, attr, post = build_user_table(ev)
    post = add_empty_lags(post)
    ev.to_parquet(args.out / "ev.parquet", index=False)
    attr.to_parquet(args.out / "attr.parquet", index=False)
    post.to_parquet(args.out / "post.parquet", index=False)
    user.to_parquet(args.out / "user.parquet", index=False)
    print("wrote", args.out)
    print("ev", len(ev), "attr", len(attr), "post", len(post), "user", len(user))


if __name__ == "__main__":
    main()
