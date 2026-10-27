#!/usr/bin/env python3
"""Annotate the exploded TencentGR event table: gap, fatigue, multi-task labels.

Does not train the three-tower. Not unique CS/CD.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from tencentgr.config import DEFAULT_CFG  # noqa: E402
from tencentgr.event_stream_features import (  # noqa: E402
    annotate_event_stream,
    exposure_rank_table,
    split_behavior_channels,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--events", default="", help="parquet/csv; default cache/events.parquet")
    p.add_argument("--cache-dir", default=DEFAULT_CFG["cache_dir"])
    p.add_argument("--out-dir", default="experiments/tencentgr/event_stream")
    p.add_argument("--max-users", type=int, default=0)
    p.add_argument("--example-user", type=int, default=1)
    return p.parse_args()


def _read_events(path: Path) -> pd.DataFrame:
    if path.suffix == ".csv":
        df = pd.read_csv(path)
        df = df.drop(columns=[c for c in df.columns if str(c).startswith("Unnamed")], errors="ignore")
        return df
    return pd.read_parquet(path)


def main() -> None:
    args = parse_args()
    cache = Path(args.cache_dir)
    events_path = Path(args.events) if args.events else cache / "events.parquet"
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    df = _read_events(events_path)
    item_path = cache / "item_feat.parquet"
    item_feat = pd.read_parquet(item_path) if item_path.exists() else None
    if args.max_users and args.max_users > 0:
        keep = df["user_id"].drop_duplicates().head(int(args.max_users))
        df = df[df["user_id"].isin(keep)].copy()

    ev = annotate_event_stream(df, item_feat=item_feat)
    channels = split_behavior_channels(ev)
    rank = exposure_rank_table(ev)

    expo = ev["is_exposure"] == 1
    summary = {
        "n_events": int(len(ev)),
        "n_users": int(ev["user_id"].nunique()),
        "n_exposure": int(expo.sum()),
        "n_click": int((ev["is_click"] == 1).sum()),
        "n_conversion": int((ev["is_conversion"] == 1).sum()),
        "mean_gap_sec": float(ev["gap_sec"].mean()),
        "median_gap_sec": float(ev["gap_sec"].median()),
        "mean_session_len": float(ev.drop_duplicates(["user_id", "session_id"])["session_len"].mean()),
        "y_ctr_on_exposure": float(rank["y_ctr"].mean()) if len(rank) else None,
        "y_cvr_on_exposure": float(rank["y_cvr"].mean()) if len(rank) else None,
        "y_next_is_click_on_exposure": float(rank["y_next_is_click"].mean()) if len(rank) else None,
        "y_cvr_on_click": float(ev.loc[ev["is_click"] == 1, "y_cvr"].mean()) if (ev["is_click"] == 1).any() else None,
        "frac_item_exp_1d_gt0": float((ev["item_exp_1d"] > 0).mean()),
        "mean_user_exp_1d": float(ev["user_exp_1d"].mean()),
        "channel_lengths": {k: int(len(v)) for k, v in channels.items()},
        "note": (
            "gap/fatigue are features (history before t). y_ctr/y_cvr are same-item "
            "future labels with attribution windows. Not CATE, not unique CS/CD."
        ),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    rank.head(5000).to_csv(out / "exposure_rank_sample.csv", index=False)
    ex = ev.loc[ev["user_id"] == int(args.example_user)].head(20)
    if ex.empty:
        ex = ev.head(20)
    cols = [
        "user_id", "item_id", "action_type", "timestamp", "time", "gap_sec",
        "session_id", "pos_in_session", "item_exp_1h", "item_exp_1d", "user_exp_1d",
        "y_ctr", "y_cvr", "y_next_is_click",
    ]
    cols = [c for c in cols if c in ex.columns]
    ex[cols].to_csv(out / "example_user.csv", index=False)
    print(json.dumps(summary, indent=2))
    print("example_user rows", len(ex), "wrote", out)


if __name__ == "__main__":
    main()
