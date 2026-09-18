"""Tests for time-window FE gap + event filter."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "tencent_gr"))
from time_window_feats import (  # noqa: E402
    DAY,
    filter_events,
    fsds_feature_columns,
    propose_two_windows,
)
import pandas as pd


def test_propose_gap_at_least_month():
    t0 = 1_700_000_000
    t1 = t0 + 200 * DAY
    w1, w2 = propose_two_windows(t0, t1, window_days=45, gap_days=30)
    assert (w2.t_start - w1.t_end) / DAY >= 30 - 1e-6
    assert w1.t_end <= w2.t_start


def test_filter_events_half_open():
    evs = [(1, 0, 100), (2, 1, 200), (3, 2, 300)]
    got = filter_events(evs, 100, 300)
    assert got == [(1, 0, 100), (2, 1, 200)]


def test_fsds_drops_leak_cols():
    g = pd.DataFrame(
        {
            "user_id": [1],
            "item_id": [2],
            "y_convert": [1],
            "e_n_exp": [3.0],
            "e_n_cnv": [1.0],
            "u_cvr": [0.2],
            "i_share_linear": [0.1],
        }
    )
    cols = fsds_feature_columns(g)
    assert "e_n_exp" in cols
    assert "i_share_linear" in cols
    assert "e_n_cnv" not in cols
    assert "u_cvr" not in cols
    assert "y_convert" not in cols


def test_w1w2_candidate_drops_within_window_ranks():
    from run_w1w2_feature_select import w1w2_candidate_columns

    cols = [
        "i_share_linear",
        "u_n_events",
        "u_user_activity_rank",
        "i_item_pop_rank",
        "i_item_credit_rank",
        "ui_pop_mismatch",
    ]
    got = w1w2_candidate_columns(cols)
    assert got == ["i_share_linear", "u_n_events", "ui_pop_mismatch"]
    assert all("rank" not in c for c in got)
