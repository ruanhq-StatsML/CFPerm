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


def test_mmd_po_cmean_helpers():
    from run_standardize_mmd_fsds import (
        fit_standardizer,
        rbf_mmd2,
        standardize,
        w1w2_candidate_columns,
    )

    rng = np.random.default_rng(0)
    raw0 = rng.normal(scale=[1.0, 100.0, 0.01, 50.0], size=(80, 4))
    raw1 = raw0 + np.array([0.5, 20.0, 0.0, -5.0])
    sc = fit_standardizer(raw0)
    X0, X1 = standardize(sc, raw0), standardize(sc, raw1)
    assert rbf_mmd2(X0, X1, max_n=40, rng=rng) > rbf_mmd2(X0, X0 + 0.01, max_n=40, rng=rng)
    assert rbf_mmd2(X0, X1, max_n=40, rng=rng) >= 0.0
    assert "f_rank" not in w1w2_candidate_columns(["f", "f_rank"])


def test_fsds_pipeline_starts_with_standardizer():
    from run_standardize_mmd_fsds import run_fsds

    rng = np.random.default_rng(1)
    n = 120
    g = pd.DataFrame(
        {
            "user_id": rng.integers(0, 20, n),
            "item_id": rng.integers(0, 10, n),
            "y_convert": np.concatenate([np.zeros(60, int), np.ones(60, int)]),
            "f_big": rng.normal(scale=1e6, size=n),
            "f_small": rng.normal(scale=1e-3, size=n),
            "f_sig": np.concatenate([rng.normal(0, 1, 60), rng.normal(2, 1, 60)]),
        }
    )
    cols = ["f_big", "f_small", "f_sig"]
    tr, te = g.iloc[:80], g.iloc[80:]
    res = run_fsds(tr, te, cols, select_k=2, seed=0)
    assert res["ok"]
    assert res["pipeline"].startswith("StandardScaler")
    assert "f_sig" in res["selected"] or res["ranking"].iloc[0]["feature"] == "f_sig"


def test_gt_subset_evaluator(tmp_path):
    from gt_subset_evaluator import evaluate_gt, evaluate_item_subset, load_gt_items

    items = tmp_path / "gt_items.csv"
    items.write_text("item_id\n10\n20\n30\n40\n")
    orders = tmp_path / "gt_orders.csv"
    orders.write_text("order_id,sku_id\nA,10\nA,99\nB,20\nC,50\n")
    assert load_gt_items(items) == {10, 20, 30, 40}
    pred = [10, 20, 99, 7]
    m = evaluate_item_subset(pred, {10, 20, 30, 40}, ks=(2, 4))
    assert m["hit@2"] == 2.0
    assert abs(m["precision@2"] - 1.0) < 1e-9
    assert abs(m["recall@2"] - 0.5) < 1e-9
    blob = evaluate_gt(pred, gt_items_path=items, gt_orders_path=orders, ks=(2, 4))
    assert blob["available"]
    assert blob["orders"]["n_orders"] == 3.0
    assert blob["orders"]["order_coverage"] > 0.0
