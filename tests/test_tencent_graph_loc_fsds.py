"""Smoke test for TencentGR graph localization → FSDS."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from scripts.tencent_gr.graph_loc_fsds import (
    DIMS,
    dim_of,
    eval_inject,
    load_board,
    run_l1,
    top_k,
)

ROOT = Path(__file__).resolve().parents[1] / "results" / "tencent_gr_fs150"


def test_dim_mapping_covers_board_feats():
    X, Y, W, names, pack = load_board()
    assert X.shape[0] == 6000
    assert set(pack["groups"]) == set(DIMS)
    assert sum(pack["meta"]["dim_sizes"].values()) == len(names)
    assert all(dim_of(n) in DIMS for n in names)


def test_inject_hit_at_least_top2_on_strong():
    X, Y, W, names, pack = load_board()
    groups = pack["groups"]
    # Strong, no confound — should be easy
    out = eval_inject(
        X, Y, W, names, groups, strengths=(2.5,), confounds=(0.0,), seed=0
    )
    assert out["summary"]["hit_at_2_graph"] >= 0.75
    assert out["summary"]["n_trials"] == len(DIMS)


def test_natural_l1_returns_s_star():
    X, Y, W, names, pack = load_board()
    l1 = run_l1(X, Y, W, pack["groups"], seed=0, heavy=False)
    assert len(l1["S_star"]) == 2
    assert set(l1["S_star"]).issubset(set(DIMS))
    assert abs(sum(l1["blended"].values()) - 1.0) < 1e-6
