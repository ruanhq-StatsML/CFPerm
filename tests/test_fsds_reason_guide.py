"""Tests for FSDS reasoning empowerment — adjust next step + ROI ledger."""
from __future__ import annotations

import numpy as np

from agod.fsds_reason_guide import (
    FEATURE_NAMES,
    GuideConfig,
    EconParams,
    fsds_attribute,
    guide_decision,
    justify_bullets,
    net_value_path,
    node_importance,
    run_fsds_guided,
    run_score_greedy,
    run_suite,
    simulate_trajectory,
    stack_xy,
)


def test_simulate_has_opt_path_and_feature_dim():
    tr = simulate_trajectory(0, seed=1, max_depth=3, branching=2)
    assert tr.opt_path
    assert all(tr.nodes[i].y_opt == 1 for i in tr.opt_path)
    n = next(iter(tr.nodes.values()))
    assert n.x.shape == (len(FEATURE_NAMES),)
    leaves = [i for i, n in tr.nodes.items() if not n.children]
    assert sum(tr.nodes[i].y_opt for i in leaves) == 1


def test_fsds_attribute_auc_and_regime():
    trajs = [simulate_trajectory(i, seed=2, max_depth=3, branching=2) for i in range(24)]
    X, y, names = stack_xy(trajs)
    attr = fsds_attribute(X, y, names, seed=0)
    assert attr["auc_y"] >= 0.55  # latent leak should beat chance
    assert attr["regime"] in ("covariate_shift", "concept_drift")
    assert len(attr["top_features"]) >= 3
    assert abs(sum(attr["fused_importance"].values()) - 1.0) < 1e-6


def test_guide_prune_and_expand():
    tr = simulate_trajectory(3, seed=0, max_depth=3, branching=2)
    leaf = next(i for i, n in tr.nodes.items() if not n.children)
    cfg = GuideConfig(imp_prune=0.9, gap_prune=10.0, max_depth=2)  # force depth prune
    d = guide_decision(tr.nodes[leaf], importance=0.1, cfg=cfg)
    assert d["action"] == "prune"
    root = tr.nodes[tr.root_id]
    child = tr.nodes[root.children[0]]
    cfg2 = GuideConfig(imp_expand=0.4, gap_expand=-10, max_depth=10)
    child.x = child.x.copy()
    child.x[list(FEATURE_NAMES).index("score_trend")] = 0.5
    child.x[list(FEATURE_NAMES).index("score_gap")] = 0.5
    d2 = guide_decision(child, importance=0.9, cfg=cfg2)
    assert d2["action"] in ("expand", "reorder")


def test_guided_vs_greedy_suite_smoke():
    suite = run_suite(n_train=20, n_test=12, seed=7, max_depth=3, branching=2)
    assert "incremental" in suite
    assert suite["auc_y"] > 0.5
    assert suite["guided"]["success_rate"] >= 0.0
    # ledger keys
    for k in ("roi", "net_incremental", "delta_success", "delta_tokens"):
        assert k in suite["incremental"]
    lines = justify_bullets(suite)
    assert any("ROI" in L or "roi" in L.lower() for L in lines)


def test_net_value_and_search_paths():
    tr = simulate_trajectory(5, seed=3, max_depth=3, branching=2)
    X, y, names = stack_xy([tr] + [simulate_trajectory(100 + i, seed=3) for i in range(15)])
    attr = fsds_attribute(X, y, names, seed=1)
    b = run_score_greedy(tr, beam=2)
    g = run_fsds_guided(
        tr,
        attr["proba_model"],
        cfg=GuideConfig(max_depth=3, beam=2),
        imp_mean=0.5,
        imp_std=0.2,
    )
    nb = net_value_path(tr, b, p=EconParams())
    ng = net_value_path(tr, g, p=EconParams())
    assert "net" in nb and "net" in ng
    assert b.nodes_visited >= 1 and g.nodes_visited >= 1
