"""Tests for smart_bandit_router_v1."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from smart_bandit_router import (  # noqa: E402
    ACTIONS,
    BayesianSubgroupBandit,
    FairContextualBandit,
    business_reward,
    go_live_gates,
    grow_context_leaves,
    leaf_support_ok,
    make_openml_stream,
    make_synthetic_credit_stream,
    run_dataset_eval,
    run_policy,
    shaped_reward,
)


def test_granularity_respects_min_support():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(2000, 4))
    y = (X[:, 0] + 0.5 * X[:, 1] > 0).astype(int)
    gran = grow_context_leaves(X, y, min_support=100, max_leaf_nodes=16, seed=0)
    assert gran["n_leaves"] >= 2
    assert leaf_support_ok(gran["leaf_ids"], gran["min_support_used"])


def test_business_reward_signs():
    labels = np.array([1, 0, 1, 0])
    scores = np.array([0.9, 0.9, 0.1, 0.1])
    r = business_reward("Keep", labels, scores)
    assert r[0] > 0  # TP
    assert r[1] < 0  # FP
    assert r[2] < 0  # FN
    assert r[3] > 0  # TN


def test_bayesian_subgroup_updates():
    b = BayesianSubgroupBandit(n_actions=3, seed=1)
    ctx = np.array([0.7, 0.2, 0.1, 0.3])
    a = b.choose_action(ctx, leaf_id=3, protected=0)
    assert a in (0, 1, 2)
    b.update(3, a, reward=1.0, protected=0)
    assert b.n_all == 1
    assert b.leaf_counts[3] == 1
    report = b.subgroup_bayes_report()
    assert "n_leaves_mature" in report


def test_fair_bandit_alias_backcompat():
    b = FairContextualBandit(n_actions=3, context_dim=4, seed=1)
    ctx = np.array([0.7, 0.2, 0.1, 0.3])
    a = b.choose_action(ctx, protected=0)
    assert a in (0, 1, 2)
    b.update(ctx, a, reward=1.0, protected=0)
    assert b.n_all == 1


def test_synthetic_bandit_beats_keep():
    stream = make_synthetic_credit_stream(n_samples=8000, n_batches=12, seed=2)
    rec = run_dataset_eval(stream, seed=2)
    band = rec["go_live"]["summary_reward"]["bandit"]
    keep = rec["go_live"]["summary_reward"]["keep_baseline"]
    assert band >= keep - 1e-6
    assert rec["policies"]["bandit"]["action_share_overall"]["Keep"] >= 0.2
    assert "lift_vs_keep" in rec["go_live"]["gates"]
    assert rec["go_live"]["verdict"] in ("READY", "CONDITIONAL", "HOLD")
    assert "bayesian_subgroup_confidence" in rec["go_live"]["gates"]


def test_breast_cancer_stream_runs():
    stream = make_openml_stream("breast_cancer", n_batches=8, min_support=30, seed=3, max_n=4000)
    assert stream.contexts.shape[1] == 4
    base = run_policy(stream, policy="baseline_keep", seed=3)
    band = run_policy(stream, policy="bandit", seed=3)
    assert np.isfinite(band["total_reward"])
    assert np.isfinite(base["total_reward"])
    assert "bayes_subgroups" in band


def test_go_live_gate_structure():
    stream = make_synthetic_credit_stream(n_samples=4000, n_batches=8, seed=4)
    band = run_policy(stream, policy="bandit", seed=4)
    base = run_policy(stream, policy="baseline_keep", seed=4)
    gates = go_live_gates(
        band, base, oracle_total=band["total_reward"] + 1000, gran_meta=stream.meta["granularity"]
    )
    assert "ready_for_online" in gates
    assert gates["verdict"] in ("READY", "CONDITIONAL", "HOLD")
    assert "evaluator_rubric" in gates
    assert set(ACTIONS).issubset(set(band["action_share_overall"]))


def test_shaped_reward_finite():
    r = shaped_reward(
        "Keep",
        np.array([1]),
        np.array([0.8]),
        np.array([0.2]),
        np.array([0.3]),
        np.array([0.1]),
        np.array([0.2]),
    )
    assert np.isfinite(r).all()
