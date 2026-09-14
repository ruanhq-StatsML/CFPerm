"""Rolling PO-learner refit and gated re-adjustment."""
from __future__ import annotations

import numpy as np

from agod.po_refit import (
    assign_hop,
    assign_pair,
    batch_contrast,
    make_batch_stream,
    mix_lambda,
    refit_po_weights,
    residual_hop_ratio,
    run_oracle_switch,
    run_resid_stream,
    run_refit_stream,
    run_switch_stream,
    run_uniform_on_same_rows,
    should_readjust,
)


def test_gate_only_fires_when_t1_risk_higher():
    assert should_readjust(1.0, gate=1.25) is False
    assert should_readjust(1.4, gate=1.25) is True


def test_hop_and_pair_masks():
    batch = np.array([0, 0, 1, 1, 2, 2, 3, 3])
    t0, t1 = assign_hop(batch, 2)
    assert t0.tolist() == [False, False, True, True, False, False, False, False]
    assert t1.tolist() == [False, False, False, False, True, True, False, False]
    t0, t1 = assign_pair(batch, 2)
    assert t0.tolist() == [True, True, False, False, False, False, False, False]
    assert t1.tolist() == [False, False, True, True, True, True, False, False]


def test_similar_stream_gate_rarely_fires():
    stream = make_batch_stream(n_batches=6, n_per=80, p=8, seed=0, cov=0.0)
    gated = run_refit_stream(stream, assign="pair", gate=1.25, always=False)
    always = run_refit_stream(stream, assign="pair", gate=1.25, always=True)
    uni = run_uniform_on_same_rows(stream, assign="pair")
    assert gated["fire_rate"] <= 0.5
    # similar batches: always-on PO should not beat uniform by much
    assert uni["online_mse"] <= always["online_mse"] + 0.08


def test_concept_hop_raises_ratio_and_gate():
    similar = make_batch_stream(n_batches=6, n_per=90, p=10, seed=1, cov=0.0)
    jumped = make_batch_stream(
        n_batches=6, n_per=90, p=10, seed=1, cov=0.0, concept_at=3, concept=1.0
    )
    g_sim = run_refit_stream(similar, assign="hop", gate=1.15, always=False)
    g_jmp = run_refit_stream(jumped, assign="hop", gate=1.15, always=False)
    assert g_jmp["fire_rate"] >= g_sim["fire_rate"]
    # hop after the concept cut should show a large ratio
    ratios = [h["ratio"] for h in g_jmp["history"] if h["t"] >= 3]
    assert max(ratios) > 1.1


def test_refit_weights_mean_one_when_fired():
    stream = make_batch_stream(
        n_batches=5, n_per=60, p=8, seed=2, concept_at=2, concept=1.0
    )
    t0, t1 = assign_hop(stream.batch, 3)
    rec = refit_po_weights(stream.X, stream.y, t0, t1, gate=1.0, always=True)
    w = rec["weights"][t1]
    assert rec["fired"]
    assert abs(float(w.mean()) - 1.0) < 0.05
    assert rec["ok"]


def test_adaptive_matches_pair_uniform_when_quiet():
    from agod.po_refit import run_adaptive_stream

    stream = make_batch_stream(n_batches=6, n_per=80, p=8, seed=3, cov=0.0)
    adp = run_adaptive_stream(stream, gate=1.25)
    uni = run_uniform_on_same_rows(stream, assign="pair")
    assert adp["fire_rate"] <= 0.5
    assert abs(adp["online_mse"] - uni["online_mse"]) < 0.08


def test_soft_lambda_zero_below_gate():
    assert mix_lambda(1.0, gate=1.25) == 0.0
    assert mix_lambda(1.25, gate=1.25) == 0.0
    assert mix_lambda(1.25 + 2.5, gate=1.25, soft_scale=2.0) == 1.0
    assert 0.0 < mix_lambda(1.6, gate=1.25) < 1.0


def test_switch_matches_pair_uniform_when_quiet():
    stream = make_batch_stream(n_batches=6, n_per=80, p=8, seed=3, cov=0.0)
    sw = run_switch_stream(stream, gate=1.25)
    uni = run_uniform_on_same_rows(stream, assign="pair")
    assert sw["fire_rate"] <= 0.5
    assert abs(sw["online_mse"] - uni["online_mse"]) < 0.08


def test_oracle_drops_old_batch_after_cut():
    jumped = make_batch_stream(
        n_batches=6, n_per=80, p=8, seed=4, concept_at=3, concept=1.0
    )
    ora = run_oracle_switch(jumped)
    fires = [h["fired"] for h in ora["history"]]
    # hops t=1,2 are before cut=3; t>=3 fire
    assert fires == [False, False, True, True]


def test_soft_refit_weights_stay_near_one_when_quiet():
    stream = make_batch_stream(n_batches=5, n_per=60, p=8, seed=5, cov=0.0)
    t0, t1 = assign_pair(stream.batch, 3)
    rec = refit_po_weights(stream.X, stream.y, t0, t1, gate=1.25, soft=True)
    w = rec["weights"][t1]
    assert abs(float(w.mean()) - 1.0) < 0.05
    if not rec["fired"]:
        assert np.allclose(w, 1.0)


def test_residual_gate_fires_on_concept_not_similar():
    similar = make_batch_stream(n_batches=6, n_per=80, p=8, seed=6, cov=0.0)
    jumped = make_batch_stream(
        n_batches=6, n_per=80, p=8, seed=6, concept_at=3, concept=1.0
    )
    sim = run_resid_stream(similar, gate=2.0)
    jmp = run_resid_stream(jumped, gate=2.0)
    assert sim["fire_rate"] <= 0.25
    assert jmp["fire_rate"] >= sim["fire_rate"]
    # hop t=3 is the first labeled post-cut batch
    fired = {h["t"]: h["fired"] for h in jmp["history"]}
    assert fired[3] is True
    t0, t1 = assign_hop(jumped.batch, 3)
    rho, _, _ = residual_hop_ratio(jumped.X, jumped.y, t0, t1)
    assert rho > 2.0
