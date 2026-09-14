"""Gated online-RF √po_risk0 and PO-tail localization."""
from __future__ import annotations

import numpy as np

from agod.online_rfperm import (
    hop_fires,
    po_tail_mask,
    run_rfperm_stream,
)
from agod.po_iptw import po_iptw_weights
from agod.po_refit import make_batch_stream, run_uniform_last_two


def test_hop_fires_skips_first_and_needs_jump():
    assert hop_fires(2.0, None, gate=1.5) is False
    assert hop_fires(1.4, 1.0, gate=1.5) is False
    assert hop_fires(1.6, 1.0, gate=1.5) is True


def test_sqrt_po_weights_mean_one():
    rng = np.random.default_rng(0)
    po = rng.uniform(0.2, 4.0, size=40)
    w = po_iptw_weights(po, mode="sqrt")
    assert abs(float(w.mean()) - 1.0) < 0.08
    assert np.all(w > 0)


def test_po_tail_keeps_high_risk_fraction():
    r = np.arange(100, dtype=float)
    mask = po_tail_mask(r, q=0.30, min_n=8)
    assert int(mask.sum()) == 30
    assert set(np.flatnonzero(mask).tolist()) == set(range(70, 100))


def test_similar_stream_rfperm_rarely_fires():
    stream = make_batch_stream(n_batches=6, n_per=80, p=8, seed=0, cov=0.0)
    rec = run_rfperm_stream(stream, gate=1.5)
    uni = run_uniform_last_two(stream)
    assert rec["fire_rate"] <= 0.25
    # quiet: last-two uniform, same rows as rfperm
    assert abs(uni["online_mse"] - rec["online_mse"]) < 1e-9


def test_concept_stream_rfperm_fires_at_cut():
    similar = make_batch_stream(n_batches=6, n_per=90, p=10, seed=1, cov=0.0)
    jumped = make_batch_stream(
        n_batches=6, n_per=90, p=10, seed=1, cov=0.0, concept_at=3, concept=1.0
    )
    sim = run_rfperm_stream(similar, gate=1.5)
    jmp = run_rfperm_stream(jumped, gate=1.5)
    assert sim["fire_rate"] <= 0.25
    assert jmp["fire_rate"] >= sim["fire_rate"]
    fired = {h["t"]: h["fired"] for h in jmp["history"]}
    assert fired[1] is False
    assert fired[3] is True
    hop3 = next(h for h in jmp["history"] if h["t"] == 3)
    assert hop3["mean_r1"] > hop3["mean_r0"] * 1.5


def test_fit_predict_single_class_does_not_crash():
    from agod.po_refit import fit_predict

    rng = np.random.default_rng(0)
    X = rng.normal(size=(20, 4))
    y = np.ones(20, dtype=int)
    for kind in ("rf", "xgb"):
        pred = fit_predict(X, y, X[:5], learner=kind, seed=0, task="acc")
        assert pred.shape == (5,)
        assert np.all(pred == 1)


def test_local_tail_trains_fewer_rows_when_fired():
    jumped = make_batch_stream(
        n_batches=6, n_per=80, p=8, seed=4, cov=0.0, concept_at=3, concept=1.0
    )
    full = run_rfperm_stream(jumped, gate=1.5, localize=False)
    loc = run_rfperm_stream(jumped, gate=1.5, localize=True, q=0.30)
    assert loc["fire_rate"] == full["fire_rate"]
    fired_hops = [h for h, g in zip(loc["history"], full["history"]) if h["fired"]]
    assert fired_hops
    for h, g in zip(loc["history"], full["history"]):
        if not h["fired"]:
            assert h["n_train"] == g["n_train"]
        else:
            assert h["n_train"] < g["n_train"]
            assert h["n_train"] <= int(np.ceil(0.30 * g["n_train"])) + 1
