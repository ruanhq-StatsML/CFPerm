"""Gated online RFPerm + PO-risk."""
from __future__ import annotations

import numpy as np

from agod.online_rfperm import (
    hop_fires,
    last_two_sqrt_weights,
    run_rfperm_stream,
)
from agod.po_iptw import po_iptw_weights
from agod.po_refit import make_batch_stream, run_uniform_last_two


def test_online_rfperm_is_the_stream_entry():
    from agod import run_online_rfperm
    from agod.online_rfperm import run_rfperm_stream

    assert run_online_rfperm is run_rfperm_stream


def test_hop_fires_skips_first_and_needs_jump():
    assert hop_fires(2.0, None, gate=1.5) is False
    assert hop_fires(1.4, 1.0, gate=1.5) is False
    assert hop_fires(1.6, 1.0, gate=1.5) is True


def test_hop_fires_rejects_vacuous_zero_denominator():
    # Occupancy empty-room: e_prev=0, e_now large, γ is meaningless.
    assert hop_fires(0.80, 0.0, gate=1.5, e_floor=0.02) is False
    assert hop_fires(0.80, 0.01, gate=1.5, e_floor=0.02) is False
    assert hop_fires(0.80, 0.40, gate=1.5, e_floor=0.02) is True


def test_sqrt_po_weights_mean_one():
    rng = np.random.default_rng(0)
    po = rng.uniform(0.2, 4.0, size=40)
    w = po_iptw_weights(po, mode="sqrt")
    assert abs(float(w.mean()) - 1.0) < 0.08
    assert np.all(w > 0)


def test_last_two_reweight_keeps_control_near_one():
    batch = np.repeat([0, 1, 2], 10)
    rng = np.random.default_rng(1)
    po1 = rng.uniform(0.2, 5.0, size=10)
    tr, w = last_two_sqrt_weights(batch, 2, po1)
    assert int(tr.sum()) == 20
    assert abs(float(w.mean()) - 1.0) < 0.08
    treated = batch[tr] == 2
    # T=0 stays flatter; T=1 carries the √PO shape
    assert float(w[treated].std()) >= float(w[~treated].std())


def test_similar_stream_rfperm_rarely_fires():
    stream = make_batch_stream(n_batches=6, n_per=80, p=8, seed=0, cov=0.0)
    rec = run_rfperm_stream(stream, gate=1.5)
    uni = run_uniform_last_two(stream)
    assert rec["fire_rate"] <= 0.25
    assert abs(uni["online_mse"] - rec["online_mse"]) < 1e-9


def test_constant_label_stretch_does_not_vacuous_fire():
    from agod.po_refit import Stream

    rng = np.random.default_rng(0)
    n_per, n_batches, p = 40, 6, 6
    X = rng.normal(size=(n_per * n_batches, p))
    y = np.zeros(n_per * n_batches, dtype=int)
    y[3 * n_per :] = 1
    batch = np.repeat(np.arange(n_batches), n_per)
    stream = Stream(X=X, y=y, batch=batch, name="flip", task="acc")
    rec = run_rfperm_stream(stream, gate=1.5)
    # The first post-flip hop has e_prev=0 on the constant stretch.
    hop = next(h for h in rec["history"] if h["t"] == 3)
    assert hop["mean_r0"] < 0.02
    assert hop["fired"] is False


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
    assert abs(hop3["mean_w"] - 1.0) < 0.08


def test_observation_weights_rank_po_on_treated():
    jumped = make_batch_stream(
        n_batches=6, n_per=80, p=8, seed=4, cov=0.0, concept_at=3, concept=1.0
    )
    rec = run_rfperm_stream(jumped, gate=1.5, detail=True)
    hop = next(h for h in rec["history"] if h["fired"])
    po = np.asarray(hop["po_risk0"], float)
    w = np.asarray(hop["w"], float)
    treated = np.asarray(hop["treated"], int) == 1
    assert po.shape == w.shape
    assert np.isfinite(po).all() and np.isfinite(w).all()
    assert np.all(w > 0)
    # T=1: √PO is rank-preserving before clip; after mean-1 still monotone
    order = np.argsort(po[treated])
    assert np.all(np.diff(w[treated][order]) >= -1e-9)
    assert hop["po_t1"]["p90"] >= hop["po_t1"]["p50"]
    assert hop["w_t1"]["p90"] >= hop["w_t1"]["p50"]


def test_quiet_observations_have_unit_weights():
    stream = make_batch_stream(n_batches=6, n_per=80, p=8, seed=0, cov=0.0)
    rec = run_rfperm_stream(stream, gate=1.5, detail=True)
    quiet = [h for h in rec["history"] if not h["fired"]]
    assert quiet
    for h in quiet:
        w = np.asarray(h["w"], float)
        assert np.allclose(w, 1.0)
        assert np.isfinite(np.asarray(h["po_risk0"], float)).all()


def test_dre_last_two_same_rows_and_ignores_concept():
    from agod.po_refit import run_dre_last_two

    jumped = make_batch_stream(
        n_batches=6, n_per=80, p=8, seed=4, cov=0.0, concept_at=3, concept=1.0
    )
    uni = run_uniform_last_two(jumped)
    dre = run_dre_last_two(jumped)
    po = run_rfperm_stream(jumped, gate=1.5)
    assert dre["n_hops"] == uni["n_hops"]
    for h, g in zip(dre["history"], uni["history"]):
        assert h["n_train"] == g["n_train"]
    # P(X) is unchanged, so DRE cannot systematically beat uniform.
    assert dre["online_mse"] >= uni["online_mse"] - 0.08
    assert po["fire_rate"] > 0.0


def test_reweight_keeps_last_two_rows_when_fired():
    jumped = make_batch_stream(
        n_batches=6, n_per=80, p=8, seed=4, cov=0.0, concept_at=3, concept=1.0
    )
    rec = run_rfperm_stream(jumped, gate=1.5)
    uni = run_uniform_last_two(jumped)
    for h, g in zip(rec["history"], uni["history"]):
        assert h["n_train"] == g["n_train"]
    fired = [h for h in rec["history"] if h["fired"]]
    assert fired
    assert all(abs(h["mean_w"] - 1.0) < 0.08 for h in fired)
    assert all("po_t1" in h and "w_t1" in h for h in rec["history"])


def test_fit_predict_single_class_does_not_crash():
    from agod.po_refit import fit_predict

    rng = np.random.default_rng(0)
    X = rng.normal(size=(20, 4))
    y = np.ones(20, dtype=int)
    for kind in ("rf", "xgb"):
        pred = fit_predict(X, y, X[:5], learner=kind, seed=0, task="acc")
        assert pred.shape == (5,)
        assert np.all(pred == 1)
