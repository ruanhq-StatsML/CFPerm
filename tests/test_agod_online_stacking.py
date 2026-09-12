"""Unit tests for modality-agnostic online stacking (no GPU / shards)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from agod.online_stacking import (
    OnlineStacker,
    honest_step,
    leaky_step,
    oracle_convex_combo,
    pi_to_lr,
    project_simplex,
    best_expert_loss,
    directional_scores,
    direction_match_weights,
    linear_gain_is_vertex,
)


def _switch_stream(t: int = 70, seed: int = 0):
    rng = np.random.default_rng(seed)
    y = np.zeros(t)
    va, vb, vc = [], [], []
    mid = t // 2
    for i in range(t):
        if i < mid:
            va.append(rng.normal(0.0, 0.12))
            vb.append(rng.normal(0.0, 1.10))
        else:
            va.append(rng.normal(0.0, 1.10))
            vb.append(rng.normal(0.0, 0.12))
        vc.append(rng.normal(0.0, 0.70))
    votes = np.column_stack([va, vb, vc])
    return votes, y, ("a", "b", "c")


def test_project_simplex_sums_to_one():
    rng = np.random.default_rng(1)
    for _ in range(8):
        v = rng.normal(size=4)
        w = project_simplex(v)
        assert abs(w.sum() - 1.0) < 1e-8
        assert np.all(w >= -1e-12)


def test_honest_hedge_tracks_switch():
    votes, y, names = _switch_stream(t=100, seed=0)
    st = OnlineStacker(names, method="hedge", eta=0.85, share=0.08)
    for t in range(len(y)):
        honest_step(st, {n: float(votes[t, i]) for i, n in enumerate(names)}, float(y[t]))
    pi = st.pi()
    assert pi["b"] > pi["a"]
    assert pi["b"] > 0.45


def test_osl_disc_picks_an_expert():
    votes, y, names = _switch_stream(t=40, seed=2)
    st = OnlineStacker(names, method="osl_disc")
    for t in range(len(y)):
        honest_step(st, {n: float(votes[t, i]) for i, n in enumerate(names)}, float(y[t]))
    pi = st.pi()
    assert abs(max(pi.values()) - 1.0) < 1e-9


def test_bg_downweights_noisy_expert():
    rng = np.random.default_rng(3)
    names = ("good", "bad")
    st = OnlineStacker(names, method="bg", ewma=0.25)
    for _ in range(40):
        votes = {"good": float(rng.normal(0, 0.1)), "bad": float(rng.normal(0, 1.2))}
        honest_step(st, votes, 0.0)
    assert st.pi()["good"] > st.pi()["bad"]


def test_gls_reduces_clone_mass():
    rng = np.random.default_rng(4)
    names = ("e0", "e1", "e2")
    st_g = OnlineStacker(names, method="gls_ewma", ewma=0.15)
    st_e = OnlineStacker(names, method="equal")
    for _ in range(80):
        common = rng.normal(0, 0.45)
        votes = {
            "e0": float(common + rng.normal(0, 0.04)),
            "e1": float(common + rng.normal(0, 0.04)),  # clone
            "e2": float(rng.normal(0, 0.28)),
        }
        honest_step(st_g, votes, 0.0)
        honest_step(st_e, votes, 0.0)
    clone_g = st_g.pi()["e0"] + st_g.pi()["e1"]
    clone_e = st_e.pi()["e0"] + st_e.pi()["e1"]
    assert clone_g < clone_e - 0.04
    assert st_g.pi()["e2"] > st_e.pi()["e2"]


def test_leaky_underestimates_preq_loss():
    votes, y, names = _switch_stream(t=90, seed=5)
    hon, leak = OnlineStacker(names, method="hedge", eta=0.4), OnlineStacker(
        names, method="hedge", eta=0.4
    )
    lh, ll = [], []
    for t in range(len(y)):
        v = {n: float(votes[t, i]) for i, n in enumerate(names)}
        lh.append(honest_step(hon, v, float(y[t]))["preq_loss"])
        ll.append(leaky_step(leak, v, float(y[t]))["preq_loss"])
    # leaky sees in-sample π; should not be worse on average
    assert float(np.mean(ll)) <= float(np.mean(lh)) + 0.02


def test_oracle_combo_beats_best_expert():
    rng = np.random.default_rng(6)
    t = 50
    y = rng.normal(size=t)
    votes = np.column_stack(
        [y + rng.normal(0, 0.4, t), y + rng.normal(0, 0.5, t), rng.normal(0, 1.2, t)]
    )
    be = best_expert_loss(votes, y)
    oc = oracle_convex_combo(votes, y)
    assert oc["loss"] <= be["loss"] + 1e-12


def test_pi_to_lr_positive():
    lr = pi_to_lr({"video": 0.7, "text": 0.2, "audio": 0.1}, ["video", "text", "audio"])
    assert all(lr[m] > 0 for m in ("video", "text", "audio"))
    assert lr["video"] > lr["audio"]


def test_linear_holdout_gain_is_a_vertex():
    scores = {"a": 0.91, "b": 0.40, "c": -0.10}
    assert linear_gain_is_vertex(scores, scores.keys()) == "a"


def test_direction_match_splits_complementary_grads():
    d = 8
    ga, gb, gc = np.zeros(d), np.zeros(d), np.zeros(d)
    ga[0] = 1.0
    gb[1] = 1.0
    gc[2] = 1.0
    hold = np.zeros(d)
    hold[0] = hold[1] = 1.0 / np.sqrt(2.0)
    packed = direction_match_weights(
        {"a": ga, "b": gb, "c": gc}, hold, ["a", "b", "c"]
    )
    assert packed["pi"]["a"] > 0.30 and packed["pi"]["b"] > 0.30
    assert packed["pi"]["c"] < 0.12
    assert packed["vertex"] in ("a", "b")
    # linear scores would pick a single axis; match must not be one-hot
    assert max(packed["pi"].values()) < 0.90


def test_directional_scores_cosine():
    g = {"a": np.array([1.0, 0.0]), "b": np.array([0.0, 1.0])}
    s = directional_scores(g, np.array([1.0, 0.0]), ["a", "b"])
    assert s["a"] > 0.99
    assert abs(s["b"]) < 0.05


if __name__ == "__main__":
    for fn in list(globals().values()):
        if callable(fn) and fn.__name__.startswith("test_"):
            fn()
            print(f"  {fn.__name__}: OK")
    print("test_agod_online_stacking: OK")
