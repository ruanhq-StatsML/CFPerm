"""Tracking, OOF features, and π variation."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from agod.stack_track import (
    RollingOSL,
    alignment_leak,
    disjoint_probe_holdout,
    oof_vs_leaky_stack,
    path_tv,
    run_switch_methods,
    sticky_mix,
    sticky_path,
    switch_delay,
    tv,
)


def _switch(t: int = 80, seed: int = 0):
    rng = np.random.default_rng(seed)
    y = np.zeros(t)
    mid = t // 2
    rows = []
    for i in range(t):
        if i < mid:
            a, b = rng.normal(0, 0.12), rng.normal(0, 1.10)
        else:
            a, b = rng.normal(0, 1.10), rng.normal(0, 0.12)
        rows.append((a, b, rng.normal(0, 0.70)))
    return np.array(rows), y, ("a", "b", "c"), mid


def test_tv_is_half_l1():
    assert abs(tv({"a": 1.0, "b": 0.0}, {"a": 0.0, "b": 1.0}) - 1.0) < 1e-12
    assert tv({"a": 0.5, "b": 0.5}, {"a": 0.5, "b": 0.5}) < 1e-12


def test_share_tracks_faster_than_share0():
    votes, y, names, mid = _switch(t=90, seed=1)
    packed = run_switch_methods(votes, y, names, t_switch=mid, eta=0.9)
    d0 = packed["hedge_share0"]["delay_from_switch"]
    ds = packed["hedge_share"]["delay_from_switch"]
    assert ds is not None
    if d0 is None:
        return
    assert ds <= d0


def test_osl_disc_lags_share_and_window():
    votes, y, names, mid = _switch(t=80, seed=2)
    packed = run_switch_methods(votes, y, names, t_switch=mid)
    b_disc = packed["osl_disc"]["final_pi"]["b"]
    assert packed["hedge_share"]["final_pi"]["b"] > b_disc + 0.05
    assert packed["osl_window"]["final_pi"]["b"] > b_disc + 0.05
    assert not packed["hedge_share"]["locked"]


def test_window_osl_forgets_first_half():
    names = ("a", "b")
    roll = RollingOSL(names, window=5)
    for _ in range(20):
        roll.update({"a": 0.01, "b": 1.0})
    assert roll.pi()["a"] > 0.99
    for _ in range(8):
        roll.update({"a": 1.0, "b": 0.01})
    assert roll.pi()["b"] > 0.99


def test_sticky_cuts_path_tv():
    votes, y, names, mid = _switch(t=70, seed=3)
    packed = run_switch_methods(votes, y, names, t_switch=mid)
    fast = packed["hedge_share"]["traj"]
    slow = sticky_path(fast, names, lam=0.12)
    assert path_tv(slow) < path_tv(fast) - 1e-6


def test_sticky_mix_convex():
    slow = {"a": 1.0, "b": 0.0}
    fast = {"a": 0.0, "b": 1.0}
    mid = sticky_mix(slow, fast, ("a", "b"), lam=0.25)
    assert abs(mid["b"] - 0.25) < 1e-9


def test_oof_less_optimistic_than_leaky():
    rng = np.random.default_rng(6)
    n, d = 120, 60
    xs = rng.normal(size=(n, 1))
    y = xs[:, 0] + 0.20 * rng.normal(size=n)
    xn = rng.normal(size=(n, d))
    out = oof_vs_leaky_stack(xs, xn, y, n_folds=4, lam=1e-3)
    assert out["pi_oof"]["good"] > 0.80
    assert out["corr_leak_noise"] > out["corr_oof_noise"] + 0.25
    assert out["pi_leak"]["noise"] > out["pi_oof"]["noise"]


def test_same_batch_alignment_is_leaky():
    rng = np.random.default_rng(5)
    g = rng.normal(size=24)
    probe, hold = g + 0.5 * rng.normal(size=24), g + 0.5 * rng.normal(size=24)
    a = alignment_leak(probe, hold, g)
    assert a["leaky_cos"] > 0.99
    assert a["honest_cos"] < 0.95
    p, h = disjoint_probe_holdout(40, probe_frac=0.5, rng=rng)
    assert len(set(p) & set(h)) == 0
    assert len(p) + len(h) == 40


def test_switch_delay_none_when_never_leads():
    traj = [{"a": 0.9, "b": 0.1} for _ in range(10)]
    assert switch_delay(traj, new_leader="b", t_switch=3, thresh=0.45) is None


if __name__ == "__main__":
    for fn in list(globals().values()):
        if callable(fn) and fn.__name__.startswith("test_"):
            fn()
            print(f"  {fn.__name__}: OK")
    print("test_agod_stack_track: OK")
