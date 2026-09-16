"""Online stacking of three-tower votes: switch, leak, snapshot GLS vs π(t)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from agod.online_portrait import OnlinePortrait
from agod.online_stacking import OnlineStacker
from agod.three_tower import CATS, NegPool
from agod.tower_stack import (
    run_tower_stack,
    stack_votes,
    switch_events,
    switch_report,
    tower_stack_step,
)


def _item(pid, cat, price=20.0, y=1.0):
    return {"parent_asin": pid, "category": cat, "price": price, "y": y}


def _pool(n=30, seed=0):
    rng = np.random.default_rng(seed)
    cats = list(CATS[:-1])
    items = []
    for k in range(n):
        items.append(
            {
                "parent_asin": f"c{k}",
                "category": cats[int(rng.integers(0, len(cats)))],
                "price": float(rng.uniform(12, 40)),
            }
        )
    return NegPool(items)


def test_stack_votes_ignore_y():
    p = OnlinePortrait()
    for _ in range(6):
        p.update("Tools", 20.0)
    snap = p.snapshot()
    from agod.three_tower import score_triplet

    trip = score_triplet(snap, _item("p", "Tools"), [_item("n", "Fashion", 80.0)])
    v = stack_votes(snap, _item("p", "Tools"), trip)
    assert v["user"] > 0.99
    v2 = stack_votes(snap, _item("q", "Fashion"), trip)
    assert v2["user"] == 0.0


def test_honest_first_event_has_empty_hist():
    p = OnlinePortrait()
    st = OnlineStacker(("user", "item", "neg"), method="hedge", share=0.08)
    seen: set[str] = set()
    rng = np.random.default_rng(0)
    out = tower_stack_step(
        p, st, _item("a", "Tools"), _pool(), seen, 1.0, rng, leak="none"
    )
    assert out["snap"]["review_cnt"] == 0
    assert out["votes"]["user"] == 0.0
    assert p.snapshot()["review_cnt"] == 1


def test_leaky_portrait_counts_current():
    p = OnlinePortrait()
    st = OnlineStacker(("user", "item", "neg"), method="hedge")
    seen: set[str] = set()
    rng = np.random.default_rng(1)
    out = tower_stack_step(
        p, st, _item("a", "Tools"), _pool(), seen, 1.0, rng, leak="portrait"
    )
    assert out["snap"]["review_cnt"] == 1
    assert out["votes"]["user"] == 1.0


def test_hedge_moves_to_item_after_cat_switch():
    events, mid = switch_events(18)
    pool = _pool(40, seed=4)
    hon = run_tower_stack(events, pool, method="hedge", leak="none", seed=4, share=0.10, eta=0.9)
    rep = switch_report(hon, mid)
    assert rep["vote_user_pre"] > 0.7
    assert rep["vote_user_post"] < rep["vote_user_pre"] - 0.3
    assert rep["pi_item_post"] > rep["pi_item_pre"]
    assert hon["path_tv"] > hon["path_tv_gls"]


def test_osl_disc_is_a_vertex():
    events, _mid = switch_events(12)
    pool = _pool(24, seed=5)
    disc = run_tower_stack(events, pool, method="osl_disc", leak="none", seed=5)
    last = disc["traj"][-1]["pi"]
    assert abs(max(last.values()) - 1.0) < 1e-9


def test_leaky_stacker_is_optimistic():
    events, _mid = switch_events(14)
    pool = _pool(24, seed=6)
    hon = run_tower_stack(events, pool, method="hedge", leak="none", seed=6)
    leak = run_tower_stack(events, pool, method="hedge", leak="stacker", seed=6)
    assert leak["mean_preq"] <= hon["mean_preq"] + 1e-12


if __name__ == "__main__":
    for fn in list(globals().values()):
        if callable(fn) and fn.__name__.startswith("test_"):
            fn()
            print(f"  {fn.__name__}: OK")
    print("test_agod_tower_stack: OK")
