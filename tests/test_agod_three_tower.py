"""Three-tower + online portrait: cat_hist, review_cnt, neg pool, AGOD."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from agod.online_portrait import OnlinePortrait, honest_portrait_step, votes_from_snapshot
from agod.online_stacking import OnlineStacker
from agod.three_tower import (
    NegPool,
    cat_match,
    count_confidence,
    honest_three_tower_step,
    item_tower,
    pool_eval,
    score_triplet,
    user_tower,
)


def _item(pid, cat, price=20.0):
    return {"parent_asin": pid, "category": cat, "price": price}


def test_portrait_emit_then_update():
    p = OnlinePortrait()
    snap0 = p.snapshot()
    assert snap0["review_cnt"] == 0
    assert snap0["cat_hist"] == {}
    p.update("Tools", 10.0)
    snap1 = p.snapshot()
    assert snap1["review_cnt"] == 1
    assert snap1["cat_hist"]["Tools"] == 1.0
    p.update("Sports", 40.0)
    snap2 = p.snapshot()
    assert snap2["review_cnt"] == 2
    assert abs(snap2["cat_hist"]["Tools"] - 0.5) < 1e-9
    assert snap2["price_mean"] == 25.0


def test_votes_use_frozen_snapshot_not_current_y():
    p = OnlinePortrait()
    p.update("Tools", 20.0)
    snap = p.snapshot()
    v = votes_from_snapshot(snap, {"category": "Tools", "price": 22.0})
    assert v["cat"] > 0.99
    v2 = votes_from_snapshot(snap, {"category": "Fashion", "price": 22.0})
    assert v2["cat"] == 0.0


def test_cat_match_and_count_shrink():
    p = OnlinePortrait()
    assert count_confidence(p.snapshot()) < 0.05
    u0 = user_tower(p.snapshot())
    assert abs(u0[:5].max() - u0[:5].min()) < 0.05  # ~uniform
    for _ in range(12):
        p.update("Tools", 15.0)
    snap = p.snapshot()
    assert count_confidence(snap) > 0.85
    assert cat_match(snap, _item("x", "Tools")) > 0.99
    assert cat_match(snap, _item("y", "Fashion")) == 0.0
    u = user_tower(snap)
    assert u[0] > 0.7  # Tools mass after shrink


def test_neg_pool_excludes_seen():
    items = [_item(f"i{k}", "Tools" if k % 2 == 0 else "Sports") for k in range(20)]
    pool = NegPool(items)
    seen = {"i0", "i2", "i4"}
    rng = np.random.default_rng(0)
    negs = pool.sample(6, seen, rng)
    assert len(negs) == 6
    assert pool.contamination(negs, seen) == 0.0


def test_hard_negs_close_the_gap():
    p = OnlinePortrait()
    for _ in range(8):
        p.update("Tools", 20.0)
    snap = p.snapshot()
    catalog = [_item(f"t{k}", "Tools", 20.0) for k in range(8)] + [
        _item(f"s{k}", "Sports", 80.0) for k in range(8)
    ]
    pool = NegPool(catalog)
    seen = {"t0"}
    rng = np.random.default_rng(1)
    ev = pool_eval(snap, _item("t0", "Tools", 20.0), pool, seen, rng, k_neg=5)
    assert ev["seen_excluded"]["contam"] == 0.0
    assert ev["same_cat_hard"]["gap"] < ev["seen_excluded"]["gap"] + 0.05
    assert ev["same_cat_hard"]["cat_match_neg"] >= ev["seen_excluded"]["cat_match_neg"] - 1e-9


def test_honest_step_does_not_count_current():
    p = OnlinePortrait()
    stacker = OnlineStacker(("cat", "count", "price"), method="hedge", share=0.05)
    item = {"category": "Tools", "price": 12.0}
    out = honest_portrait_step(p, stacker, item, y=1.0)
    assert out["snap"]["review_cnt"] == 0
    assert p.snapshot()["review_cnt"] == 1


def test_cold_start_damps_lr_to_one():
    p = OnlinePortrait()
    catalog = [_item(f"a{k}", "Sports", 50.0) for k in range(10)]
    pool = NegPool(catalog)
    seen: set[str] = set()
    rng = np.random.default_rng(2)
    out = honest_three_tower_step(
        p, _item("pos", "Tools", 18.0), pool, seen, y=1.0, rng=rng, k_neg=3
    )
    assert out["review_cnt"] == 0
    assert out["trip"]["conf"] < 0.05
    for k in ("user", "item", "neg"):
        assert abs(out["lr"][k] - 1.0) < 1e-9


def test_user_score_tracks_cat_hist_match():
    """⟨cat_hist, pos.cat⟩ is the user-tower vote; it rises once hist matches."""
    from agod.three_tower import agod_from_triplet

    p = OnlinePortrait()
    cold = score_triplet(p.snapshot(), _item("p", "Tools", 20.0), [_item("n", "Fashion", 90.0)])
    for _ in range(10):
        p.update("Tools", 20.0)
    warm = score_triplet(p.snapshot(), _item("p", "Tools", 20.0), [_item("n", "Fashion", 90.0)])
    a_cold = agod_from_triplet(cold, 1.0)
    a_warm = agod_from_triplet(warm, 1.0)
    assert a_warm["scores"]["user"] > a_cold["scores"]["user"] + 0.1
    assert a_warm["pi"]["neg"] < 0.05
    assert a_cold["pi"]["neg"] < 0.05


def test_three_tower_agod_lr_positive():
    p = OnlinePortrait()
    for _ in range(6):
        p.update("Tools", 18.0)
    catalog = [_item(f"a{k}", "Sports", 50.0) for k in range(10)]
    pool = NegPool(catalog)
    seen: set[str] = set()
    rng = np.random.default_rng(2)
    out = honest_three_tower_step(
        p, _item("pos", "Tools", 18.0), pool, seen, y=1.0, rng=rng, k_neg=3
    )
    assert out["review_cnt"] == 6
    assert all(v > 0 for v in out["lr"].values())
    assert abs(sum(out["pi"].values()) - 1.0) < 1e-6
    assert "pos" in seen


def test_score_gap_positive_when_hist_matches():
    p = OnlinePortrait()
    for _ in range(10):
        p.update("Tools", 20.0)
    snap = p.snapshot()
    pos = _item("p", "Tools", 20.0)
    negs = [_item("n", "Fashion", 90.0)]
    trip = score_triplet(snap, pos, negs)
    assert trip["gap"] > 0.05
    assert trip["cat_match_pos"] > trip["cat_match_neg"]


if __name__ == "__main__":
    for fn in list(globals().values()):
        if callable(fn) and fn.__name__.startswith("test_"):
            fn()
            print(f"  {fn.__name__}: OK")
    print("test_agod_three_tower: OK")
