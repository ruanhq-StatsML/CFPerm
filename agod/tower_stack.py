"""Online stacking of three-tower votes. Snapshot GLS is not π(t).

The three-tower prototype GLS-packs (cat_match, s_pos, s_neg) independently
each step. With y≡1 the oracle-axis snapshot barely moves. This module is
the missing prequential loop:

    freeze portrait AND freeze π_t
    votes from (snap, item, pool)     # no y
    honest_step(stacker, votes, y)    # π_{t+1}
    sticky mix → conf-damped LR       # FWD on
    THEN update portrait / seen

Category switch (Tools → Sports, same price band) is the recsys analogue
of the half-stream expert switch: cat_hist lags, the user vote drops,
Hedge+Fixed-Share spends path TV leaving user; discrete OSL is a vertex;
the per-step GLS snapshot barely moves. A leaky stacker scores in-sample
and looks optimistic. Leaky portrait counts the current event (cold start
cat_match 0→1).

MMoE on (U,I,N) is π(x) — contrast only.
"""
from __future__ import annotations

from typing import Mapping

import numpy as np

from .online_portrait import OnlinePortrait, votes_from_snapshot
from .online_stacking import OnlineStacker, _simplex, honest_step, leaky_step, pi_to_lr
from .stack_track import path_tv, sticky_mix, switch_delay, tv
from .three_tower import (
    TOWERS,
    NegPool,
    agod_from_triplet,
    cat_match,
    count_confidence,
    score_triplet,
)

LEAK_MODES = ("none", "portrait", "stacker", "both")


def stack_votes(snap: Mapping, pos: Mapping, trip: Mapping) -> dict[str, float]:
    """Scalar votes predicting y. Frozen snapshot only; y does not enter.

    user  — cat_match = ⟨cat_hist, onehot(pos.cat)⟩  (dies on a cat switch)
    item  — price affinity, category-agnostic item-tower coordinate
    neg   — 1 − mean cos(U, N): easy pool looks good, hard pool does not
    """
    price = float(votes_from_snapshot(snap, pos)["price"])
    s_neg = float(trip.get("s_neg") or 0.0)
    return {
        "user": float(np.clip(cat_match(snap, pos), 0.0, 1.0)),
        "item": float(np.clip(price, 0.0, 1.0)),
        "neg": float(np.clip(1.0 - s_neg, 0.0, 1.0)),
    }


def mmoe_towers(trip: Mapping, weight: np.ndarray | None = None) -> dict[str, float]:
    """Their head on the three embeddings: π(x)=softmax(W mean(U,I,N))."""
    e = np.stack(
        [np.asarray(trip["u"], float), np.asarray(trip["i"], float), np.asarray(trip["n"], float)]
    )
    n, d = e.shape
    w = np.asarray(weight, float) if weight is not None else np.random.default_rng(0).normal(
        scale=0.15, size=(n, d)
    )
    logits = w @ e.mean(axis=0)
    logits = logits - float(logits.max())
    pi = _simplex(np.exp(logits))
    return {TOWERS[i]: float(pi[i]) for i in range(3)}


def _sample_negs(pool: NegPool, snap: Mapping, seen: set[str], rng, k_neg: int, hard: bool) -> list:
    top_cat = None
    if hard and snap.get("cat_hist"):
        top_cat = max(snap["cat_hist"], key=snap["cat_hist"].get)
    negs = pool.sample(k_neg, seen, rng, hard_cat=top_cat if hard else None)
    if hard and len(negs) < k_neg:
        negs = negs + pool.sample(k_neg - len(negs), seen, rng)
    return negs


def tower_stack_step(
    portrait: OnlinePortrait,
    stacker: OnlineStacker,
    pos: Mapping,
    pool: NegPool,
    seen: set[str],
    y: float,
    rng: np.random.Generator,
    *,
    leak: str = "none",
    k_neg: int = 4,
    hard: bool = False,
    sticky_pi: dict[str, float] | None = None,
    lam: float = 0.25,
) -> dict:
    """One step. leak='none' is emit-then-update + honest π. Else optimistic."""
    if leak not in LEAK_MODES:
        raise ValueError(f"leak must be one of {LEAK_MODES}")
    cat = str(pos.get("category") or pos.get("main_category") or "UNK")
    pid = str(pos.get("parent_asin") or pos.get("id") or "")

    if leak in ("portrait", "both"):
        portrait.update(cat, pos.get("price"))
        snap = portrait.snapshot()
    else:
        snap = portrait.snapshot()

    negs = _sample_negs(pool, snap, seen, rng, k_neg, hard)
    trip = score_triplet(snap, pos, negs)
    votes = stack_votes(snap, pos, trip)
    y = float(y)
    if leak in ("stacker", "both"):
        scored = leaky_step(stacker, votes, y)
    else:
        scored = honest_step(stacker, votes, y)

    gls = agod_from_triplet(trip, y)
    conf = count_confidence(snap)
    pi_fast = scored["pi"]
    if sticky_pi is None:
        sticky_pi = dict(pi_fast)
    pi_act = sticky_mix(sticky_pi, pi_fast, TOWERS, lam=lam)
    lr = pi_to_lr(pi_act, TOWERS, beta=0.10)
    lr_damped = {k: 1.0 + conf * (v - 1.0) for k, v in lr.items()}
    moe = mmoe_towers(trip)

    if leak not in ("portrait", "both"):
        portrait.update(cat, pos.get("price"))
    if pid:
        seen.add(pid)

    return {
        "snap": snap,
        "trip": {k: trip[k] for k in ("s_pos", "s_neg", "gap", "cat_match_pos", "cat_match_neg", "conf")},
        "votes": votes,
        "pi": pi_fast,
        "pi_act": pi_act,
        "pi_gls": gls["pi"],
        "pi_mmoe": moe,
        "lr": lr_damped,
        "preq_loss": scored.get("preq_loss"),
        "n_eff": scored.get("n_eff"),
        "review_cnt": snap["review_cnt"],
        "leak": leak,
        "category": cat,
    }


def switch_events(n_each: int = 16, *, price: float = 20.0) -> tuple[list[dict], int]:
    """Tools positives, then Sports positives at the same price band."""
    n_each = int(n_each)
    hist = [
        {"parent_asin": f"t{i}", "category": "Tools", "price": price + 0.3 * i, "y": 1.0}
        for i in range(n_each)
    ]
    hist += [
        {"parent_asin": f"s{i}", "category": "Sports", "price": price + 0.3 * i, "y": 1.0}
        for i in range(n_each)
    ]
    return hist, n_each


def run_tower_stack(
    events: list[dict],
    pool: NegPool,
    *,
    method: str = "hedge",
    leak: str = "none",
    seed: int = 0,
    share: float = 0.08,
    eta: float = 0.85,
    lam: float = 0.25,
) -> dict:
    rng = np.random.default_rng(seed)
    portrait = OnlinePortrait()
    stacker = OnlineStacker(TOWERS, method=method, share=share, eta=eta)
    seen: set[str] = set()
    sticky = {e: 1.0 / len(TOWERS) for e in TOWERS}
    traj = []
    for ev in events:
        pos = {k: ev[k] for k in ("parent_asin", "category", "price")}
        step = tower_stack_step(
            portrait,
            stacker,
            pos,
            pool,
            seen,
            float(ev["y"]),
            rng,
            leak=leak,
            sticky_pi=sticky,
            lam=lam,
        )
        sticky = step["pi_act"]
        traj.append(step)
    return {
        "traj": traj,
        "final_snap": portrait.snapshot(),
        "method": method,
        "leak": leak,
        "seen": len(seen),
        "path_tv": path_tv([r["pi"] for r in traj]),
        "path_tv_act": path_tv([r["pi_act"] for r in traj]),
        "path_tv_gls": path_tv([r["pi_gls"] for r in traj]),
        "mean_preq": float(np.mean([r["preq_loss"] for r in traj])),
    }


def switch_report(payload: dict, t_switch: int) -> dict:
    """How π moves after the Tools → Sports cut."""
    pis = [r["pi"] for r in payload["traj"]]
    delay = switch_delay(pis, new_leader="item", t_switch=t_switch, thresh=0.40)
    pre = payload["traj"][:t_switch]
    post = payload["traj"][t_switch:]
    def _mean_pi(rows, k):
        if not rows:
            return float("nan")
        return float(np.mean([r["pi"][k] for r in rows]))

    def _mean_vote(rows, k):
        if not rows:
            return float("nan")
        return float(np.mean([r["votes"][k] for r in rows]))

    return {
        "delay_item": delay,
        "pi_user_pre": _mean_pi(pre, "user"),
        "pi_user_post": _mean_pi(post, "user"),
        "pi_item_pre": _mean_pi(pre, "item"),
        "pi_item_post": _mean_pi(post, "item"),
        "vote_user_pre": _mean_vote(pre, "user"),
        "vote_user_post": _mean_vote(post, "user"),
        "vote_item_pre": _mean_vote(pre, "item"),
        "vote_item_post": _mean_vote(post, "item"),
        "preq_pre": float(np.mean([r["preq_loss"] for r in pre])) if pre else float("nan"),
        "preq_post": float(np.mean([r["preq_loss"] for r in post])) if post else float("nan"),
        "path_tv": payload["path_tv"],
        "path_tv_gls": payload["path_tv_gls"],
        "tv_at_switch": tv(pis[t_switch - 1], pis[t_switch]) if t_switch < len(pis) and t_switch > 0 else 0.0,
    }
