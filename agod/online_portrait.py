"""Online causal portrait → inference alignment. No feature store.

Running state is four Welford / count statistics, nothing else:

    cat_hist     past category frequencies
    review_cnt   number of past events
    price_mean   past price μ
    price_std    past price σ

Protocol (same honesty as stacking, same causality as the recsys Beam DoFn):

    snap = portrait.snapshot()          # freeze
    votes / vecs from (snap, item)      # item catalog only; not y_t
    π = direction_match(vecs, g_hold)   # inference alignment
    LR = pi_to_lr(π)
    portrait.update(item)               # THEN fold the event in
    stacker.update(...)

MMoE is shown only as a contrast: π(x)=softmax(Wx) on expert vecs.
Alignment π is a function of holdout scores, not of x.
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from .online_stacking import (
    EPS,
    OnlineStacker,
    _simplex,
    direction_match_weights,
    honest_step,
    pi_to_lr,
)

EXPERTS = ("cat", "count", "price")


def _welford_update(n: int, mean: float, m2: float, x: float) -> tuple[int, float, float]:
    n1 = n + 1
    delta = x - mean
    mean1 = mean + delta / n1
    return n1, mean1, m2 + delta * (x - mean1)


def _welford_std(n: int, m2: float) -> float:
    if n <= 1:
        return 0.0
    return float((m2 / (n - 1)) ** 0.5)


class OnlinePortrait:
    """Per-user running cat / count / price moments. Emit-then-update."""

    def __init__(self):
        self._cat_cnt: dict[str, int] = {}
        self._n = 0
        self._price_n = 0
        self._price_mean = 0.0
        self._price_m2 = 0.0

    def snapshot(self) -> dict:
        total = self._n
        cat_hist = (
            {k: round(v / total, 6) for k, v in self._cat_cnt.items()} if total else {}
        )
        return {
            "cat_hist": dict(cat_hist),
            "review_cnt": int(self._n),
            "price_mean": None if self._price_n == 0 else float(self._price_mean),
            "price_std": _welford_std(self._price_n, self._price_m2),
        }

    def update(self, category: str, price: float | None) -> dict:
        cat = category or "UNK"
        self._cat_cnt[cat] = self._cat_cnt.get(cat, 0) + 1
        self._n += 1
        if price is not None:
            try:
                px = float(price)
            except (TypeError, ValueError):
                px = None
            if px is not None:
                self._price_n, self._price_mean, self._price_m2 = _welford_update(
                    self._price_n, self._price_mean, self._price_m2, px
                )
        return self.snapshot()


def votes_from_snapshot(snap: Mapping, item: Mapping) -> dict[str, float]:
    """Scalar expert votes from frozen portrait + current item (not y).

    cat    — past frequency of this item's category
    count  — 1 − exp(−n/5), how much history we have
    price  — Gaussian affinity of item price to past μ/σ
    """
    cat = str(item.get("category") or item.get("main_category") or "UNK")
    hist = snap.get("cat_hist") or {}
    n = float(snap.get("review_cnt") or 0)
    mu = snap.get("price_mean")
    sd = float(snap.get("price_std") or 0.0)
    px = item.get("price")
    try:
        px = float(px) if px is not None else None
    except (TypeError, ValueError):
        px = None

    cat_v = float(hist.get(cat, 0.0))
    count_v = float(1.0 - np.exp(-n / 5.0))
    if mu is None or px is None:
        price_v = 0.5
    else:
        z = abs(px - float(mu)) / max(sd, 1.0)
        price_v = float(np.exp(-0.5 * z * z))
    return {"cat": cat_v, "count": count_v, "price": price_v}


def vote_vectors(snap: Mapping, item: Mapping, *, dim: int = 8) -> dict[str, np.ndarray]:
    """Pack the same three votes as directions in R^d for GLS alignment."""
    v = votes_from_snapshot(snap, item)
    out = {}
    for i, name in enumerate(EXPERTS):
        vec = np.zeros(dim, float)
        vec[i] = 1.0
        vec[i + 3] = float(v[name])
        out[name] = vec
    return out


def holdout_target(item: Mapping, y: float, *, dim: int = 8) -> np.ndarray:
    """Holdout direction: which expert *should* have fired, from y and item.

    This is the Super Learner target, not a feature. y and the current item
    may enter here; they must not enter ``votes_from_snapshot``.
    """
    g = np.zeros(dim, float)
    y = float(y)
    v = votes_from_snapshot(
        {
            "cat_hist": {str(item.get("category") or "UNK"): 1.0},
            "review_cnt": 5,
            "price_mean": item.get("price"),
            "price_std": 1.0,
        },
        item,
    )
    # y>0.5 → want the experts that matched the item; else the orthogonal leftover
    if y >= 0.5:
        g[0] = 1.0
        g[3] = v["cat"]
        g[1] = 0.35
        g[4] = v["count"]
    else:
        g[2] = 1.0
        g[5] = v["price"]
    n = float(np.linalg.norm(g))
    return g if n < EPS else g / n


def mmoe_gate(expert_vecs: np.ndarray, weight: np.ndarray | None = None) -> dict:
    """Their head, for contrast only: π(x)=softmax(W mean(experts)).

    ``expert_vecs`` is [N, D]. Query is the mean over experts (as in
    ``TwoTaskMMoE``). This depends on x, not on a holdout target.
    """
    e = np.asarray(expert_vecs, float)
    if e.ndim == 1:
        e = e.reshape(1, -1)
    n, d = e.shape
    query = e.mean(axis=0)
    w = np.zeros((n, d)) if weight is None else np.asarray(weight, float)
    if w.shape != (n, d):
        rng = np.random.default_rng(0)
        w = rng.normal(scale=0.15, size=(n, d))
    logits = w @ query
    logits = logits - logits.max()
    pi = _simplex(np.exp(logits))
    fused = pi @ e
    return {"pi": {EXPERTS[i]: float(pi[i]) for i in range(min(n, len(EXPERTS)))}, "fused": fused}


def align_from_snapshot(
    snap: Mapping,
    item: Mapping,
    y: float,
) -> dict:
    """Inference alignment: GLS / direction-match of frozen votes onto holdout."""
    vecs = vote_vectors(snap, item)
    g_hold = holdout_target(item, y)
    packed = direction_match_weights(vecs, g_hold, EXPERTS)
    lr = pi_to_lr(packed["pi"], EXPERTS, beta=0.10)
    return {
        "votes": votes_from_snapshot(snap, item),
        "pi": packed["pi"],
        "scores": packed["scores"],
        "match_mse": packed["match_mse"],
        "vertex": packed["vertex"],
        "lr": lr,
        "g_hold": g_hold,
    }


def honest_portrait_step(
    portrait: OnlinePortrait,
    stacker: OnlineStacker,
    item: Mapping,
    y: float,
) -> dict:
    """One landed step: freeze portrait → align π → LR → then update both."""
    snap = portrait.snapshot()
    aligned = align_from_snapshot(snap, item, y)
    votes = aligned["votes"]
    scored = honest_step(stacker, votes, float(y))
    portrait.update(str(item.get("category") or "UNK"), item.get("price"))
    return {
        "snap": snap,
        "votes": votes,
        "align_pi": aligned["pi"],
        "align_lr": aligned["lr"],
        "match_mse": aligned["match_mse"],
        "stack_pi": scored["pi"],
        "preq_loss": scored.get("preq_loss"),
        "n_eff": scored.get("n_eff"),
    }


def leaky_portrait_step(
    portrait: OnlinePortrait,
    stacker: OnlineStacker,
    item: Mapping,
    y: float,
) -> dict:
    """Wrong order: update portrait with the current event, then score."""
    portrait.update(str(item.get("category") or "UNK"), item.get("price"))
    snap = portrait.snapshot()
    aligned = align_from_snapshot(snap, item, y)
    from .online_stacking import leaky_step

    scored = leaky_step(stacker, aligned["votes"], float(y))
    return {
        "snap": snap,
        "votes": aligned["votes"],
        "align_pi": aligned["pi"],
        "preq_loss": scored.get("preq_loss"),
    }
