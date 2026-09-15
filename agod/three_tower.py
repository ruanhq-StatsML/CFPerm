"""Three-tower prototype: causal portrait, item, negative pool → AGOD align.

Towers (not MMoE, not a feature store):

    U  user   — frozen cat_hist + review_cnt (+ price μ)
    I  item   — current catalog category + price
    N  neg    — same item tower on a seen-excluded pool (their CausalPosNeg)

Score is cosine(U, I) vs cosine(U, N). AGOD does *not* GLS the raw
embeddings — that vertices on item. Votes are packed scalars
(cat_match, s_pos, s_neg) on orthogonal axes; holdout is y-only;
then π → conf-damped LR. Portrait is emit-then-update.

How the two stats are characterized (this is the point):

    cat_hist   user direction on the category simplex;
               match = ⟨cat_hist, onehot(item.cat)⟩
    review_cnt confidence / shrink: n→0 mixes U toward uniform,
               so AGOD must not crank LR on a cold portrait
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from .online_portrait import OnlinePortrait
from .online_stacking import EPS, direction_match_weights, pi_to_lr, _unit

CATS = ("Tools", "Sports", "Fashion", "Electronics", "UNK")
TOWERS = ("user", "item", "neg")
DIM = 8  # 5 cats + n_feat + price + pad


def _cat_index(name: str) -> int:
    n = str(name or "UNK")
    return CATS.index(n) if n in CATS else CATS.index("UNK")


def cat_match(snap: Mapping, item: Mapping) -> float:
    """⟨cat_hist, onehot(item.cat)⟩ — past frequency of this category."""
    cat = str(item.get("category") or item.get("main_category") or "UNK")
    return float((snap.get("cat_hist") or {}).get(cat, 0.0))


def count_confidence(snap: Mapping, *, n0: float = 5.0) -> float:
    """review_cnt → [0, 1] shrink weight. 0 = uniform prior, 1 = trust hist."""
    n = float(snap.get("review_cnt") or 0)
    return float(1.0 - np.exp(-n / max(n0, 1e-6)))


def user_tower(snap: Mapping, *, n0: float = 5.0) -> np.ndarray:
    """U: shrunk category simplex + log-count + price μ. Frozen snapshot only."""
    conf = count_confidence(snap, n0=n0)
    hist = snap.get("cat_hist") or {}
    h = np.array([float(hist.get(c, 0.0)) for c in CATS], float)
    if h.sum() <= EPS:
        h = np.full(len(CATS), 1.0 / len(CATS))
    else:
        h = h / h.sum()
    unif = np.full(len(CATS), 1.0 / len(CATS))
    h = (1.0 - conf) * unif + conf * h
    vec = np.zeros(DIM, float)
    vec[: len(CATS)] = h
    vec[len(CATS)] = conf
    mu = snap.get("price_mean")
    vec[len(CATS) + 1] = 0.0 if mu is None else float(mu) / 100.0
    return vec


def item_tower(item: Mapping) -> np.ndarray:
    """I: one-hot category + price. Catalog, not a rating."""
    vec = np.zeros(DIM, float)
    vec[_cat_index(item.get("category") or item.get("main_category"))] = 1.0
    try:
        px = float(item.get("price")) if item.get("price") is not None else 0.0
    except (TypeError, ValueError):
        px = 0.0
    vec[len(CATS) + 1] = px / 100.0
    return vec


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.dot(_unit(a), _unit(b)))


class NegPool:
    """Seen-excluded item pool. Global catalog is the hole they still have."""

    def __init__(self, items: Sequence[Mapping]):
        self.items = [dict(x) for x in items]

    def sample(
        self,
        k: int,
        seen: set[str],
        rng: np.random.Generator,
        *,
        hard_cat: str | None = None,
    ) -> list[dict]:
        k = int(max(k, 0))
        out: list[dict] = []
        idxs = rng.permutation(len(self.items))
        for j in idxs:
            it = self.items[int(j)]
            pid = str(it.get("parent_asin") or it.get("id") or "")
            if pid and pid in seen:
                continue
            if hard_cat is not None:
                cat = str(it.get("category") or it.get("main_category") or "")
                if cat != hard_cat:
                    continue
            out.append(it)
            if len(out) >= k:
                break
        return out

    def contamination(self, picked: Sequence[Mapping], seen: set[str]) -> float:
        if not picked:
            return 0.0
        hits = 0
        for it in picked:
            pid = str(it.get("parent_asin") or it.get("id") or "")
            hits += int(bool(pid) and pid in seen)
        return float(hits / len(picked))


def score_triplet(
    snap: Mapping,
    pos: Mapping,
    negs: Sequence[Mapping],
    *,
    n0: float = 5.0,
) -> dict:
    """cos(U,I) vs mean cos(U,N_k). Portrait frozen."""
    u = user_tower(snap, n0=n0)
    i = item_tower(pos)
    s_pos = cosine(u, i)
    nvecs = [item_tower(n) for n in negs] if negs else [np.zeros(DIM)]
    s_negs = [cosine(u, nv) for nv in nvecs]
    s_neg = float(np.mean(s_negs))
    n_bar = np.mean(np.stack(nvecs, axis=0), axis=0)
    return {
        "u": u,
        "i": i,
        "n": n_bar,
        "s_pos": s_pos,
        "s_neg": s_neg,
        "gap": float(s_pos - s_neg),
        "cat_match_pos": cat_match(snap, pos),
        "cat_match_neg": float(np.mean([cat_match(snap, n) for n in negs])) if negs else 0.0,
        "conf": count_confidence(snap, n0=n0),
    }


def tower_vote_vectors(trip: Mapping) -> dict[str, np.ndarray]:
    """Pack tower scalars onto orthogonal axes so GLS can mix them.

    Raw embeddings {U,I,N} cannot be the GLS library if ĝ_hold lives in
    their span (item becomes the vertex — linear gain). Votes here are
    the *characterization* of each tower, and do not see y:

        user  cat_match_pos = ⟨cat_hist, onehot(pos.cat)⟩
        item  s_pos         = cos(U, I)
        neg   s_neg         = mean cos(U, N_k)
    """
    vals = {
        "user": float(trip.get("cat_match_pos") or 0.0),
        "item": float(trip.get("s_pos") or 0.0),
        "neg": float(trip.get("s_neg") or 0.0),
    }
    out = {}
    for i, name in enumerate(TOWERS):
        vec = np.zeros(DIM, float)
        vec[i] = 1.0
        vec[i + 3] = vals[name]
        out[name] = vec
    return out


def holdout_from_triplet(trip: Mapping, y: float) -> np.ndarray:
    """Oracle holdout direction: which tower *should* have fired.

    Super Learner target — may see y. Must not be a reconstruction of
    {U, I, N} (that vertices on item) and must not copy the vote scalars
    (that is alignment leak). y=1 wants user/item axes; y=0 wants neg.
    """
    del trip  # votes live on the expert side; target is y-only
    g = np.zeros(DIM, float)
    if float(y) >= 0.5:
        g[0] = 1.0
        g[3] = 1.0
        g[1] = 0.35
        g[4] = 0.70
    else:
        g[2] = 1.0
        g[5] = 1.0
    return _unit(g)


def agod_from_triplet(trip: Mapping, y: float) -> dict:
    """Packed tower votes → GLS onto holdout → π → conf-damped LR."""
    vecs = tower_vote_vectors(trip)
    g_hold = holdout_from_triplet(trip, y)
    packed = direction_match_weights(vecs, g_hold, TOWERS)
    lr = pi_to_lr(packed["pi"], TOWERS, beta=0.10)
    conf = float(trip.get("conf") or 0.0)
    lr_damped = {k: 1.0 + conf * (v - 1.0) for k, v in lr.items()}
    return {
        "pi": packed["pi"],
        "scores": packed["scores"],
        "match_mse": packed["match_mse"],
        "vertex": packed["vertex"],
        "lr": lr,
        "lr_damped": lr_damped,
        "conf": conf,
    }


def honest_three_tower_step(
    portrait: OnlinePortrait,
    pos: Mapping,
    pool: NegPool,
    seen: set[str],
    y: float,
    rng: np.random.Generator,
    *,
    k_neg: int = 4,
    hard: bool = False,
) -> dict:
    """Freeze portrait → score pos/neg → AGOD π/LR → then update portrait."""
    snap = portrait.snapshot()
    top_cat = None
    if hard and snap.get("cat_hist"):
        top_cat = max(snap["cat_hist"], key=snap["cat_hist"].get)
    negs = pool.sample(k_neg, seen, rng, hard_cat=top_cat if hard else None)
    if hard and len(negs) < k_neg:
        negs = negs + pool.sample(k_neg - len(negs), seen, rng)
    trip = score_triplet(snap, pos, negs)
    agod = agod_from_triplet(trip, y)
    portrait.update(str(pos.get("category") or pos.get("main_category") or "UNK"), pos.get("price"))
    pid = str(pos.get("parent_asin") or pos.get("id") or "")
    if pid:
        seen.add(pid)
    return {
        "snap": snap,
        "trip": {k: trip[k] for k in ("s_pos", "s_neg", "gap", "cat_match_pos", "cat_match_neg", "conf")},
        "pi": agod["pi"],
        "scores": agod["scores"],
        "lr": agod["lr_damped"],
        "match_mse": agod["match_mse"],
        "contam": pool.contamination(negs, seen - {pid} if pid else seen),
        "n_neg": len(negs),
        "review_cnt": snap["review_cnt"],
    }


def pool_eval(
    snap: Mapping,
    pos: Mapping,
    pool: NegPool,
    seen: set[str],
    rng: np.random.Generator,
    *,
    k_neg: int = 8,
) -> dict:
    """Random vs seen-excluded vs same-cat hard negatives."""
    top_cat = None
    if snap.get("cat_hist"):
        top_cat = max(snap["cat_hist"], key=snap["cat_hist"].get)
    random_all = pool.sample(k_neg, set(), rng)
    clean = pool.sample(k_neg, seen, rng)
    hard = pool.sample(k_neg, seen, rng, hard_cat=top_cat)
    rows = {}
    for name, negs, seen_set in (
        ("random_global", random_all, set()),
        ("seen_excluded", clean, seen),
        ("same_cat_hard", hard, seen),
    ):
        trip = score_triplet(snap, pos, negs)
        rows[name] = {
            "gap": trip["gap"],
            "s_pos": trip["s_pos"],
            "s_neg": trip["s_neg"],
            "cat_match_pos": trip["cat_match_pos"],
            "cat_match_neg": trip["cat_match_neg"],
            "contam": pool.contamination(negs, seen),
            "n": len(negs),
        }
    return rows
