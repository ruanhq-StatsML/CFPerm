"""Causal recsys features — the Amazon pipeline's OOF analogue.

Source: JingxiangQU/mmoe-multimodal-rec (``data4moe_beam.py``,
``data4model.py``), the WebDataset behind
``jingxiang11111/amazon_reviews_for_rec``.

The feature contract is time, not a fold index:

    for each user, sort events by timestamp
        user_feat ← aggregates of *past* reviews only
        emit (user_feat, item_meta, image, labels)
        *then* fold the current review into the running state

Current review title/text never enters the feature blob (it would make
``label_good`` a sentiment-of-this-review task). Item side is catalog
meta, not a target encoding of ratings. Split is by event date
(train ≤ 2023-06-30 < valid), not a random row split.

That is the recsys version of Super Learner OOF: the meta-features at t
must not have been computed from y_t. Random-split target encoding is
the leaky stacking analogue.

The *model* on top is MMoE: π(x)=softmax(Wx) over user/item/image
experts. Features are causal; the gate is MoE. Do not mix the two.

AGOD's Amazon smoke currently concatenates ``item\\nuser`` into one
hashing bag — it throws this contract away.
"""
from __future__ import annotations

from collections import deque
from datetime import date, datetime
from typing import Sequence

import numpy as np


TRAIN_END = date(2023, 6, 30)
VALID_END = date(2023, 9, 30)


def _try_float(x):
    if x in (None, ""):
        return None
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _welford_update(n, mean, m2, x):
    n1 = n + 1
    delta = x - mean
    mean1 = mean + delta / n1
    m2_1 = m2 + delta * (x - mean1)
    return n1, mean1, m2_1


def _welford_std(n, m2):
    if n <= 1:
        return 0.0
    return float((m2 / (n - 1)) ** 0.5)


def split_by_date(event_date: date, *, train_end: date = TRAIN_END, valid_end: date = VALID_END) -> str:
    """Their ``SplitByDate``: train ≤ t_end < valid ≤ v_end < test."""
    if event_date <= train_end:
        return "train"
    if event_date <= valid_end:
        return "valid"
    return "test"


def causal_user_feat(
    past: Sequence[dict],
    *,
    hist_len: int = 3,
) -> dict:
    """Running user portrait from reviews *strictly before* the current event."""
    cat_cnt: dict[str, int] = {}
    price_n, price_mean, price_m2 = 0, 0.0, 0.0
    hist: deque[dict] = deque(maxlen=int(hist_len))
    for r in past:
        cat = r.get("main_category") or "UNK"
        cat_cnt[cat] = cat_cnt.get(cat, 0) + 1
        px = _try_float(r.get("price"))
        if px is not None:
            price_n, price_mean, price_m2 = _welford_update(price_n, price_mean, price_m2, px)
        hist.append(
            {
                "title": str(r.get("review_title") or r.get("title") or ""),
                "text": str(r.get("review_text") or r.get("text") or ""),
            }
        )
    total = len(past)
    return {
        "cat_hist": {k: round(v / total, 4) for k, v in cat_cnt.items()} if total else {},
        "review_cnt": total,
        "price_mean": round(price_mean, 4) if price_n > 0 else None,
        "price_std": round(_welford_std(price_n, price_m2), 4) if price_n > 1 else 0.0,
        "history": list(hist),
    }


def leaky_user_feat(
    past_and_current: Sequence[dict],
    *,
    hist_len: int = 3,
) -> dict:
    """Include the current review in the portrait — recsys target leakage."""
    return causal_user_feat(past_and_current, hist_len=hist_len)


def current_review_leaks_into_feat(feat: dict, current: dict) -> bool:
    """True if this event's own review text sits in user_feat.history."""
    cur = (str(current.get("review_text") or current.get("text") or "")).strip()
    if not cur:
        return False
    for h in feat.get("history") or []:
        if cur and cur in str(h.get("text") or ""):
            return True
    return False


def emit_causal_rows(user_events: Sequence[dict], *, hist_len: int = 3) -> list[dict]:
    """Per-user time order: feat from past → emit → then update. Their Beam DoFn."""

    def _ts(r):
        ts = r.get("sort_timestamp")
        return ts if isinstance(ts, (int, float)) else -1

    rows = sorted(user_events, key=_ts)
    past: list[dict] = []
    out = []
    for r in rows:
        feat = causal_user_feat(past, hist_len=hist_len)
        rec = dict(r)
        rec["user_feat"] = feat
        rec["split"] = split_by_date(_as_date(r))
        rec["leaks_current"] = int(feat.get("review_cnt") or 0) != len(past)
        rec["text_in_history"] = current_review_leaks_into_feat(feat, r)
        out.append(rec)
        past.append(r)
    return out


def _as_date(r: dict) -> date:
    if r.get("event_date"):
        d = r["event_date"]
        if isinstance(d, date):
            return d
        return date.fromisoformat(str(d)[:10])
    ts = r.get("sort_timestamp")
    if isinstance(ts, (int, float)):
        return datetime.utcfromtimestamp(ts / 1000.0).date()
    return TRAIN_END


def build_user_text(feat: dict) -> str:
    """``data4model.build_user_text`` — structured portrait flattened to a prompt."""
    cat_hist = {k: v for k, v in (feat.get("cat_hist") or {}).items() if v}
    if cat_hist:
        cat_hist_str = "; ".join(f"{cat}: {cnt * 100:.0f}%" for cat, cnt in cat_hist.items())
    else:
        cat_hist_str = "No browsing history"
    n = feat.get("review_cnt") or 0
    pm = feat.get("price_mean")
    ps = feat.get("price_std") or 0.0
    parts = []
    for h in feat.get("history") or []:
        piece = (h.get("text") or h.get("title") or "").strip()
        if piece:
            parts.append(piece)
    hist = " ".join(f"Review{i + 1}: {p}" for i, p in enumerate(parts)) if parts else "No review history."
    return (
        f"Category history: {cat_hist_str}. "
        f"Total reviews: {n if n else 'No reviews'}. "
        f"Avg price: {pm if pm is not None else 'N/A'}. "
        f"Price std: {ps if ps else 'No price variation'}. "
        f"Review history: {hist}"
    )


def build_item_text(rec: dict) -> str:
    cat = rec.get("main_category") or "Unknown category"
    title = rec.get("product_title") or rec.get("title") or "No title"
    price = rec.get("price")
    price_str = f"{float(price):.2f}" if _try_float(price) is not None else "N/A"
    return (
        f"Item category: {cat}. Item title: {title}. Item price: {price_str}."
    )


def agod_concat_blob(rec: dict) -> str:
    """What the Amazon AGOD smoke actually hashes: item and user in one bag."""
    feat = rec.get("user_feat") or {}
    return f"{build_item_text(rec)}\n{build_user_text(feat)}"


def leak_report(rows: Sequence[dict]) -> dict:
    """How often current-review text lands in features; split counts."""
    n = len(rows)
    leaks = sum(1 for r in rows if r.get("leaks_current"))
    text_hits = sum(1 for r in rows if r.get("text_in_history"))
    splits = {"train": 0, "valid": 0, "test": 0}
    for r in rows:
        splits[r.get("split", "test")] = splits.get(r.get("split", "test"), 0) + 1
    return {
        "n": n,
        "n_leak_current": leaks,
        "leak_rate": float(leaks / n) if n else 0.0,
        "n_text_overlap": text_hits,
        "text_overlap_rate": float(text_hits / n) if n else 0.0,
        "splits": splits,
    }


def sentiment_label_shortcut(text: str) -> float:
    """Toy polarity: current-review text predicting its own star label."""
    t = (text or "").lower()
    pos = sum(w in t for w in ("great", "love", "excellent", "perfect", "good"))
    neg = sum(w in t for w in ("bad", "hate", "broken", "terrible", "poor"))
    return float(pos - neg)


def shortcut_gap(user_events: Sequence[dict]) -> dict:
    """Causal portrait vs leaking the current review into the blob.

    If the current review text is a feature, corr(text, label) is a
    sentiment task. Causal user_feat should be much weaker.
    """
    by_user: dict[str, list[dict]] = {}
    for r in user_events:
        by_user.setdefault(str(r.get("user_id", "")), []).append(r)
    emitted: list[dict] = []
    for rows in by_user.values():
        emitted.extend(emit_causal_rows(rows))
    y = np.array([int(r.get("label_good", 0)) for r in emitted], float)
    leak_s = np.array(
        [sentiment_label_shortcut(str(r.get("review_text") or "")) for r in emitted],
        float,
    )
    causal_s = np.array(
        [float((r["user_feat"].get("review_cnt") or 0)) for r in emitted],
        float,
    )

    def _corr(a, b):
        if a.std() < 1e-12 or b.std() < 1e-12:
            return 0.0
        return float(np.corrcoef(a, b)[0, 1])

    return {
        "corr_current_review_text": _corr(leak_s, y),
        "corr_causal_review_cnt": _corr(causal_s, y),
        "leak_rate": leak_report(emitted)["leak_rate"],
    }
