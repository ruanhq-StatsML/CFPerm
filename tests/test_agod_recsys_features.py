"""Causal recsys features vs leaky current-review / random-split encoding."""
from __future__ import annotations

import sys
from datetime import date, datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from agod.recsys_features import (
    TRAIN_END,
    agod_concat_blob,
    build_user_text,
    current_review_leaks_into_feat,
    emit_causal_rows,
    leak_report,
    leaky_user_feat,
    shortcut_gap,
    split_by_date,
)


def _ts(d: date) -> int:
    return int(datetime(d.year, d.month, d.day, tzinfo=timezone.utc).timestamp() * 1000)


def _user_stream():
    return [
        {
            "user_id": "u0",
            "parent_asin": "A",
            "main_category": "Tools",
            "price": 20.0,
            "product_title": "drill",
            "review_title": "great drill",
            "review_text": "love this excellent tool",
            "sort_timestamp": _ts(date(2023, 1, 10)),
            "event_date": date(2023, 1, 10),
            "label_good": 1,
        },
        {
            "user_id": "u0",
            "parent_asin": "B",
            "main_category": "Sports",
            "price": 40.0,
            "product_title": "tent",
            "review_title": "bad tent",
            "review_text": "terrible broken poles",
            "sort_timestamp": _ts(date(2023, 4, 1)),
            "event_date": date(2023, 4, 1),
            "label_good": 0,
        },
        {
            "user_id": "u0",
            "parent_asin": "C",
            "main_category": "Tools",
            "price": 15.0,
            "product_title": "bit",
            "review_title": "perfect",
            "review_text": "good perfect fit",
            "sort_timestamp": _ts(date(2023, 8, 1)),
            "event_date": date(2023, 8, 1),
            "label_good": 1,
        },
    ]


def test_first_event_has_empty_history():
    rows = emit_causal_rows(_user_stream())
    assert rows[0]["user_feat"]["review_cnt"] == 0
    assert rows[0]["user_feat"]["history"] == []
    assert not rows[0]["leaks_current"]
    assert not rows[0]["text_in_history"]


def test_second_event_sees_only_the_first():
    rows = emit_causal_rows(_user_stream())
    feat = rows[1]["user_feat"]
    assert feat["review_cnt"] == 1
    assert "Tools" in feat["cat_hist"]
    assert feat["history"][0]["text"] == "love this excellent tool"
    assert not current_review_leaks_into_feat(feat, rows[1])


def test_leaky_feat_contains_current_review():
    stream = _user_stream()
    leak = leaky_user_feat(stream[:2])
    assert current_review_leaks_into_feat(leak, stream[1])


def test_temporal_split_matches_their_cut():
    assert split_by_date(date(2023, 6, 30)) == "train"
    assert split_by_date(date(2023, 7, 1)) == "valid"
    assert split_by_date(date(2023, 10, 1)) == "test"
    rows = emit_causal_rows(_user_stream())
    assert rows[0]["split"] == "train"
    assert rows[1]["split"] == "train"
    assert rows[2]["split"] == "valid"
    assert rows[2]["event_date"] > TRAIN_END


def test_causal_stream_has_zero_current_leak():
    rep = leak_report(emit_causal_rows(_user_stream()))
    assert rep["leak_rate"] == 0.0
    assert rep["splits"]["train"] == 2
    assert rep["splits"]["valid"] == 1


def test_current_review_text_is_a_label_shortcut():
    gap = shortcut_gap(_user_stream())
    assert gap["corr_current_review_text"] > 0.9
    assert gap["corr_current_review_text"] > gap["corr_causal_review_cnt"]


def test_user_text_does_not_embed_current_review():
    rows = emit_causal_rows(_user_stream())
    blob = build_user_text(rows[1]["user_feat"])
    assert "terrible broken poles" not in blob
    assert "love this excellent tool" in blob
    concat = agod_concat_blob(rows[1])
    assert "tent" in concat
    assert "terrible broken poles" not in concat


if __name__ == "__main__":
    for fn in list(globals().values()):
        if callable(fn) and fn.__name__.startswith("test_"):
            fn()
            print(f"  {fn.__name__}: OK")
    print("test_agod_recsys_features: OK")
