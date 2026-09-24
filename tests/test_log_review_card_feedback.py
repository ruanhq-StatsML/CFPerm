"""Tests for review-card feedback logger."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "tencent_gr"))
from log_review_card_feedback import build_feedback  # noqa: E402


def test_build_feedback_useful():
    card = {
        "source_summary": "s.json",
        "direction": {"sign_Dy": "pos"},
        "review_hint": {"queue_bucket": "末跳嫌疑"},
        "tips": ["i_share_last", "u_span_sec"],
        "ticket_custom_fields": {
            "graph_shift_sign_dy": "pos",
            "graph_shift_queue_bucket": "末跳嫌疑",
            "graph_shift_tip_top3": "i_share_last,u_span_sec",
        },
    }
    row = build_feedback(card, label="useful", note="hit", reviewer="t")
    assert row["label"] == "useful"
    assert row["sign_Dy"] == "pos"
    assert row["tips_top3"][0] == "i_share_last"
