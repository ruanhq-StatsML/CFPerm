"""Smoke test for review-agent card wave diff."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "tencent_gr"))
from diff_review_agent_cards import diff_cards  # noqa: E402


def test_diff_detects_flips():
    prev = {
        "card_type": "review_agent_context",
        "source_summary": "w1",
        "tips": ["u_span_sec", "i_share_last"],
        "direction": {"sign_Dy": "pos", "tip_signs": {"u_span_sec": "+", "i_share_last": "+"}},
        "review_hint": {"queue_bucket": "活跃跨度异常（短刷/长挂）"},
    }
    cur = {
        "card_type": "review_agent_context",
        "source_summary": "w2",
        "tips": ["u_span_sec", "i_n_users"],
        "direction": {"sign_Dy": "neg", "tip_signs": {"u_span_sec": "-", "i_n_users": "+"}},
        "review_hint": {"queue_bucket": "触达用户数变动"},
    }
    d = diff_cards(prev, cur)
    assert d["sign_Dy"]["flipped"] is True
    assert d["queue_bucket"]["flipped"] is True
    assert "i_n_users" in d["tips_added"]
    assert "i_share_last" in d["tips_removed"]
    assert any(x["feature"] == "u_span_sec" for x in d["tip_sign_flips"])
    assert "翻转" in d["paste_for_agent"]


def test_diff_quiet_when_same():
    card = {
        "card_type": "review_agent_context",
        "source_summary": "x",
        "tips": ["a"],
        "direction": {"sign_Dy": "flat", "tip_signs": {"a": "0"}},
        "review_hint": {"queue_bucket": "q"},
    }
    d = diff_cards(card, card)
    assert d["sign_Dy"]["flipped"] is False
    assert "无显著差分" in d["paste_for_agent"]
