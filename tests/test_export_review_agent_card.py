"""Smoke test for review-agent card export."""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "tencent_gr"))
from export_review_agent_card import build_card  # noqa: E402


def test_build_card_minimal():
    blob = {
        "localize_k": 10,
        "localized_items_head": [1, 2, 3],
        "n_localized_edges": {"W1_train": 5, "W2": 9},
        "fsds_W1_holdout": {"top_features": ["u_span_sec", "i_share_last"]},
        "direction": {
            "Dy": 0.01,
            "sign_Dy": "pos",
            "tip_signs": {"u_span_sec": "+", "i_share_last": "+"},
            "report": "[Direction] Dy=0.01 (pos)",
        },
    }
    card = build_card(blob, source="test")
    assert card["card_type"] == "review_agent_context"
    assert card["direction"]["sign_Dy"] == "pos"
    assert "u_span_sec" in card["tips"]
    assert "paste_for_agent" in card and "sign_Dy=pos" in card["paste_for_agent"]
    assert "非定罪" in card["disclaimer"]
    assert "ticket_custom_fields" in card
    assert card["ticket_custom_fields"]["graph_shift_sign_dy"] == "pos"
    assert "u_span_sec" in card["ticket_custom_fields"]["graph_shift_tip_top3"]
    assert card["gray_flags"]["allow"] is True


def test_kill_switch_disables_card():
    blob = {
        "localize_k": 1,
        "fsds_W1_holdout": {"top_features": ["u_span_sec"]},
        "direction": {"sign_Dy": "pos", "tip_signs": {"u_span_sec": "+"}},
    }
    card = build_card(
        blob, source="x", flags={"enabled": True, "sample_rate": 1.0, "kill_switch": True}
    )
    assert card["gray_flags"]["allow"] is False
    assert card["review_hint"]["suggested_action_level"] == "L0_observe"
    assert card["ticket_custom_fields"]["graph_shift_gray_allow"] is False
