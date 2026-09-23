"""Tests for ticket field dry-run / POST helper."""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "tencent_gr"))
from post_review_card_ticket_fields import build_payload, load_ticket_fields  # noqa: E402


def test_load_from_card(tmp_path: Path):
    card = {
        "card_type": "review_agent_context",
        "ticket_custom_fields": {
            "graph_shift_sign_dy": "pos",
            "graph_shift_sla_level": "urgent",
        },
    }
    p = tmp_path / "card.json"
    p.write_text(json.dumps(card))
    fields = load_ticket_fields(p)
    assert fields["graph_shift_sign_dy"] == "pos"
    payload = build_payload(fields, source=str(p))
    assert payload["event"] == "graph_shift_review_card"
    assert payload["custom_fields"]["graph_shift_sla_level"] == "urgent"
    assert payload["disclaimer"] == "clue_not_conviction"


def test_load_flat_fields(tmp_path: Path):
    p = tmp_path / "tf.json"
    p.write_text(json.dumps({"graph_shift_queue_bucket": "q"}))
    assert load_ticket_fields(p)["graph_shift_queue_bucket"] == "q"
