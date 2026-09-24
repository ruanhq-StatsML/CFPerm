"""Smoke test for batch review-agent card export."""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "tencent_gr"))
from export_review_agent_cards_batch import export_one, find_summaries  # noqa: E402


def test_find_and_export(tmp_path: Path):
    root = tmp_path / "results"
    run = root / "wave_a"
    run.mkdir(parents=True)
    blob = {
        "localize_k": 3,
        "fsds_W1_holdout": {"top_features": ["u_span_sec"]},
        "direction": {"sign_Dy": "pos", "tip_signs": {"u_span_sec": "+"}},
    }
    summary = run / "summary.json"
    summary.write_text(json.dumps(blob))
    found = find_summaries(root)
    assert found == [summary]
    out = tmp_path / "cards"
    dest = export_one(
        summary,
        out_parent=out,
        root=root,
        flags={"enabled": True, "sample_rate": 1.0, "kill_switch": False},
    )
    assert dest.name == "wave_a"
    card = json.loads((dest / "review_agent_card.json").read_text())
    assert card["gray_flags"]["allow"] is True
    assert (dest / "ticket_custom_fields.json").exists()


def test_slug_when_summary_at_root(tmp_path: Path):
    root = tmp_path / "one_run"
    root.mkdir()
    summary = root / "summary.json"
    summary.write_text(
        json.dumps(
            {
                "localize_k": 1,
                "fsds_W1_holdout": {"top_features": ["u_span_sec"]},
                "direction": {"sign_Dy": "flat", "tip_signs": {}},
            }
        )
    )
    dest = export_one(
        summary,
        out_parent=tmp_path / "out",
        root=root,
        flags={"enabled": True, "sample_rate": 1.0, "kill_switch": False},
    )
    assert dest.name == "one_run"
    assert dest != tmp_path / "out"
