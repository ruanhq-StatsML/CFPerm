"""Tests for 审出加速 human-use gate."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "tencent_gr"))
from check_review_accel_human_gate import evaluate  # noqa: E402


def test_not_ready_with_only_smoke():
    rep = evaluate(
        [{"label": "useful", "reviewer": "smoke"}, {"label": "useful", "reviewer": "ci"}],
        min_human=1,
    )
    assert rep["ready"] is False
    assert rep["status"] == "NOT_READY"
    assert rep["feature_freeze"] is True


def test_ready_with_human():
    rep = evaluate(
        [
            {"label": "useful", "reviewer": "smoke"},
            {"label": "not_useful", "reviewer": "alice"},
        ],
        min_human=1,
    )
    assert rep["ready"] is True
    assert rep["status"] == "READY"
    assert "alice" in rep["human_reviewers"]
