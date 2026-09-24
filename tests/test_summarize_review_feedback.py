"""Tests for review feedback rollup."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "tencent_gr"))
from summarize_review_feedback import summarize  # noqa: E402


def test_summarize_by_bucket():
    rows = [
        {"label": "useful", "queue_bucket": "末跳"},
        {"label": "useful", "queue_bucket": "末跳"},
        {"label": "not_useful", "queue_bucket": "末跳"},
        {"label": "not_useful", "queue_bucket": "其它"},
    ]
    s = summarize(rows)
    assert s["overall"]["useful"] == 2
    assert s["overall"]["useful_rate"] == 0.5
    assert s["by_queue_bucket"][0]["queue_bucket"] == "末跳"
    assert s["by_queue_bucket"][0]["useful_rate"] == round(2 / 3, 4)
