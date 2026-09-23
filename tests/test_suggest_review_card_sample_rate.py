"""Tests for useful_rate → gray sample_rate suggestion."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "tencent_gr"))
from suggest_review_card_sample_rate import suggest_flags  # noqa: E402


def test_hold_when_insufficient_n():
    rep = suggest_flags(
        {"overall": {"useful": 1, "not_useful": 0, "useful_rate": 1.0}},
        {"enabled": True, "sample_rate": 1.0, "kill_switch": False},
        min_n=5,
    )
    assert rep["action"] == "hold"
    assert rep["suggested_flags"]["sample_rate"] == 1.0


def test_decrease_on_low_useful_rate():
    rep = suggest_flags(
        {"overall": {"useful": 2, "not_useful": 8, "useful_rate": 0.2}},
        {"enabled": True, "sample_rate": 1.0, "kill_switch": False},
        min_n=5,
        kill_rate=0.15,
        low_rate=0.4,
    )
    # 0.2 is above kill_rate 0.15 but below low_rate → decrease
    assert rep["action"] == "decrease_sample_rate"
    assert rep["suggested_flags"]["sample_rate"] == 0.5


def test_kill_switch_on_very_low_rate():
    rep = suggest_flags(
        {"overall": {"useful": 1, "not_useful": 19, "useful_rate": 0.05}},
        {"enabled": True, "sample_rate": 0.8, "kill_switch": False},
        min_n=5,
    )
    assert rep["action"] == "kill_switch_on"
    assert rep["suggested_flags"]["kill_switch"] is True


def test_increase_on_high_rate():
    rep = suggest_flags(
        {"overall": {"useful": 8, "not_useful": 2, "useful_rate": 0.8}},
        {"enabled": True, "sample_rate": 0.4, "kill_switch": True},
        min_n=5,
    )
    assert rep["action"] == "increase_sample_rate"
    assert rep["suggested_flags"]["kill_switch"] is False
    assert rep["suggested_flags"]["sample_rate"] > 0.4
