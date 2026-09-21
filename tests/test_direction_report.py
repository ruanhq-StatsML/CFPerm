"""Tests for signed direction JSON payload."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "tencent_gr"))
from direction_report import build_direction_dict  # noqa: E402


def test_build_direction_dict_pos_and_tip_signs():
    g0 = pd.DataFrame(
        {
            "y_convert": [0, 0, 1, 0],
            "a": [0.0, 0.0, 0.0, 0.0],
            "b": [1.0, 1.0, 1.0, 1.0],
        }
    )
    g1 = pd.DataFrame(
        {
            "y_convert": [1, 1, 1, 0],
            "a": [2.0, 2.0, 2.0, 2.0],
            "b": [0.0, 0.0, 0.0, 0.0],
        }
    )
    feat = pd.DataFrame(
        {
            "feature": ["a", "b"],
            "mean_W1": [0.0, 1.0],
            "mean_W2": [2.0, 0.0],
            "cmean_abs": [2.0, 1.0],
        }
    )
    d = build_direction_dict(g0, g1, ["a", "b"], feat_diag=feat)
    assert d["sign_Dy"] == "pos"
    assert d["Dy"] is not None and d["Dy"] > 0
    assert d["tip_signs"]["a"] == "+"
    assert d["tip_signs"]["b"] == "-"
    assert "report" in d and "tip_signs" in d["report"]
