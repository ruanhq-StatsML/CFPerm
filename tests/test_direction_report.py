"""Tests for direction ensure + S1/S2/S3 scenario mapping."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "tencent_gr"))
from direction_report import (  # noqa: E402
    build_direction_dict,
    ensure_direction,
    scenario_from_direction,
)
from export_review_agent_card import build_card  # noqa: E402


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


def test_ensure_direction_fills_tips_without_inventing_dy(tmp_path: Path):
    feat = pd.DataFrame(
        {
            "feature": ["i_share_linear", "i_credit_last"],
            "mean_W1": [0.0, 0.1],
            "mean_W2": [1.0, 0.9],
        }
    )
    diag = tmp_path / "feature_shift_diagnostics.csv"
    feat.to_csv(diag, index=False)
    blob = {
        "fsds_W1_holdout": {"top_features": ["i_share_linear", "i_credit_last"]},
        "W1_meta": {"pos_rate": 0.001},
        "W2_meta": {"pos_rate": 0.002},  # must NOT invent Dy from this
    }
    d = ensure_direction(blob, feat_diag_path=diag)
    assert d["sign_Dy"] == "flat"
    assert d["dy_missing"] is True
    assert d["Dy"] is None
    assert d["tip_signs"]["i_share_linear"] == "+"
    assert d["tip_signs"]["i_credit_last"] == "+"


def test_scenario_s3_when_flat_linear():
    d = {
        "sign_Dy": "flat",
        "Dy": None,
        "dy_missing": True,
        "tip_signs": {"i_share_linear": "+", "i_credit_linear": "+"},
    }
    sc = scenario_from_direction(d)
    assert sc["family_code"] == "S3_drift"
    assert sc["auto_ban"] is False


def test_scenario_s1_last_hop():
    d = {"sign_Dy": "pos", "Dy": 0.01, "tip_signs": {"i_share_last": "+"}}
    sc = scenario_from_direction(d)
    assert sc["family_code"] == "S1_brush"
    assert sc["sub_code"] == "last_hop"


def test_build_card_auto_ensures_tip_signs(tmp_path: Path):
    feat = pd.DataFrame(
        {
            "feature": ["i_share_linear", "i_log1p_n_users"],
            "mean_W1": [0.0, 0.0],
            "mean_W2": [1.0, -0.5],
        }
    )
    summary = tmp_path / "summary.json"
    feat.to_csv(tmp_path / "feature_shift_diagnostics.csv", index=False)
    blob = {
        "localize_k": 10,
        "n_localized_edges": {"W1_train": 5, "W2": 8},
        "timeline": {"gap_days": 30},
        "fsds_W1_holdout": {
            "top_features": ["i_share_linear", "i_log1p_n_users"]
        },
    }
    summary.write_text(json.dumps(blob))
    card = build_card(blob, source=str(summary))
    assert card["direction"]["sign_Dy"] == "flat"
    assert card["direction"]["tip_signs"]["i_share_linear"] == "+"
    assert card["direction"]["tip_signs"]["i_log1p_n_users"] == "-"
    assert card["scenario"]["family_code"] == "S3_drift"
    assert "i_share_linear (+)" in card["paste_for_agent"]
    assert card["ticket_custom_fields"]["graph_shift_scenario_family"] == "S3_drift"
