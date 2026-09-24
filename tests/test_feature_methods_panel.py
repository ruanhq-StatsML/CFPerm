"""Tests for multi-method feature panel (cmean / MMD / PO / FSDS)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "tencent_gr"))
from export_review_agent_card import build_card  # noqa: E402
from feature_methods_panel import build_feature_methods_panel  # noqa: E402


def _toy_tables():
    diag = pd.DataFrame(
        {
            "feature": ["a", "b", "c", "d", "e"],
            "cmean_abs": [5.0, 4.0, 1.0, 0.5, 3.0],
            "mean_W1": [0, 0, 0, 0, 0],
            "mean_W2": [1, 1, 0, 0, 1],
            "mmd_loco": [0.9, 0.1, 0.8, 0.05, 0.7],
            "mmd_full": [1.0] * 5,
            "po_vimp": [0.01, 0.5, 0.4, 0.02, 0.3],
        }
    )
    rank = pd.DataFrame(
        {
            "feature": ["a", "b", "c", "d", "e"],
            "f_score": [2.0, 1.5, 0.2, 3.0, 0.1],
            "rank": [2, 3, 4, 1, 5],
            "selected": [1, 1, 0, 1, 0],
        }
    )
    return diag, rank


def test_panel_separates_fsds_only_and_shift_consensus():
    diag, rank = _toy_tables()
    panel = build_feature_methods_panel(diag, rank, tip_features=["a", "b", "d"], top_k=3)
    assert "a" in panel["tops"]["cmean"]
    assert "d" in panel["tops"]["fsds"]  # high F but low shift → fsds_only candidate
    assert "d" in panel["fsds_only"] or "d" in panel["consensus_top"]
    assert panel["note"]
    assert "cmean" in panel["methods"]
    assert panel["methods"]["fsds"]["role"] == "supervised_y"
    assert panel["methods"]["po_vimp"]["role"] == "shift_po"


def test_card_includes_feature_methods(tmp_path: Path):
    diag, rank = _toy_tables()
    diag.to_csv(tmp_path / "feature_shift_diagnostics.csv", index=False)
    rank.to_csv(tmp_path / "fsds_feature_ranking.csv", index=False)
    blob = {
        "localize_k": 5,
        "timeline": {"gap_days": 30},
        "fsds_W1_holdout": {"top_features": ["a", "b", "d"]},
    }
    summary = tmp_path / "summary.json"
    summary.write_text(json.dumps(blob))
    card = build_card(blob, source=str(summary))
    fm = card["feature_methods"]
    assert fm["tops"]["cmean"]
    assert fm["tops"]["fsds"]
    assert "特征多法" in card["paste_for_agent"]
    assert "方法Top" in card["paste_for_agent"]
    assert "graph_shift_feat_consensus_top3" in card["ticket_custom_fields"]
    md = __import__("export_review_agent_card", fromlist=["card_to_md"]).card_to_md(card)
    assert "特征多法 Top" in md
