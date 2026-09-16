"""Tests for style-drift iteration prototype."""
from __future__ import annotations

import json
from pathlib import Path

from scripts.agod.style_drift_iter_proto import (
    STYLE_NAMES,
    _synth_hh,
    action_step,
    main,
    run_iteration,
)

OUT = Path(__file__).resolve().parents[1] / "results" / "agod" / "style_drift_iter"
DOCS = Path(__file__).resolve().parents[1] / "docs" / "biz"


def test_style_iter_before_after_improves_on_synth():
    summary = run_iteration(_synth_hh(360), seed=0)
    b = summary["before_after"]
    assert b["auc_pre_post_before"] > 0.7
    assert b["fired_before"] is True
    assert b["top_family_before"] in ("length", "punct", "register")
    assert b["top_metric_before"] in STYLE_NAMES
    assert b["auc_pre_post_after"] < b["auc_pre_post_before"] - 0.02
    assert b["improved"] is True
    assert b["fired_after"] is False
    assert "另账" in summary["external_one_liner_cn"]
    ticket = summary["round0_before"]["action"]["ticket"]
    assert ticket == "style_only_creative_decoding"
    dont = " ".join(summary["round0_before"]["action"]["dont"])
    assert "偏好" in dont and "客服" in dont


def test_action_quiet_when_not_fired():
    detect = {"fired": False, "ledger": "另账"}
    attr = {"top_family": "register", "top_metric": "formal"}
    act = action_step(detect, attr)
    assert act["ticket"] == "none"


def test_main_synth_writes_artifacts():
    assert main(["--synth", "--seed", "0", "--n", "600"]) == 0
    summary = json.loads((OUT / "summary.json").read_text())
    assert (OUT / "REPORT.md").exists()
    assert (DOCS / "STYLE_DRIFT_ITER_PROTOTYPE.md").exists()
    assert summary["before_after"]["improved"]
    assert "iteration_flow" in summary
