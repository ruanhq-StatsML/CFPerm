"""Smoke test for 审出加速包 e2e runner."""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "tencent_gr"))
from smoke_review_accel_pack import run_smoke  # noqa: E402


def test_run_smoke(tmp_path: Path):
    summary = tmp_path / "summary.json"
    summary.write_text(
        json.dumps(
            {
                "localize_k": 3,
                "gap_days": 30,
                "fsds_W1_holdout": {"top_features": ["u_span_sec"]},
                "direction": {"sign_Dy": "pos", "tip_signs": {"u_span_sec": "+"}},
            }
        )
    )
    flags = tmp_path / "flags.json"
    flags.write_text(
        json.dumps({"enabled": True, "sample_rate": 1.0, "kill_switch": False})
    )
    out = tmp_path / "smoke"
    report = run_smoke(
        summary,
        out_dir=out,
        flags_path=flags,
        tip_overlay_path=None,
        label="useful",
    )
    assert report["steps_ok"][-1] == "dry_run_ticket_post"
    assert (out / "smoke_report.md").exists()
    assert (out / "review_agent_card.json").exists()
    assert (out / "ticket_post_receipt.json").exists()
    assert report["eta_soft_hint"] == "lower_confidence"
