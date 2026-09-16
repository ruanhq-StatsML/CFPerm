"""Smoke tests for OnlineRFPerm live-infer use-case (mock backend, no GPU)."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "agod" / "online_rfperm_live_infer.py"
OUT = ROOT / "results" / "agod" / "online_rfperm_live_infer_test"


def test_mock_live_infer_fires_after_cut():
    out = OUT
    if out.exists():
        for p in out.glob("*"):
            p.unlink()
    else:
        out.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable,
        str(SCRIPT),
        "--backend",
        "mock",
        "--n-per",
        "16",
        "--n-batches",
        "6",
        "--cut-batch",
        "3",
        "--out",
        str(out),
    ]
    env = {**dict(**{k: v for k, v in __import__("os").environ.items()}), "PYTHONPATH": str(ROOT)}
    r = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stdout + "\n" + r.stderr
    summary = json.loads((out / "summary.json").read_text())
    assert summary["n"] == 96
    assert summary["cut_batch"] == 3
    assert summary["first_fire_batch"] is not None
    assert summary["first_fire_batch"] >= summary["cut_batch"]
    assert summary["detection_delay_batch"] is not None
    assert 0 <= summary["detection_delay_batch"] <= 2
    assert summary["mean_y_bad_hop"] > summary["mean_y_bad_quiet"]
    quiet_fires = [h for h in summary["hops"] if h["regime"] == "quiet" and h["fired"]]
    assert not quiet_fires, "should not fire in quiet regime for mock hop"
    assert (out / "REPORT.md").exists()
    fired_actions = [h["action"] for h in summary["hops"] if h["fired"]]
    assert fired_actions, "expected at least one fire"
    assert any(a != "quiet_pass" for a in fired_actions)
