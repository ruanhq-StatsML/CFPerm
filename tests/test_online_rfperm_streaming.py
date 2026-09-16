"""Smoke tests for continuous streaming OnlineRFPerm."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "agod" / "online_rfperm_streaming_test.py"
OUT = ROOT / "results" / "agod" / "online_rfperm_streaming_test_pytest"
STREAM = ROOT / "results" / "agod" / "online_rfperm_multi_datasets" / "halueval" / "stream.jsonl"


def test_streaming_catches_hop_or_skips_if_no_source():
    if not STREAM.exists():
        return  # multi streams not materialized in this checkout
    if OUT.exists():
        for p in OUT.glob("*"):
            if p.is_file():
                p.unlink()
    else:
        OUT.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(SCRIPT),
        "--datasets",
        "halueval",
        "--win",
        "20",
        "--gate",
        "1.25",
        "--out",
        str(OUT),
    ]
    env = {**dict(**{k: v for k, v in __import__("os").environ.items()}), "PYTHONPATH": str(ROOT)}
    r = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, text=True, timeout=180)
    assert r.returncode == 0, r.stdout + "\n" + r.stderr
    summary = json.loads((OUT / "summary.json").read_text())
    ds = summary["datasets"][0]
    assert ds["dataset"] == "halueval"
    # injected hop should be caught by at least one continuous mode
    caught = (
        ds["orf_freeze"]["first1"] is not None or ds["orf_slide"]["first1"] is not None
    )
    assert caught, ds
    # smooth control should stay quiet
    assert ds["smooth_control"]["freeze_first1"] is None
    assert ds["smooth_control"]["slide_first1"] is None
