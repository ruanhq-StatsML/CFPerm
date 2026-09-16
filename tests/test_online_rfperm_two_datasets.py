"""Two-dataset OnlineRFPerm live-infer smoke (mock backend)."""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "agod" / "online_rfperm_two_datasets.py"
OUT = ROOT / "results" / "agod" / "online_rfperm_two_datasets_test"


def test_two_datasets_mock_emits_latex():
    OUT.mkdir(parents=True, exist_ok=True)
    # ensure squad cache exists via mock path: script will download if missing
    cmd = [
        sys.executable,
        str(SCRIPT),
        "--backend",
        "mock",
        "--n-per",
        "10",
        "--n-batches",
        "6",
        "--cut-batch",
        "3",
        "--out",
        str(OUT),
    ]
    env = {**os.environ, "PYTHONPATH": str(ROOT)}
    r = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, text=True, timeout=180)
    assert r.returncode == 0, r.stdout + "\n" + r.stderr
    summary = json.loads((OUT / "summary.json").read_text())
    assert len(summary["datasets"]) == 2
    names = {d["dataset"] for d in summary["datasets"]}
    assert names == {"halueval", "squad"}
    for d in summary["datasets"]:
        assert d["first_fire_batch"] is not None
        assert d["first_fire_batch"] >= d["cut_batch"]
        assert d["detection_delay_batch"] is not None
        assert d["detection_delay_batch"] >= 0
    tex = (OUT / "two_datasets.tex").read_text()
    assert "HaluEval" in tex or "halueval" in tex
    assert "SQuAD" in tex or "squad" in tex
    assert "OnlineRFPerm" in tex
    assert "\\begin{tabular}" in tex
