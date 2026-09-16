"""Smoke: OOB-baseline bench emits LaTeX on existing streams."""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "agod" / "bench_infer_detectors.py"
STREAM = ROOT / "results" / "agod" / "online_rfperm_two_datasets" / "halueval" / "stream.jsonl"
OUT = ROOT / "results" / "agod" / "bench_infer_detectors_test"


def test_bench_oob_baselines_latex():
    assert STREAM.exists(), "need prior two-dataset streams"
    OUT.mkdir(parents=True, exist_ok=True)
    r = subprocess.run(
        [sys.executable, str(SCRIPT), "--out", str(OUT)],
        cwd=ROOT,
        env={**os.environ, "PYTHONPATH": str(ROOT)},
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert r.returncode == 0, r.stdout + "\n" + r.stderr
    summary = json.loads((OUT / "summary.json").read_text())
    assert "halueval" in summary["results"] and "squad" in summary["results"]
    methods = {m["method"] for m in summary["results"]["halueval"]}
    assert "OnlineRFPerm" in methods and "BOCPD" in methods and "ADWIN" in methods
    tex = (OUT / "bench_detectors.tex").read_text()
    assert "BOCPD" in tex and "OnlineRFPerm" in tex
