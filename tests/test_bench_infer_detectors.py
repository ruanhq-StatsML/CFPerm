"""Smoke: manuscript first_k trail-batch stats + dataset export."""
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
EXPORT = ROOT / "data" / "hf_cache" / "infer_bench_export"


def test_first_k_consecutive_logic():
    sys.path.insert(0, str(ROOT / "scripts" / "agod"))
    from bench_infer_detectors import first_k_consecutive

    det = [False, True, True, False, True, True, True]
    assert first_k_consecutive(det, 1) == 1
    assert first_k_consecutive(det, 2) == 1  # 头一回连续两个 fire 从 trail batch 1
    assert first_k_consecutive(det, 3) == 4
    assert first_k_consecutive([False, False], 1) is None


def test_bench_manuscript_tables():
    assert STREAM.exists()
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
    assert "halueval" in summary["results"]
    row = next(m for m in summary["results"]["halueval"]["batch"] if m["method"] == "OnlineRFPerm")
    assert "first1" in row and "first2" in row and "first3" in row
    # trail-relative: first1 should be >= 0 if fired, and < n_trail
    if row["first1"] is not None:
        assert 0 <= row["first1"] < summary["results"]["halueval"]["n_trail_batches"]
    tex = (OUT / "bench_detectors.tex").read_text()
    assert "first$_k$" in tex or "first" in tex
    assert (OUT / "first_k_logic.py").exists()
    assert (EXPORT / "halueval_preview.csv").exists()
    assert (EXPORT / "squad_preview.csv").exists()
