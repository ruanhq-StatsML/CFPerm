"""Smoke: mid-stream fracture is caught; smooth quiet stays quiet."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "agod" / "online_rfperm_fracture_perturb.py"
OUT = ROOT / "results" / "agod" / "online_rfperm_fracture_pytest"
CARD = ROOT / "data" / "hf_cache" / "infer_bench_export" / "halueval.jsonl"


def test_invent_fracture_caught_on_halueval():
    if not CARD.exists():
        return
    OUT.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(SCRIPT),
        "--datasets",
        "halueval",
        "--perturbations",
        "invent_fracture",
        "label_flip",
        "--n",
        "200",
        "--n-ref",
        "100",
        "--win",
        "20",
        "--out",
        str(OUT),
    ]
    env = {**dict(**{k: v for k, v in __import__("os").environ.items()}), "PYTHONPATH": str(ROOT)}
    r = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, text=True, timeout=180)
    assert r.returncode == 0, r.stdout + "\n" + r.stderr
    summary = json.loads((OUT / "summary.json").read_text())
    assert summary["n_caught"] == summary["n_cases"] == 2
    for c in summary["cases"]:
        assert c["caught"]
        assert c["smooth_quiet"]
        assert c["fracture"]["orf_freeze"]["detection_delay_obs"] == 0
