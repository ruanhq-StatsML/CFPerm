"""HotpotQA faithfulness prototype: supporting_facts + answer-precision labels."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "agod" / "online_rfperm_hotpot_proto.py"
OUT = ROOT / "results" / "agod" / "online_rfperm_hotpot_proto_test"


def test_hotpot_proto_recommended_fires():
    out = OUT
    out.mkdir(parents=True, exist_ok=True)
    for p in out.glob("*"):
        if p.is_file():
            p.unlink()

    cmd = [
        sys.executable,
        str(SCRIPT),
        "--backend",
        "mock",
        "--n-per",
        "20",
        "--n-batches",
        "10",
        "--cut-batch",
        "5",
        "--out",
        str(out),
    ]
    env = {**dict(**{k: v for k, v in __import__("os").environ.items()}), "PYTHONPATH": str(ROOT)}
    r = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stdout + "\n" + r.stderr

    summary = json.loads((out / "summary.json").read_text())
    assert summary["n_ref"] == 100
    contrast = summary["label_contrast"]

    # Legacy Jaccard-vs-full should remain near-saturated (the bug we document).
    assert contrast["y_jaccard_full"]["quiet_mean"] >= 0.9

    # Recommended label opens a quiet→hop gap.
    rec = summary["recommended"]
    assert rec["mean_y_bad_hop"] > rec["mean_y_bad_quiet"] + 0.3
    assert rec["first_fire_batch"] is not None
    assert rec["detection_delay_batch"] is not None
    assert rec["detection_delay_batch"] <= 1

    assert (out / "REPORT.md").exists()
    assert (out / "stream_table.parquet").exists()
