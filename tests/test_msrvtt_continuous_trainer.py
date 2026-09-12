"""Tests for the continuous trainer prototype."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from msrvtt_continuous_trainer import run_continuous_trainer  # noqa: E402
from msrvtt_multimodal_attribution import GROUP_NAMES, make_synthetic_bundle  # noqa: E402


def test_continuous_trainer_logs_round_to_round_pi():
    b = make_synthetic_bundle(
        n_videos=6, n_windows=20, seed=8, video_shift=1.8, audio_shift=0.2, text_shift=0.02
    )
    summary, probe = run_continuous_trainer(b, n_batches=10, eta0=0.04, steps_per_batch=4, n_estimators=25, seed=8)
    assert summary["n_batches"] == 10
    assert len(summary["history"]) == 10
    assert summary["history"][0]["phase"] == "warmup"
    assert summary["history"][0]["active"] == list(GROUP_NAMES)
    rot = [h for h in summary["history"] if h["phase"] == "rotate"]
    assert len(rot) == 9
    assert rot[0]["active"][:3] == list(GROUP_NAMES)
    for row in rot:
        s = sum(row["pi"][g] for g in GROUP_NAMES)
        assert abs(s - 1.0) < 1e-6
        assert set(row["lr"]) == set(GROUP_NAMES)
        assert "delta_pi" in row
    assert summary["mean_pi"]["video"] > summary["mean_pi"]["text"]
    assert probe.W["video"].shape[1] == 6


def test_continuous_trainer_plot(tmp_path):
    from msrvtt_continuous_trainer import plot_continuous_trainer

    b = make_synthetic_bundle(n_videos=4, n_windows=20, seed=3, video_shift=1.2)
    summary, _ = run_continuous_trainer(b, n_batches=5, steps_per_batch=2, n_estimators=15, seed=3)
    out = plot_continuous_trainer(summary, tmp_path / "cont.png")
    assert out.exists()
    assert out.stat().st_size > 2000
