"""Typed shift stepsize: identification and signed-LR tests (no zip required)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from msrvtt_continuous_trainer import SeparateHeadProbe  # noqa: E402
from msrvtt_multimodal_attribution import GROUP_NAMES  # noqa: E402
from typed_shift_stepsize import (  # noqa: E402
    concept_intensity,
    covariate_intensity,
    make_typed_stream,
    run_method,
    tss_lr,
)


def _split(X):
    from msrvtt_multimodal_attribution import GROUPS

    return {n: X[:, sl] for n, sl in GROUPS.items()}


def test_covariate_intensity_ranks_video_block():
    stream = make_typed_stream(
        n_batches=6, n_per=32, seed=3, cov={"video": 0.7, "audio": 0.05, "text": 0.0}
    )
    i0 = np.flatnonzero(stream.batch == 0)
    i1 = np.flatnonzero(stream.batch == 5)
    c, share, _ = covariate_intensity(stream.X[i0], stream.X[i1])
    assert c["video"] > c["audio"]
    assert c["video"] > c["text"]
    assert share["video"] > 0.5


def test_concept_intensity_jumps_after_rotation():
    stream = make_typed_stream(
        n_batches=8,
        n_per=36,
        seed=5,
        cov={"video": 0.0, "audio": 0.0, "text": 0.0},
        concept={"video": 1.0},
        concept_at=4,
    )
    probe = SeparateHeadProbe(n_classes=4, seed=5)
    i0 = np.flatnonzero(stream.batch == 0)
    Xs0 = _split(stream.X[i0])
    y0 = stream.y[i0]
    for _ in range(12):
        probe.step(Xs0, y0, {g: 0.08 for g in GROUP_NAMES})
    pre = np.flatnonzero(stream.batch == 3)
    post = np.flatnonzero(stream.batch == 5)
    d_pre, _ = concept_intensity(probe, Xs0, y0, _split(stream.X[pre]), stream.y[pre])
    d_post, _ = concept_intensity(probe, Xs0, y0, _split(stream.X[post]), stream.y[post])
    assert d_post["video"] > d_pre["video"] + 0.05
    assert d_post["video"] > d_post["text"]


def test_tss_lr_signs_and_quiet_freeze():
    quiet = tss_lr({"video": 0.02, "audio": 0.01, "text": 0.0}, {"video": 0.0, "audio": 0.0, "text": 0.01})
    assert quiet["video"] == 0.0
    assert quiet["text"] == 0.0
    cov = tss_lr({"video": 1.2, "audio": 0.1, "text": 0.0}, {"video": 0.0, "audio": 0.0, "text": 0.0}, eta0=0.1)
    concept = tss_lr({"video": 0.0, "audio": 0.0, "text": 0.0}, {"video": 0.8, "audio": 0.0, "text": 0.0}, eta0=0.1)
    assert concept["video"] > cov["video"]
    assert concept["video"] > concept["text"]
    # inverting π would give text the largest step; TSS must not
    pi = tss_lr({"video": 1.2, "audio": 0.3, "text": 0.0}, {"video": 0.05, "audio": 0.02, "text": 0.0}, eta0=0.1)
    assert pi["text"] == 0.0


def test_run_method_tss_on_cov_only_shrinks_video_lr():
    stream = make_typed_stream(n_batches=6, n_per=28, seed=9, cov={"video": 0.6})
    tss, _ = run_method(stream, method="tss", eta0=0.08, steps_per_batch=3, seed=9)
    const, _ = run_method(stream, method="constant", eta0=0.08, steps_per_batch=3, seed=9)
    pi, _ = run_method(stream, method="fsds_pi", eta0=0.08, steps_per_batch=3, seed=9)
    assert tss["mean_c"]["video"] > tss["mean_c"]["text"]
    assert tss["mean_lr"]["video"] < const["mean_lr"]["video"]
    assert pi["mean_lr"]["video"] > pi["mean_lr"]["text"]


def test_run_method_tss_raises_lr_at_concept():
    stream = make_typed_stream(
        n_batches=8, n_per=28, seed=11, concept={"video": 1.0}, concept_at=4
    )
    tss, _ = run_method(stream, method="tss", eta0=0.08, steps_per_batch=3, seed=11)
    cosine, _ = run_method(stream, method="cosine", eta0=0.08, steps_per_batch=3, seed=11)
    pre = [h["lr"]["video"] for h in tss["history"] if h["phase"] == "adapt" and h["round"] < 4]
    post = [h["lr"]["video"] for h in tss["history"] if h["phase"] == "adapt" and h["round"] >= 4]
    assert np.mean(post) > np.mean(pre)
    last_t = tss["history"][-1]["lr"]["video"]
    last_c = cosine["history"][-1]["lr"]["video"]
    assert last_t > last_c
    assert tss["n_classes"] == 4
