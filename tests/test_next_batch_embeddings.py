"""Identification tests for Wu InstDisc, Yu SDC, and typed GPM (no zip)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from gpm_fsds import (  # noqa: E402
    gpm_gate,
    gpm_identification,
    project_grad,
    representation_bases,
    residual_energy,
    run_gpm_method,
)
from instance_discrimination import (  # noqa: E402
    LinearProjector,
    MemoryBank,
    nonparametric_softmax_step,
    run_instdisc_stream,
)
from prototype_drift import (  # noqa: E402
    class_prototypes,
    ncm_accuracy,
    prototype_cosine,
    run_prototype_stream,
    sdc_compensate,
)
from typed_shift_stepsize import make_typed_stream  # noqa: E402


def test_instdisc_retrieves_own_slot():
    rng = np.random.default_rng(0)
    n, d, k = 24, 16, 8
    X = rng.normal(size=(n, d))
    X += np.eye(n, d) * 3.0
    proj = LinearProjector(d, k, seed=0)
    bank = MemoryBank(n, k, seed=1)
    acc = 0.0
    for _ in range(12):
        _, acc = nonparametric_softmax_step(proj, bank, X, np.arange(n), eta=0.2, tau=0.1)
    assert acc > 0.6


def test_instdisc_collision_rises_on_video_covariate():
    stream = make_typed_stream(
        n_batches=6, n_per=36, seed=3, cov={"video": 0.55, "audio": 0.0, "text": 0.0}
    )
    rec, _ = run_instdisc_stream(stream, dim=16, epochs=3, eta=0.12, seed=3)
    assert rec["mean_collision"]["video"] > rec["mean_collision"]["text"] + 0.2
    assert rec["mean_acc"]["video"] < rec["mean_acc"]["text"]


def test_sdc_recovers_location_shift_prototypes():
    rng = np.random.default_rng(4)
    y = np.repeat(np.arange(4), 12)
    mu = rng.normal(size=(4, 8))
    X0 = mu[y] + 0.05 * rng.normal(size=(len(y), 8))
    shift = np.linspace(0.4, 1.1, 8)
    X1 = X0 + shift
    mu0, _ = class_prototypes(X0, y, n_classes=4)
    mu1, _ = class_prototypes(X1, y, n_classes=4)
    stale = ncm_accuracy(X1, y, mu0)
    mu_sdc = sdc_compensate(mu0, X1 - shift, X1)
    sdc_acc = ncm_accuracy(X1, y, mu_sdc)
    assert np.linalg.norm(mu_sdc - mu1) < 0.5 * np.linalg.norm(mu0 - mu1)
    assert prototype_cosine(mu_sdc, mu1) > prototype_cosine(mu0, mu1) + 0.05
    assert sdc_acc >= stale


def test_prototype_stream_sdc_beats_stale_under_covariate():
    stream = make_typed_stream(
        n_batches=6, n_per=40, seed=5, cov={"video": 0.7, "audio": 0.0, "text": 0.0}
    )
    rec = run_prototype_stream(stream)
    assert rec["mean_proto_cos_sdc"]["video"] > rec["mean_proto_cos_stale"]["video"] + 0.1
    assert rec["mean_ncm_sdc"]["video"] >= rec["mean_ncm_stale"]["video"] - 1e-9


def test_gpm_residual_and_orthogonal_projection():
    rng = np.random.default_rng(6)
    M, _, k = representation_bases(rng.normal(size=(40, 10)), thresh=0.8, max_k=3)
    assert k >= 1
    dW = rng.normal(size=(10, 4))
    dWp = project_grad(dW, M)
    assert np.max(np.abs(M.T @ dWp)) < 1e-6
    Xin = rng.normal(size=(20, M.shape[1])) @ M.T
    assert residual_energy(Xin, M) < 1e-6


def test_gpm_gate_matches_tss_signs():
    assert gpm_gate(1.0, 0.0) is True
    assert gpm_gate(0.0, 0.8) is False
    assert gpm_gate(0.02, 0.01) is False
    assert gpm_gate(1.0, 0.8) is False


def test_typed_gpm_releases_on_concept_always_does_not():
    stream = make_typed_stream(
        n_batches=8, n_per=36, seed=11, concept={"video": 1.0}, concept_at=4
    )
    typed = run_gpm_method(stream, mode="typed", steps_per_batch=3, warmup_steps=2, seed=11)
    always = run_gpm_method(stream, mode="always", steps_per_batch=3, warmup_steps=2, seed=11)
    none = run_gpm_method(stream, mode="none", steps_per_batch=3, warmup_steps=2, seed=11)
    post = [h for h in typed["history"] if h["round"] >= 4]
    assert np.mean([h["gate"]["video"] for h in post]) < 0.5
    assert typed["post_acc"] >= always["post_acc"] - 1e-9
    assert none["post_acc"] >= always["post_acc"] - 1e-9


def test_gpm_identification_video_has_c_on_covariate():
    stream = make_typed_stream(
        n_batches=5, n_per=30, seed=2, cov={"video": 0.4, "audio": 0.0, "text": 0.0}
    )
    rec = gpm_identification(stream, thresh=0.9, max_k=6)
    assert rec["mean_c"]["video"] > rec["mean_c"]["text"]
    assert rec["mean_rank"]["video"] >= 1
