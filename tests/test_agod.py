"""AGOD: MSG decomposition, clever-covariate routing, and online distillation."""

from __future__ import annotations

import numpy as np

from agod.chronoberg import (
    AUDIO_DRIFT_YEAR,
    CHRONOBERG_WINDOWS,
    ChronoBergStream,
    SyntheticChronoBergConfig,
)
from agod.distill import Student, Teacher, agod_loss, sgd_step
from agod.experiment import run_baseline_comparison
from agod.metrics import attribution_consistency, recall_at_k
from agod.msg import compute_msg_state, rf_domain_auc_vimp
from agod.online import AGODConfig, run_stream
from agod.routing import gate_weights, softmax_weights


def _tiny_stream(seed: int = 7) -> ChronoBergStream:
    return ChronoBergStream(
        SyntheticChronoBergConfig(
            n_per_window=72,
            dims={"audio": 12, "image": 12, "text": 12},
            seed=seed,
        )
    )


def test_chronoberg_calendar_and_shift_flags():
    stream = _tiny_stream()
    assert stream.years == CHRONOBERG_WINDOWS
    ref = stream.reference()
    assert ref.year == 1750
    assert ref.drifted_audio is False
    late = stream.make_window(len(stream.years) - 1)
    assert late.year == 1950
    assert late.year >= AUDIO_DRIFT_YEAR
    assert late.drifted_audio is True
    assert late.X["audio"].shape == (72, 12)
    assert late.Y.shape == (72,)


def test_rf_domain_auc_near_half_under_no_shift():
    rng = np.random.default_rng(0)
    x0 = rng.normal(size=(80, 8))
    x1 = rng.normal(size=(80, 8))
    auc, vimp, _, overlap, _ = rf_domain_auc_vimp(x0, x1, seed=0)
    assert 0.35 < auc < 0.65
    assert vimp < 1.5
    assert overlap > 0.7


def test_rf_domain_detects_covariate_shift():
    rng = np.random.default_rng(1)
    x0 = rng.normal(size=(90, 8))
    x1 = rng.normal(size=(90, 8)) + np.array([2.2] + [0] * 7)
    auc, vimp, _, overlap, top = rf_domain_auc_vimp(x0, x1, seed=1)
    assert auc > 0.8
    assert vimp > 1.0
    assert overlap < 0.45
    assert top[0] == 0


def test_msg_separates_concept_drift_from_covariate_shift():
    stream = _tiny_stream(seed=11)
    ref = stream.reference()
    late = stream.make_window(len(stream.years) - 1)
    state = compute_msg_state(ref, late, gamma=1.5, seed=11)
    # Text should look like covariate shift; audio like concept drift; image quiet.
    assert state.details["text"].auc > state.details["image"].auc + 0.08
    assert state.details["audio"].auc > state.details["image"].auc
    assert state.gaps["audio"] > state.gaps["image"]
    assert state.gaps["audio"] >= state.gaps["text"]


def test_softmax_and_gate():
    alpha = softmax_weights(np.array([0.1, 0.2, 2.0]), tau=0.25)
    assert np.isclose(alpha.sum(), 1.0)
    assert alpha[2] > 0.75
    gated, loc = gate_weights(np.array([0.2, 0.2, 0.6]), np.array([0.51, 0.90, 0.70]))
    assert np.isclose(gated.sum(), 1.0)
    assert gated[0] < gated[2]
    assert loc == (1,)


def test_sgd_step_improves_local_alignment():
    rng = np.random.default_rng(3)
    stream = _tiny_stream(seed=3)
    batch = stream.reference()
    teacher = Teacher.from_dims(stream.dims, out_dim=8, rng=rng)
    student = Student.from_teacher(teacher, rng, scale=0.4)
    t_mod = teacher.embed_batch(batch)
    alpha = {"audio": 1.0, "image": 0.0, "text": 0.0}
    before = agod_loss(student.embed_batch(batch), t_mod, alpha).local["audio"]
    for _ in range(25):
        sgd_step(student, batch, t_mod, alpha, lr=0.15)
    after = agod_loss(student.embed_batch(batch), t_mod, alpha).local["audio"]
    assert after < before * 0.85


def test_recall_at_k_perfect_on_identity():
    z = np.eye(6)
    assert recall_at_k(z, z, k=1) == 1.0


def test_attribution_consistency_monotone():
    est = {"audio": 0.9, "image": 0.1, "text": 0.4}
    tru = {"audio": 2.0, "image": 0.0, "text": 1.0}
    assert attribution_consistency(est, tru) > 0.9


def test_b2_prefers_text_b3_prefers_audio_on_drift_window():
    stream = _tiny_stream(seed=21)
    cfg = AGODConfig(
        seed=21,
        pretrain_steps=8,
        steps_per_window=4,
        teacher_dim=8,
        tau=0.30,
        gamma=1.8,
        momentum=0.05,
    )
    b2 = run_stream(stream, baseline="B2", config=cfg)
    stream = _tiny_stream(seed=21)
    b3 = run_stream(stream, baseline="B3", config=cfg)
    drift2 = [lg for lg in b2.logs if lg.drift_recall is not None and lg.year >= 1850]
    drift3 = [lg for lg in b3.logs if lg.year >= 1850]
    assert drift2 and drift3
    mean_text_b2 = np.mean([lg.alpha["text"] for lg in drift2])
    mean_audio_b2 = np.mean([lg.alpha["audio"] for lg in drift2])
    mean_audio_b3 = np.mean([lg.alpha["audio"] for lg in drift3])
    mean_text_b3 = np.mean([lg.alpha["text"] for lg in drift3])
    assert mean_text_b2 > mean_audio_b2
    assert mean_audio_b3 > mean_text_b3 or mean_audio_b3 > mean_audio_b2


def test_baseline_comparison_smoke():
    cfg = SyntheticChronoBergConfig(
        n_per_window=48,
        dims={"audio": 10, "image": 10, "text": 10},
        seed=5,
    )
    tcfg = AGODConfig(
        seed=5,
        pretrain_steps=6,
        steps_per_window=3,
        teacher_dim=8,
        tau=0.35,
        gamma=1.5,
        momentum=0.1,
    )
    result = run_baseline_comparison(config=cfg, trainer_config=tcfg)
    assert set(result.metrics) == {"B1", "B2", "B3"}
    assert result.ranking
    for name in result.metrics:
        assert 0.0 <= result.metrics[name]["drift_subgroup_recall"] <= 1.0
        assert result.metrics[name]["steps"] == 4.0


def test_newsgroups_and_amazon_streams():
    from agod.suites import make_suite

    news = make_suite("newsgroups", n_per_window=24, seed=1, dim=10)
    amz = make_suite("amazon", n_per_window=24, seed=1, dim=10)
    assert news.years[0] == 1990
    assert amz.make_window(3).drifted_audio is True
    ref = news.reference()
    late = news.make_window(len(news.years) - 1)
    assert ref.X["text"].shape[1] == 10
    assert late.drifted_audio is True


def test_budget_agod_uses_fewer_flops_than_static():
    cfg = SyntheticChronoBergConfig(
        n_per_window=40,
        dims={"audio": 8, "image": 8, "text": 8},
        seed=4,
    )
    tcfg = AGODConfig(seed=4, pretrain_steps=4, steps_per_window=4, teacher_dim=6)
    result = run_baseline_comparison(config=cfg, trainer_config=tcfg, baselines=("B1", "B5"))
    assert result.metrics["B5"]["rel_flops"] < 0.55
    assert result.metrics["B1"]["rel_flops"] == 1.0


def test_llm_boost_beats_always_student_on_ood():
    from agod.llm_boost import run_llm_boost

    policies, _ = run_llm_boost(seed=3, n_id=48, n_ood=48, dim=32)
    assert policies["agod-boost"].ood_accuracy >= policies["always-student"].ood_accuracy - 1e-9
    assert policies["agod-boost"].latency < policies["always-teacher"].latency
    assert policies["always-student"].teacher_frac == 0.0

