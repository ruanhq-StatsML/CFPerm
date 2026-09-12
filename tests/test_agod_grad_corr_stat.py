"""Unit tests for modality-grad correlation → next LR (no GPU / shards)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from agod.grad_corr_stat import (
    characterize_corr,
    fisher_z_test,
    gls_weights,
    shrink_correlation,
    correlation_from_pairs,
    effective_ensemble_size,
    partial_uniqueness,
    stat_corr_lr,
    temporal_gain,
    variance_stabilizing_scale,
)
from agod.lr_controller import alpha_to_lr, schedule_modality_lr

MODS = ("video", "text", "audio")


def _iso(c: float, mods=MODS) -> dict[str, float]:
    return {
        f"{a}|{b}": float(c)
        for i, a in enumerate(mods)
        for b in mods[i + 1 :]
    }


def test_independent_reduces_to_soft():
    mods = list(MODS)
    alpha = {m: 1.0 / 3.0 for m in mods}
    pair = _iso(0.0)
    packed = stat_corr_lr(alpha, pair, mods, use_temporal=False, use_conflict=False)
    soft = alpha_to_lr(alpha, mods, beta=0.10)
    st = packed["stats"]
    assert st["n_eff"]["n_eff"] > 2.4
    assert packed["eta_global"] > 0.85
    for m in mods:
        assert abs(packed["pi"][m] - alpha[m]) < 0.08
        assert abs(packed["lr"][m] - soft[m]) < 0.25


def test_equicorr_shrinks_neff_and_eta():
    mods = list(MODS)
    alpha = {m: 1.0 / 3.0 for m in mods}
    pair = _iso(0.90)
    packed = stat_corr_lr(alpha, pair, mods, use_temporal=False, use_conflict=False)
    neff = packed["stats"]["n_eff"]["n_eff"]
    # M / (1+(M-1)ρ) = 3/2.8 ≈ 1.07
    assert 0.9 < neff < 1.8
    assert packed["eta_global"] < 0.85
    # exchangeable voters: allocation stays flat
    pis = [packed["pi"][m] for m in mods]
    assert max(pis) - min(pis) < 0.08


def test_unequal_alpha_high_corr_concentrates_on_leader():
    mods = list(MODS)
    alpha = {"video": 0.70, "text": 0.20, "audio": 0.10}
    pair = _iso(0.88)
    packed = stat_corr_lr(alpha, pair, mods, use_temporal=False, use_conflict=False)
    assert packed["pi"]["video"] > 0.55
    assert packed["pi"]["video"] > packed["pi"]["text"] > packed["pi"]["audio"] - 1e-9
    assert packed["lr"]["video"] > packed["lr"]["audio"]


def test_gls_snr_beats_equal_under_hetero_alpha():
    mods = list(MODS)
    r = correlation_from_pairs(_iso(0.7), mods)
    r = shrink_correlation(r, lam=0.1)["R_shrink"]
    alpha = {"video": 0.65, "text": 0.25, "audio": 0.10}
    gls = gls_weights(r, alpha, mods)
    pi = np.array([gls["pi"][m] for m in mods])
    a = np.array([alpha[m] for m in mods])
    eq = np.full(3, 1.0 / 3.0)
    snr_gls = float((pi @ a) ** 2 / (pi @ r @ pi))
    snr_eq = float((eq @ a) ** 2 / (eq @ r @ eq))
    assert snr_gls >= snr_eq - 1e-9


def test_partial_uniqueness_drops_with_correlation():
    mods = list(MODS)
    u0 = partial_uniqueness(correlation_from_pairs(_iso(0.0), mods), mods)["uniqueness"]
    u1 = partial_uniqueness(
        shrink_correlation(correlation_from_pairs(_iso(0.85), mods), lam=0.05)["R_shrink"],
        mods,
    )["uniqueness"]
    assert np.mean(list(u0.values())) > np.mean(list(u1.values()))


def test_conflict_damps_opposing_shared():
    mods = list(MODS)
    alpha = {m: 1.0 / 3.0 for m in mods}
    pair = _iso(0.2)
    packed_ok = stat_corr_lr(
        alpha, pair, mods, use_temporal=False, cos_to_shared={m: 0.8 for m in mods}
    )
    packed_bad = stat_corr_lr(
        alpha,
        pair,
        mods,
        use_temporal=False,
        cos_to_shared={"video": -0.9, "text": 0.8, "audio": 0.8},
    )
    assert packed_bad["lr"]["video"] < packed_ok["lr"]["video"]
    assert packed_bad["conflict_gain"]["video"] < 0.7


def test_temporal_gain_monotone():
    assert temporal_gain(0.95) > temporal_gain(0.0) > temporal_gain(-0.8)


def test_fisher_z_detects_high_rho():
    hi = fisher_z_test(0.85, 40.0)
    lo = fisher_z_test(0.05, 8.0)
    assert hi["significant_pos"]
    assert hi["ci95"][0] > 0.6
    assert not lo["significant_pos"]


def test_variance_scale_collinear_vs_independent():
    n = 3
    pi = np.full(n, 1.0 / n)
    eta_i = variance_stabilizing_scale(pi, np.eye(n))
    r_hi = 0.05 * np.eye(n) + 0.95 * np.ones((n, n))
    np.fill_diagonal(r_hi, 1.0)
    eta_c = variance_stabilizing_scale(pi, r_hi)
    assert abs(eta_i - 1.0) < 0.05
    assert eta_c < eta_i


def test_schedule_dispatch_gls_and_stat():
    mods = list(MODS)
    alpha = {"video": 0.5, "text": 0.3, "audio": 0.2}
    pair = _iso(0.75)
    lr_gls = schedule_modality_lr("soft_gls", alpha, mods, pair_cos=pair)
    lr_stat = schedule_modality_lr("soft_stat", alpha, mods, pair_cos=pair)
    lr_soft = schedule_modality_lr("soft", alpha, mods)
    for m in mods:
        assert np.isfinite(lr_gls[m]) and lr_gls[m] > 0
        assert np.isfinite(lr_stat[m]) and lr_stat[m] > 0
    # empty geometry → soft fallback
    lr_fb = schedule_modality_lr("soft_stat", alpha, mods, pair_cos={})
    for m in mods:
        assert abs(lr_fb[m] - lr_soft[m]) < 1e-9


def test_characterize_corr_keys():
    info = characterize_corr(_iso(0.4), list(MODS), n_obs=20)
    assert "n_eff" in info and "gls" in info and "fisher_z" in info
    assert 1.0 <= info["n_eff"]["n_eff"] <= 3.0 + 1e-6


def test_blend_stays_near_alpha_when_independent():
    mods = list(MODS)
    alpha = {"video": 0.70, "text": 0.20, "audio": 0.10}
    packed = stat_corr_lr(alpha, _iso(0.0), mods, use_temporal=False, use_conflict=False)
    assert packed["gamma_gls"] < 0.15
    assert abs(packed["pi"]["video"] - 0.70) < 0.08


def test_blend_moves_toward_gls_when_collinear():
    mods = list(MODS)
    alpha = {"video": 0.70, "text": 0.20, "audio": 0.10}
    packed = stat_corr_lr(alpha, _iso(0.90), mods, use_temporal=False, use_conflict=False)
    assert packed["gamma_gls"] > 0.55
    assert packed["pi"]["video"] >= 0.70 - 1e-6
    assert packed["pi"]["audio"] <= 0.10 + 1e-6
    assert packed["pi"]["video"] > packed["pi_gls"]["video"] - 1e-6 or packed["pi"]["video"] > 0.80


def test_kish_matches_equicorr_formula():
    mods = list(MODS)
    rho = 0.4
    r = correlation_from_pairs(_iso(rho), mods)
    out = effective_ensemble_size(r)
    closed = 3.0 / (1.0 + 2.0 * rho)
    assert abs(out["n_eff_kish"] - closed) < 1e-6
    assert abs(out["equicorr_neff"] - closed) < 1e-6


if __name__ == "__main__":
    for fn in list(globals().values()):
        if callable(fn) and fn.__name__.startswith("test_"):
            fn()
            print(f"  {fn.__name__}: OK")
    print("test_agod_grad_corr_stat: OK")
