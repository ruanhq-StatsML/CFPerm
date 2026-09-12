"""Unit tests for finite-step / entropy / MGDA stacking logics (not MoE)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from agod.stack_logics import (
    cauchy_step,
    entropy_regularized,
    frank_wolfe_linear,
    mgda_min_norm,
    mixloss_weights,
    pseudo_bma_weights,
    quadratic_finite_step,
    rayleigh_gls,
    taylor_remainder_scale,
)


EXP = ("a", "b", "c")
SCORES = {"a": 0.85, "b": 0.45, "c": 0.10}


def _R():
    r = np.eye(3)
    r[0, 1] = r[1, 0] = 0.80
    r[0, 2] = r[2, 0] = 0.10
    r[1, 2] = r[2, 1] = 0.10
    return r


def test_eta_zero_is_vertex():
    pi = quadratic_finite_step(SCORES, _R(), EXP, eta=0.0)
    assert pi["a"] > 0.99


def test_large_eta_spreads_mass():
    small = quadratic_finite_step(SCORES, _R(), EXP, eta=0.05)
    large = quadratic_finite_step(SCORES, _R(), EXP, eta=8.0)
    assert small["a"] > large["a"]
    assert large["c"] >= small["c"] - 1e-9


def test_entropy_tau_path_vertex_to_equal():
    sharp = entropy_regularized(SCORES, EXP, tau=0.02)
    flat = entropy_regularized(SCORES, EXP, tau=8.0)
    assert sharp["a"] > 0.90
    assert max(flat.values()) - min(flat.values()) < 0.12


def test_mgda_ignores_scores_uses_gram():
    # high corr a-b → min-norm should not put 2/3 on {a,b}
    pi = mgda_min_norm(_R(), EXP)
    assert pi["c"] > 0.28
    assert pi["a"] + pi["b"] < 0.80


def test_mgda_not_equal_to_finite_step():
    stack = quadratic_finite_step(SCORES, _R(), EXP, eta=1.0)
    mgda = mgda_min_norm(_R(), EXP)
    # stacking still prefers high-s expert a; MGDA does not
    assert stack["a"] > mgda["a"]


def test_pseudo_bma_n_collapses():
    """η=1 on per-sample loss is mild; η=n (ELPD) collapses — Yao et al. 2018."""
    L = {"a": 0.10, "b": 0.40, "c": 0.55}
    mild = pseudo_bma_weights(L, EXP, n_eff=1.0)
    hard = pseudo_bma_weights(L, EXP, n_eff=50.0)
    assert mild["a"] < 0.55
    assert hard["a"] > 0.99
    mx = mixloss_weights(L, EXP, eta=4.0)
    assert mx["a"] > mild["a"]


def test_typical_lr_is_linear_regime():
    t = taylor_remainder_scale(3e-3)
    assert t["linear_regime"]
    t2 = taylor_remainder_scale(1.0)
    assert not t2["linear_regime"]


def test_entropy_has_no_input_gate():
    """π = softmax(s/τ) depends only on scores, not on a feature x — not MoE."""
    pi1 = entropy_regularized(SCORES, EXP, tau=0.5)
    pi2 = entropy_regularized(SCORES, EXP, tau=0.5)
    assert pi1 == pi2


def test_frank_wolfe_is_vertex():
    pi = frank_wolfe_linear(SCORES, EXP)
    assert pi["a"] > 0.99
    assert pi["b"] < 1e-9


def test_train_lr_quadratic_stays_at_vertex():
    pi = quadratic_finite_step(SCORES, _R(), EXP, eta=3e-3)
    assert pi["a"] > 0.95


def test_joint_eta_pi_is_gls_not_vertex():
    r = _R()
    joint = rayleigh_gls(SCORES, r, EXP)
    vertex = quadratic_finite_step(SCORES, r, EXP, eta=0.0)
    seq = quadratic_finite_step(SCORES, r, EXP, eta=3e-3)
    assert vertex["a"] > 0.99
    assert abs(seq["a"] - vertex["a"]) < 0.05
    # joint SNR combo downweights the clone of a
    assert joint["pi"]["a"] < vertex["a"]
    assert joint["eta_star"] > 0.0


def test_complementary_unit_grads_need_quadratic():
    """s = (1/√2, 1/√2, 0), R = I: linear picks a vertex; joint splits a/b."""
    s = {"a": 0.5**0.5, "b": 0.5**0.5, "c": 0.0}
    r = np.eye(3)
    vertex = quadratic_finite_step(s, r, EXP, eta=0.0)
    assert vertex["a"] > 0.99 or vertex["b"] > 0.99
    joint = rayleigh_gls(s, r, EXP)
    assert abs(joint["pi"]["a"] - 0.5) < 0.08
    assert abs(joint["pi"]["b"] - 0.5) < 0.08
    assert joint["pi"]["c"] < 0.05


def test_cauchy_on_equal_pi():
    pi = {"a": 1 / 3, "b": 1 / 3, "c": 1 / 3}
    eta = cauchy_step(pi, SCORES, _R(), EXP)
    assert eta > 0.0


def test_entropy_ignores_gram():
    """Entropy has no R knob; large-η quadratic with ρ_ab=0.8 funds unique c."""
    ent = entropy_regularized(SCORES, EXP, tau=1.0)
    indep = quadratic_finite_step(SCORES, np.eye(3), EXP, eta=8.0)
    corr = quadratic_finite_step(SCORES, _R(), EXP, eta=8.0)
    assert abs(sum(ent.values()) - 1.0) < 1e-9
    assert corr["c"] > indep["c"]


if __name__ == "__main__":
    for fn in list(globals().values()):
        if callable(fn) and fn.__name__.startswith("test_"):
            fn()
            print(f"  {fn.__name__}: OK")
    print("test_agod_stack_logics: OK")
