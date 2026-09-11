"""Gap-guided LLM routing: synthetic DGP, no vendor LLM."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from clever_covariate_gap import make_synthetic_shift  # noqa: E402
from gap_guided_inference import (  # noqa: E402
    RouterConfig,
    blend_shares,
    budgeted_inference_eval,
    catalog_for_spec,
    critic_check,
    fit_frozen_gap,
    make_diffuse_equal_shift,
    openai_tool_schema,
    pack_blocks,
    pack_user_context,
    population_pack,
    predict_gap_scores,
    query_update,
    render_system_prompt,
    route_from_shares,
    simplex_entropy,
    worked_example_packet,
)


def test_blend_and_entropy_shapes():
    pi = np.array([0.1, 0.6, 0.2, 0.1])
    s = np.array([[0.05, 0.8, 0.1, 0.05], [0.25, 0.25, 0.25, 0.25]])
    r = blend_shares(pi, s, lam=0.4)
    assert r.shape == (2, 4)
    assert np.allclose(r.sum(axis=1), 1.0, atol=1e-6)
    r_pop = blend_shares(pi, s, lam=1.0)
    assert np.allclose(r_pop[0], pi)
    h_flat = simplex_entropy(np.full(4, 0.25))
    h_peak = simplex_entropy(np.array([0.85, 0.05, 0.05, 0.05]))
    assert h_flat > h_peak


def test_concentrated_route_enables_gt_tool_only():
    names = ["text", "valence", "arousal", "dominance"]
    spec_names = names
    from clever_covariate_gap import ModalitySpec

    spec = ModalitySpec(names=list(spec_names), slices=[slice(0, 1)] * 4)
    cat = catalog_for_spec(spec)
    pi = np.array([0.12, 0.70, 0.10, 0.08])
    s = np.array([0.08, 0.78, 0.09, 0.05])
    d = route_from_shares(names, pi, s, cat, RouterConfig(tau_hi=0.42, k_max=2))
    assert d.critic_must_cite == "valence"
    assert d.tools_enabled == ["vad_valence"]
    assert d.packed_blocks == ["valence"]
    assert d.pack_mode == "same_expert"
    assert "read_tokens" in d.tools_disabled
    assert d.read_order[0] == "valence"
    assert d.abstain is False
    prompt = render_system_prompt(d)
    assert "vad_valence" in prompt
    assert "read_tokens" in prompt
    assert "W" not in json.dumps(d.packet)
    schema = openai_tool_schema(d, cat)
    assert [t["function"]["name"] for t in schema] == ["vad_valence"]


def test_flat_blend_abstains():
    from clever_covariate_gap import ModalitySpec

    names = ["text", "valence", "arousal", "dominance"]
    spec = ModalitySpec(names=list(names), slices=[slice(0, 1)] * 4)
    cat = catalog_for_spec(spec)
    pi = np.full(4, 0.25)
    s = np.full(4, 0.25)
    d = route_from_shares(names, pi, s, cat, RouterConfig(lam=0.5, tau_lo=0.18))
    assert d.abstain is True
    assert d.tools_enabled == []
    assert d.packed_blocks == []
    assert d.pack_mode == "hitl_abstain"
    assert d.ask_missing
    assert d.critic_must_cite is None


def test_critic_rejects_wrong_channel():
    from clever_covariate_gap import ModalitySpec

    names = ["text", "valence", "arousal", "dominance"]
    spec = ModalitySpec(names=list(names), slices=[slice(0, 1)] * 4)
    cat = catalog_for_spec(spec)
    d = route_from_shares(
        names, np.array([0.1, 0.7, 0.1, 0.1]), np.array([0.1, 0.7, 0.1, 0.1]), cat,
    )
    bad = critic_check(d, ["text"])
    good = critic_check(d, ["valence"])
    assert bad.mismatch and not bad.ok
    assert "valence" in bad.reask_prompt
    assert good.ok and not good.mismatch


def test_frozen_scores_sum_and_no_w_needed():
    X, W, Y, spec = make_synthetic_shift(
        n=280, d_text=16, d_vad=4, gt="valence", mean_shift=1.3, seed=4,
    )
    models = fit_frozen_gap(X, W, spec, seed=4, n_estimators=30, light=False)
    scores = predict_gap_scores(models, X[:40])
    assert scores["instance_share"].shape == (40, 4)
    assert np.allclose(scores["instance_share"].sum(axis=1), 1.0, atol=1e-5)
    assert np.all(np.isfinite(scores["z"]))
    top = spec.names[int(np.argmax(models.pi))]
    assert top == "valence"


def test_budgeted_pi_recovers_valence_when_vimp_likes_text():
    """Wide text block + concentrated valence shift: π routes to valence."""
    X, W, Y, spec = make_synthetic_shift(
        n=480, d_text=28, d_vad=4, gt="valence", mean_shift=1.25, seed=9,
    )
    row = budgeted_inference_eval(
        X, W, spec, gt="valence", seed=9, n_estimators=40, test_size=0.35,
    )
    assert row["selected_block"]["pi_top1"] == "valence"
    assert row["hit_gt"]["pi_top1"] == 1.0
    assert row["hit_gt"]["blend_top1"] >= 0.55
    assert row["auc"]["pi_top1"] >= row["auc"]["random_row"] - 0.02
    assert row["auc"]["oracle_gt"] >= row["auc"]["random_row"]
    assert row["auc"]["blend_top1"] >= 0.70


def test_worked_packet_is_copy_pasteable():
    X, W, Y, spec = make_synthetic_shift(
        n=220, d_text=12, d_vad=3, gt="valence", mean_shift=1.2, seed=2,
    )
    pkt = worked_example_packet(X, W, spec, seed=2, n_estimators=25)
    assert "system_prompt" in pkt
    assert "Population gap shares" in pkt["system_prompt"]
    assert pkt["openai_tools"]
    assert "cited" in pkt["system_prompt"]
    dec = pkt["decision"]
    assert "tools_enabled" in dec
    json.dumps(pkt["decision"]["packet"])
    assert "channel: valence" in pkt["user_message"]
    assert "channel: text" not in pkt["user_message"]


def test_population_pack_same_expert_vs_pack_all():
    names = ["text", "valence", "arousal", "dominance"]
    cfg = RouterConfig(tau_hi=0.42, tau_lo=0.18)
    E, mode = population_pack(names, np.array([0.08, 0.75, 0.10, 0.07]), config=cfg)
    assert mode == "same_expert"
    assert E == ["valence"]
    E2, mode2 = population_pack(names, np.full(4, 0.25), config=cfg)
    assert mode2 == "pack_all"
    assert E2 == names


def test_context_pack_omits_disabled_and_concat_order():
    from clever_covariate_gap import ModalitySpec

    names = ["text", "valence", "arousal", "dominance"]
    spec = ModalitySpec(names=list(names), slices=[slice(0, 1)] * 4)
    cat = catalog_for_spec(spec)
    d = route_from_shares(
        names, np.array([0.1, 0.7, 0.1, 0.1]), np.array([0.08, 0.8, 0.07, 0.05]), cat,
    )
    assert d.pack_mode == "same_expert"
    assert d.packed_blocks == ["valence"]
    blob = pack_user_context(
        {"text": "TOKEN DUMP", "valence": "VAD=0.9", "arousal": "A", "dominance": "D"},
        d,
    )
    assert "VAD=0.9" in blob
    assert "TOKEN DUMP" not in blob
    X = np.arange(12, dtype=float).reshape(3, 4)
    spec2 = ModalitySpec(names=names, slices=[slice(0, 1), slice(1, 2), slice(2, 3), slice(3, 4)])
    packed = pack_blocks(X, spec2, ["valence", "text"])
    assert packed.shape == (3, 2)
    np.testing.assert_array_equal(packed[:, 0], X[:, 1])
    np.testing.assert_array_equal(packed[:, 1], X[:, 0])


def test_packed_reader_pi_beats_vimp_on_inject():
    from clever_covariate_gap import inject_block_shift

    X, W, Y, spec = make_synthetic_shift(
        n=480, d_text=28, d_vad=4, gt="text", mean_shift=0.12, seed=9,
    )
    X = inject_block_shift(X, W, spec, "valence", alpha=1.15)
    row = budgeted_inference_eval(
        X, W, spec, gt="valence", seed=9, n_estimators=40, test_size=0.35,
    )
    assert row["population_pack"]["mode"] == "same_expert"
    assert row["pack_auc"]["same_expert_pi"] >= row["pack_auc"]["same_expert_vimp"] - 0.02
    assert row["pack_auc"]["same_expert_pi"] >= 0.68
    assert row["pack_auc"]["abstain_policy"] == row["pack_auc"]["same_expert_pi"]


def test_diffuse_shift_pack_all_beats_one_expert():
    """Equal per-block signal: dropping to one specialist loses AUC. Packing all is the abstain analogue."""
    X, W, Y, spec = make_diffuse_equal_shift(n=420, d_text=20, d_vad=4, mean_shift=0.95, seed=3)
    row = budgeted_inference_eval(
        X, W, spec, gt=None, seed=3, n_estimators=35, test_size=0.35,
    )
    assert row["pack_auc"]["pack_all"] >= row["pack_auc"]["same_expert_pi"] - 0.005
    assert row["pack_auc"]["pack_all"] >= 0.70
    # Gate on a flat constructed π (not the DGP's π): dropping is unjustified.
    E, mode = population_pack(spec.names, np.full(4, 0.25), config=RouterConfig())
    assert mode == "pack_all" and E == spec.names


def test_query_update_changes_E_only_when_s_overrides_pi():
    names = ["text", "valence", "arousal", "dominance"]
    pi = np.array([0.10, 0.70, 0.12, 0.08])
    s_agree = np.array([0.08, 0.80, 0.07, 0.05])
    s_fight = np.array([0.82, 0.08, 0.05, 0.05])
    keep = query_update(names, pi, s_agree, lam=0.4)
    assert keep["E_prior"] == ["valence"]
    assert keep["switched"][0] == False
    fight0 = query_update(names, pi, s_fight, lam=0.0)
    assert fight0["switched"][0] == True
    assert fight0["E_query"][0][0] == "text"
    frozen = query_update(names, pi, s_fight, lam=1.0)
    assert frozen["switched"][0] == False
    flat = query_update(names, np.full(4, 0.25), np.full(4, 0.25), lam=0.4)
    assert flat["mode_prior"] == "pack_all"
    assert flat["switched"][0] == False


def test_block_then_adapt_is_two_stage():
    X, W, Y, spec = make_synthetic_shift(
        n=360, d_text=20, d_vad=4, gt="valence", mean_shift=1.2, seed=4,
    )
    row = budgeted_inference_eval(X, W, spec, gt="valence", seed=4, n_estimators=35)
    a = row["stage_A_block"]
    b = row["stage_B_adapt"]
    assert a["block_auc"]["valence"] >= 0.75
    assert a["block_auc"]["valence"] >= a["block_auc"]["text"] - 0.02
    assert b["adapt_prior_hit"] == 1.0
    assert "adapt_query_only_auc" in b
    assert 0.0 <= b["switch_rate"] <= 1.0
    # Noisy query-only gate should not beat the frozen prior on this DGP.
    assert b["adapt_prior_auc"] >= b["adapt_query_only_auc"] - 0.03
