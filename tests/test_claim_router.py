"""Tests for global claim router (overclaim check + fallback)."""
from __future__ import annotations

from agod.claim_router import (
    ClaimLevel,
    Evidence,
    harden_reading,
    max_supported_level,
    route_claim,
    route_ml_suite_claims,
    route_pack_interpretation,
    scan_overclaim_text,
)


def test_scan_flags_causal_and_skill():
    flags = scan_overclaim_text("persistent drivers show transferable skill")
    tags = {f["tag"] for f in flags}
    assert "causal_driver" in tags
    assert "skill_word" in tags


def test_l4_always_fallback():
    r = route_claim(
        "promote model to production — great ROI",
        Evidence(excess_auc=0.3, burn=True, rel_vs_uniform=0.9),
    )
    assert r["ok"] is False
    assert r["fallback_applied"] is True
    assert "causal / ship / ROI" in r["text"] or "diagnostic only" in r["text"]


def test_skill_word_blocked_without_excess():
    r = route_claim(
        "model has transferable skill",
        Evidence(excess_auc=0.01),
    )
    assert r["fallback_applied"] is True
    assert "skill" not in r["text"].lower() or "avoid" in r["text"].lower() or "near null" in r["text"]


def test_skill_allowed_with_excess_and_partial():
    r = route_claim(
        "excess indicates transferable skill under null",
        Evidence(excess_auc=0.12, excess_partial=0.10),
        requested_level=ClaimLevel.L2_EVIDENTIAL,
    )
    assert r["ok"] is True
    assert r["fallback_applied"] is False


def test_partial_collapse_blocks_skill():
    r = route_claim(
        "transferable skill",
        Evidence(excess_auc=0.20, excess_partial=0.0, delta_excess=0.20),
    )
    assert r["fallback_applied"] is True
    assert "partial" in r["text"].lower() or "Δexcess" in r["text"] or "delta" in r["text"].lower() or "skill" in r["text"].lower()


def test_pack_interpretation_no_persistent_drivers():
    r = route_pack_interpretation(mean_auc=0.9, mean_jaccard=0.7, mean_excess=0.15)
    assert "persistent drivers" not in r["text"]
    assert "correlates" in r["text"] or "transfer" in r["text"]


def test_harden_reading_short_ph_series():
    text = harden_reading(
        "alarm at index 3",
        Evidence(ph_alarm=True, ph_series_len=3),
    )
    assert "too short" in text or "statistic only" in text


def test_max_supported_levels():
    assert max_supported_level(Evidence()) == ClaimLevel.L0_DESCRIPTIVE
    assert max_supported_level(Evidence(mean_auc=0.7)) == ClaimLevel.L1_DIAGNOSTIC
    assert (
        max_supported_level(Evidence(excess_auc=0.1, excess_partial=0.08))
        == ClaimLevel.L2_EVIDENTIAL
    )


def test_route_ml_suite_claims_smoke():
    suite = {
        "excess_auc": 0.02,
        "partial_excess": {
            "delta_excess": 0.2,
            "excess_raw": 0.22,
            "excess_partial": 0.02,
        },
        "learning_curve": {"points": [{}, {}, {}, {}], "reading": "excess softens / noisy with n"},
        "page_hinkley": {"alarm": True, "n": 12, "reading": "skill-drop alarm at 6"},
    }
    out = route_ml_suite_claims(suite)
    assert "headline" in out
    assert out["headline"]["fallback_applied"] in (True, False)
    assert "skill" not in out["headline"]["text"].lower() or "not claim" in out["headline"]["text"].lower() or "do not claim" in out["headline"]["text"].lower()
