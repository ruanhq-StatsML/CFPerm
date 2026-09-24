"""Tests for board → reason-code generation."""
from __future__ import annotations

from scripts.generate_board_reason_codes import CATALOG, generate_reason_codes


def _mini_summary():
    return {
        "ship_gate": {
            "promote_HGB_to_production": False,
            "reason": "probe only",
        },
        "ops_content_gap": [
            {
                "chunk": 1000,
                "ops_auc": 0.99,
                "content_auc": 0.80,
                "auc_gap_ops_minus_content": 0.19,
                "ops_top0": ["e_log1p_exp", "e_n_exp"],
                "content_top0": ["i_credit_last", "i_share_last"],
                "reading": "volume explains most transfer",
            }
        ],
        "cross_pack": [
            {
                "dataset": "diffusiondb",
                "chunk": 1000,
                "mean_auc": 0.60,
                "mean_logreg_auc": 0.61,
                "mean_delta_Y": 0.01,
                "mean_jaccard": 0.2,
                "mean_fsds_cmean_jaccard": 0.0,
                "reading": "weak transfer",
            },
            {
                "dataset": "tencent_gr",
                "chunk": 1000,
                "mean_auc": 0.99,
                "mean_logreg_auc": 0.83,
                "mean_delta_Y": -0.01,
                "mean_jaccard": 0.7,
                "mean_fsds_cmean_jaccard": 0.05,
                "reading": "strong transfer + stable",
            },
            {
                "dataset": "tencent_gr_content",
                "chunk": 1000,
                "mean_auc": 0.84,
                "mean_logreg_auc": 0.84,
                "mean_delta_Y": -0.01,
                "mean_jaccard": 0.55,
                "mean_fsds_cmean_jaccard": 0.3,
                "reading": "transfer holds",
            },
            {
                "dataset": "metro_interstate",
                "chunk": 1000,
                "mean_auc": 0.97,
                "mean_logreg_auc": 0.60,
                "mean_delta_Y": 10.0,
                "mean_jaccard": 0.22,
                "mean_fsds_cmean_jaccard": 0.25,
                "reading": "shifting drivers",
            },
        ],
        "datasets": {
            "diffusiondb": {
                "ok": True,
                "by_chunk": {
                    "1000": {
                        "rows": [
                            {"top_fsds": ["mucha", "alphonse"], "hgb_auc": 0.6}
                        ]
                    }
                },
            },
            "tencent_gr": {
                "ok": True,
                "by_chunk": {
                    "1000": {
                        "rows": [
                            {
                                "top_fsds": ["e_log1p_exp", "e_n_exp", "i_log1p_n_exp"],
                                "hgb_auc": 0.99,
                            }
                        ]
                    }
                },
            },
            "tencent_gr_content": {
                "ok": True,
                "by_chunk": {
                    "1000": {
                        "rows": [
                            {
                                "top_fsds": [
                                    "i_credit_last",
                                    "i_share_last",
                                    "i_item_credit_rank",
                                ],
                                "hgb_auc": 0.84,
                            }
                        ]
                    }
                },
            },
            "metro_interstate": {
                "ok": True,
                "by_chunk": {
                    "1000": {"rows": [{"top_fsds": ["m0", "m1"], "hgb_auc": 0.97}]}
                },
            },
        },
    }


def test_catalog_codes_are_stable():
    assert "RC_BOARD_NOT_SHIP" in CATALOG
    assert "RC_CONTENT_CREDIT_CANDIDATE" in CATALOG


def test_generate_emits_core_codes():
    payload = generate_reason_codes(_mini_summary())
    codes = {c["code"] for c in payload["codes"]}
    assert "RC_BOARD_NOT_SHIP" in codes
    assert "RC_SHORTLIST_NOT_DRIVER" in codes
    assert "RC_DIFFDB_TOKEN_NO_TRAVEL" in codes
    assert "RC_INTENSITY_BASELINE" in codes
    assert "RC_CONTENT_CREDIT_CANDIDATE" in codes
    assert "RC_OPS_CONTENT_GAP" in codes
    assert "RC_SHIFTING_DRIVERS" in codes
    assert "RC_PROBE_DISAGREE" in codes  # metro HGB vs LogReg
    assert all(c["allows_ship_model"] is False for c in payload["codes"])
    assert "RC_BOARD_NOT_SHIP" in payload["paste_for_agent"]


def test_block_claims_sorted_first():
    payload = generate_reason_codes(_mini_summary())
    sevs = [c["severity"] for c in payload["codes"]]
    assert sevs[0] == "block_claim"
    assert sevs[1] == "block_claim"


def test_ad_scenario_from_tencent_neg_dy():
    s = _mini_summary()
    # mean_delta_Y already -0.01 on tencent in mini → neg → S2
    payload = generate_reason_codes(s)
    ad = payload["ad_scenario"]
    assert ad["ok"] is True
    assert ad["family_code"] == "S2_inject"
    assert ad["sign_Dy_board"] == "neg"
    assert ad["allows_ship_model"] is False
    codes = {c["code"] for c in payload["codes"]}
    assert "RC_AD_FUNNEL_CONTEXT" in codes
    assert "RC_AD_BUY_INTENSITY" in codes
    assert "RC_AD_LAST_TOUCH_CANDIDATE" in codes
    assert "RC_AD_CONVERT_DIP" in codes
    assert "广告漏斗" in payload["paste_for_agent"] or "广告" in payload["paste_for_agent"]
