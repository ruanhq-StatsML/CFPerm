"""Balance / typed correlation / effective-rank: identification tests, no zip."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from typed_aux_losses import (  # noqa: E402
    amazon_heatmap_ranks,
    balance_loss,
    corr_with_batch,
    covariance_erank,
    effective_rank,
    heatmap_erank,
    pair_report,
    rank_modalities,
    report_typed_streams,
    spectrum_psd,
    two_by_two_gram_erank,
)
from typed_shift_stepsize import make_typed_stream  # noqa: E402


def test_effective_rank_rank1_vs_full():
    assert abs(effective_rank(np.array([1.0, 0.0, 0.0])) - 1.0) < 1e-8
    ones = np.ones((6, 6))
    assert heatmap_erank(ones) < 1.05
    eye = np.eye(5)
    assert abs(heatmap_erank(eye) - 5.0) < 1e-6
    rng = np.random.default_rng(0)
    X = rng.normal(size=(40, 8))
    assert covariance_erank(X) > 4.0
    X1 = np.outer(rng.normal(size=40), rng.normal(size=8))
    assert covariance_erank(X1) < 2.2


def test_spectrum_psd_symmetric():
    A = np.array([[2.0, 0.5], [0.1, 1.0]])
    w = spectrum_psd(A)
    assert w.min() >= -1e-10
    assert w.size == 2


def test_covariance_erank_location_invariant():
    rng = np.random.default_rng(1)
    X = rng.normal(size=(50, 16))
    assert abs(covariance_erank(X) - covariance_erank(X + 4.0)) < 1e-8


def test_balance_inv_upweights_lagging_head():
    uni = {"video": 0.4, "audio": 0.5, "text": 2.2}
    inv = balance_loss(uni, kind="inv")
    var = balance_loss(uni, kind="var")
    assert inv["weights"]["text"] > inv["weights"]["video"]
    assert var["gap"]["text"] > 0
    assert var["gap"]["video"] < 0


def test_corr_with_W_ranks_video_on_covariate_stream():
    stream = make_typed_stream(n_batches=6, n_per=40, seed=4, cov={"video": 0.15, "audio": 0.03, "text": 0.0})
    i0 = np.flatnonzero(stream.batch == 0)
    i1 = np.flatnonzero(stream.batch == 5)
    rho = corr_with_batch(stream.X[i0], stream.X[i1])
    assert abs(rho["video"]) > abs(rho["text"])
    assert abs(rho["video"]) > abs(rho["audio"])


def test_inv_ce_inverts_c_ranking_on_cov():
    recs = report_typed_streams(n_per=48, n_batches=8, seed=2026)
    cov = recs["cov_only"]
    assert cov["rank_c"][0] == "video"
    assert cov["rank_rho_W"][0] == "video"
    assert cov["rank_inv_weight"] != cov["rank_c"]
    assert cov["balance_inv"]["weights"]["text"] >= cov["balance_inv"]["weights"]["video"]


def test_typed_corr_fires_on_covariate_not_concept():
    recs = report_typed_streams(n_per=48, n_batches=8, seed=2026)
    cov_l = recs["cov_only"]["typed_corr"]["value"]
    con_l = recs["concept_only"]["typed_corr"]["value"]
    assert cov_l > 0.15
    assert con_l < 0.05
    assert cov_l > 10.0 * con_l
    assert recs["cov_only"]["typed_corr"]["terms"]["video>audio"] > recs["cov_only"]["typed_corr"]["terms"]["text>video"]


def test_two_by_two_gram_erank_monotone_in_cosine():
    assert two_by_two_gram_erank(1.0) < 1.05
    assert abs(two_by_two_gram_erank(0.0) - 2.0) < 1e-6
    assert two_by_two_gram_erank(0.2) > two_by_two_gram_erank(0.9)


def test_amazon_heatmap_ranks_gift_and_hops():
    labels = ["Gift", "Music", "Beauty", "Software", "SubBox", "Instr", "Mag", "Handmade", "Fashion"]
    R = np.eye(9)
    R[0, 1:] = 0.52
    R[1:, 0] = 0.52
    R[1:, 1:] = 0.85
    np.fill_diagonal(R, 1.0)
    rec = amazon_heatmap_ranks(R, labels=labels)
    assert rec["erank"] > 1.2
    assert rec["erank"] < 6.0
    assert rec["categories_by_loo_erank"][0]["label"] == "Gift"
    assert rec["hops_by_c"][0]["from"] == "Gift"
    by_c = [h["from"] + h["to"] for h in rec["hops_by_c"]]
    by_er = [h["from"] + h["to"] for h in sorted(rec["hops_by_c"], key=lambda r: -r["gram_erank"])]
    assert by_c == by_er


def test_rank_modalities_stable():
    assert rank_modalities({"video": 0.2, "audio": 0.2, "text": 0.1})[2] == "text"
    assert pair_report is not None
