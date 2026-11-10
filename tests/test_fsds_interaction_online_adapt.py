"""Unit tests for FSDS × Interaction Value online adaptation."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from fsds_interaction_online_adapt import (  # noqa: E402
    P_X,
    SHARED_DIM,
    adapt_scores,
    clever_covariate_W,
    fsds_lambda,
    interaction_at,
    interaction_map,
    interaction_value_residual,
    interpret_hop,
    make_trimodal_stream,
    rf_hop_po_risk,
    run_method,
    run_suite,
    to_shared_embedding,
    token_index,
    triangle_inconsistency,
)


def test_stream_shapes_and_align():
    stream = make_trimodal_stream(scenario="healthy", n_batches=4, n_per=40, seed=1)
    assert stream.X.shape == (160, P_X)
    assert stream.y.shape == (160,)
    assert set(stream.align) == {"TI", "TA", "IA"}
    assert stream.align["TI"].shape == (160,)


def test_interaction_value_positive_when_joint_helps():
    err_a = np.array([1.0, 1.0, 1.0])
    err_b = np.array([1.0, 1.0, 1.0])
    err_j = np.array([0.2, 0.3, 0.1])
    phi = interaction_value_residual(err_a, err_b, err_j)
    assert float(phi.mean()) > 1.0


def test_shared_embedding_and_token_index():
    stream = make_trimodal_stream(scenario="healthy", n_batches=2, n_per=20, seed=2)
    Z, proj = to_shared_embedding(stream.X)
    assert Z.shape == (40, 3 * SHARED_DIM)
    assert token_index("text", 0) == 0
    assert token_index("video", 0) == SHARED_DIM
    assert token_index("audio", SHARED_DIM - 1) == 3 * SHARED_DIM - 1
    Z2, _ = to_shared_embedding(stream.X[:5], projections=proj)
    assert Z2.shape[1] == Z.shape[1]


def test_interaction_at_position_api():
    """Input = dim/token indices in shared space; output = Φ at that position."""
    rng = np.random.default_rng(0)
    n, d = 30, SHARED_DIM
    # build Z where score = z_text[0]*z_video[0] + noise  ⇒ Φ at (text0, video0) > 0
    Zt = rng.normal(size=(n, d))
    Zv = rng.normal(size=(n, d))
    Za = rng.normal(size=(n, d)) * 0.1
    Z = np.concatenate([Zt, Zv, Za], axis=1)
    coef_t0 = token_index("text", 0)
    coef_v0 = token_index("video", 0)

    def predict_fn(Zx):
        return Zx[:, coef_t0] * Zx[:, coef_v0]

    phi = interaction_at(predict_fn, Z, coef_t0, coef_v0)
    # for product f=x_i x_j: Φ = x_i x_j − 0 − 0 + 0 = x_i x_j
    assert np.allclose(phi, Z[:, coef_t0] * Z[:, coef_v0])
    # off-diagonal unrelated dims ≈ 0 for this f
    phi_noise = interaction_at(predict_fn, Z, token_index("text", 3), token_index("audio", 3))
    assert float(np.mean(np.abs(phi_noise))) < 1e-10
    amap = interaction_map(predict_fn, Z, "text", "video")
    assert amap.shape == (n, d)
    assert np.allclose(amap[:, 0], phi)


def test_triangle_flags_mismatched_audio():
    healthy = triangle_inconsistency({"TI": 0.8, "TA": 0.75, "IA": 0.7})
    attacked = triangle_inconsistency({"TI": 0.85, "TA": 0.05, "IA": 0.02})
    assert attacked > healthy
    assert attacked > 0.5


def test_alignment_attack_raises_triangle():
    h = make_trimodal_stream(scenario="healthy", n_batches=6, n_per=50, seed=3)
    a = make_trimodal_stream(scenario="alignment_attack", n_batches=6, n_per=50, seed=3)
    # late batches only
    late = a.batch >= a.meta["break_start"]
    tri_h = np.array(
        [
            triangle_inconsistency({k: h.align[k][i] for k in ("TI", "TA", "IA")})
            for i in np.flatnonzero(h.batch >= h.meta["break_start"])
        ]
    )
    tri_a = np.array(
        [
            triangle_inconsistency({k: a.align[k][i] for k in ("TI", "TA", "IA")})
            for i in np.flatnonzero(late)
        ]
    )
    assert float(tri_a.mean()) > float(tri_h.mean())


def test_adapt_scores_shrinks_under_shift():
    raw = np.array([2.0, -1.5, 0.8])
    W = np.array([1.0, 2.0, 0.5])
    adapted, pen = adapt_scores(raw, W, lam=1.0, is_shift=True)
    assert np.all(np.abs(adapted) <= np.abs(raw) + 1e-12)
    assert float(pen.mean()) > 0
    same, _ = adapt_scores(raw, W, lam=1.0, is_shift=False)
    assert np.allclose(same, raw)


def test_fsds_lambda_reacts_to_concept():
    # high delta → larger λ than high c / quiet delta
    lam_concept = fsds_lambda(c=0.1, delta=0.5, pi_change=0.0, triangle_mean=0.2, phi_drop=0.1)
    lam_cov = fsds_lambda(c=0.9, delta=0.02, pi_change=0.0, triangle_mean=0.05, phi_drop=0.0)
    assert lam_concept > lam_cov


def test_clever_W_nonnegative():
    W = clever_covariate_W(
        po_risks=np.array([0.2, 0.5, 0.9]),
        triangle_gap=np.array([0.1, 0.4, 0.8]),
        phi_mean=np.array([0.5, 0.2, 0.0]),
        ref_po=[0.2, 0.25],
        ref_tri=[0.15],
        ref_phi=[0.6, 0.55],
    )
    assert W.shape == (3,)
    assert np.all(W >= 0)


def test_batch_phi_nonzero_with_shared_products():
    from fsds_interaction_online_adapt import compute_batch_interaction

    stream = make_trimodal_stream(scenario="healthy", n_batches=3, n_per=60, seed=5)
    tr = stream.batch == 0
    te = stream.batch == 1
    _, _, phi, _, phi_maps, _ = compute_batch_interaction(
        stream.X[te], stream.y[te], stream.X[tr], stream.y[tr]
    )
    assert phi_maps["TI"].shape[1] == SHARED_DIM
    # products make Φ identifiable — not identically zero on complementary stream
    assert float(np.mean(np.abs(phi["TI"]))) > 1e-6


def test_rf_hop_po_risk_outcome_model():
    stream = make_trimodal_stream(scenario="concept_drift", n_batches=4, n_per=50, seed=8)
    Z, _ = to_shared_embedding(stream.X)
    b0 = stream.batch == 0
    b1 = stream.batch == 1
    pack = rf_hop_po_risk(Z[b0], stream.y[b0], Z[b1], stream.y[b1], seed=8, n_estimators=20)
    assert pack["po_abs"].shape == (int(b1.sum()),)
    assert pack["r_po"] >= 0.0
    assert abs(sum(pack["pi_po"].values()) - 1.0) < 1e-6
    assert set(pack["pi_po"]) == {"text", "video", "audio"}


def test_interpret_hop_has_read():
    phi_maps = {
        "TI": np.zeros((10, SHARED_DIM)),
        "TA": np.zeros((10, SHARED_DIM)),
        "IA": np.zeros((10, SHARED_DIM)),
    }
    phi_maps["TI"][:, 2] = 0.4
    pack = {"r_po": 0.1, "mean_abs_po": 0.2, "pi_po": {"text": 0.2, "video": 0.3, "audio": 0.5}}
    out = interpret_hop(
        phi_maps,
        {"text": 0.2, "video": 0.5, "audio": 0.3},
        pack["pi_po"],
        np.linspace(0.6, 0.9, 10),
        pack,
        {"c": 0.1, "delta": 0.2},
        phi_drop=0.0,
        r_po_ratio=1.4,
    )
    assert out["dominant_po"] == "audio"
    assert "top_tokens" in out and out["top_tokens"]["TI"][0]["token"] == 2
    assert "PO" in out["read"] or "concept" in out["read"]


def test_fsds_clever_beats_localize_on_compound():
    stream = make_trimodal_stream(scenario="compound", n_batches=6, n_per=80, seed=11)
    loc = run_method(stream, method="localize_only", seed=11, n_estimators=20)
    fsds = run_method(stream, method="fsds_clever_cov", seed=11, n_estimators=20)
    # treatment should not be worse than diagnosis-only by a large margin
    assert fsds["online_mse"] <= loc["online_mse"] + 0.08
    assert "interpretations" in fsds and len(fsds["interpretations"]) == 5
    assert "r_po" in fsds["history"][0]


def test_suite_quick_keys():
    suite = run_suite(quick=True, seeds=[2026])
    assert "compound" in suite["summary"]
    assert "fsds_clever_cov" in suite["summary"]["compound"]
    assert "interpret_board" in suite
    assert suite["interpret_board"].get("compound")
    assert suite["lifts"]["compound"]["fsds_clever_cov"]["vs_localize"] == suite["lifts"][
        "compound"
    ]["fsds_clever_cov"]["vs_localize"]
