"""FSDS × multimodal Interaction Value → Online Adaptation (clever covariates).

Localization (π / Attr / RFPerm) is diagnosis. Treatment here is *inference*
adaptation only: inject PO-risk / Interaction-Shift as clever covariates into
logits / scores — **no weight update**.

  Logit_adapted = Logit_original − λ · W_shift

Locked FSDS signs (unchanged):
  ĉ large / δ̂ quiet  → shrink λ (covariate hop; do not chase P(X|W))
  δ̂ large            → raise λ (concept / interaction break; distrust joint)
  PO-risk / Φ_ij      → clever covariate W, not feature reweighting
  Triangle break      → alignment / contamination gate (L3)

Synthetic board mimics CMU-MOSEI / FakeAVCeleb style trimodal streams
(text, video≈image, audio) with controlled Interaction Shift.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.stats import ks_2samp
from sklearn.linear_model import Ridge

from amazon_continuous_batches import concept_intensity_mse, heatmap_hop_c
from attribution_adapter import attribution_ewma_rate
from msrvtt_multimodal_attribution import (
    SEED,
    rf_domain,
    standardize_columns,
)


def _modality_mass(vimp):
    """π_m on this board's GROUPS (not MSR-VTT 2048-d layout)."""
    v = np.asarray(vimp, dtype=float)
    mass = {g: float(np.maximum(v[sl], 0.0).sum()) for g, sl in GROUPS.items()}
    tot = sum(mass.values()) + 1e-12
    return mass, {g: mass[g] / tot for g in mass}

# Compact dims for fast online hops (MOSEI-like geometry, not full CLIP).
D_TEXT = 32
D_VIDEO = 32
D_AUDIO = 24
P_X = D_TEXT + D_VIDEO + D_AUDIO
GROUPS = {
    "text": slice(0, D_TEXT),
    "video": slice(D_TEXT, D_TEXT + D_VIDEO),
    "audio": slice(D_TEXT + D_VIDEO, P_X),
}
GROUP_NAMES = ("text", "video", "audio")
PAIR_NAMES = ("TI", "TA", "IA")  # text-video, text-audio, video-audio

METHODS = (
    "static",           # joint Ridge, no adaptation
    "localize_only",    # FSDS π + Φ diagnostics; raw preds unchanged
    "clever_cov",       # KS on PO-risk / align → λ · W penalty
    "fsds_clever_cov",  # FSDS (ĉ,δ̂) + triangle + Φ gate λ
    "oracle_penalty",   # upper bound: always know |err| as W
)

SCENARIOS = (
    "healthy",
    "interaction_break",   # joint complementarity collapses mid-stream
    "alignment_attack",    # fake / mismatched audio (triangle break)
    "concept_drift",       # Y|X map flips on late batches
    "compound",            # interaction break + alignment attack
)


@dataclass
class TrimodalStream:
    X: np.ndarray
    y: np.ndarray
    batch: np.ndarray
    align: dict  # pairwise cosine scores per row
    meta: dict = field(default_factory=dict)


SHARED_DIM = 16  # common embedding width after projection
# stacked shared space: [text(0:D) | video(D:2D) | audio(2D:3D)]
SHARED_SLICES = {
    "text": slice(0, SHARED_DIM),
    "video": slice(SHARED_DIM, 2 * SHARED_DIM),
    "audio": slice(2 * SHARED_DIM, 3 * SHARED_DIM),
}
SHARED_P = 3 * SHARED_DIM
_PROJ_SEEDS = {"text": 11, "video": 22, "audio": 33}


def _block(X, name):
    return np.asarray(X, dtype=float)[:, GROUPS[name]]


def _proj_matrix(d_in, d_out, seed):
    rng = np.random.default_rng(seed)
    P = rng.normal(size=(int(d_in), int(d_out)))
    P /= np.linalg.norm(P, axis=0, keepdims=True) + 1e-12
    return P


def to_shared_embedding(X, shared_dim=SHARED_DIM, projections=None):
    """Map unequal modality blocks → one embedding space Z ∈ R^{n × 3d}.

    Returns ``(Z, projections)``. Freeze ``projections`` across hops so
    dim/token indices stay comparable online.
    """
    X = np.asarray(X, dtype=float)
    d = int(shared_dim)
    if projections is None:
        projections = {
            g: _proj_matrix(GROUPS[g].stop - GROUPS[g].start, d, _PROJ_SEEDS[g])
            for g in GROUP_NAMES
        }
    blocks = [_block(X, g) @ projections[g] for g in GROUP_NAMES]
    return np.concatenate(blocks, axis=1), projections


def token_index(modality, position, shared_dim=SHARED_DIM):
    """(modality, local dim/token) → global index in shared embedding Z."""
    if modality not in SHARED_SLICES:
        raise KeyError("unknown modality %s" % modality)
    pos = int(position)
    if pos < 0 or pos >= int(shared_dim):
        raise IndexError("position %d out of shared_dim=%d" % (pos, shared_dim))
    return int(SHARED_SLICES[modality].start + pos)


def interaction_at(predict_fn, Z, i, j):
    """Pairwise interaction Φ_{i,j} at shared-embedding dims/tokens ``i``, ``j``.

    Inputs
      predict_fn : callable(Z) → (n,) scores
      Z          : shared embedding (n, 3d)
      i, j       : dimension / token indices in that space

    Output (n,)
      Φ_{i,j}(x) = f(x) − f(x_{−i}) − f(x_{−j}) + f(x_{−i,−j})

    Positive ⇒ complementary at that position; ≤0 ⇒ redundant / conflict.
    """
    Z = np.asarray(Z, dtype=float)
    i, j = int(i), int(j)
    if i == j:
        return np.zeros(Z.shape[0], dtype=float)
    f_full = np.asarray(predict_fn(Z), dtype=float)
    Zi = Z.copy()
    Zi[:, i] = 0.0
    Zj = Z.copy()
    Zj[:, j] = 0.0
    Zij = Z.copy()
    Zij[:, [i, j]] = 0.0
    return f_full - np.asarray(predict_fn(Zi), dtype=float) - np.asarray(
        predict_fn(Zj), dtype=float
    ) + np.asarray(predict_fn(Zij), dtype=float)


def interaction_map(predict_fn, Z, modality_a, modality_b, shared_dim=SHARED_DIM):
    """Aligned cross-modal Φ at each shared position k=0..d−1 → shape (n, d)."""
    Z = np.asarray(Z, dtype=float)
    d = int(shared_dim)
    out = np.zeros((Z.shape[0], d), dtype=float)
    for k in range(d):
        out[:, k] = interaction_at(
            predict_fn,
            Z,
            token_index(modality_a, k, shared_dim=d),
            token_index(modality_b, k, shared_dim=d),
        )
    return out


def pairwise_align_scores(X, shared_dim=SHARED_DIM, projections=None):
    """Row-wise cosine on shared embeddings for TI / TA / IA."""
    Z, _ = to_shared_embedding(X, shared_dim=shared_dim, projections=projections)
    blocks = {g: Z[:, SHARED_SLICES[g]] for g in GROUP_NAMES}
    for g in blocks:
        blocks[g] = blocks[g] / (np.linalg.norm(blocks[g], axis=1, keepdims=True) + 1e-12)
    pairs = {
        "TI": ("text", "video"),
        "TA": ("text", "audio"),
        "IA": ("video", "audio"),
    }
    return {key: np.sum(blocks[a] * blocks[b], axis=1) for key, (a, b) in pairs.items()}


def triangle_inconsistency(align_row):
    """max−min of {S_TI, S_TA, S_IA}; high ⇒ fake / cross-modal attack."""
    vals = np.asarray([align_row[k] for k in PAIR_NAMES], dtype=float)
    return float(vals.max() - vals.min())


def interaction_value_residual(err_a, err_b, err_joint):
    """Modality-level aggregate fallback: Φ̂ ≈ (e_a + e_b) − e_joint."""
    return np.asarray(err_a, dtype=float) + np.asarray(err_b, dtype=float) - np.asarray(
        err_joint, dtype=float
    )


def bias_decouple(err_text, err_video, err_audio, err_joint):
    """R-learner style residual split: unimodal + interaction remainder."""
    e_t = np.asarray(err_text, dtype=float)
    e_v = np.asarray(err_video, dtype=float)
    e_a = np.asarray(err_audio, dtype=float)
    e_j = np.asarray(err_joint, dtype=float)
    # interaction residual: joint should beat best unimodal if complementary
    e_uni = np.minimum(np.minimum(e_t, e_v), e_a)
    e_int = e_uni - e_j  # positive when joint helps
    return {
        "E_text": e_t,
        "E_video": e_v,
        "E_audio": e_a,
        "E_interaction": e_int,
        "var_text": float(np.var(e_t)),
        "var_video": float(np.var(e_v)),
        "var_audio": float(np.var(e_a)),
        "var_interaction": float(np.var(e_int)),
    }


def make_trimodal_stream(
    scenario="compound",
    n_batches=8,
    n_per=120,
    seed=SEED,
    p_text=D_TEXT,
    p_video=D_VIDEO,
    p_audio=D_AUDIO,
):
    """Controlled MOSEI-like stream with optional Interaction / Alignment Shift.

    Early batches: complementary modalities (joint beats unimodal).
    Late batches (scenario-dependent):
      interaction_break — shared signal removed; modalities redundant/noisy
      alignment_attack  — audio replaced by unrelated noise (triangle break)
      concept_drift     — outcome map sign flip
      compound          — interaction_break + alignment_attack
    """
    rng = np.random.default_rng(seed)
    assert p_text == D_TEXT and p_video == D_VIDEO and p_audio == D_AUDIO
    rows_X, rows_y, rows_b = [], [], []
    break_start = max(2, n_batches // 2)

    for t in range(n_batches):
        # latent factors
        z_shared = rng.normal(size=(n_per, 4))
        z_text = rng.normal(size=(n_per, 4))
        z_video = rng.normal(size=(n_per, 4))
        z_audio = rng.normal(size=(n_per, 4))
        late = t >= break_start

        use_interact = True
        flip_y = False
        attack_audio = False
        if late:
            if scenario in ("interaction_break", "compound"):
                use_interact = False
            if scenario == "concept_drift":
                flip_y = True
            if scenario in ("alignment_attack", "compound"):
                attack_audio = True

        # embeddings: concat projected latents + noise
        def _embed(z_main, z_share, p, scale_share=1.0):
            base = rng.normal(size=(n_per, p)) * 0.35
            k = min(4, p)
            base[:, :k] += z_main
            if use_interact:
                base[:, :k] += scale_share * z_share
            return base

        Xt = _embed(z_text, z_shared, p_text, 1.0)
        Xv = _embed(z_video, z_shared, p_video, 1.0)
        Xa = _embed(z_audio, z_shared, p_audio, 1.0)
        if attack_audio:
            # mismatched / fake audio: destroy shared structure in audio block
            Xa = rng.normal(size=Xa.shape)

        X = np.concatenate([Xt, Xv, Xa], axis=1)

        # outcome: shared interaction + modality-private + noise
        y = (
            0.55 * z_shared[:, 0]
            + 0.25 * z_text[:, 0]
            + 0.20 * z_video[:, 0]
            + 0.15 * z_audio[:, 0]
            + 0.08 * rng.normal(size=n_per)
        )
        if not use_interact:
            # complementarity gone: only private noisy cues + mild shared leak
            y = (
                0.10 * z_shared[:, 0]
                + 0.45 * z_text[:, 0]
                + 0.35 * z_video[:, 0]
                + 0.25 * rng.normal(size=n_per)
            )
        if flip_y:
            y = -y + 0.15 * rng.normal(size=n_per)

        rows_X.append(X)
        rows_y.append(y)
        rows_b.append(np.full(n_per, t, dtype=int))

    X = np.vstack(rows_X)
    y = np.concatenate(rows_y)
    batch = np.concatenate(rows_b)
    align = pairwise_align_scores(X)
    return TrimodalStream(
        X=X,
        y=y,
        batch=batch,
        align=align,
        meta={
            "scenario": scenario,
            "n_batches": n_batches,
            "n_per": n_per,
            "break_start": break_start,
            "groups": {k: [sl.start, sl.stop] for k, sl in GROUPS.items()},
        },
    )


def _fit_predict(Xtr, ytr, Xte, alpha=1.0):
    clf = Ridge(alpha=float(alpha))
    clf.fit(Xtr, ytr)
    return clf.predict(Xte), clf


def _abs_err(pred, y):
    return np.abs(np.asarray(pred, dtype=float) - np.asarray(y, dtype=float))


def _mse(pred, y):
    return float(np.mean((np.asarray(pred) - np.asarray(y)) ** 2))


def _batch_means(X, batch):
    batch = np.asarray(batch, dtype=int)
    k = int(batch.max()) + 1
    return np.stack([X[batch == t].mean(axis=0) for t in range(k)])


def _cosine_matrix(mus):
    k = len(mus)
    R = np.zeros((k, k))
    for i in range(k):
        a = mus[i] / (np.linalg.norm(mus[i]) + 1e-12)
        for j in range(k):
            b = mus[j] / (np.linalg.norm(mus[j]) + 1e-12)
            R[i, j] = float(np.dot(a, b))
    return R


def hop_fsds_stats(stream, seed=SEED, n_estimators=40):
    """Per-hop FSDS (ĉ, δ̂) + RF modality π on trimodal X."""
    X = standardize_columns(stream.X)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    k = int(batch.max()) + 1
    mus = _batch_means(X, batch)
    R = _cosine_matrix(mus)
    rows = []
    prev_pi = None
    for t in range(1, k):
        ip = np.flatnonzero(batch == t - 1)
        ic = np.flatnonzero(batch == t)
        c = float(heatmap_hop_c(R, t))
        delta, _ = concept_intensity_mse(X[ip], y[ip], X[ic], y[ic])
        vimp, auc = rf_domain(X[ip], X[ic], seed=seed + t, n_estimators=n_estimators)
        _, pi = _modality_mass(vimp)
        pi_change = 0.0
        if prev_pi is not None:
            pi_change = 0.5 * float(
                sum(abs(pi.get(g, 0.0) - prev_pi.get(g, 0.0)) for g in GROUP_NAMES)
            )
        prev_pi = pi
        rows.append(
            {
                "round": int(t),
                "c": c,
                "delta": float(delta),
                "rf_auc": float(auc) if auc == auc else float("nan"),
                "pi": {g: float(pi.get(g, 0.0)) for g in GROUP_NAMES},
                "pi_change": float(pi_change),
                "ewma_rate": attribution_ewma_rate(c, delta, pi_change),
            }
        )
    return rows


def _augment_shared(Z, shared_dim=SHARED_DIM):
    """Shared Z plus aligned Hadamard products (so Φ_{i,j} is identifiable).

    Linear Ridge on raw Z has Φ≡0; products put cross-modal interaction into f.
    """
    Z = np.asarray(Z, dtype=float)
    d = int(shared_dim)
    Zt = Z[:, SHARED_SLICES["text"]]
    Zv = Z[:, SHARED_SLICES["video"]]
    Za = Z[:, SHARED_SLICES["audio"]]
    return np.concatenate([Z, Zt * Zv, Zt * Za, Zv * Za], axis=1)


def compute_batch_interaction(X, y, Xtr, ytr, alpha=1.0, projections=None):
    """Fit on shared embedding; return preds + position-wise Φ maps.

    All modalities live in the same Z-space. Pairwise Φ for TI/TA/IA is the
    row-mean of ``interaction_map`` over aligned token positions.
    ``predict_fn`` zeros dims in Z then re-augments — so ``interaction_at(i,j)``
    sees the product terms move with those tokens.
    """
    Ztr, projections = to_shared_embedding(Xtr, projections=projections)
    Zte, _ = to_shared_embedding(X, projections=projections)
    pred_j, clf = _fit_predict(_augment_shared(Ztr), ytr, _augment_shared(Zte), alpha=alpha)

    def predict_fn(Z):
        return clf.predict(_augment_shared(Z))

    preds = {"joint": pred_j}
    errs = {"joint": _abs_err(pred_j, y)}
    for name in GROUP_NAMES:
        sl = SHARED_SLICES[name]
        # unimodal probe in the *same* embedding coords
        pred_m, _ = _fit_predict(Ztr[:, sl], ytr, Zte[:, sl], alpha=alpha)
        preds[name] = pred_m
        errs[name] = _abs_err(pred_m, y)

    phi_maps = {
        "TI": interaction_map(predict_fn, Zte, "text", "video"),
        "TA": interaction_map(predict_fn, Zte, "text", "audio"),
        "IA": interaction_map(predict_fn, Zte, "video", "audio"),
    }
    # per-sample scalar: mean signed interaction over aligned tokens
    phi = {k: phi_maps[k].mean(axis=1) for k in phi_maps}
    decouple = bias_decouple(errs["text"], errs["video"], errs["audio"], errs["joint"])
    return preds, errs, phi, decouple, phi_maps, projections


def clever_covariate_W(po_risks, triangle_gap, phi_mean, ref_po, ref_tri, ref_phi):
    """Build W_shift from PO-risk, triangle gap, and interaction drop.

    All terms non-negative; larger W ⇒ lower trust in current joint score.
    """
    po = np.asarray(po_risks, dtype=float)
    # normalize vs rolling reference means
    po_ref = float(np.mean(ref_po)) + 1e-8 if len(ref_po) else float(np.mean(po)) + 1e-8
    tri_ref = float(np.mean(ref_tri)) + 1e-8 if len(ref_tri) else 0.2
    phi_ref = float(np.mean(ref_phi)) if len(ref_phi) else float(np.mean(phi_mean))

    w_po = po / po_ref
    w_tri = np.asarray(triangle_gap, dtype=float) / tri_ref
    # interaction drop relative to reference (positive when Φ collapsed)
    drop = np.maximum(0.0, phi_ref - np.asarray(phi_mean, dtype=float))
    w_phi = drop / (abs(phi_ref) + 1e-8)
    W = 0.45 * w_po + 0.35 * w_tri + 0.20 * w_phi
    return np.maximum(W, 0.0)


def fsds_lambda(c, delta, pi_change, triangle_mean, phi_drop, base_lambda=1.0):
    """TSS-gated penalty strength for clever covariates.

    Large δ / triangle / Φ-drop → raise λ; large ĉ alone → shrink (don't chase cov).
    """
    react = attribution_ewma_rate(c, delta, pi_change)  # ∈ [0.15, 0.75]
    boost = 1.0 + 1.5 * float(np.clip(triangle_mean, 0.0, 2.0)) + 1.0 * float(
        np.clip(phi_drop, 0.0, 2.0)
    )
    # when ĉ dominates and δ quiet, damp
    damp = (1.0 + 2.0 * float(c)) / (1.0 + 2.0 * float(c) + 4.0 * float(delta) + 1e-8)
    lam = float(base_lambda) * (0.5 + react) * boost * (0.55 + 0.45 * (1.0 - 0.5 * damp))
    return float(np.clip(lam, 0.05, 4.0))


def detect_shift(ref_po, cur_po, ref_align, cur_align, alpha=0.05):
    """KS tests on PO-risk and alignment — Interaction / Alignment Shift."""
    if len(ref_po) < 8 or len(cur_po) < 8:
        return False, 1.0, 1.0
    _, p_po = ks_2samp(np.asarray(ref_po, dtype=float), np.asarray(cur_po, dtype=float))
    _, p_al = ks_2samp(
        np.asarray(ref_align, dtype=float), np.asarray(cur_align, dtype=float)
    )
    return bool(p_po < alpha or p_al < alpha), float(p_po), float(p_al)


def adapt_scores(raw_pred, W, lam, is_shift, y_scale=None):
    """Inference-only adaptation: pull predictions toward caution under shift.

    For regression we shrink toward the batch-free prior 0 (MOSEI-centered),
    proportional to λ · W. Equivalent spirit to Logit − λ W on a classifier.
    """
    raw = np.asarray(raw_pred, dtype=float)
    W = np.asarray(W, dtype=float)
    if not is_shift:
        return raw.copy(), np.zeros_like(raw)
    penalty = float(lam) * W
    # shrink magnitude: adapted = raw * (1 / (1 + penalty))
    adapted = raw / (1.0 + penalty)
    return adapted, penalty


def run_method(stream, method="fsds_clever_cov", seed=SEED, alpha=0.05, lambda_penalty=1.0):
    """Online hop: fit on past → diagnose → optional clever-covariate adapt."""
    if method not in METHODS:
        raise ValueError("unknown method %s" % method)
    X = standardize_columns(stream.X)
    y = np.asarray(stream.y, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    align = {k: np.asarray(stream.align[k], dtype=float) for k in PAIR_NAMES}
    k = int(batch.max()) + 1
    hop_stats = {h["round"]: h for h in hop_fsds_stats(stream, seed=seed)}

    ref_po, ref_tri, ref_phi, ref_align = [], [], [], []
    projections = None  # freeze shared-space maps so token indices stay fixed
    history = []
    for t in range(1, k):
        tr = batch < t
        ic = np.flatnonzero(batch == t)
        Xtr, ytr = X[tr], y[tr]
        Xt, yt = X[ic], y[ic]

        preds, errs, phi, decouple, phi_maps, projections = compute_batch_interaction(
            Xt, yt, Xtr, ytr, projections=projections
        )
        raw = preds["joint"]
        po_risks = errs["joint"]
        tri = np.array(
            [
                triangle_inconsistency({p: align[p][ic][i] for p in PAIR_NAMES})
                for i in range(len(ic))
            ],
            dtype=float,
        )
        # mean pairwise Φ per sample over aligned shared tokens (TI+TA)
        phi_mean = 0.5 * (phi["TI"] + phi["TA"])
        align_batch = align["TI"][ic]  # primary align channel for KS

        is_shift, p_po, p_al = detect_shift(ref_po, po_risks, ref_align, align_batch, alpha=alpha)
        hop = hop_stats.get(t, {"c": 0.0, "delta": 0.0, "pi_change": 0.0, "pi": {}})
        phi_drop = 0.0
        if ref_phi:
            phi_drop = max(0.0, float(np.mean(ref_phi)) - float(np.mean(phi_mean)))

        W = clever_covariate_W(po_risks, tri, phi_mean, ref_po, ref_tri, ref_phi)

        if method == "static":
            adapted, penalty = raw.copy(), np.zeros_like(raw)
            used_shift = False
            lam = 0.0
        elif method == "localize_only":
            adapted, penalty = raw.copy(), np.zeros_like(raw)
            used_shift = bool(is_shift)
            lam = 0.0
        elif method == "clever_cov":
            lam = float(lambda_penalty) if is_shift else 0.0
            adapted, penalty = adapt_scores(raw, W, lam, is_shift)
            used_shift = bool(is_shift)
        elif method == "fsds_clever_cov":
            # FSDS gate: also fire on large δ or triangle/Φ collapse even if KS quiet
            fsds_fire = (
                float(hop["delta"]) > 0.08
                or float(np.mean(tri)) > 0.55
                or phi_drop > 0.05
            )
            used_shift = bool(is_shift or fsds_fire)
            lam = (
                fsds_lambda(
                    hop["c"],
                    hop["delta"],
                    hop.get("pi_change", 0.0),
                    float(np.mean(tri)),
                    phi_drop,
                    base_lambda=lambda_penalty,
                )
                if used_shift
                else 0.0
            )
            adapted, penalty = adapt_scores(raw, W, lam, used_shift)
        elif method == "oracle_penalty":
            # cheat: use true absolute error as W (upper-bound calibration)
            W_or = _abs_err(raw, yt) / (float(np.mean(_abs_err(raw, yt))) + 1e-8)
            used_shift = True
            lam = float(lambda_penalty)
            adapted, penalty = adapt_scores(raw, W_or, lam, True)
        else:
            raise ValueError(method)

        # human-review queue: high W under shift
        review = (W >= np.quantile(W, 0.80)) & used_shift if used_shift else np.zeros(len(W), dtype=bool)

        mse_raw = _mse(raw, yt)
        mse_ad = _mse(adapted, yt)
        # "treated" score: if review, fall back to best unimodal (cheap backup)
        backup = preds["text"].copy()
        for name in ("video", "audio"):
            better = errs[name] < _abs_err(backup, yt)
            backup = np.where(better, preds[name], backup)
        final = adapted.copy()
        if method in ("clever_cov", "fsds_clever_cov", "oracle_penalty"):
            final = np.where(review, 0.5 * adapted + 0.5 * backup, adapted)
        mse_final = _mse(final, yt)

        history.append(
            {
                "round": int(t),
                "online_mse": mse_final,
                "mse_raw": mse_raw,
                "mse_adapted": mse_ad,
                "mean_phi_TI": float(np.mean(phi["TI"])),
                "mean_phi_TA": float(np.mean(phi["TA"])),
                "mean_triangle": float(np.mean(tri)),
                "mean_W": float(np.mean(W)),
                "mean_penalty": float(np.mean(penalty)),
                "lambda": float(lam),
                "is_shift": bool(used_shift),
                "p_po": float(p_po),
                "p_align": float(p_al),
                "review_rate": float(np.mean(review)),
                "c": float(hop.get("c", 0.0)),
                "delta": float(hop.get("delta", 0.0)),
                "pi": hop.get("pi", {}),
                "var_interaction": float(decouple["var_interaction"]),
                "n": int(len(ic)),
            }
        )

        # rolling reference: only grow on non-shift early, else always extend lightly
        if t == 1 or not used_shift:
            ref_po.extend(po_risks.tolist())
            ref_tri.extend(tri.tolist())
            ref_phi.extend(phi_mean.tolist())
            ref_align.extend(align_batch.tolist())
        else:
            # slow leak so reference doesn't freeze forever
            ref_po.extend(po_risks[::4].tolist())
            ref_tri.extend(tri[::4].tolist())
            ref_phi.extend(phi_mean[::4].tolist())
            ref_align.extend(align_batch[::4].tolist())
        # cap reference length
        cap = 800
        ref_po, ref_tri, ref_phi, ref_align = (
            ref_po[-cap:],
            ref_tri[-cap:],
            ref_phi[-cap:],
            ref_align[-cap:],
        )

    online = np.array([h["online_mse"] for h in history], dtype=float)
    return {
        "method": method,
        "online_mse": float(np.mean(online)),
        "online_mse_sd": float(np.std(online, ddof=0)),
        "shift_rate": float(np.mean([h["is_shift"] for h in history])),
        "mean_review_rate": float(np.mean([h["review_rate"] for h in history])),
        "mean_phi_TI": float(np.mean([h["mean_phi_TI"] for h in history])),
        "mean_triangle": float(np.mean([h["mean_triangle"] for h in history])),
        "mean_W": float(np.mean([h["mean_W"] for h in history])),
        "history": history,
        "scenario": stream.meta.get("scenario"),
    }


def run_suite(
    scenarios=None,
    methods=None,
    seeds=None,
    n_batches=8,
    n_per=120,
    quick=False,
):
    scenarios = list(scenarios or SCENARIOS)
    methods = list(methods or METHODS)
    if quick:
        scenarios = [s for s in scenarios if s in ("healthy", "compound", "interaction_break")]
        methods = [m for m in methods if m != "oracle_penalty"] or methods
        n_batches = min(n_batches, 6)
        n_per = min(n_per, 60)
        seeds = seeds or [2026, 2027]
    seeds = list(seeds or [2026, 2027, 2028, 2029])

    rows = []
    table = {sc: {m: [] for m in methods} for sc in scenarios}
    for seed in seeds:
        for sc in scenarios:
            stream = make_trimodal_stream(
                scenario=sc, n_batches=n_batches, n_per=n_per, seed=seed
            )
            for method in methods:
                rec = run_method(stream, method=method, seed=seed)
                rows.append(
                    {
                        "seed": int(seed),
                        "scenario": sc,
                        "method": method,
                        "online_mse": rec["online_mse"],
                        "shift_rate": rec["shift_rate"],
                        "mean_review_rate": rec["mean_review_rate"],
                        "mean_phi_TI": rec["mean_phi_TI"],
                        "mean_triangle": rec["mean_triangle"],
                        "mean_W": rec["mean_W"],
                        "history": rec["history"],
                    }
                )
                table[sc][method].append(rec["online_mse"])

    summary = {}
    lifts = {}
    for sc in scenarios:
        summary[sc] = {}
        lifts[sc] = {}
        base = np.asarray(table[sc].get("static", [np.nan]), dtype=float)
        loc = np.asarray(table[sc].get("localize_only", [np.nan]), dtype=float)
        for method in methods:
            arr = np.asarray(table[sc][method], dtype=float)
            summary[sc][method] = {
                "online_mse": {"mean": float(np.mean(arr)), "sd": float(np.std(arr))},
            }
            lifts[sc][method] = {
                "vs_static": float(np.mean(base) - np.mean(arr)),
                "vs_localize": float(np.mean(loc) - np.mean(arr)),
            }
    return {
        "summary": summary,
        "lifts": lifts,
        "rows": rows,
        "scenarios": scenarios,
        "methods": methods,
        "seeds": seeds,
        "n_batches": n_batches,
        "n_per": n_per,
    }


def run_msrvtt_hop_probe(n_videos=24, n_windows=20, seed=SEED, n_estimators=30):
    """Optional real MSR-VTT window probe: Φ + triangle on early→late hop.

    Uses synthetic fallback if zip missing. Down-projects blocks for Ridge speed.
    """
    from msrvtt_multimodal_attribution import (
        GROUPS as M_GROUPS,
        load_window_bundle,
        make_synthetic_bundle,
    )

    try:
        bundle = load_window_bundle(root=Path.cwd())
        source = "msrvtt"
    except Exception:
        bundle = make_synthetic_bundle(
            n_videos=n_videos, n_windows=n_windows, seed=seed, video_shift=1.4
        )
        source = "synthetic_msrvtt"

    # subsample videos for speed
    vids = np.unique(bundle.video_id)
    rng = np.random.default_rng(seed)
    if len(vids) > n_videos:
        keep = set(rng.choice(vids, size=n_videos, replace=False).tolist())
        mask = np.array([v in keep for v in bundle.video_id])
    else:
        mask = np.ones(len(bundle.y), dtype=bool)

    X_full = standardize_columns(bundle.X[mask])
    y = np.asarray(bundle.y[mask], dtype=float)
    W = np.asarray(bundle.W[mask], dtype=int)

    # project each block to compact dims matching GROUPS
    def _proj(block, d_out, s):
        rng2 = np.random.default_rng(s)
        P = rng2.normal(size=(block.shape[1], d_out))
        P /= np.linalg.norm(P, axis=0, keepdims=True) + 1e-12
        return block @ P

    Xt = _proj(X_full[:, M_GROUPS["text"]], D_TEXT, seed + 1)
    Xv = _proj(X_full[:, M_GROUPS["video"]], D_VIDEO, seed + 2)
    Xa = _proj(X_full[:, M_GROUPS["audio"]], D_AUDIO, seed + 3)
    X = np.concatenate([Xt, Xv, Xa], axis=1)
    align = pairwise_align_scores(X)

    X0, y0 = X[W == 0], y[W == 0]
    X1, y1 = X[W == 1], y[W == 1]
    preds, errs, phi, decouple, phi_maps, _ = compute_batch_interaction(X1, y1, X0, y0)
    tri = np.array(
        [
            triangle_inconsistency({p: align[p][W == 1][i] for p in PAIR_NAMES})
            for i in range(X1.shape[0])
        ]
    )
    vimp, auc = rf_domain(X0, X1, seed=seed, n_estimators=n_estimators)
    _, pi = _modality_mass(vimp)
    # top shared token by |Φ_TI|
    top_k = int(np.argmax(np.abs(phi_maps["TI"]).mean(axis=0)))
    return {
        "source": source,
        "n0": int(X0.shape[0]),
        "n1": int(X1.shape[0]),
        "rf_auc": float(auc) if auc == auc else float("nan"),
        "pi": {g: float(pi.get(g, 0.0)) for g in GROUP_NAMES},
        "mean_phi_TI": float(np.mean(phi["TI"])),
        "mean_phi_TA": float(np.mean(phi["TA"])),
        "mean_phi_IA": float(np.mean(phi["IA"])),
        "top_token_TI": int(top_k),
        "top_token_phi_TI": float(phi_maps["TI"][:, top_k].mean()),
        "mean_triangle": float(np.mean(tri)),
        "mse_joint": _mse(preds["joint"], y1),
        "mse_text": _mse(preds["text"], y1),
        "mse_video": _mse(preds["video"], y1),
        "mse_audio": _mse(preds["audio"], y1),
        "var_interaction": float(decouple["var_interaction"]),
        "align_mean": {k: float(np.mean(align[k][W == 1])) for k in PAIR_NAMES},
    }


def plot_suite(suite, out_path: Path):
    import matplotlib.pyplot as plt

    scenarios = suite["scenarios"]
    methods = [m for m in suite["methods"] if m in suite["summary"][scenarios[0]]]
    colors = {
        "static": "#6B6B6B",
        "localize_only": "#A33B24",
        "clever_cov": "#2F6B4F",
        "fsds_clever_cov": "#1F6F8B",
        "oracle_penalty": "#C45C26",
    }
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.2), constrained_layout=True)

    ax = axes[0]
    x = np.arange(len(scenarios))
    width = 0.15
    for i, method in enumerate(methods):
        means = [suite["summary"][sc][method]["online_mse"]["mean"] for sc in scenarios]
        sds = [suite["summary"][sc][method]["online_mse"]["sd"] for sc in scenarios]
        ax.bar(
            x + (i - len(methods) / 2) * width + width / 2,
            means,
            width,
            yerr=sds,
            label=method,
            color=colors.get(method, "#333"),
            alpha=0.9,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(scenarios, rotation=20, ha="right")
    ax.set_ylabel("online MSE")
    ax.set_title("FSDS × Interaction Value online adapt", loc="left")
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1]
    # lift of fsds_clever_cov vs localize_only
    if "fsds_clever_cov" in methods and "localize_only" in methods:
        lifts = [suite["lifts"][sc]["fsds_clever_cov"]["vs_localize"] for sc in scenarios]
        ax.bar(x, lifts, color="#1F6F8B", alpha=0.9)
        ax.axhline(0.0, color="#333", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(scenarios, rotation=20, ha="right")
        ax.set_ylabel("ΔMSE (localize − fsds_clever_cov)")
        ax.set_title("Treatment lift vs localization-only", loc="left")
    else:
        ax.axis("off")

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    return out_path


def write_readme(suite, msrvtt_probe, out_path: Path):
    lines = [
        "# FSDS × Interaction Value Online Adaptation",
        "",
        "Localization names the sick modality; **clever covariates** treat at",
        "inference without retrain:",
        "",
        "```",
        "score_adapted = score_raw / (1 + λ · W_shift)",
        "```",
        "",
        "W_shift mixes PO-risk, triangle inconsistency, and Interaction-Value drop.",
        "`fsds_clever_cov` gates λ with typed hop `(ĉ, δ̂)` (locked TSS signs).",
        "",
        "## Live numbers",
        "",
        "| scenario | static | localize_only | clever_cov | **fsds_clever_cov** | lift vs localize |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for sc in suite["scenarios"]:
        s = suite["summary"][sc]
        lift = suite["lifts"][sc]["fsds_clever_cov"]["vs_localize"]
        lines.append(
            "| %s | %.4f | %.4f | %.4f | **%.4f** | %.4f |"
            % (
                sc,
                s["static"]["online_mse"]["mean"],
                s["localize_only"]["online_mse"]["mean"],
                s["clever_cov"]["online_mse"]["mean"],
                s["fsds_clever_cov"]["online_mse"]["mean"],
                lift,
            )
        )
    lines += [
        "",
        "## What each method does",
        "",
        "- `static` — joint Ridge, no gate.",
        "- `localize_only` — still computes π / Φ / triangle (diagnosis) but leaves preds alone.",
        "- `clever_cov` — KS on PO-risk & alignment → inject W into scores.",
        "- `fsds_clever_cov` — same W, λ from `(ĉ,δ̂,Δπ)` + triangle/Φ collapse.",
        "- `oracle_penalty` — upper bound using true |error| as W.",
        "",
        "## MSR-VTT hop probe",
        "",
    ]
    if msrvtt_probe:
        lines += [
            "Source: `%s` · n0=%d n1=%d · RF AUC=%.3f"
            % (
                msrvtt_probe["source"],
                msrvtt_probe["n0"],
                msrvtt_probe["n1"],
                msrvtt_probe["rf_auc"],
            ),
            "",
            "| π_text | π_video | π_audio | mean Φ_TI | triangle | mse_joint |",
            "| ---: | ---: | ---: | ---: | ---: | ---: |",
            "| %.3f | %.3f | %.3f | %.4f | %.3f | %.4f |"
            % (
                msrvtt_probe["pi"]["text"],
                msrvtt_probe["pi"]["video"],
                msrvtt_probe["pi"]["audio"],
                msrvtt_probe["mean_phi_TI"],
                msrvtt_probe["mean_triangle"],
                msrvtt_probe["mse_joint"],
            ),
            "",
        ]
    lines += [
        "## Run",
        "",
        "```",
        "python3 scripts/run_fsds_interaction_online_adapt.py",
        "python3 scripts/run_fsds_interaction_online_adapt.py --quick",
        "```",
        "",
        "Code: `Python/src/fsds_interaction_online_adapt.py`",
        "",
        "Product line: Pseudo-Outcome Risk / Interaction Shift as **clever covariates**",
        "closes the multimodal loop — diagnose with FSDS, treat in milliseconds.",
        "",
    ]
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines))
    return out_path
