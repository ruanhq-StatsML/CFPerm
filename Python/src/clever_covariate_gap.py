"""Modality-specific gap decomposition as clever covariates.

Stage 1 — decompose a batch shift W into per-modality contributions:
  * block RF-domain AUC / VIMP mass / block MMD
  * leave-one-modality-out (LOMO) drop in domain AUC
  * instance-level gap  ê(X) − ê(X_{-m})

Stage 2 — turn those contributions into features:
  * z_share_m(x)  = instance relative contribution (X-only, no W leak)
  * z_logit_e_m   = logit ê_m(X_m)  (stacking summary of modality m)
  * H_m           = π_m · (W − ê_m) / (ê_m (1−ê_m))   TMLE clever covariate

H_m is for outcome targeting / PO-risk, not for re-predicting W.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Sequence, Tuple

import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

EPS = 1e-12
CLIP_E = 0.01


def _as_1d(v) -> np.ndarray:
    return np.asarray(v, dtype=float).reshape(-1)


def relu_normalize(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float).reshape(-1)
    v = np.maximum(v, 0.0)
    s = float(v.sum())
    if s <= EPS:
        return np.full_like(v, 1.0 / max(len(v), 1))
    return v / s


def logit(p: np.ndarray, clip: float = CLIP_E) -> np.ndarray:
    p = np.clip(np.asarray(p, dtype=float), clip, 1.0 - clip)
    return np.log(p / (1.0 - p))


def _rf_clf(p: int, n: int, seed: int, n_estimators: int = 100):
    # Fixed capacity: sqrt(n) leaves underfit at large n and jitter small-n AUCs.
    return RandomForestClassifier(
        n_estimators=n_estimators,
        max_depth=8,
        min_samples_leaf=5,
        max_features="sqrt",
        n_jobs=-1,
        random_state=seed,
    )


def _rf_reg(seed: int, n_estimators: int = 80):
    return RandomForestRegressor(
        n_estimators=n_estimators,
        max_depth=8,
        min_samples_leaf=5,
        n_jobs=-1,
        random_state=seed,
    )


def _median_bandwidth(Z: np.ndarray) -> float:
    if Z.shape[0] < 2:
        return 1.0
    # subsample for the median heuristic
    rng = np.random.default_rng(0)
    if Z.shape[0] > 80:
        Z = Z[rng.choice(Z.shape[0], 80, replace=False)]
    diff = Z[:, None, :] - Z[None, :, :]
    dist = np.sqrt(np.sum(diff * diff, axis=2))
    med = np.median(dist[dist > 0])
    return float(med) if med > 0 else 1.0


def rbf_mmd2(X: np.ndarray, Y: np.ndarray, max_n: int = 120, seed: int = 0) -> float:
    rng = np.random.default_rng(seed)
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    if X.shape[0] > max_n:
        X = X[rng.choice(X.shape[0], max_n, replace=False)]
    if Y.shape[0] > max_n:
        Y = Y[rng.choice(Y.shape[0], max_n, replace=False)]
    if X.shape[0] < 2 or Y.shape[0] < 2:
        return 0.0
    pooled = np.vstack([X, Y])
    med = _median_bandwidth(pooled)
    gamma = 1.0 / (2.0 * med * med + EPS)

    def sq(A, B):
        A2 = np.sum(A * A, axis=1, keepdims=True)
        B2 = np.sum(B * B, axis=1, keepdims=True)
        return A2 + B2.T - 2.0 * (A @ B.T)

    n, m = X.shape[0], Y.shape[0]
    Kxx = np.exp(-gamma * sq(X, X))
    Kyy = np.exp(-gamma * sq(Y, Y))
    Kxy = np.exp(-gamma * sq(X, Y))
    np.fill_diagonal(Kxx, 0.0)
    np.fill_diagonal(Kyy, 0.0)
    return float(max(Kxx.sum() / (n * (n - 1)) + Kyy.sum() / (m * (m - 1)) - 2.0 * Kxy.mean(), 0.0))


def crossfit_propensity(
    X: np.ndarray,
    W: np.ndarray,
    *,
    n_splits: int = 3,
    seed: int = 0,
    n_estimators: int = 80,
    clip: float = CLIP_E,
) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    W = _as_1d(W).astype(int)
    n, p = X.shape
    e = np.zeros(n, dtype=float)
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for k, (tr, te) in enumerate(cv.split(X, W)):
        clf = _rf_clf(p, len(tr), seed + k, n_estimators=n_estimators)
        clf.fit(X[tr], W[tr])
        e[te] = clf.predict_proba(X[te])[:, 1]
    return np.clip(e, clip, 1.0 - clip)


def fit_vimp(X: np.ndarray, W: np.ndarray, *, seed: int = 0, n_estimators: int = 100) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    W = _as_1d(W).astype(int)
    clf = _rf_clf(X.shape[1], len(W), seed, n_estimators=n_estimators)
    clf.fit(X, W)
    return clf.feature_importances_.astype(float)


def cv_domain_auc(
    X: np.ndarray,
    W: np.ndarray,
    *,
    seed: int = 0,
    n_estimators: int = 100,
    n_splits: int = 5,
) -> Tuple[float, float, np.ndarray]:
    """K-fold domain AUC (mean, sd) plus mean OOF-fold VIMP."""
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    W = _as_1d(W).astype(int)
    n, p = X.shape
    n_splits = int(min(n_splits, int((W == 0).sum()), int((W == 1).sum())))
    n_splits = max(n_splits, 2)
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    aucs, vimps = [], []
    for k, (tr, te) in enumerate(cv.split(X, W)):
        clf = _rf_clf(p, len(tr), seed + k, n_estimators=n_estimators)
        clf.fit(X[tr], W[tr])
        aucs.append(float(roc_auc_score(W[te], clf.predict_proba(X[te])[:, 1])))
        vimps.append(clf.feature_importances_.astype(float))
    return float(np.mean(aucs)), float(np.std(aucs, ddof=1) if len(aucs) > 1 else 0.0), np.mean(np.stack(vimps), axis=0)


def holdout_domain_auc(
    X: np.ndarray,
    W: np.ndarray,
    *,
    seed: int = 0,
    n_estimators: int = 100,
    test_size: float = 0.25,
) -> Tuple[float, np.ndarray]:
    """Back-compat wrapper: CV AUC mean + full-data VIMP."""
    auc, _sd, _ = cv_domain_auc(X, W, seed=seed, n_estimators=n_estimators)
    vimp = fit_vimp(X, W, seed=seed, n_estimators=n_estimators)
    return auc, vimp


def po_risk_fit(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    *,
    seed: int = 0,
    n_splits: int = 3,
    clever_H: np.ndarray | None = None,
    n_estimators: int = 60,
) -> Tuple[np.ndarray, float, np.ndarray]:
    """DR / TMLE-style pseudo-outcome, then RF CATE VIMP.

    Standard PO: (Y − m̂)(W − ê)
    Clever PO:   (Y − m̂) · H   with H = Σ_m π_m (W − ê_m)/(ê_m(1−ê_m))
    """
    X = np.asarray(X, dtype=float)
    Y = _as_1d(Y)
    W = _as_1d(W).astype(int)
    n, p = X.shape
    m_hat = np.zeros(n)
    e_hat = np.zeros(n)
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for k, (tr, te) in enumerate(cv.split(X, W)):
        m = _rf_reg(seed + k, n_estimators=n_estimators)
        e = _rf_clf(p, len(tr), seed + 40 + k, n_estimators=n_estimators)
        m.fit(X[tr], Y[tr])
        e.fit(X[tr], W[tr])
        m_hat[te] = m.predict(X[te])
        e_hat[te] = e.predict_proba(X[te])[:, 1]
    e_hat = np.clip(e_hat, CLIP_E, 1.0 - CLIP_E)
    residual = Y - m_hat
    if clever_H is None:
        po = residual * (W - e_hat)
    else:
        H = np.asarray(clever_H, dtype=float)
        if H.ndim == 2:
            H = H.sum(axis=1)
        po = residual * _as_1d(H)
    tau = _rf_reg(seed + 7, n_estimators=max(80, n_estimators))
    tau.fit(X, po)
    pred = tau.predict(X)
    return tau.feature_importances_.astype(float), float(np.mean(pred ** 2)), po


@dataclass
class ModalitySpec:
    names: List[str]
    slices: List[slice]

    def split(self, X: np.ndarray) -> Dict[str, np.ndarray]:
        return {name: np.asarray(X[:, sl], dtype=float) for name, sl in zip(self.names, self.slices)}

    @property
    def n_mod(self) -> int:
        return len(self.names)


def vimp_mass_share(vimp: np.ndarray, spec: ModalitySpec) -> Dict[str, float]:
    vimp = np.asarray(vimp, dtype=float)
    mass = {name: float(np.sum(vimp[sl])) for name, sl in zip(spec.names, spec.slices)}
    tot = sum(mass.values()) + EPS
    return {k: v / tot for k, v in mass.items()}


def selection_auc_block(scores: np.ndarray, spec: ModalitySpec, gt: str) -> float:
    scores = np.asarray(scores, dtype=float)
    pos = np.zeros(len(scores), dtype=int)
    sl = spec.slices[spec.names.index(gt)]
    pos[sl] = 1
    if pos.min() == pos.max():
        return float("nan")
    return float(roc_auc_score(pos, scores))


@dataclass
class GapDecomposition:
    names: List[str]
    block_auc: Dict[str, float]
    block_mmd: Dict[str, float]
    lomo_auc_drop: Dict[str, float]
    vimp_share: Dict[str, float]
    pi_auc: np.ndarray
    pi_mmd: np.ndarray
    pi_lomo: np.ndarray
    pi_vimp: np.ndarray
    pi_instance_mean: np.ndarray
    pi_consensus: np.ndarray
    e_full: np.ndarray
    e_block: np.ndarray  # (n, M)
    instance_gap: np.ndarray  # (n, M)
    instance_share: np.ndarray  # (n, M)
    z_logit_e: np.ndarray  # (n, M)
    H_clever: np.ndarray  # (n, M)
    full_auc: float
    notes: Dict[str, str] = field(default_factory=dict)

    def as_dict(self) -> dict:
        names = self.names

        def pack(pi):
            return {n: round(float(v), 4) for n, v in zip(names, pi)}

        return {
            "full_auc": round(float(self.full_auc), 4),
            "block_auc": {k: round(float(v), 4) for k, v in self.block_auc.items()},
            "block_mmd": {k: round(float(v), 6) for k, v in self.block_mmd.items()},
            "lomo_auc_drop": {k: round(float(v), 4) for k, v in self.lomo_auc_drop.items()},
            "vimp_share": {k: round(float(v), 4) for k, v in self.vimp_share.items()},
            "pi_auc": pack(self.pi_auc),
            "pi_mmd": pack(self.pi_mmd),
            "pi_lomo": pack(self.pi_lomo),
            "pi_vimp": pack(self.pi_vimp),
            "pi_instance_mean": pack(self.pi_instance_mean),
            "pi_consensus": pack(self.pi_consensus),
        }


def decompose_modality_gap(
    X: np.ndarray,
    W: np.ndarray,
    spec: ModalitySpec,
    *,
    seed: int = 0,
    n_splits: int = 3,
    n_estimators: int = 100,
    mmd_max_n: int = 120,
    light: bool = False,
) -> GapDecomposition:
    """Estimate each modality's relative contribution to P(X) shift (W).

    `light=True` skips LOMO (used for multi-seed efficiency / GT repeats).
    Block AUCs come from OOF propensities, not a second holdout RF.
    """
    X = np.asarray(X, dtype=float)
    W = _as_1d(W).astype(int)
    blocks = spec.split(X)
    names = spec.names
    M = len(names)
    n = len(W)

    full_vimp = fit_vimp(X, W, seed=seed, n_estimators=n_estimators)
    vimp_share = vimp_mass_share(full_vimp, spec)

    e_full = crossfit_propensity(X, W, n_splits=n_splits, seed=seed, n_estimators=n_estimators)
    full_auc = float(roc_auc_score(W, e_full))
    e_block = np.zeros((n, M), dtype=float)
    block_auc: Dict[str, float] = {}
    block_mmd: Dict[str, float] = {}
    lomo_auc_drop: Dict[str, float] = {}
    instance_gap = np.zeros((n, M), dtype=float)

    X0, X1 = X[W == 0], X[W == 1]
    for j, name in enumerate(names):
        Xj = blocks[name]
        e_block[:, j] = crossfit_propensity(
            Xj, W, n_splits=n_splits, seed=seed + 10 * (j + 1), n_estimators=n_estimators
        )
        block_auc[name] = float(roc_auc_score(W, e_block[:, j]))
        block_mmd[name] = rbf_mmd2(X0[:, spec.slices[j]], X1[:, spec.slices[j]], max_n=mmd_max_n, seed=seed + j)

        if light:
            lomo_auc_drop[name] = 0.0
            continue
        keep = [i for i in range(X.shape[1]) if i not in range(spec.slices[j].start, spec.slices[j].stop)]
        X_lomo = X[:, keep] if keep else np.ones((n, 1))
        e_lomo = crossfit_propensity(
            X_lomo, W, n_splits=n_splits, seed=seed + 40 * (j + 1), n_estimators=n_estimators
        )
        auc_lomo = float(roc_auc_score(W, e_lomo))
        lomo_auc_drop[name] = max(float(full_auc - auc_lomo), 0.0)
        instance_gap[:, j] = e_full - e_lomo

    pi_auc = relu_normalize(np.array([max(block_auc[n] - 0.5, 0.0) for n in names]))
    pi_mmd = relu_normalize(np.array([block_mmd[n] for n in names]))
    pi_vimp = relu_normalize(np.array([vimp_share[n] for n in names]))
    abs_gap = np.abs(instance_gap)
    instance_share = abs_gap / (abs_gap.sum(axis=1, keepdims=True) + EPS)
    if light:
        pi_lomo = np.full(M, 1.0 / M)
        pi_instance_mean = np.full(M, 1.0 / M)
        pi_consensus = relu_normalize(pi_auc + pi_mmd + pi_vimp)
    else:
        pi_lomo = relu_normalize(np.array([lomo_auc_drop[n] for n in names]))
        pi_instance_mean = relu_normalize(abs_gap.mean(axis=0))
        pi_consensus = relu_normalize(pi_auc + pi_mmd + pi_lomo + pi_vimp + pi_instance_mean)

    z_logit_e = logit(e_block)
    H = np.zeros((n, M), dtype=float)
    for j in range(M):
        ej = e_block[:, j]
        H[:, j] = pi_consensus[j] * (W - ej) / (ej * (1.0 - ej) + EPS)

    return GapDecomposition(
        names=list(names),
        block_auc=block_auc,
        block_mmd=block_mmd,
        lomo_auc_drop=lomo_auc_drop,
        vimp_share=vimp_share,
        pi_auc=pi_auc,
        pi_mmd=pi_mmd,
        pi_lomo=pi_lomo,
        pi_vimp=pi_vimp,
        pi_instance_mean=pi_instance_mean,
        pi_consensus=pi_consensus,
        e_full=e_full,
        e_block=e_block,
        instance_gap=instance_gap,
        instance_share=instance_share,
        z_logit_e=z_logit_e,
        H_clever=H,
        full_auc=float(full_auc),
    )


def clever_z(gap: GapDecomposition) -> np.ndarray:
    """Contribution-weighted logit propensity: Z_m = π_m · logit ê_m(X_m)."""
    pi = np.asarray(gap.pi_consensus, dtype=float).reshape(1, -1)
    return gap.z_logit_e * pi


def clever_design(X: np.ndarray, gap: GapDecomposition) -> np.ndarray:
    """Raw X stacked with clever-Z."""
    return np.hstack([X, clever_z(gap)])


def compare_raw_vs_clever(
    X: np.ndarray,
    W: np.ndarray,
    Y: np.ndarray | None,
    spec: ModalitySpec,
    gap: GapDecomposition,
    *,
    seed: int = 0,
    gt: str | None = None,
    with_po: bool = False,
    n_estimators: int = 100,
    n_splits: int = 5,
) -> dict:
    W = _as_1d(W).astype(int)
    Z = clever_z(gap)
    XZ = clever_design(X, gap)
    auc_raw, sd_raw, _ = cv_domain_auc(X, W, seed=seed, n_estimators=n_estimators, n_splits=n_splits)
    auc_z, sd_z, vimp_z = cv_domain_auc(Z, W, seed=seed + 1, n_estimators=n_estimators, n_splits=n_splits)
    auc_xz, sd_xz, _ = cv_domain_auc(XZ, W, seed=seed + 2, n_estimators=n_estimators, n_splits=n_splits)
    vimp_raw = fit_vimp(X, W, seed=seed, n_estimators=n_estimators)
    out = {
        "domain_auc_raw": round(float(auc_raw), 4),
        "domain_auc_raw_sd": round(float(sd_raw), 4),
        "domain_auc_clever": round(float(auc_z), 4),
        "domain_auc_clever_sd": round(float(sd_z), 4),
        "domain_auc_stack_xz": round(float(auc_xz), 4),
        "domain_auc_stack_sd": round(float(sd_xz), 4),
        "domain_auc_delta": round(float(auc_z - auc_raw), 4),
        "domain_auc_delta_stack": round(float(auc_xz - auc_raw), 4),
        "vimp_share_raw": {k: round(v, 4) for k, v in vimp_mass_share(vimp_raw, spec).items()},
        "p_raw": int(X.shape[1]),
        "p_clever": int(Z.shape[1]),
    }
    tot_z = float(np.sum(np.maximum(vimp_z, 0.0))) + EPS
    out["z_vimp"] = {n: round(float(max(vimp_z[i], 0.0) / tot_z), 4) for i, n in enumerate(spec.names)}
    out["z_logit_vimp"] = out["z_vimp"]
    out["z_share_vimp"] = out["z_vimp"]

    if with_po and Y is not None:
        Y = _as_1d(Y)
        _, po_raw, _ = po_risk_fit(X, Y, W, seed=seed, n_estimators=n_estimators)
        _, po_z, _ = po_risk_fit(Z, Y, W, seed=seed + 2, n_estimators=n_estimators)
        _, po_h, _ = po_risk_fit(X, Y, W, seed=seed + 3, clever_H=gap.H_clever, n_estimators=n_estimators)
        out["po_risk_raw"] = round(float(po_raw), 6)
        out["po_risk_clever_Z"] = round(float(po_z), 6)
        out["po_risk_tmle_H"] = round(float(po_h), 6)

    if gt is not None and gt in spec.names:
        out["selection_auc_raw"] = round(selection_auc_block(vimp_raw, spec, gt), 4)
        out["mass_on_gt_raw"] = round(out["vimp_share_raw"][gt], 4)
        out["z_share_on_gt"] = out["z_vimp"][gt]
        out["z_logit_on_gt"] = out["z_vimp"][gt]
        out["pi_consensus_on_gt"] = round(float(gap.pi_consensus[spec.names.index(gt)]), 4)
    return out


def sample_efficiency_curve(
    X: np.ndarray,
    W: np.ndarray,
    spec: ModalitySpec,
    *,
    ns: Sequence[int] = (200, 400, 800, 1200),
    seed: int = 0,
    n_estimators: int = 80,
    n_repeats: int = 4,
) -> List[dict]:
    """Subsample n over several seeds; light gap (no LOMO) + 5-fold CV AUC."""
    X = np.asarray(X, dtype=float)
    W = _as_1d(W).astype(int)
    i0 = np.where(W == 0)[0]
    i1 = np.where(W == 1)[0]
    rows = []
    for n in ns:
        n0 = min(len(i0), n // 2)
        n1 = min(len(i1), n // 2)
        if n0 < 50 or n1 < 50:
            continue
        raws, zs, xzs = [], [], []
        for r in range(n_repeats):
            rng = np.random.default_rng(seed + 17 * n + r)
            sel = np.concatenate([
                rng.choice(i0, n0, replace=False),
                rng.choice(i1, n1, replace=False),
            ])
            Xs, Ws = X[sel], W[sel]
            gap = decompose_modality_gap(
                Xs, Ws, spec, seed=seed + n + r, n_splits=3,
                n_estimators=n_estimators, light=True,
            )
            row = compare_raw_vs_clever(
                Xs, Ws, None, spec, gap, seed=seed + n + r,
                with_po=False, n_estimators=n_estimators, n_splits=4,
            )
            raws.append(row["domain_auc_raw"])
            zs.append(row["domain_auc_clever"])
            xzs.append(row["domain_auc_stack_xz"])
        rows.append({
            "n": int(n0 + n1),
            "n_repeats": n_repeats,
            "auc_raw": round(float(np.mean(raws)), 4),
            "auc_raw_sd": round(float(np.std(raws, ddof=1)), 4),
            "auc_clever": round(float(np.mean(zs)), 4),
            "auc_clever_sd": round(float(np.std(zs, ddof=1)), 4),
            "auc_stack": round(float(np.mean(xzs)), 4),
            "auc_stack_sd": round(float(np.std(xzs, ddof=1)), 4),
            "delta": round(float(np.mean(zs) - np.mean(raws)), 4),
            "delta_stack": round(float(np.mean(xzs) - np.mean(raws)), 4),
            "p_raw": int(X.shape[1]),
            "p_clever": spec.n_mod,
        })
    return rows


def make_synthetic_shift(
    *,
    n: int = 500,
    d_text: int = 24,
    d_vad: int = 4,
    gt: str = "valence",
    mean_shift: float = 1.15,
    seed: int = 0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, ModalitySpec]:
    """Four-block DGP: text / valence / arousal / dominance. Shift lives in `gt`."""
    rng = np.random.default_rng(seed)
    W = rng.integers(0, 2, size=n)
    names = ["text", "valence", "arousal", "dominance"]
    dims = [d_text, d_vad, d_vad, d_vad]
    blocks = []
    for name, d in zip(names, dims):
        Xj = rng.normal(0.0, 1.0, size=(n, d))
        if name == gt:
            Xj[W == 1] += mean_shift
        blocks.append(Xj)
    X = np.hstack(blocks)
    # Concept-drift outcome driven by the GT block on the late batch
    sl_start = int(np.cumsum([0] + dims)[names.index(gt)])
    y_signal = X[:, sl_start]
    Y = 0.4 * y_signal + 1.1 * W * y_signal + rng.normal(0, 0.7, size=n)
    slices, start = [], 0
    for d in dims:
        slices.append(slice(start, start + d))
        start += d
    spec = ModalitySpec(names=names, slices=slices)
    return X, W, Y, spec


def inject_block_shift(
    X: np.ndarray,
    W: np.ndarray,
    spec: ModalitySpec,
    gt: str,
    *,
    alpha: float = 0.85,
) -> np.ndarray:
    """Add a mean shift to the late batch on one modality (controlled GT)."""
    Xc = np.array(X, dtype=float, copy=True)
    sl = spec.slices[spec.names.index(gt)]
    width = sl.stop - sl.start
    u = np.ones(width, dtype=float) / np.sqrt(width)
    Xc[W == 1, sl] = Xc[W == 1, sl] + alpha * u
    return Xc


def shuffle_block(
    X: np.ndarray,
    spec: ModalitySpec,
    name: str,
    *,
    seed: int = 0,
) -> np.ndarray:
    """Break association between one modality and W (null control)."""
    Xc = np.array(X, dtype=float, copy=True)
    sl = spec.slices[spec.names.index(name)]
    rng = np.random.default_rng(seed)
    Xc[:, sl] = Xc[rng.permutation(len(Xc)), sl]
    return Xc
