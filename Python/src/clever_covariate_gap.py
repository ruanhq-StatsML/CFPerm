"""Modality-specific gap decomposition, then block-aware Stage-2 training.

Stage 0 — the analyst partitions coordinates into blocks B_1,...,B_M
  (modalities / feature groups). This is the prior; it is not learned.

Stage 1 — decompose a batch shift W into per-block contributions π:
  * block RF-domain AUC / VIMP mass / block MMD
  * leave-one-modality-out (LOMO) drop in domain AUC
  * instance-level gap  ê(X) − ê(X_{-m})
  Consensus π is a simplex weight over blocks.

Stage 2 — freeze π and use it to *guide the next training step*, not only
as extra columns. Trees are monotone-scale invariant, so “multiply column j
by √π_m” does not change RF splits; the weights have to change the *procedure*:

  * π-opinion pool:  ê_π = Σ_m π_m ê_m(X_m)
      linear opinion pool of block-wise OOF RFs (Genest–Zidek).
  * block-aware forest (BAWF): each tree samples feature j with
      p_j ∝ π_{m(j)} / |B_m|   (biased random subspace / informed bagging).
  * adaptive group-logit: after standardizing, scale column j by
      √(π_m / |B_m|) and fit L2 logistic (group-ridge / adaptive-lasso analogue).

Feature view of the same weights (previous prototype):
  * Z_m = π_m · logit ê_m(X_m)   then RF on Z or on [X, Z]
  * H_m = π_m · (W − ê_m) / (ê_m (1−ê_m))   TMLE clever covariate for PO-risk

H_m is for outcome targeting, not for re-predicting W.
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


def smoothed_pi(pi: np.ndarray, floor: float = 0.05) -> np.ndarray:
    """Dirichlet floor so a block with π=0 is not dropped from Stage-2."""
    pi = relu_normalize(np.asarray(pi, dtype=float).reshape(-1))
    m = len(pi)
    if m == 0:
        return pi
    floor = float(np.clip(floor, 0.0, 1.0))
    out = (1.0 - floor) * pi + floor / m
    s = float(out.sum())
    return out / s if s > EPS else np.full(m, 1.0 / m)


def block_column_probs(
    spec: ModalitySpec,
    pi: np.ndarray,
    *,
    n_features: int | None = None,
    floor: float = 0.05,
) -> np.ndarray:
    """Per-coordinate sampling / ridge weights: p_j ∝ π̃_{m(j)} / |B_{m(j)}|.

    Splitting the block mass equally inside the block keeps a wide text
    block from dominating a 5-d VAD block with the same π_m.
    """
    pi = smoothed_pi(pi, floor=floor)
    pdim = int(n_features if n_features is not None else max(sl.stop for sl in spec.slices))
    p = np.zeros(pdim, dtype=float)
    for m, sl in enumerate(spec.slices):
        width = int(sl.stop - sl.start)
        if width <= 0:
            continue
        p[sl] = pi[m] / width
    s = float(p.sum())
    if s <= EPS:
        p[:] = 1.0 / max(pdim, 1)
        return p
    p /= s
    return p


def pi_column_replicates(
    spec: ModalitySpec,
    pi: np.ndarray,
    n_features: int,
    *,
    floor: float = 0.05,
) -> np.ndarray:
    """Repeat column j  ≈  d · p_j  times so uniform max_features ≈ sampling p_j.

    This is the RF-fair way to feed π into sklearn: same trees, same per-split
    uniform draw, but the *multiset* of columns is π-weighted. Random subspace
    (BAWF) is a different estimator class and should not be compared to RF AUC.
    """
    p = block_column_probs(spec, pi, n_features=n_features, floor=floor)
    reps = np.maximum(1, np.rint(p * n_features).astype(int))
    return np.repeat(np.arange(int(n_features)), reps)


def expand_by_pi(
    X: np.ndarray,
    spec: ModalitySpec,
    pi: np.ndarray,
    *,
    floor: float = 0.05,
) -> Tuple[np.ndarray, np.ndarray]:
    X = np.asarray(X, dtype=float)
    idx = pi_column_replicates(spec, pi, X.shape[1], floor=floor)
    return X[:, idx], idx


def collapse_expanded_vimp(vimp: np.ndarray, idx: np.ndarray, spec: ModalitySpec) -> Dict[str, float]:
    d = int(np.max(idx)) + 1
    raw = np.zeros(d, dtype=float)
    vimp = np.asarray(vimp, dtype=float)
    for imp, j in zip(vimp, idx):
        raw[int(j)] += float(imp)
    return vimp_mass_share(raw, spec)


def cv_pi_weighted_rf_auc(
    X: np.ndarray,
    W: np.ndarray,
    spec: ModalitySpec,
    pi: np.ndarray,
    *,
    seed: int = 0,
    n_estimators: int = 100,
    n_splits: int = 5,
    floor: float = 0.05,
) -> Tuple[float, float, Dict[str, float]]:
    """K-fold AUC of sklearn RF on π-replicated columns (fair RF analogue)."""
    Xw, idx = expand_by_pi(X, spec, pi, floor=floor)
    auc, sd, vimp = cv_domain_auc(Xw, W, seed=seed, n_estimators=n_estimators, n_splits=n_splits)
    return auc, sd, collapse_expanded_vimp(vimp, idx, spec)


def opinion_pool_scores(gap: GapDecomposition, *, log_pool: bool = False) -> np.ndarray:
    """Linear (or log) opinion pool of block-wise OOF propensities, weights = π."""
    pi = np.asarray(gap.pi_consensus, dtype=float).reshape(1, -1)
    if log_pool:
        return (gap.z_logit_e * pi).sum(axis=1)
    return (gap.e_block * pi).sum(axis=1)


def _tree_pos_proba(tree, Xc: np.ndarray) -> np.ndarray:
    proba = tree.predict_proba(Xc)
    classes = list(tree.classes_)
    if 1 not in classes:
        return np.zeros(Xc.shape[0], dtype=float)
    return proba[:, classes.index(1)].astype(float)


def _fit_subspace_tree(
    X: np.ndarray,
    W: np.ndarray,
    cols: np.ndarray,
    boot: np.ndarray,
    max_depth: int,
    min_samples_leaf: int,
    rs: int,
):
    from sklearn.tree import DecisionTreeClassifier

    tree = DecisionTreeClassifier(
        max_depth=max_depth,
        min_samples_leaf=min_samples_leaf,
        max_features=None,
        random_state=rs,
    )
    tree.fit(X[boot][:, cols], W[boot])
    return tree, np.asarray(cols, dtype=int)


def fit_block_aware_forest(
    X: np.ndarray,
    W: np.ndarray,
    spec: ModalitySpec,
    pi: np.ndarray,
    *,
    n_estimators: int = 100,
    max_depth: int = 8,
    min_samples_leaf: int = 5,
    seed: int = 0,
    floor: float = 0.05,
    n_jobs: int = -1,
) -> dict:
    """Bagged trees with π-biased random-subspace feature sampling.

    sklearn RF uses *uniform* max_features at each split, which is invariant
    to column scaling. Here each tree draws k=√d columns with
    p_j ∝ π_m / |B_m|, then splits on that subset (Ho 1998, weighted).
    """
    from joblib import Parallel, delayed

    X = np.asarray(X, dtype=float)
    W = _as_1d(W).astype(int)
    n, d = X.shape
    k = max(1, min(d, int(np.sqrt(d))))
    probs = block_column_probs(spec, pi, n_features=d, floor=floor)
    rng = np.random.default_rng(seed)
    jobs = []
    for t in range(int(n_estimators)):
        boot = rng.integers(0, n, n)
        cols = rng.choice(d, size=k, replace=False, p=probs)
        jobs.append((boot, cols, int(seed + t)))
    fitted = Parallel(n_jobs=n_jobs, prefer="threads")(
        delayed(_fit_subspace_tree)(X, W, cols, boot, max_depth, min_samples_leaf, rs)
        for boot, cols, rs in jobs
    )
    return {
        "trees": [t for t, _ in fitted],
        "cols": [c for _, c in fitted],
        "probs": probs,
    }


def predict_block_aware_forest(model: dict, X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    trees, cols = model["trees"], model["cols"]
    if not trees:
        return np.full(X.shape[0], 0.5)
    acc = np.zeros(X.shape[0], dtype=float)
    for tree, c in zip(trees, cols):
        acc += _tree_pos_proba(tree, X[:, c])
    return acc / len(trees)


def block_aware_importance(model: dict, spec: ModalitySpec) -> Dict[str, float]:
    """Mean tree impurity mass mapped back to original columns, then blocks.

    Not Horvitz–Thompson-corrected: this is 'what the guided forest used'.
    """
    d = int(len(model["probs"]))
    imp = np.zeros(d, dtype=float)
    n_t = len(model["trees"])
    for tree, cols in zip(model["trees"], model["cols"]):
        fi = np.asarray(tree.feature_importances_, dtype=float)
        imp[cols] += fi
    if n_t:
        imp /= n_t
    return vimp_mass_share(imp, spec)


def cv_block_aware_auc(
    X: np.ndarray,
    W: np.ndarray,
    spec: ModalitySpec,
    pi: np.ndarray,
    *,
    seed: int = 0,
    n_estimators: int = 100,
    n_splits: int = 5,
    floor: float = 0.05,
) -> Tuple[float, float, Dict[str, float]]:
    """K-fold AUC of the π-guided forest; importance from a full-sample fit."""
    X = np.asarray(X, dtype=float)
    W = _as_1d(W).astype(int)
    n_splits = int(min(n_splits, int((W == 0).sum()), int((W == 1).sum())))
    n_splits = max(n_splits, 2)
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    aucs = []
    masses = []
    for k, (tr, te) in enumerate(cv.split(X, W)):
        model = fit_block_aware_forest(
            X[tr], W[tr], spec, pi, n_estimators=n_estimators,
            seed=seed + k, floor=floor,
        )
        pred = predict_block_aware_forest(model, X[te])
        aucs.append(float(roc_auc_score(W[te], pred)))
        masses.append(block_aware_importance(model, spec))
    mass = {name: float(np.mean([m[name] for m in masses])) for name in spec.names}
    return (
        float(np.mean(aucs)),
        float(np.std(aucs, ddof=1) if len(aucs) > 1 else 0.0),
        mass,
    )


def adaptive_group_scale(
    X: np.ndarray,
    spec: ModalitySpec,
    pi: np.ndarray,
    *,
    floor: float = 0.05,
) -> np.ndarray:
    """Column scales √(d · p_j) with p_j ∝ π_m / |B_m| (mean scale ≈ 1)."""
    d = int(np.asarray(X).shape[1])
    p = block_column_probs(spec, pi, n_features=d, floor=floor)
    return np.sqrt(np.maximum(p * d, EPS))


def cv_adaptive_group_logit_auc(
    X: np.ndarray,
    W: np.ndarray,
    spec: ModalitySpec,
    pi: np.ndarray,
    *,
    seed: int = 0,
    n_splits: int = 5,
    C: float = 1.0,
    floor: float = 0.05,
) -> Tuple[float, float]:
    """K-fold AUC of L2 logistic with π-adaptive group scaling.

    Scaling X_j by √π_m is a no-op for RF splits; it *is* the adaptive-ridge
    reparameterization for a linear model (Zou 2006, group analogue).
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    X = np.asarray(X, dtype=float)
    W = _as_1d(W).astype(int)
    scales = adaptive_group_scale(X, spec, pi, floor=floor)
    n_splits = int(min(n_splits, int((W == 0).sum()), int((W == 1).sum())))
    n_splits = max(n_splits, 2)
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    aucs = []
    for k, (tr, te) in enumerate(cv.split(X, W)):
        scaler = StandardScaler()
        Xtr = scaler.fit_transform(X[tr]) * scales.reshape(1, -1)
        Xte = scaler.transform(X[te]) * scales.reshape(1, -1)
        clf = LogisticRegression(
            C=C, solver="lbfgs", max_iter=500, random_state=seed + k,
        )
        clf.fit(Xtr, W[tr])
        aucs.append(float(roc_auc_score(W[te], clf.predict_proba(Xte)[:, 1])))
    return (
        float(np.mean(aucs)),
        float(np.std(aucs, ddof=1) if len(aucs) > 1 else 0.0),
    )


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
    with_block_train: bool = True,
) -> dict:
    W = _as_1d(W).astype(int)
    Z = clever_z(gap)
    XZ = clever_design(X, gap)
    pi = np.asarray(gap.pi_consensus, dtype=float)
    auc_raw, sd_raw, _ = cv_domain_auc(X, W, seed=seed, n_estimators=n_estimators, n_splits=n_splits)
    auc_z, sd_z, vimp_z = cv_domain_auc(Z, W, seed=seed + 1, n_estimators=n_estimators, n_splits=n_splits)
    auc_xz, sd_xz, _ = cv_domain_auc(XZ, W, seed=seed + 2, n_estimators=n_estimators, n_splits=n_splits)
    vimp_raw = fit_vimp(X, W, seed=seed, n_estimators=n_estimators)
    pool = opinion_pool_scores(gap)
    auc_pool = float(roc_auc_score(W, pool))
    auc_logpool = float(roc_auc_score(W, opinion_pool_scores(gap, log_pool=True)))
    out = {
        "domain_auc_raw": round(float(auc_raw), 4),
        "domain_auc_raw_sd": round(float(sd_raw), 4),
        "domain_auc_clever": round(float(auc_z), 4),
        "domain_auc_clever_sd": round(float(sd_z), 4),
        "domain_auc_stack_xz": round(float(auc_xz), 4),
        "domain_auc_stack_sd": round(float(sd_xz), 4),
        "domain_auc_pool": round(float(auc_pool), 4),
        "domain_auc_logpool": round(float(auc_logpool), 4),
        "domain_auc_delta": round(float(auc_z - auc_raw), 4),
        "domain_auc_delta_stack": round(float(auc_xz - auc_raw), 4),
        "domain_auc_delta_pool": round(float(auc_pool - auc_raw), 4),
        "vimp_share_raw": {k: round(v, 4) for k, v in vimp_mass_share(vimp_raw, spec).items()},
        "p_raw": int(X.shape[1]),
        "p_clever": int(Z.shape[1]),
    }
    if with_block_train:
        auc_bawf, sd_bawf, mass_bawf = cv_block_aware_auc(
            X, W, spec, pi, seed=seed + 3, n_estimators=n_estimators, n_splits=n_splits,
        )
        pi_unif = np.full(len(pi), 1.0 / max(len(pi), 1))
        auc_sub, sd_sub, _ = cv_block_aware_auc(
            X, W, spec, pi_unif, seed=seed + 3, n_estimators=n_estimators, n_splits=n_splits,
        )
        auc_adapt, sd_adapt = cv_adaptive_group_logit_auc(
            X, W, spec, pi, seed=seed + 4, n_splits=n_splits,
        )
        auc_logit_unif, sd_logit_unif = cv_adaptive_group_logit_auc(
            X, W, spec, pi_unif, seed=seed + 4, n_splits=n_splits,
        )
        out["domain_auc_bawf"] = round(float(auc_bawf), 4)
        out["domain_auc_bawf_sd"] = round(float(sd_bawf), 4)
        out["domain_auc_subspace"] = round(float(auc_sub), 4)
        out["domain_auc_subspace_sd"] = round(float(sd_sub), 4)
        out["domain_auc_adapt"] = round(float(auc_adapt), 4)
        out["domain_auc_adapt_sd"] = round(float(sd_adapt), 4)
        out["domain_auc_logit_unif"] = round(float(auc_logit_unif), 4)
        out["domain_auc_logit_unif_sd"] = round(float(sd_logit_unif), 4)
        out["domain_auc_delta_bawf"] = round(float(auc_bawf - auc_raw), 4)
        out["domain_auc_delta_bawf_vs_sub"] = round(float(auc_bawf - auc_sub), 4)
        out["domain_auc_delta_adapt"] = round(float(auc_adapt - auc_raw), 4)
        out["domain_auc_delta_adapt_vs_unif"] = round(float(auc_adapt - auc_logit_unif), 4)
        out["bawf_vimp"] = {k: round(v, 4) for k, v in mass_bawf.items()}
        auc_pirf, sd_pirf, mass_pirf = cv_pi_weighted_rf_auc(
            X, W, spec, pi, seed=seed + 5, n_estimators=n_estimators, n_splits=n_splits,
        )
        out["domain_auc_pi_rf"] = round(float(auc_pirf), 4)
        out["domain_auc_pi_rf_sd"] = round(float(sd_pirf), 4)
        out["domain_auc_delta_pi_rf"] = round(float(auc_pirf - auc_raw), 4)
        out["pi_rf_vimp"] = {k: round(v, 4) for k, v in mass_pirf.items()}
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
        if "bawf_vimp" in out:
            out["bawf_on_gt"] = round(float(out["bawf_vimp"][gt]), 4)
        if "pi_rf_vimp" in out:
            out["pi_rf_on_gt"] = round(float(out["pi_rf_vimp"][gt]), 4)
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
        raws, zs, xzs, pools, bawfs, adapts, subs, pirfs = [], [], [], [], [], [], [], []
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
            pools.append(row["domain_auc_pool"])
            bawfs.append(row.get("domain_auc_bawf", row["domain_auc_raw"]))
            adapts.append(row.get("domain_auc_adapt", row["domain_auc_raw"]))
            subs.append(row.get("domain_auc_subspace", row["domain_auc_raw"]))
            pirfs.append(row.get("domain_auc_pi_rf", row["domain_auc_raw"]))
        rows.append({
            "n": int(n0 + n1),
            "n_repeats": n_repeats,
            "auc_raw": round(float(np.mean(raws)), 4),
            "auc_raw_sd": round(float(np.std(raws, ddof=1)), 4),
            "auc_clever": round(float(np.mean(zs)), 4),
            "auc_clever_sd": round(float(np.std(zs, ddof=1)), 4),
            "auc_stack": round(float(np.mean(xzs)), 4),
            "auc_stack_sd": round(float(np.std(xzs, ddof=1)), 4),
            "auc_pool": round(float(np.mean(pools)), 4),
            "auc_pool_sd": round(float(np.std(pools, ddof=1)), 4),
            "auc_bawf": round(float(np.mean(bawfs)), 4),
            "auc_bawf_sd": round(float(np.std(bawfs, ddof=1)), 4),
            "auc_subspace": round(float(np.mean(subs)), 4),
            "auc_subspace_sd": round(float(np.std(subs, ddof=1)), 4),
            "auc_adapt": round(float(np.mean(adapts)), 4),
            "auc_adapt_sd": round(float(np.std(adapts, ddof=1)), 4),
            "delta": round(float(np.mean(zs) - np.mean(raws)), 4),
            "delta_stack": round(float(np.mean(xzs) - np.mean(raws)), 4),
            "delta_pool": round(float(np.mean(pools) - np.mean(raws)), 4),
            "delta_bawf": round(float(np.mean(bawfs) - np.mean(raws)), 4),
            "delta_bawf_vs_sub": round(float(np.mean(bawfs) - np.mean(subs)), 4),
            "auc_pi_rf": round(float(np.mean(pirfs)), 4),
            "auc_pi_rf_sd": round(float(np.std(pirfs, ddof=1)), 4),
            "delta_pi_rf": round(float(np.mean(pirfs) - np.mean(raws)), 4),
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
