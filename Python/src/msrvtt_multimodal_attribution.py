"""MSR-VTT multimodal FSDS attribution with resampling inference.

Layout of ``s`` (last column = label):
    [768 video | 512 audio | 768 text | 1 label]

Batch assignment W is the within-video early vs late window split.
"""
from __future__ import annotations

import json
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy import stats
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict

D_VIDEO = 768
D_AUDIO = 512
D_TEXT = 768
P_X = D_VIDEO + D_AUDIO + D_TEXT
GROUPS = {
    "video": slice(0, D_VIDEO),
    "audio": slice(D_VIDEO, D_VIDEO + D_AUDIO),
    "text": slice(D_VIDEO + D_AUDIO, P_X),
}
GROUP_NAMES = ("video", "audio", "text")
SEED = 2026
ALPHA_CLIP = 0.01
ZIP_NAMES = (
    "feature_video_audio.zip",
    "data/msrvtt/feature_video_audio.zip",
    "data/feature_video_audio.zip",
    "datasets/feature_video_audio.zip",
    "msr_vtt/feature_video_audio.zip",
    "msr_vtt/features_window/feature_video_audio.zip",
)


def _as_1d(x):
    return np.asarray(x, dtype=float).reshape(-1)


def standardize_columns(X):
    X = np.asarray(X, dtype=float)
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd = np.where(sd < 1e-8, 1.0, sd)
    return (X - mu) / sd


def drop_group(X, name):
    mask = np.ones(X.shape[1], dtype=bool)
    mask[GROUPS[name]] = False
    return X[:, mask]


def modality_mass(vimp, positive=True):
    v = np.asarray(vimp, dtype=float)
    out = {}
    for name, sl in GROUPS.items():
        block = v[sl]
        out[name] = float(np.clip(block, 0, None).sum() if positive else block.sum())
    tot = sum(abs(x) for x in out.values()) + 1e-12
    share = {k: out[k] / tot for k in out}
    return out, share


def _rf_clf(p, n, seed, n_estimators=150):
    return RandomForestClassifier(
        n_estimators=n_estimators,
        max_depth=max(2, int(round(np.sqrt(p)))),
        min_samples_leaf=max(1, int(round(np.sqrt(max(n, 2)) // 2))),
        max_features="sqrt",
        n_jobs=-1,
        random_state=seed,
    )


def _rf_reg(p, n, seed, n_estimators=150):
    return RandomForestRegressor(
        n_estimators=n_estimators,
        max_depth=max(2, int(round(np.sqrt(p)))),
        min_samples_leaf=max(1, int(round(np.sqrt(max(n, 2)) // 2))),
        max_features=max(1, int(round(p / 3.0))),
        n_jobs=-1,
        random_state=seed,
    )


def oof_auc(X, W, seed=SEED, n_splits=5, n_estimators=150):
    W = np.asarray(W, dtype=int)
    if len(np.unique(W)) < 2:
        return float("nan")
    n_splits = max(2, min(n_splits, int(np.min(np.bincount(W)))))
    if n_splits < 2:
        return float("nan")
    clf = _rf_clf(X.shape[1], len(W), seed, n_estimators=n_estimators)
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    try:
        proba = cross_val_predict(clf, X, W, cv=cv, method="predict_proba")
        return float(roc_auc_score(W, proba[:, 1]))
    except Exception:
        return float("nan")


def rf_domain(X0, X1, seed=SEED, n_estimators=150):
    X = np.vstack([X0, X1])
    W = np.concatenate([np.zeros(len(X0), dtype=int), np.ones(len(X1), dtype=int)])
    clf = _rf_clf(X.shape[1], len(W), seed, n_estimators=n_estimators)
    clf.fit(X, W)
    vimp = clf.feature_importances_.astype(float)
    auc = oof_auc(X, W, seed=seed + 1, n_estimators=max(60, n_estimators // 2))
    return vimp, auc


def _import_mmd():
    import sys

    root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(root / "vendor" / "fsds"))
    from VIMP_mmd_benchmark import MMD, _rbf_mmd2_unbiased, _median_bandwidth, _subsample_batch

    return MMD, _rbf_mmd2_unbiased, _median_bandwidth, _subsample_batch


def _mmd_auto(A, B, median_bw, rbf, seed):
    pooled = np.vstack([A, B])
    med = median_bw(pooled, seed=seed)
    gamma = 1.0 / (2.0 * med * med + 1e-12)
    return rbf(A, B, gamma)


def group_mmd_loco(X0, X1, max_n=80, seed=SEED):
    """Group MMD: block-wise MMD (own bandwidth) + leave-one-group-out delta."""
    MMD, rbf, median_bw, subsample = _import_mmd()
    rng = np.random.default_rng(seed)
    A = subsample(X0, max_n, rng)
    B = subsample(X1, max_n, rng)
    full = _mmd_auto(A, B, median_bw, rbf, seed)
    loco = {}
    block = {}
    for name, sl in GROUPS.items():
        block[name] = float(_mmd_auto(A[:, sl], B[:, sl], median_bw, rbf, seed + 1))
        g = _mmd_auto(drop_group(A, name), drop_group(B, name), median_bw, rbf, seed + 2)
        loco[name] = float(full - g)
    return float(full), loco, block


def coord_mmd_vimp(X0, X1, max_n=80, seed=SEED):
    MMD, rbf, median_bw, subsample = _import_mmd()
    rng = np.random.default_rng(seed)
    A = subsample(X0, max_n, rng)
    B = subsample(X1, max_n, rng)
    p = A.shape[1]
    vimp = np.zeros(p)
    for j in range(p):
        Aj, Bj = A[:, [j]], B[:, [j]]
        med = np.median(np.abs(np.concatenate([Aj, Bj]) - np.median(np.concatenate([Aj, Bj]))))
        med = float(med) if med > 0 else 1.0
        gamma = 1.0 / (2.0 * med * med + 1e-12)
        vimp[j] = rbf(Aj, Bj, gamma)
    return vimp


def feature_mmd_loco(X0, X1, max_n=50, seed=SEED, max_features=None):
    """Leave-one-coordinate-out MMD VIMP (shared bandwidth)."""
    MMD, rbf, median_bw, subsample = _import_mmd()
    rng = np.random.default_rng(seed)
    A = subsample(X0, max_n, rng)
    B = subsample(X1, max_n, rng)
    p = A.shape[1]
    gamma = 1.0 / (2.0 * median_bw(np.vstack([A, B]), seed=seed) ** 2 + 1e-12)
    full = rbf(A, B, gamma)
    vimp = np.zeros(p)
    idx = np.arange(p)
    if max_features is not None and max_features < p:
        idx = rng.choice(p, max_features, replace=False)
        vimp[:] = np.nan
    for j in idx:
        mask = np.ones(p, dtype=bool)
        mask[j] = False
        vimp[j] = full - rbf(A[:, mask], B[:, mask], gamma)
    return vimp


def crossfit_po(X, Y, W, seed=SEED, n_splits=5, n_estimators=80):
    n = len(Y)
    Yf = _as_1d(Y)
    W = np.asarray(W, dtype=int)
    m_hat = np.zeros(n)
    e_hat = np.zeros(n)
    n_splits = max(2, min(n_splits, int(np.min(np.bincount(W)))))
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for fold, (tr, te) in enumerate(cv.split(X, W)):
        m = _rf_reg(X.shape[1], len(tr), seed + fold, n_estimators=n_estimators)
        e = _rf_clf(X.shape[1], len(tr), seed + 40 + fold, n_estimators=n_estimators)
        m.fit(X[tr], Yf[tr])
        e.fit(X[tr], W[tr])
        m_hat[te] = m.predict(X[te])
        e_hat[te] = e.predict_proba(X[te])[:, 1]
    e_hat = np.clip(e_hat, ALPHA_CLIP, 1.0 - ALPHA_CLIP)
    po = (Yf - m_hat) * (W - e_hat)
    return po, m_hat, e_hat


def po_tau_vimp(X, po, seed=SEED, n_estimators=150):
    tau = _rf_reg(X.shape[1], len(po), seed, n_estimators=n_estimators)
    tau.fit(X, po)
    hat = tau.predict(X)
    return float(np.mean(hat ** 2)), tau.feature_importances_.astype(float), tau


def group_permute_po_risk(X, po, tau, seed=SEED, n_repeat=5):
    """Leave-one-group-permuted PO-risk: shuffle a modality block, keep tau fixed."""
    rng = np.random.default_rng(seed)
    r_full = float(np.mean(tau.predict(X) ** 2))
    contrib = {}
    for name, sl in GROUPS.items():
        deltas = []
        for k in range(n_repeat):
            Xp = X.copy()
            Xp[:, sl] = rng.permutation(Xp[:, sl])
            deltas.append(float(np.mean(tau.predict(Xp) ** 2)) - r_full)
        contrib[name] = float(np.mean(deltas))
    return contrib


def logo_po_risk(X, po, seed=SEED, n_estimators=120):
    r_full, vimp, tau = po_tau_vimp(X, po, seed=seed, n_estimators=n_estimators)
    contrib = {}
    for name in GROUP_NAMES:
        Xg = drop_group(X, name)
        tau_g = _rf_reg(Xg.shape[1], len(po), seed + 11 + GROUP_NAMES.index(name), n_estimators=n_estimators)
        tau_g.fit(Xg, po)
        r_g = float(np.mean(tau_g.predict(Xg) ** 2))
        contrib[name] = r_g - r_full
    perm = group_permute_po_risk(X, po, tau, seed=seed + 21)
    return r_full, contrib, vimp, perm


def holm_adjust(pvals):
    pvals = np.asarray(pvals, dtype=float)
    m = len(pvals)
    order = np.argsort(pvals)
    adj = np.empty(m)
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, (m - rank) * pvals[idx])
        adj[idx] = min(1.0, running)
    return adj


def percentile_ci(samples, level=0.95):
    a = np.asarray(samples, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return (float("nan"), float("nan"))
    lo = float(np.quantile(a, (1.0 - level) / 2.0))
    hi = float(np.quantile(a, 1.0 - (1.0 - level) / 2.0))
    return lo, hi


def bootstrap_pvalue(diffs):
    """Two-sided percentile bootstrap p-value for H0: mean(diff)=0."""
    d = np.asarray(diffs, dtype=float)
    d = d[np.isfinite(d)]
    if d.size == 0:
        return float("nan")
    p = 2.0 * min(np.mean(d <= 0.0), np.mean(d >= 0.0))
    return float(min(1.0, max(p, 1.0 / (d.size + 1.0))))


def cluster_indices(video_ids, rng):
    videos = np.unique(video_ids)
    draw = rng.choice(videos, size=len(videos), replace=True)
    return np.concatenate([np.flatnonzero(video_ids == v) for v in draw])


def group_partition_stat(vimp):
    _, share = modality_mass(vimp)
    vals = np.array([share[g] for g in GROUP_NAMES])
    return float(vals.max() - vals.min())


def permute_group_labels(vimp, rng):
    """Random partition into 768/512/768; tests whether named blocks are special."""
    p = np.asarray(vimp, dtype=float).copy()
    rng.shuffle(p)
    return group_partition_stat(p)


@dataclass
class WindowBundle:
    X: np.ndarray
    y: np.ndarray
    W: np.ndarray
    video_id: np.ndarray
    window_idx: np.ndarray
    y_align: np.ndarray
    meta: dict = field(default_factory=dict)

    @property
    def s(self):
        return np.hstack([self.X, self.y.reshape(-1, 1)])


def _early_late(window_idx, video_id):
    W = np.zeros(len(window_idx), dtype=int)
    for v in np.unique(video_id):
        idx = np.flatnonzero(video_id == v)
        order = np.argsort(window_idx[idx])
        half = len(idx) // 2
        late = idx[order[half:]]
        W[late] = 1
    return W


def _align_score(X):
    v = X[:, GROUPS["video"]]
    t = X[:, GROUPS["text"]]
    vn = v / (np.linalg.norm(v, axis=1, keepdims=True) + 1e-12)
    tn = t / (np.linalg.norm(t, axis=1, keepdims=True) + 1e-12)
    return np.sum(vn * tn, axis=1)


def _load_npy_map(root: Path):
    out = {}
    for p in root.rglob("*"):
        if p.suffix.lower() == ".npy":
            out[p.name] = p
        elif p.suffix.lower() == ".npz":
            out[p.name] = p
    return out


def _maybe_unzip(zip_path: Path, dest: Path):
    dest.mkdir(parents=True, exist_ok=True)
    marker = dest / ".unzipped_from"
    if marker.exists() and marker.read_text().strip() == str(zip_path.resolve()):
        return dest
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(dest)
    marker.write_text(str(zip_path.resolve()))
    return dest


def find_feature_zip(root: Path | None = None) -> Path | None:
    root = Path(root) if root is not None else Path.cwd()
    for rel in ZIP_NAMES:
        p = (root / rel).resolve()
        if p.exists():
            return p
    hits = list(root.rglob("feature_video_audio.zip"))
    return hits[0] if hits else None


def load_window_bundle(path: Path | None = None, root: Path | None = None) -> WindowBundle:
    """Load ``s`` from feature_video_audio.zip or an extracted directory."""
    root = Path(root) if root is not None else Path.cwd()
    zip_path = None
    extract_dir = root / "data" / "msrvtt" / "features_window"
    if path is not None:
        path = Path(path)
        if path.is_file() and path.suffix == ".zip":
            zip_path = path
        elif path.is_file() and path.suffix == ".npy":
            s = np.load(path)
            return bundle_from_s(s)
        elif path.is_dir():
            extract_dir = path
    else:
        zip_path = find_feature_zip(root)

    files = {}
    if zip_path is not None:
        extract_dir = _maybe_unzip(zip_path, extract_dir)
        files = _load_npy_map(extract_dir)
        # also allow members named s.npy inside nested folders
    else:
        if extract_dir.exists():
            files = _load_npy_map(extract_dir)
        files.update(_load_npy_map(root))

    if "s.npy" in files:
        s = np.load(files["s.npy"])
        extra = {}
        if "video_labels.npy" in files:
            extra["video_id"] = np.load(files["video_labels.npy"]).reshape(-1)
        if "window_idx.npy" in files:
            extra["window_idx"] = np.load(files["window_idx.npy"]).reshape(-1)
        return bundle_from_s(s, **extra)

    # 3d window tensors
    if "video_feat_3d.npy" in files:
        vf = np.load(files["video_feat_3d.npy"])
        af = np.load(files["audio_feat_3d.npy"])
        tf = np.load(files["text_feat_3d.npy"])
        n_v, n_w = vf.shape[0], vf.shape[1]
        X = np.concatenate(
            [vf.reshape(n_v * n_w, -1), af.reshape(n_v * n_w, -1), tf.reshape(n_v * n_w, -1)],
            axis=1,
        )
        video_id = np.repeat(np.arange(n_v), n_w)
        window_idx = np.tile(np.arange(n_w), n_v)
        y = video_id.astype(float)
        if "video_labels.npy" in files:
            y = np.load(files["video_labels.npy"]).reshape(-1)[: X.shape[0]].astype(float)
        return _bundle(X, y, video_id, window_idx, {"source": "3d"})

    if "video_feat.npy" in files and "audio_feat.npy" in files and "text_feat.npy" in files:
        vf = np.load(files["video_feat.npy"])
        af = np.load(files["audio_feat.npy"])
        tf = np.load(files["text_feat.npy"])
        X = np.concatenate([vf, af, tf], axis=1)
        n = X.shape[0]
        if "video_labels.npy" in files:
            video_id = np.load(files["video_labels.npy"]).reshape(-1)[:n]
        else:
            video_id = np.zeros(n, dtype=int)
        if "window_idx.npy" in files:
            window_idx = np.load(files["window_idx.npy"]).reshape(-1)[:n]
        else:
            # reconstruct windows if labels are constant on blocks
            window_idx = np.zeros(n, dtype=int)
            for v in np.unique(video_id):
                idx = np.flatnonzero(video_id == v)
                window_idx[idx] = np.arange(len(idx))
        y = video_id.astype(float)
        return _bundle(X, y, video_id, window_idx, {"source": "flat_feats"})

    searched = zip_path or extract_dir
    raise FileNotFoundError(
        "Could not load MSR-VTT window features. Expected feature_video_audio.zip "
        "(or s.npy with 768+512+768+1 columns) under %s" % searched
    )


def _bundle(X, y, video_id, window_idx, meta):
    X = np.asarray(X, dtype=float)
    if X.shape[1] == P_X + 1:
        y = X[:, -1]
        X = X[:, :P_X]
    if X.shape[1] != P_X:
        raise ValueError("expected X with %d columns, got %s" % (P_X, X.shape))
    y = _as_1d(y)
    video_id = np.asarray(video_id).reshape(-1)
    window_idx = np.asarray(window_idx).reshape(-1)
    W = _early_late(window_idx, video_id)
    return WindowBundle(
        X=X,
        y=y,
        W=W,
        video_id=video_id,
        window_idx=window_idx,
        y_align=_align_score(X),
        meta=meta,
    )


def bundle_from_s(s, video_id=None, window_idx=None):
    s = np.asarray(s, dtype=float)
    meta = {"source": "s", "s_shape": list(s.shape)}
    if s.ndim == 3:
        n_v, n_w, d = s.shape
        s2 = s.reshape(n_v * n_w, d)
        vid = np.repeat(np.arange(n_v), n_w) if video_id is None else np.asarray(video_id).reshape(-1)
        widx = np.tile(np.arange(n_w), n_v) if window_idx is None else np.asarray(window_idx).reshape(-1)
        y = s2[:, -1] if d == P_X + 1 else vid.astype(float)
        X = s2[:, :P_X] if d >= P_X else s2
        return _bundle(X, y, vid, widx, meta)
    if s.ndim != 2:
        raise ValueError("s must be 2d or 3d, got %s" % (s.shape,))
    n, d = s.shape
    if d == P_X + 1:
        X, y = s[:, :P_X], s[:, -1]
    elif d == P_X:
        X, y = s, np.zeros(n)
    else:
        raise ValueError("s has %d columns; expected %d or %d" % (d, P_X, P_X + 1))
    if video_id is None:
        # treat integer label as video id when cardinality is small
        uniq = np.unique(y)
        if 1 < uniq.size <= n // 4 + 1:
            video_id = y.astype(int)
        else:
            video_id = np.zeros(n, dtype=int)
    video_id = np.asarray(video_id).reshape(-1)
    if window_idx is None:
        window_idx = np.zeros(n, dtype=int)
        for v in np.unique(video_id):
            idx = np.flatnonzero(video_id == v)
            window_idx[idx] = np.arange(len(idx))
    return _bundle(X, y, video_id, window_idx, meta)


def make_synthetic_bundle(
    n_videos=8,
    n_windows=40,
    seed=SEED,
    video_shift=1.25,
    audio_shift=0.15,
    text_shift=0.05,
):
    """Late windows carry a mean shift concentrated in the video block."""
    rng = np.random.default_rng(seed)
    rows = []
    y = []
    vid = []
    widx = []
    for v in range(n_videos):
        X = rng.normal(size=(n_windows, P_X))
        half = n_windows // 2
        X[half:, GROUPS["video"]] += video_shift
        X[half:, GROUPS["audio"]] += audio_shift
        X[half:, GROUPS["text"]] += text_shift
        # label = video id plus a small concept-drift in late windows
        lab = np.full(n_windows, float(v))
        lab[half:] += 0.35 * X[half:, 0]
        rows.append(X)
        y.append(lab)
        vid.append(np.full(n_windows, v))
        widx.append(np.arange(n_windows))
    X = np.vstack(rows)
    return _bundle(X, np.concatenate(y), np.concatenate(vid), np.concatenate(widx), {"synthetic": True})


def per_video_table(bundle: WindowBundle, seed=SEED, n_estimators=120):
    rows = []
    for v in np.unique(bundle.video_id):
        idx = np.flatnonzero(bundle.video_id == v)
        X = standardize_columns(bundle.X[idx])
        W = bundle.W[idx]
        if W.sum() < 4 or (len(W) - W.sum()) < 4:
            continue
        vimp, auc = rf_domain(X[W == 0], X[W == 1], seed=seed + int(v), n_estimators=n_estimators)
        mass, share = modality_mass(vimp)
        rows.append(
            {
                "video_id": int(v) if float(v).is_integer() else str(v),
                "n": int(len(idx)),
                "n0": int((W == 0).sum()),
                "n1": int((W == 1).sum()),
                "auc": auc,
                "share": share,
                "dominant": max(share, key=share.get),
            }
        )
    return rows


def pooled_point_estimates(bundle: WindowBundle, seed=SEED, n_estimators=120, mmd_n=80):
    X = standardize_columns(bundle.X)
    X0, X1 = X[bundle.W == 0], X[bundle.W == 1]
    rf_vimp, rf_auc = rf_domain(X0, X1, seed=seed, n_estimators=n_estimators)
    rf_mass, rf_share = modality_mass(rf_vimp)
    mmd_full, mmd_loco, mmd_block = group_mmd_loco(X0, X1, max_n=mmd_n, seed=seed)
    coord = coord_mmd_vimp(X0, X1, max_n=mmd_n, seed=seed)
    _, coord_share = modality_mass(coord)
    mmd_share = coord_share  # FSDS-style: coordinate MMD aggregated to modality
    po, _, _ = crossfit_po(X, bundle.y, bundle.W, seed=seed, n_estimators=max(40, n_estimators // 2))
    r_full, logo, po_vimp, po_perm = logo_po_risk(X, po, seed=seed + 4, n_estimators=n_estimators)
    logo_share = _share_from_contrib(logo)
    perm_share = _share_from_contrib(po_perm)
    _, po_feat_share = modality_mass(po_vimp)
    return {
        "rf_auc": rf_auc,
        "rf_vimp": rf_vimp,
        "rf_share": rf_share,
        "rf_feature_indices": top_local_indices(rf_vimp),
        "mmd_full": mmd_full,
        "mmd_logo": mmd_loco,
        "mmd_block": mmd_block,
        "mmd_block_share": _share_from_contrib(mmd_block),
        "mmd_share": mmd_share,
        "coord_mmd_vimp": coord,
        "coord_mmd_share": coord_share,
        "coord_feature_indices": top_local_indices(coord),
        "po_risk": r_full,
        "po_logo": logo,
        "po_logo_share": logo_share,
        "po_perm": po_perm,
        "po_perm_share": perm_share,
        "po_vimp": po_vimp,
        "po_feat_share": po_feat_share,
        "po_feature_indices": top_local_indices(po_vimp),
        "po": po,
    }


def top_local_indices(vimp, k=20):
    vimp = np.asarray(vimp, dtype=float)
    out = {}
    for name, sl in GROUPS.items():
        block = np.nan_to_num(vimp[sl], nan=-np.inf)
        k_use = min(k, block.size)
        out[name] = [int(i) for i in np.argsort(-block)[:k_use]]
    return out


def _contrib_to_vimp(contrib):
    v = np.zeros(P_X)
    for name, sl in GROUPS.items():
        v[sl] = contrib[name] / max((sl.stop - sl.start), 1)
    return v


def _share_from_contrib(contrib):
    # keep sign; normalize by L1 so shares compare across methods
    vals = {k: float(contrib[k]) for k in GROUP_NAMES}
    tot = sum(abs(v) for v in vals.values()) + 1e-12
    return {k: vals[k] / tot for k in GROUP_NAMES}


def _shares_from_split(X, y, W, seed, n_estimators, mmd_n):
    X = standardize_columns(X)
    X0, X1 = X[W == 0], X[W == 1]
    if min(len(X0), len(X1)) < 8:
        return None
    rf_vimp, auc = rf_domain(X0, X1, seed=seed, n_estimators=n_estimators)
    _, rf_share = modality_mass(rf_vimp)
    coord = coord_mmd_vimp(X0, X1, max_n=mmd_n, seed=seed)
    _, mmd_share = modality_mass(coord)
    po, _, _ = crossfit_po(X, y, W, seed=seed, n_splits=3, n_estimators=max(30, n_estimators // 2))
    r_full, logo, po_vimp, po_perm = logo_po_risk(X, po, seed=seed + 4, n_estimators=n_estimators)
    _, po_share = modality_mass(po_vimp)
    return {
        "auc": auc,
        "po_risk": r_full,
        "rf": rf_share,
        "mmd": mmd_share,
        "po": po_share,
        "po_logo": logo,
        "po_perm": po_perm,
    }


def cluster_bootstrap(bundle, B=10, seed=SEED, n_estimators=60, mmd_n=60):
    rng = np.random.default_rng(seed)
    recs = []
    for b in range(B):
        idx = cluster_indices(bundle.video_id, rng)
        got = _shares_from_split(
            bundle.X[idx], bundle.y[idx], bundle.W[idx], seed + b, n_estimators, mmd_n
        )
        if got is not None:
            recs.append(got)
    return recs


def within_video_permute_auc(bundle, B=10, seed=SEED, n_estimators=80):
    rng = np.random.default_rng(seed)
    X = standardize_columns(bundle.X)
    _, obs_auc = rf_domain(
        X[bundle.W == 0], X[bundle.W == 1], seed=seed, n_estimators=n_estimators
    )
    null = []
    for b in range(B):
        Wp = bundle.W.copy()
        for v in np.unique(bundle.video_id):
            idx = np.flatnonzero(bundle.video_id == v)
            Wp[idx] = rng.permutation(Wp[idx])
        if Wp.sum() < 8 or (len(Wp) - Wp.sum()) < 8:
            continue
        _, auc = rf_domain(X[Wp == 0], X[Wp == 1], seed=seed + b, n_estimators=n_estimators)
        null.append(auc)
    null = np.asarray(null, dtype=float)
    p = (1.0 + np.sum(null >= obs_auc)) / (len(null) + 1.0) if len(null) else float("nan")
    return float(obs_auc), null, float(p)


def pairwise_from_bootstrap(recs, method):
    pairs = (("video", "audio"), ("video", "text"), ("audio", "text"))
    out = {}
    raw_p = []
    names = []
    for a, b in pairs:
        diffs = np.array([r[method][a] - r[method][b] for r in recs], dtype=float)
        mean = float(np.mean(diffs))
        p = bootstrap_pvalue(diffs)
        sd = float(np.std(diffs, ddof=1)) if diffs.size > 1 else float("nan")
        var = float(np.var(diffs, ddof=1)) if diffs.size > 1 else float("nan")
        # B=10: interval from bootstrap SD rather than unstable percentiles.
        if np.isfinite(sd):
            lo, hi = mean - 1.96 * sd, mean + 1.96 * sd
        else:
            lo, hi = percentile_ci(diffs)
        out["%s-%s" % (a, b)] = {
            "mean_diff": mean,
            "sd": sd,
            "var": var,
            "ci95": [float(lo), float(hi)],
            "p_bootstrap": p,
            "excludes_zero": bool(lo > 0 or hi < 0),
        }
        raw_p.append(p)
        names.append("%s-%s" % (a, b))
    adj = holm_adjust(raw_p)
    for i, name in enumerate(names):
        out[name]["p_holm"] = float(adj[i])
    return out


def friedman_wilcoxon(per_video, method_share_key="share"):
    if len(per_video) < 3:
        return {"friedman_p": float("nan"), "wilcoxon": {}}
    mat = np.array([[r[method_share_key][g] for g in GROUP_NAMES] for r in per_video], dtype=float)
    try:
        stat, p = stats.friedmanchisquare(mat[:, 0], mat[:, 1], mat[:, 2])
        friedman = {"stat": float(stat), "p": float(p)}
    except Exception:
        friedman = {"stat": float("nan"), "p": float("nan")}
    wx = {}
    pairs = (("video", "audio"), ("video", "text"), ("audio", "text"))
    p_raw = []
    keys = []
    for i, (a, b) in enumerate(pairs):
        ia, ib = GROUP_NAMES.index(a), GROUP_NAMES.index(b)
        try:
            res = stats.wilcoxon(mat[:, ia], mat[:, ib], zero_method="pratt", alternative="two-sided")
            wx["%s-%s" % (a, b)] = {"stat": float(res.statistic), "p": float(res.pvalue)}
            p_raw.append(float(res.pvalue))
        except Exception:
            wx["%s-%s" % (a, b)] = {"stat": float("nan"), "p": float("nan")}
            p_raw.append(1.0)
        keys.append("%s-%s" % (a, b))
    adj = holm_adjust(p_raw)
    for i, k in enumerate(keys):
        wx[k]["p_holm"] = float(adj[i])
    return {"friedman": friedman, "wilcoxon": wx, "n_videos": int(len(per_video))}


def group_label_permutation_test(vimp, B=99, seed=SEED):
    rng = np.random.default_rng(seed)
    obs = group_partition_stat(vimp)
    null = np.array([permute_group_labels(vimp, rng) for _ in range(B)])
    p = (1.0 + np.sum(null >= obs)) / (B + 1.0)
    return {"obs_gap": float(obs), "p": float(p), "null_mean": float(null.mean())}


def _mean_sd_ci(samples):
    a = np.asarray(samples, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return {
            "mean": float("nan"),
            "sd": float("nan"),
            "var": float("nan"),
            "ci95": [float("nan"), float("nan")],
        }
    mean = float(np.mean(a))
    if a.size < 2:
        return {"mean": mean, "sd": float("nan"), "var": float("nan"), "ci95": [mean, mean]}
    sd = float(np.std(a, ddof=1))
    var = float(np.var(a, ddof=1))
    # B is small (default 10): report mean/variance and a normal interval from the bootstrap SE.
    return {
        "mean": mean,
        "sd": sd,
        "var": var,
        "ci95": [mean - 1.96 * sd, mean + 1.96 * sd],
    }


def summarize_bootstrap(recs):
    summary = {"n_valid": len(recs), "B": len(recs)}
    for method in ("rf", "mmd", "po"):
        summary[method] = {}
        for g in GROUP_NAMES:
            samples = np.array([r[method][g] for r in recs], dtype=float)
            summary[method][g] = _mean_sd_ci(samples)
        summary[method]["pairwise"] = pairwise_from_bootstrap(recs, method)
    aucs = np.array([r["auc"] for r in recs], dtype=float)
    summary["auc"] = _mean_sd_ci(aucs)
    return summary


def run_attribution(
    bundle: WindowBundle,
    B=10,
    B_perm=10,
    seed=SEED,
    n_estimators=120,
    bootstrap_estimators=60,
    mmd_n=80,
):
    point = pooled_point_estimates(bundle, seed=seed, n_estimators=n_estimators, mmd_n=mmd_n)
    videos = per_video_table(bundle, seed=seed, n_estimators=max(80, n_estimators // 2))
    boots = cluster_bootstrap(
        bundle, B=B, seed=seed + 17, n_estimators=bootstrap_estimators, mmd_n=min(mmd_n, 60)
    )
    boot_sum = summarize_bootstrap(boots)
    tests = {
        "per_video_rf_shares": friedman_wilcoxon(videos, "share"),
        "rf_group_label_perm": group_label_permutation_test(point["rf_vimp"], B=99, seed=seed),
        "po_group_label_perm": group_label_permutation_test(point["po_vimp"], B=99, seed=seed + 3),
        "coord_mmd_group_label_perm": group_label_permutation_test(
            point["coord_mmd_vimp"], B=99, seed=seed + 5
        ),
    }
    obs_auc, null_auc, p_auc = within_video_permute_auc(
        bundle, B=B_perm, seed=seed + 9, n_estimators=bootstrap_estimators
    )
    tests["within_video_auc_perm"] = {
        "obs_auc": obs_auc,
        "p": p_auc,
        "null_mean": float(np.mean(null_auc)) if len(null_auc) else float("nan"),
    }
    n_vid = len(np.unique(bundle.video_id))
    n = len(bundle.y)
    return {
        "n": int(n),
        "n_videos": int(n_vid),
        "n0": int((bundle.W == 0).sum()),
        "n1": int((bundle.W == 1).sum()),
        "layout": {"video": D_VIDEO, "audio": D_AUDIO, "text": D_TEXT, "label": 1},
        "point": {
            "rf_auc": point["rf_auc"],
            "rf_share": point["rf_share"],
            "rf_feature_indices": point["rf_feature_indices"],
            "mmd_full": point["mmd_full"],
            "mmd_logo": point["mmd_logo"],
            "mmd_block": point["mmd_block"],
            "mmd_block_share": point.get("mmd_block_share"),
            "mmd_share": point["mmd_share"],
            "coord_mmd_share": point["coord_mmd_share"],
            "coord_feature_indices": point["coord_feature_indices"],
            "po_risk": point["po_risk"],
            "po_logo": point["po_logo"],
            "po_logo_share": point["po_logo_share"],
            "po_perm": point["po_perm"],
            "po_perm_share": point["po_perm_share"],
            "po_feat_share": point["po_feat_share"],
            "po_feature_indices": point["po_feature_indices"],
        },
        "per_video": videos,
        "bootstrap": boot_sum,
        "tests": tests,
        "arrays": {
            "rf_vimp": point["rf_vimp"],
            "coord_mmd_vimp": point["coord_mmd_vimp"],
            "po_vimp": point["po_vimp"],
        },
        "meta": bundle.meta,
    }


def json_ready(obj):
    if isinstance(obj, dict):
        return {str(k): json_ready(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [json_ready(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if obj is None or isinstance(obj, (str, int, float, bool)):
        return obj
    return str(obj)


def write_json(path: Path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_ready(payload), indent=2))
