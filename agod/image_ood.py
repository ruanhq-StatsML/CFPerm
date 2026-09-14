"""Image-OOD benchmark: frozen embeddings + pseudo-outcome PO-risk.

1. Embedding x = frozen ViT (or cached CLIP img feat).
2. Fit pseudo-outcome μ on ID-train only.
3. Score PO-risk on ID-test ∪ OOD; report AUROC / FPR95 / AUPR + hard-rank.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

# pack -> (subdir, id_value, ood_value, domain_col)
CLIP_OOD_PACKS: Dict[str, Tuple[str, str, str, str]] = {
    "coco_outdoor_indoor": ("coco_outdoor_indoor", "outdoor", "indoor", "domain"),
    "coco_time_order": ("coco_time_order", "early_id", "late_id", "domain"),
    "indiana_cxr": ("indiana_cxr", "Frontal", "Lateral", "projection"),
    "fashion_iq": ("fashion_iq", "train", "test", "split"),
}
IMAGE_OOD_PACKS = CLIP_OOD_PACKS


@dataclass
class ImageOODSplit:
    X_id_train: np.ndarray
    y_id_train: np.ndarray
    X_id_test: np.ndarray
    y_id_test: np.ndarray
    X_ood: np.ndarray
    y_ood: np.ndarray
    pack: str
    id_domain: str
    ood_domain: str
    meta: dict


def coarse_labels(y: np.ndarray, top_k: int = 8) -> np.ndarray:
    y = np.asarray(y).ravel().astype(int)
    u, c = np.unique(y, return_counts=True)
    keep = set(u[np.argsort(-c)[:top_k]].tolist())
    remap = {lab: i for i, lab in enumerate(sorted(keep))}
    other = top_k
    return np.asarray([remap.get(int(v), other) for v in y], dtype=int)


def _load_meta(d: Path) -> pd.DataFrame:
    for name in ("df_metadata.csv", "df_metadata.parquet"):
        p = d / name
        if p.is_file():
            return pd.read_csv(p) if p.suffix == ".csv" else pd.read_parquet(p)
    raise FileNotFoundError(f"no metadata under {d}")


def _find_feat(d: Path) -> Path:
    for name in ("vit_feats.npy", "img_feats.npy", "clip_embedding.npy"):
        p = d / name
        if p.is_file():
            return p
    raise FileNotFoundError(f"no feature npy under {d}")


def _pca_id_only(X_tr, X_te, X_ood, pca_d: int, seed: int):
    if pca_d and X_tr.shape[1] > pca_d:
        pca = PCA(n_components=pca_d, random_state=seed)
        X_tr = pca.fit_transform(X_tr).astype(np.float32)
        X_te = pca.transform(X_te).astype(np.float32)
        X_ood = pca.transform(X_ood).astype(np.float32)
    return X_tr, X_te, X_ood


def _make_split(
    X, y, id_mask, ood_mask, *, pack, id_domain, ood_domain, pca_d, id_test_frac, seed, max_n, extra_meta=None
) -> ImageOODSplit:
    X_id, y_id = X[id_mask], y[id_mask]
    X_ood, y_ood = X[ood_mask], y[ood_mask]
    rng = np.random.default_rng(seed)
    if max_n is not None:
        if len(X_id) > max_n // 2:
            idx = rng.choice(len(X_id), max_n // 2, replace=False)
            X_id, y_id = X_id[idx], y_id[idx]
        if len(X_ood) > max_n // 2:
            idx = rng.choice(len(X_ood), max_n // 2, replace=False)
            X_ood, y_ood = X_ood[idx], y_ood[idx]
    strat = y_id if len(np.unique(y_id)) > 1 else None
    try:
        X_tr, X_te, y_tr, y_te = train_test_split(
            X_id, y_id, test_size=id_test_frac, random_state=seed, stratify=strat
        )
    except ValueError:
        X_tr, X_te, y_tr, y_te = train_test_split(
            X_id, y_id, test_size=id_test_frac, random_state=seed, stratify=None
        )
    X_tr, X_te, X_ood = _pca_id_only(X_tr, X_te, X_ood, pca_d, seed)
    meta = {
        "pca_d": int(X_tr.shape[1]),
        "n_id_train": int(len(X_tr)),
        "n_id_test": int(len(X_te)),
        "n_ood": int(len(X_ood)),
        "n_classes_id": int(len(np.unique(y_tr))),
        **(extra_meta or {}),
    }
    return ImageOODSplit(
        X_id_train=X_tr,
        y_id_train=np.asarray(y_tr, int),
        X_id_test=X_te,
        y_id_test=np.asarray(y_te, int),
        X_ood=X_ood,
        y_ood=np.asarray(y_ood, int),
        pack=pack,
        id_domain=id_domain,
        ood_domain=ood_domain,
        meta=meta,
    )


def load_clip_ood_pack(
    root: Path,
    pack: str,
    *,
    pca_d: int = 64,
    top_k_labels: int = 8,
    id_test_frac: float = 0.25,
    seed: int = 0,
    max_n: Optional[int] = None,
) -> ImageOODSplit:
    if pack not in CLIP_OOD_PACKS:
        raise KeyError(pack)
    sub, id_dom, ood_dom, dom_col = CLIP_OOD_PACKS[pack]
    d = Path(root) / "data/img_txt" / sub
    X = np.load(_find_feat(d)).astype(np.float32)
    y = coarse_labels(np.load(d / "labels.npy"), top_k=top_k_labels)
    meta_df = _load_meta(d)
    if dom_col not in meta_df.columns:
        raise KeyError(f"{pack}: missing {dom_col}")
    dom = meta_df[dom_col].astype(str).to_numpy()
    uniq = {v.lower(): v for v in np.unique(dom)}
    id_key = uniq.get(str(id_dom).lower(), str(id_dom))
    ood_key = uniq.get(str(ood_dom).lower(), str(ood_dom))
    return _make_split(
        X, y, dom == id_key, dom == ood_key,
        pack=pack, id_domain=str(id_key), ood_domain=str(ood_key),
        pca_d=pca_d, id_test_frac=id_test_frac, seed=seed, max_n=max_n,
        extra_meta={"backbone": "clip_cached", "feat": _find_feat(d).name},
    )


def load_food101_vit_pack(
    root: Path,
    *,
    pca_d: int = 64,
    id_frac_classes: float = 0.7,
    id_test_frac: float = 0.25,
    seed: int = 0,
    max_n: Optional[int] = None,
    model_name: str = "vit_tiny_patch16_224",
    extract_max_n: int = 8000,
) -> ImageOODSplit:
    from agod.vit_embed import embed_food101_parquets, make_class_holdout_domains

    root = Path(root)
    out_dir = root / "data/img_txt/food101_vit"
    pq_dir = root / "data/raw/food101/data"
    embed_food101_parquets(
        pq_dir, out_dir, model_name=model_name, max_n=extract_max_n, seed=seed
    )
    X = np.load(out_dir / "vit_feats.npy").astype(np.float32)
    y_raw = np.load(out_dir / "labels.npy").astype(int)
    meta_df = _load_meta(out_dir).copy()
    if "domain" in meta_df.columns and set(meta_df["domain"].astype(str)) >= {"id", "ood"}:
        dom = meta_df["domain"].astype(str).to_numpy()
    else:
        dom = make_class_holdout_domains(y_raw, id_frac_classes=id_frac_classes, seed=seed)
        meta_df["domain"] = dom
        meta_df.to_csv(out_dir / "df_metadata.csv", index=False)

    id_mask = dom == "id"
    ood_mask = dom == "ood"
    id_classes = sorted(np.unique(y_raw[id_mask]).tolist())
    remap = {c: i for i, c in enumerate(id_classes)}
    # OOD / held-out classes → −1 so nll/resid use label-free fallbacks
    y = np.asarray([remap.get(int(yi), -1) for yi in y_raw], dtype=int)
    return _make_split(
        X, y, id_mask, ood_mask,
        pack="food101_vit", id_domain="id_classes", ood_domain="heldout_classes",
        pca_d=pca_d, id_test_frac=id_test_frac, seed=seed, max_n=max_n,
        extra_meta={
            "backbone": model_name,
            "protocol": "class_holdout",
            "id_frac_classes": id_frac_classes,
            "n_id_classes": len(id_classes),
            "n_ood_classes": int(len(np.unique(y_raw[ood_mask]))),
        },
    )


def load_image_ood_pack(root: Path, pack: str, **kwargs) -> ImageOODSplit:
    kwargs = dict(kwargs)
    kwargs.pop("feat", None)
    if pack == "food101_vit":
        keys = ("pca_d", "id_frac_classes", "id_test_frac", "seed", "max_n", "model_name", "extract_max_n")
        return load_food101_vit_pack(root, **{k: kwargs[k] for k in keys if k in kwargs})
    keys = ("pca_d", "top_k_labels", "id_test_frac", "seed", "max_n")
    return load_clip_ood_pack(root, pack, **{k: kwargs[k] for k in keys if k in kwargs})


def fit_po_classifier(X, y, seed: int = 0) -> RandomForestClassifier:
    clf = RandomForestClassifier(
        n_estimators=100, max_depth=14, min_samples_leaf=2, random_state=seed, n_jobs=1
    )
    clf.fit(X, np.asarray(y, int))
    return clf


def fit_po_regressor(X, y, seed: int = 0) -> RandomForestRegressor:
    rf = RandomForestRegressor(
        n_estimators=80, max_depth=12, min_samples_leaf=2, random_state=seed, n_jobs=1
    )
    rf.fit(X, np.asarray(y, float))
    return rf


def score_po_msp(clf, X) -> np.ndarray:
    return 1.0 - clf.predict_proba(X).max(axis=1)


def score_po_nll(clf, X, y) -> np.ndarray:
    """1 − p(y|x). Unseen / sentinel labels (−1) fall back to 1 − max_c p(c|x)."""
    p = clf.predict_proba(X)
    classes = list(clf.classes_)
    y = np.asarray(y).ravel().astype(int)
    out = np.empty(len(y), float)
    for i, yi in enumerate(y):
        if yi in classes:
            out[i] = 1.0 - float(p[i, classes.index(yi)])
        else:
            # class-holdout OOD: no true ID label — use label-free uncertainty
            out[i] = 1.0 - float(p[i].max())
    return out


def score_po_resid(reg, X, y) -> np.ndarray:
    """|y − μ(x)|. Sentinel y < 0 → fractional distance of μ to nearest class index."""
    pred = np.asarray(reg.predict(X), float).ravel()
    y = np.asarray(y, float).ravel()
    out = np.abs(y - pred)
    unk = y < 0
    if np.any(unk):
        out[unk] = np.abs(pred[unk] - np.round(pred[unk]))
    return out


def score_energy(clf, X, T: float = 1.0) -> np.ndarray:
    """RF has no logits; use Shannon entropy as an energy / uncertainty proxy.

    Classical −T logsumexp(logit/T) is degenerate on already-normalized probs
    (sums to 1 → energy ≈ 0). Higher entropy ⇒ more OOD-like.
    """
    del T  # kept for API compatibility with logit-energy callers
    p = np.clip(clf.predict_proba(X), 1e-8, 1.0)
    return -np.sum(p * np.log(p), axis=1)


def score_dre_domain(X_id, X_eval, *, seed: int = 0) -> np.ndarray:
    Xs = StandardScaler().fit_transform(np.vstack([X_id, X_eval]))
    y = np.concatenate([np.zeros(len(X_id)), np.ones(len(X_eval))])
    clf = LogisticRegression(max_iter=500, random_state=seed)
    clf.fit(Xs, y)
    return clf.predict_proba(Xs[len(X_id):])[:, 1]


def fpr_at_tpr(y_true, score, *, tpr_level: float = 0.95) -> float:
    y = np.asarray(y_true).astype(int).ravel()
    s = np.asarray(score, float).ravel()
    pos, neg = s[y == 1], s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    thr = np.quantile(pos, 1.0 - tpr_level)
    return float(np.mean(neg >= thr))


def binary_ood_metrics(y_ood, score) -> Dict[str, float]:
    y = np.asarray(y_ood).astype(int).ravel()
    s = np.asarray(score, float).ravel()
    out = {
        "auroc": float("nan"),
        "aupr": float("nan"),
        "fpr95": float("nan"),
        "score_mean_id": float(np.mean(s[y == 0])) if np.any(y == 0) else float("nan"),
        "score_mean_ood": float(np.mean(s[y == 1])) if np.any(y == 1) else float("nan"),
    }
    if len(np.unique(y)) < 2:
        return out
    out["auroc"] = float(roc_auc_score(y, s))
    out["aupr"] = float(average_precision_score(y, s))
    out["fpr95"] = fpr_at_tpr(y, s, tpr_level=0.95)
    return out


def eval_scores_on_split(split: ImageOODSplit, scores: Dict[str, np.ndarray]):
    n_id, n_ood = len(split.X_id_test), len(split.X_ood)
    y = np.concatenate([np.zeros(n_id, int), np.ones(n_ood, int)])
    return {name: binary_ood_metrics(y, sc) for name, sc in scores.items()}


def list_available_packs(root: Path) -> List[str]:
    root = Path(root)
    out: List[str] = []
    if (root / "data/raw/food101/data").is_dir() or (
        root / "data/img_txt/food101_vit/vit_feats.npy"
    ).is_file():
        out.append("food101_vit")
    for pack, (sub, _, _, _) in CLIP_OOD_PACKS.items():
        d = root / "data/img_txt" / sub
        try:
            _find_feat(d)
            if (d / "labels.npy").is_file():
                out.append(pack)
        except FileNotFoundError:
            continue
    return out
