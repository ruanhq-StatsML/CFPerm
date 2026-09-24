"""Sandbox: DiffusionDB prompt theme clusters (unsupervised, not FSDS/PO).

Scenario
--------
Sample prompts from data/diffusiondb/metadata.parquet.
TF-IDF → MiniBatchKMeans → silhouette / top terms.

Deliberately far from AGOD temporal FSDS / transfer boards.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from sklearn.cluster import MiniBatchKMeans
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics import silhouette_score

ROOT = Path(__file__).resolve().parents[1]


def load_prompts(max_n: int = 4000, seed: int = 0) -> Optional[List[str]]:
    p = ROOT / "data/diffusiondb/metadata.parquet"
    if not p.is_file():
        return None
    df = pd.read_parquet(p)
    col = "prompt" if "prompt" in df.columns else None
    if col is None:
        for c in df.columns:
            if "prompt" in c.lower() or "text" in c.lower():
                col = c
                break
    if col is None:
        return None
    s = df[col].astype(str).fillna("")
    s = s[s.str.len() > 8]
    if len(s) == 0:
        return None
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(s), size=min(max_n, len(s)), replace=False)
    return s.iloc[idx].tolist()


def theme_clusters(
    prompts: List[str],
    *,
    k: int = 8,
    max_features: int = 4000,
    seed: int = 0,
) -> Dict[str, Any]:
    vec = TfidfVectorizer(
        max_features=max_features,
        ngram_range=(1, 2),
        min_df=3,
        max_df=0.9,
        stop_words="english",
    )
    X = vec.fit_transform(prompts)
    if X.shape[0] < k * 5:
        return {"ok": False, "reason": "too_few_docs"}
    km = MiniBatchKMeans(n_clusters=k, random_state=seed, batch_size=512, n_init=3)
    labels = km.fit_predict(X)
    # silhouette on a subsample (dense-ish)
    n_sil = min(1500, X.shape[0])
    rng = np.random.default_rng(seed)
    sil_idx = rng.choice(X.shape[0], size=n_sil, replace=False)
    try:
        sil = float(silhouette_score(X[sil_idx], labels[sil_idx], metric="cosine"))
    except Exception:
        sil = float("nan")

    terms = np.array(vec.get_feature_names_out())
    centers = km.cluster_centers_
    top = []
    sizes = []
    for j in range(k):
        order = np.argsort(centers[j])[::-1][:8]
        top.append(terms[order].tolist())
        sizes.append(int(np.sum(labels == j)))
    return {
        "ok": True,
        "n_docs": int(X.shape[0]),
        "k": k,
        "silhouette_cosine": sil,
        "cluster_sizes": sizes,
        "top_terms": top,
        "inertia": float(km.inertia_),
        "note": "TF-IDF MiniBatchKMeans — unsupervised themes, not FSDS/PO",
    }


def run_theme_job(
    *,
    max_n: int = 4000,
    k_grid: tuple = (6, 8, 10),
    seed: int = 0,
) -> Dict[str, Any]:
    prompts = load_prompts(max_n=max_n, seed=seed)
    if not prompts:
        return {"ok": False, "reason": "no_diffusiondb_prompts"}
    runs = [theme_clusters(prompts, k=k, seed=seed) for k in k_grid]
    ok = [r for r in runs if r.get("ok")]
    best = max(ok, key=lambda r: r.get("silhouette_cosine", float("-inf"))) if ok else None
    return {
        "ok": bool(ok),
        "n_prompts": len(prompts),
        "runs": runs,
        "best_k": None if best is None else best["k"],
        "best_silhouette": None if best is None else best["silhouette_cosine"],
        "best_top_terms_head": None if best is None else (best["top_terms"][:3]),
        "scenario": "DiffusionDB prompt theme clusters",
        "method": "TF-IDF + MiniBatchKMeans + silhouette",
        "distance_from_agod": "no temporal FSDS, no transfer AUC, no PO",
    }
