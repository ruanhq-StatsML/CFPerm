"""Text-only multiclass OVA on first-n new categories (no tabular).

Head is still CE: SGDClassifier(loss='log_loss').
Feature-level FSDS VIMP localizes; sample-level weights enter partial_fit.

Boards
------
  news20_comp   : known comp.* → novel sci.*
  news20_rec    : known rec.*  → novel talk.politics.*
  news20_mix    : known 8 topics → novel 8 topics
  cnae9         : company-doc BoW multiclass (text)
  amazon_cat    : Amazon review TF-IDF → product-category multiclass
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.preprocessing import LabelEncoder, StandardScaler

from first_category_performance import (
    LOSS_NAME,
    SEED,
    _sample_board_n_per_novel,
    evaluate_board_ova_n_first,
)
from fsds_cfperm_weight_benchmark import _openml_xy


# ---------------------------------------------------------------------------
# Text board builders
# ---------------------------------------------------------------------------


def _tfidf_xy(
    texts: Sequence[str],
    y: np.ndarray,
    *,
    max_features: int = 256,
    seed: int = SEED,
) -> Tuple[np.ndarray, np.ndarray]:
    vec = TfidfVectorizer(
        max_features=int(max_features),
        stop_words="english",
        min_df=2,
        ngram_range=(1, 2),
    )
    X = vec.fit_transform(list(texts)).toarray()
    X = StandardScaler(with_mean=False).fit_transform(X)  # sparse-friendly scale
    # densify already; use mean-scale for SGD stability
    X = StandardScaler().fit_transform(np.asarray(X, float))
    return X, np.asarray(y).ravel()


def _news20_board(
    *,
    name: str,
    known_names: Sequence[str],
    novel_names: Sequence[str],
    n_first: int,
    seed: int,
    max_features: int = 256,
) -> Optional[Dict]:
    cats = list(known_names) + list(novel_names)
    bun = fetch_20newsgroups(
        subset="all",
        categories=cats,
        remove=("headers", "footers", "quotes"),
        shuffle=True,
        random_state=int(seed),
    )
    # map to 0..K-1 in cats order
    name_to_i = {n: i for i, n in enumerate(cats)}
    y = np.array([name_to_i[bun.target_names[t]] for t in bun.target], dtype=int)
    X, y = _tfidf_xy(bun.data, y, max_features=max_features, seed=seed)
    known = list(range(len(known_names)))
    novel = list(range(len(known_names), len(cats)))
    return _sample_board_n_per_novel(
        X, y,
        known=known, novel=novel, n_first=n_first,
        n_ref=min(1200, max(200, 80 * len(known))),
        n_live_k=min(400, max(100, 40 * len(known))),
        seed=seed,
        name=name,
        note=f"20newsgroups TF-IDF; known={list(known_names)} → novel={list(novel_names)}",
    )


def board_news20_comp(*, seed: int = SEED, n_first: int = 100) -> Optional[Dict]:
    known = [
        "comp.graphics", "comp.os.ms-windows.misc",
        "comp.sys.ibm.pc.hardware", "comp.sys.mac.hardware", "comp.windows.x",
    ]
    novel = ["sci.crypt", "sci.electronics", "sci.med", "sci.space"]
    return _news20_board(
        name="news20_comp", known_names=known, novel_names=novel,
        n_first=n_first, seed=seed,
    )


def board_news20_rec(*, seed: int = SEED, n_first: int = 100) -> Optional[Dict]:
    known = ["rec.autos", "rec.motorcycles", "rec.sport.baseball", "rec.sport.hockey"]
    novel = ["talk.politics.guns", "talk.politics.mideast", "talk.politics.misc"]
    return _news20_board(
        name="news20_rec", known_names=known, novel_names=novel,
        n_first=n_first, seed=seed,
    )


def board_news20_mix(*, seed: int = SEED, n_first: int = 100) -> Optional[Dict]:
    known = [
        "alt.atheism", "comp.graphics", "misc.forsale", "rec.autos",
        "sci.med", "soc.religion.christian", "talk.politics.guns", "rec.sport.hockey",
    ]
    novel = [
        "comp.sys.mac.hardware", "rec.motorcycles", "sci.space", "sci.crypt",
        "talk.politics.mideast", "talk.religion.misc", "comp.windows.x", "rec.sport.baseball",
    ]
    return _news20_board(
        name="news20_mix", known_names=known, novel_names=novel,
        n_first=n_first, seed=seed, max_features=320,
    )


def board_cnae9_text(*, seed: int = SEED, n_first: int = 100) -> Optional[Dict]:
    """CNAE-9 Brazilian company documents — text/BoW multiclass."""
    X, y = _openml_xy("cnae-9", version=1)
    if X.shape[1] > 256:
        keep = np.argsort(-np.asarray(X, float).var(0))[:256]
        X = np.asarray(X, float)[:, keep]
    X = StandardScaler().fit_transform(np.asarray(X, float))
    known, novel = [0, 1, 2, 3, 4], [5, 6, 7, 8]
    return _sample_board_n_per_novel(
        X, y, known=known, novel=novel, n_first=n_first,
        n_ref=400, n_live_k=200, seed=seed,
        name="cnae9_text",
        note="CNAE-9 company-doc BoW multiclass (text); known 1–5 → novel 6–9",
    )


def board_amazon_category(
    *,
    seed: int = SEED,
    n_first: int = 100,
    max_features: int = 256,
    max_per_cat: int = 400,
) -> Optional[Dict]:
    """Amazon review text → product-category multiclass (category-shift cold-start)."""
    cache = Path(__file__).resolve().parents[2] / "data" / "amazon" / "heads"
    cats = [
        "Gift_Cards", "Digital_Music", "All_Beauty", "Software",
        "Subscription_Boxes", "Musical_Instruments", "Magazine_Subscriptions",
        "Handmade_Products", "Amazon_Fashion",
    ]
    texts: List[str] = []
    labels: List[str] = []
    rng = np.random.default_rng(int(seed))
    available = []
    for cat in cats:
        path = cache / f"{cat}.jsonl"
        if not path.exists():
            continue
        rows = []
        with path.open() as f:
            for line in f:
                try:
                    obj = json.loads(line)
                except Exception:
                    continue
                title = str(obj.get("title") or "")
                text = str(obj.get("text") or obj.get("body") or "")
                s = (title + " " + text).strip()
                if len(s) < 8:
                    continue
                rows.append(s)
        if len(rows) < int(n_first):
            continue
        if len(rows) > max_per_cat:
            pick = rng.choice(len(rows), size=max_per_cat, replace=False)
            rows = [rows[i] for i in pick]
        available.append(cat)
        texts.extend(rows)
        labels.extend([cat] * len(rows))

    if len(available) < 4:
        print("[text-ova] amazon_cat skip: need ≥4 cached categories with n≥n_first", flush=True)
        return None

    name_to_new = {available[i]: i for i in range(len(available))}
    y = np.array([name_to_new[c] for c in labels], dtype=int)
    X, y = _tfidf_xy(texts, y, max_features=max_features, seed=seed)
    # known = first half of available categories, novel = rest
    n_known = max(2, len(available) // 2)
    known = list(range(n_known))
    novel = list(range(n_known, len(available)))
    b = _sample_board_n_per_novel(
        X, y, known=known, novel=novel, n_first=n_first,
        n_ref=min(800, 100 * n_known), n_live_k=min(300, 50 * n_known),
        seed=seed,
        name="amazon_cat",
        note=f"Amazon review TF-IDF → category multiclass; known={available[:n_known]} → novel={available[n_known:]}",
    )
    if b:
        b["category_names"] = available
    return b


def load_text_boards_n_first(*, seed: int = SEED, n_first: int = 100) -> List[Dict]:
    """Text multiclass boards only — no tabular HAR/mice/jvowel/covtype/etc."""
    boards: List[Dict] = []
    for fn in (
        lambda: board_news20_comp(seed=seed, n_first=n_first),
        lambda: board_news20_rec(seed=seed, n_first=n_first),
        lambda: board_news20_mix(seed=seed, n_first=n_first),
        lambda: board_cnae9_text(seed=seed, n_first=n_first),
        lambda: board_amazon_category(seed=seed, n_first=n_first),
    ):
        try:
            b = fn()
            if b is not None:
                boards.append(b)
                print("[text-ova] loaded", b["name"], flush=True)
            else:
                print("[text-ova] skip (insufficient n_first):", fn.__name__ if hasattr(fn, "__name__") else fn, flush=True)
        except Exception as exc:  # noqa: BLE001
            print("[text-ova] failed:", exc, flush=True)
    return boards


def run_text_ova_nfirst(
    *,
    seed: int = SEED,
    n_first: int = 100,
    schemes: Sequence[str] = ("uniform", "weighted", "both_fused"),
) -> Dict:
    boards = load_text_boards_n_first(seed=seed, n_first=n_first)
    rows = []
    for b in boards:
        print("[text-ova] eval", b["name"], "n_first", n_first, "loss", LOSS_NAME, flush=True)
        rows.append(
            evaluate_board_ova_n_first(
                b, seed=seed, n_first=n_first, schemes=schemes,
            )
        )

    def _mean_delta(key, sch):
        vals = [r[key].get(sch) for r in rows if r.get(key, {}).get(sch) is not None]
        return float(np.mean(vals)) if vals else None

    return {
        "domain": "text_multiclass_only",
        "tabular_excluded": True,
        "loss": LOSS_NAME,
        "n_first": int(n_first),
        "seed": int(seed),
        "boards": rows,
        "n_boards": len(rows),
        "mean_delta_ova_auc_weighted": _mean_delta("delta_ova_auc_vs_uniform", "weighted"),
        "mean_delta_ova_auc_fused": _mean_delta("delta_ova_auc_vs_uniform", "both_fused"),
        "mean_delta_acc_first_weighted": _mean_delta("delta_acc_on_first_vs_uniform", "weighted"),
        "mean_delta_acc_first_fused": _mean_delta("delta_acc_on_first_vs_uniform", "both_fused"),
        "takeaway": (
            "Text multiclass only (no tabular). CE head; sample-level weights vs uniform; "
            "KPI = OVA AUC / acc on first n_first of each new category."
        ),
    }
