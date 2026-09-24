"""FW-theme: DiffDB theme sticker for flywheel step JSON.

Opportunity
-----------
Unsupervised prompt clusters (TF-IDF + k-means) as a *context sticker*
other agents can read — not a forecast, not AGOD FSDS.

Usage
-----
1. Fit ``ThemeSticker`` on DiffusionDB prompts (or any text list).
2. ``stick(text)`` → ``{cluster_id, top_terms, theme_tag}``.
3. ``annotate_steps(steps, texts)`` writes ``theme_*`` onto each step JSON
   and appends ``theme:<id>`` into ``decision_path``.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence

import numpy as np
from sklearn.cluster import MiniBatchKMeans
from sklearn.feature_extraction.text import TfidfVectorizer

from sandbox.prompt_themes import load_prompts, theme_clusters


class ThemeSticker:
    def __init__(
        self,
        *,
        k: int = 8,
        max_features: int = 3000,
        seed: int = 0,
    ):
        self.k = k
        self.max_features = max_features
        self.seed = seed
        self.vec: Optional[TfidfVectorizer] = None
        self.km: Optional[MiniBatchKMeans] = None
        self.top_terms: List[List[str]] = []
        self.fitted = False

    def fit(self, prompts: Sequence[str]) -> "ThemeSticker":
        prompts = [str(p) for p in prompts if str(p).strip()]
        if len(prompts) < self.k * 5:
            raise ValueError("need more prompts to fit ThemeSticker")
        self.vec = TfidfVectorizer(
            max_features=self.max_features,
            ngram_range=(1, 2),
            min_df=2,
            max_df=0.95,
            stop_words="english",
        )
        X = self.vec.fit_transform(prompts)
        self.km = MiniBatchKMeans(
            n_clusters=self.k, random_state=self.seed, batch_size=512, n_init=3
        )
        self.km.fit(X)
        terms = np.array(self.vec.get_feature_names_out())
        self.top_terms = []
        for j in range(self.k):
            order = np.argsort(self.km.cluster_centers_[j])[::-1][:6]
            self.top_terms.append(terms[order].tolist())
        self.fitted = True
        return self

    def stick(self, text: str) -> Dict[str, Any]:
        if not self.fitted or self.vec is None or self.km is None:
            raise RuntimeError("ThemeSticker not fitted")
        X = self.vec.transform([str(text)])
        cid = int(self.km.predict(X)[0])
        return {
            "theme_cluster": cid,
            "theme_tag": f"theme:{cid}",
            "theme_terms": list(self.top_terms[cid]),
        }

    def annotate_steps(
        self,
        steps: Sequence[Dict[str, Any]],
        texts: Sequence[str],
    ) -> List[Dict[str, Any]]:
        """Attach theme sticker to each step; pad/cycle texts if shorter."""
        if not texts:
            return [dict(s) for s in steps]
        out = []
        for i, s in enumerate(steps):
            text = texts[i % len(texts)]
            st = self.stick(text)
            sp = dict(s)
            sp.update(st)
            path = list(sp.get("decision_path") or [])
            # insert theme after observe if present
            if path and path[0] == "observe":
                path = [path[0], st["theme_tag"]] + path[1:]
            else:
                path = [st["theme_tag"]] + path
            sp["decision_path"] = path
            sp["theme_text_head"] = str(text)[:80]
            out.append(sp)
        return out


def fit_sticker_from_diffusiondb(
    *,
    max_n: int = 3000,
    k: int = 8,
    seed: int = 0,
) -> Optional[ThemeSticker]:
    prompts = load_prompts(max_n=max_n, seed=seed)
    if not prompts:
        return None
    return ThemeSticker(k=k, seed=seed).fit(prompts)


def theme_sticker_scorecard(
    *,
    max_n: int = 3000,
    k: int = 8,
    n_stick: int = 200,
    seed: int = 0,
) -> Dict[str, Any]:
    """Fit themes, stick onto synthetic flywheel-like step shells, summarize."""
    prompts = load_prompts(max_n=max_n, seed=seed)
    if not prompts:
        return {"ok": False, "reason": "no_prompts"}
    # pick k via quick silhouette helper if available
    meta = theme_clusters(prompts, k=k, seed=seed)
    sticker = ThemeSticker(k=k, seed=seed).fit(prompts)
    # Stratified sample across clusters so the scorecard isn't one-bucket.
    labels = []
    # predict in chunks
    chunk = 500
    for i in range(0, len(prompts), chunk):
        part = prompts[i : i + chunk]
        X = sticker.vec.transform(part)
        labels.extend(sticker.km.predict(X).tolist())
    labels = np.asarray(labels, dtype=int)
    rng = np.random.default_rng(seed)
    texts: List[str] = []
    per = max(1, n_stick // k)
    for cid in range(k):
        members = np.where(labels == cid)[0]
        if members.size == 0:
            continue
        take = rng.choice(members, size=min(per, members.size), replace=False)
        texts.extend(prompts[j] for j in take)
    if len(texts) < n_stick:
        extra = rng.choice(len(prompts), size=n_stick - len(texts), replace=True)
        texts.extend(prompts[j] for j in extra)
    texts = texts[:n_stick]
    rng.shuffle(texts)
    steps = []
    for i in range(n_stick):
        steps.append(
            {
                "step": i,
                "action": "idle",
                "model": "hgb",
                "residual": float(rng.normal()),
                "surprise": float(abs(rng.normal())),
                "decision_path": ["observe", "forecast:hgb", "decide:idle", "log"],
            }
        )
    annotated = sticker.annotate_steps(steps, texts)
    counts = np.bincount(
        [int(s["theme_cluster"]) for s in annotated], minlength=k
    ).tolist()
    # path must contain theme tag
    n_tagged = sum(1 for s in annotated if any(str(x).startswith("theme:") for x in s["decision_path"]))
    return {
        "ok": True,
        "k": k,
        "n_prompts_fit": len(prompts),
        "n_steps_stuck": n_stick,
        "n_paths_with_theme": n_tagged,
        "cluster_counts": counts,
        "top_terms": sticker.top_terms,
        "silhouette_ref": meta.get("silhouette_cosine"),
        "example": annotated[0] if annotated else None,
        "headline": (
            f"FW-theme: stuck theme:* on {n_tagged}/{n_stick} decision_paths; "
            f"k={k}, silhouette≈{meta.get('silhouette_cosine')}"
        ),
        "opportunity": "FW-theme",
        "note": "context sticker only — not a forecast claim",
    }
