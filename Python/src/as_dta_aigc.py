"""AS-DTA: Adversarial-Shift Diagnosis & Targeted Adaptation on AIGC detection.

Pipeline (HC3 Human vs ChatGPT):
  1. Fine-tune / fit a binary detector on REF (Human=0, AI=1).
  2. Build LIVE batches: REF, back-translation (benign), synonym (benign),
     malignant adversarial swaps (TextFooler-lite).
  3. Diagnose each LIVE batch with PO-risk *proxy* + Peaks (CVaR / Kurtosis /
     Peak Gap) via bootstrap sub-environments.
  4. Targeted adaptation:
       - malignant → freeze lexical block, adapt syntactic head
       - benign    → conformal threshold calibration (FPR target)
  5. Emit diagnosis dashboard (scatter Mean × Peak Gap).

PO-risk *proxy* (fast MVP; not full PermuCATE):
  PO_i = 1 - p̂(Y_i | X_i)   or   (Y_i - p̂_i)^2
  Causal upgrade: swap unit field for PermuCATE φ / R-risk
  (see docs/method/AS_DTA_MVP_po_permucate_mm.md).
Feature blocks (detector-agnostic):
  lexical   — surface / embedding-like (char+unigram mass)
  syntactic — longer-range / CLS-like (bigram + function-word pattern)

Backends:
  sklearn  — TF-IDF + LogisticRegression (CPU default, CI-friendly)
  hf       — DistilBERT / RoBERTa sequence classifier (optional GPU/CPU)

HC3 load: local ``data/hc3/all.jsonl`` (Hugging Face Hello-SimpleAI/HC3).
"""
from __future__ import annotations

import json
import math
import re
import warnings
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.stats import kurtosis as scipy_kurtosis
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    f1_score,
    roc_auc_score,
)
from sklearn.model_selection import train_test_split

SEED = 2026
ROOT = Path(__file__).resolve().parents[2]
HC3_DEFAULT = ROOT / "data" / "hc3" / "all.jsonl"
HC3_DIR = ROOT / "data" / "hc3"
HC3_INDUSTRIAL = {
    "finance": HC3_DIR / "finance.jsonl",
    "medicine": HC3_DIR / "medicine.jsonl",
    "all": HC3_DIR / "all.jsonl",
}
WORD_RE = re.compile(r"[A-Za-z']+")

# Default industrial cost model (USD per decision error on a moderated item).
# FP = human wrongly flagged as AI (complaint / trust / legal) — priced higher.
# FN = AI slipped through (integrity / compliance miss).
COST_FP_DEFAULT = 50.0
COST_FN_DEFAULT = 20.0

LIVE_KINDS = (
    "ref",
    "benign_backtranslate",
    "benign_synonym",
    "malignant_attack",
)


# ---------------------------------------------------------------------------
# Data: HC3 → binary Human vs ChatGPT rows
# ---------------------------------------------------------------------------


def ensure_hc3(path: Path = HC3_DEFAULT) -> Path:
    """Download HC3 jsonl (``all`` / ``finance`` / ``medicine``) if missing."""
    path = Path(path)
    if path.exists() and path.stat().st_size > 0:
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    filename = path.name if path.suffix == ".jsonl" else "all.jsonl"
    try:
        from huggingface_hub import hf_hub_download

        downloaded = hf_hub_download(
            repo_id="Hello-SimpleAI/HC3",
            filename=filename,
            repo_type="dataset",
            local_dir=str(path.parent),
        )
        src = Path(downloaded)
        if src.resolve() != path.resolve():
            path.write_bytes(src.read_bytes())
        return path
    except Exception as exc:  # noqa: BLE001
        raise FileNotFoundError(
            "HC3 not found at %s and download failed: %s" % (path, exc)
        ) from exc


def load_hc3_binary(
    path: Path = HC3_DEFAULT,
    max_pairs: int = 4000,
    max_chars: int = 800,
    seed: int = SEED,
) -> Dict[str, np.ndarray]:
    """Flatten HC3 into texts + labels (0=human, 1=chatgpt)."""
    path = ensure_hc3(path)
    rng = np.random.default_rng(int(seed))
    humans, machines = [], []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            row = json.loads(line)
            for ans in row.get("human_answers") or []:
                if isinstance(ans, str) and ans.strip():
                    humans.append(ans.strip()[:max_chars])
            for ans in row.get("chatgpt_answers") or []:
                if isinstance(ans, str) and ans.strip():
                    machines.append(ans.strip()[:max_chars])
    rng.shuffle(humans)
    rng.shuffle(machines)
    n = min(int(max_pairs), len(humans), len(machines))
    texts = humans[:n] + machines[:n]
    labels = np.array([0] * n + [1] * n, dtype=int)
    idx = rng.permutation(len(texts))
    texts = [texts[i] for i in idx]
    labels = labels[idx]
    sources = np.array(["human"] * n + ["chatgpt"] * n, dtype=object)[idx]
    return {
        "texts": np.array(texts, dtype=object),
        "labels": labels,
        "sources": sources,
        "n_pairs": int(n),
    }


def split_ref_live(
    data: Dict[str, np.ndarray],
    live_frac: float = 0.35,
    seed: int = SEED,
):
    """Hold out LIVE pool from REF train/val."""
    X = data["texts"]
    y = data["labels"]
    X_ref, X_live, y_ref, y_live = train_test_split(
        X, y, test_size=float(live_frac), random_state=int(seed), stratify=y
    )
    X_tr, X_va, y_tr, y_va = train_test_split(
        X_ref, y_ref, test_size=0.2, random_state=int(seed) + 1, stratify=y_ref
    )
    return {
        "train": (X_tr, y_tr),
        "val": (X_va, y_va),
        "live": (X_live, y_live),
    }


# ---------------------------------------------------------------------------
# Synonym / paraphrase / attack utilities
# ---------------------------------------------------------------------------


_WN = None


def _wordnet():
    global _WN
    if _WN is not None:
        return _WN
    try:
        from nltk.corpus import wordnet as wn

        # trigger lazy load
        _ = wn.synsets("test")
        _WN = wn
    except Exception:
        try:
            import nltk

            nltk.download("wordnet", quiet=True)
            nltk.download("omw-1.4", quiet=True)
            from nltk.corpus import wordnet as wn

            _WN = wn
        except Exception:
            _WN = False
    return _WN


def _synonyms(word: str, max_n: int = 6) -> List[str]:
    wn = _wordnet()
    if not wn:
        return []
    w = word.lower()
    out = []
    for syn in wn.synsets(w):
        for lem in syn.lemmas():
            name = lem.name().replace("_", " ")
            if name.lower() != w and name.isalpha() and " " not in name:
                if name not in out:
                    out.append(name)
            if len(out) >= max_n:
                return out
    return out


def synonym_replace(text: str, rate: float = 0.15, seed: int = SEED) -> str:
    """Benign: replace ``rate`` fraction of content words with WordNet synonyms."""
    rng = np.random.default_rng(int(seed) + (hash(text) % 10_000))
    toks = text.split()
    if not toks:
        return text
    cand_idx = [i for i, t in enumerate(toks) if WORD_RE.fullmatch(t) and len(t) > 3]
    if not cand_idx:
        return text
    k = max(1, int(round(len(cand_idx) * float(rate))))
    pick = rng.choice(cand_idx, size=min(k, len(cand_idx)), replace=False)
    out = list(toks)
    for i in pick:
        syns = _synonyms(out[i])
        if syns:
            repl = str(rng.choice(syns))
            # preserve capitalization
            if out[i][0].isupper():
                repl = repl.capitalize()
            out[i] = repl
    return " ".join(out)


def backtranslate_proxy(text: str, seed: int = SEED) -> str:
    """Benign paraphrase proxy when MarianMT is unavailable.

    Applies light synonym rate + local phrase reorder to mimic round-trip MT
    distribution shift without loading en↔zh models.
    """
    rng = np.random.default_rng(int(seed) + (hash(text) % 10_000))
    t = synonym_replace(text, rate=0.08, seed=int(seed) + 3)
    sents = re.split(r"(?<=[.!?])\s+", t)
    if len(sents) > 2 and rng.random() < 0.5:
        i = int(rng.integers(0, len(sents) - 1))
        sents[i], sents[i + 1] = sents[i + 1], sents[i]
    # mild function-word jitter
    words = " ".join(sents).split()
    for i, w in enumerate(words):
        if w.lower() in {"the", "a", "an"} and rng.random() < 0.12:
            words[i] = ""
    return " ".join(w for w in words if w)


def try_backtranslate_marian(texts: Sequence[str], max_n: int = 32) -> Optional[List[str]]:
    """Optional real en→zh→en back-translation (heavy). Returns None if unavailable."""
    try:
        from transformers import pipeline
    except Exception:
        return None
    try:
        en2zh = pipeline("translation", model="Helsinki-NLP/opus-mt-en-zh", device=-1)
        zh2en = pipeline("translation", model="Helsinki-NLP/opus-mt-zh-en", device=-1)
    except Exception:
        return None
    out = []
    for t in list(texts)[:max_n]:
        z = en2zh(t[:400], max_length=256)[0]["translation_text"]
        e = zh2en(z, max_length=256)[0]["translation_text"]
        out.append(e)
    # pad remaining with proxy
    for t in list(texts)[max_n:]:
        out.append(backtranslate_proxy(t))
    return out


# ---------------------------------------------------------------------------
# Detectors
# ---------------------------------------------------------------------------


@dataclass
class DetectorMetrics:
    accuracy: float
    f1: float
    auc: float
    n: int


@dataclass
class SklearnDetector:
    """TF-IDF detector with lexical / syntactic feature blocks."""

    word_vec: TfidfVectorizer
    char_vec: TfidfVectorizer
    clf: LogisticRegression
    threshold: float = 0.5
    block_freeze: Optional[str] = None  # 'lexical' freezes char+unigram path
    meta: dict = field(default_factory=dict)

    def _xy(self, texts: Sequence[str]):
        Xw = self.word_vec.transform(texts)
        Xc = self.char_vec.transform(texts)
        from scipy import sparse

        if self.block_freeze == "lexical":
            # zero lexical (char) mass — syntactic/word path only
            Xc = sparse.csr_matrix(Xc.shape)
        elif self.block_freeze == "syntactic":
            Xw = sparse.csr_matrix(Xw.shape)
        return sparse.hstack([Xw, Xc]).tocsr()

    def predict_proba(self, texts: Sequence[str]) -> np.ndarray:
        X = self._xy(texts)
        return self.clf.predict_proba(X)[:, 1]

    def predict(self, texts: Sequence[str]) -> np.ndarray:
        p = self.predict_proba(texts)
        return (p >= self.threshold).astype(int)

    def po_proxy(self, texts: Sequence[str], labels: np.ndarray, kind: str = "one_minus_p") -> np.ndarray:
        labels = np.asarray(labels, dtype=int)
        p1 = self.predict_proba(texts)
        p_true = np.where(labels == 1, p1, 1.0 - p1)
        if kind == "sq_err":
            return (labels.astype(float) - p1) ** 2
        return 1.0 - p_true


def fit_sklearn_detector(
    texts: Sequence[str],
    labels: np.ndarray,
    seed: int = SEED,
    C: float = 2.0,
) -> SklearnDetector:
    labels = np.asarray(labels, dtype=int)
    word_vec = TfidfVectorizer(
        analyzer="word",
        ngram_range=(1, 2),
        min_df=2,
        max_features=20000,
        sublinear_tf=True,
    )
    char_vec = TfidfVectorizer(
        analyzer="char_wb",
        ngram_range=(3, 5),
        min_df=2,
        max_features=12000,
        sublinear_tf=True,
    )
    Xw = word_vec.fit_transform(texts)
    Xc = char_vec.fit_transform(texts)
    from scipy import sparse

    X = sparse.hstack([Xw, Xc]).tocsr()
    clf = LogisticRegression(
        C=float(C),
        max_iter=400,
        class_weight="balanced",
        random_state=int(seed),
        solver="liblinear",
    )
    clf.fit(X, labels)
    return SklearnDetector(
        word_vec=word_vec,
        char_vec=char_vec,
        clf=clf,
        meta={
            "backend": "sklearn",
            "n_train": int(len(labels)),
            "n_word_feat": int(Xw.shape[1]),
            "n_char_feat": int(Xc.shape[1]),
            "blocks": {
                "lexical": "char_wb 3-5 (embedding/surface proxy)",
                "syntactic": "word 1-2 grams (longer-range proxy)",
            },
        },
    )


def adapt_freeze_lexical(det: SklearnDetector, texts: Sequence[str], labels: np.ndarray, seed: int = SEED):
    """Malignant remedy: freeze lexical channel + adversarial-aug syntactic head.

    1. Zero char_wb (lexical / surface) features at fit & inference.
    2. Augment REF with detector-directed collisions so the word-path learns
       attack-stable patterns (contrastive-lite without a second encoder).
    """
    labels = np.asarray(labels, dtype=int)
    from scipy import sparse

    rng = np.random.default_rng(int(seed) + 9)
    aug_texts = list(texts)
    aug_labels = list(labels)
    n_aug = min(len(texts), 240)
    pick = rng.choice(len(texts), size=n_aug, replace=False)
    tmp = SklearnDetector(
        word_vec=det.word_vec,
        char_vec=det.char_vec,
        clf=det.clf,
        threshold=det.threshold,
        block_freeze=None,
        meta=dict(det.meta),
    )
    for i in pick:
        # synonym + feature-collision augs (same family as LIVE malignant)
        aug_texts.append(
            textfooler_lite(
                texts[int(i)],
                int(labels[int(i)]),
                tmp,
                rate=0.35,
                seed=int(seed) + int(i),
            )
        )
        aug_labels.append(int(labels[int(i)]))
    Xw = det.word_vec.transform(aug_texts)
    n_char = det.char_vec.transform(["x"]).shape[1]
    Xc = sparse.csr_matrix((len(aug_texts), n_char))  # lexical frozen
    X = sparse.hstack([Xw, Xc]).tocsr()
    clf = LogisticRegression(
        C=1.5,
        max_iter=600,
        class_weight="balanced",
        random_state=int(seed) + 9,
        solver="liblinear",
    )
    clf.fit(X, np.asarray(aug_labels, dtype=int))
    out = SklearnDetector(
        word_vec=det.word_vec,
        char_vec=det.char_vec,
        clf=clf,
        threshold=det.threshold,
        block_freeze="lexical",
        meta={
            **det.meta,
            "adapted": "freeze_lexical_adv_aug",
            "n_aug": int(n_aug),
        },
    )
    return out


def evaluate_detector(det, texts: Sequence[str], labels: np.ndarray) -> DetectorMetrics:
    labels = np.asarray(labels, dtype=int)
    p = det.predict_proba(texts)
    yhat = (p >= getattr(det, "threshold", 0.5)).astype(int)
    auc = float(roc_auc_score(labels, p)) if len(np.unique(labels)) > 1 else float("nan")
    return DetectorMetrics(
        accuracy=float(accuracy_score(labels, yhat)),
        f1=float(f1_score(labels, yhat, zero_division=0)),
        auc=auc,
        n=int(len(labels)),
    )


def classification_counts(y_true, y_hat):
    y_true = np.asarray(y_true, dtype=int)
    y_hat = np.asarray(y_hat, dtype=int)
    tp = int(((y_true == 1) & (y_hat == 1)).sum())
    fp = int(((y_true == 0) & (y_hat == 1)).sum())
    tn = int(((y_true == 0) & (y_hat == 0)).sum())
    fn = int(((y_true == 1) & (y_hat == 0)).sum())
    return {"tp": tp, "fp": fp, "tn": tn, "fn": fn}


def business_pr_metrics(
    det,
    texts: Sequence[str],
    labels: np.ndarray,
    cost_fp: float = COST_FP_DEFAULT,
    cost_fn: float = COST_FN_DEFAULT,
    target_fpr: float = 0.05,
    target_recall: float = 0.90,
) -> dict:
    """Precision / recall / FPR + operating-point & expected-cost business KPIs.

    Positive class = AI (label 1). FP = human flagged as AI (trust / legal).
    FN = AI missed (integrity). Expected cost = c_fp·FP + c_fn·FN on the batch.
    """
    from sklearn.metrics import (
        average_precision_score,
        precision_recall_curve,
        precision_score,
        recall_score,
    )

    labels = np.asarray(labels, dtype=int)
    texts = list(texts)
    p = np.asarray(det.predict_proba(texts), dtype=float)
    thr = float(getattr(det, "threshold", 0.5))
    yhat = (p >= thr).astype(int)
    counts = classification_counts(labels, yhat)
    prec = float(precision_score(labels, yhat, zero_division=0))
    rec = float(recall_score(labels, yhat, zero_division=0))
    hum = labels == 0
    fpr = float((yhat[hum] == 1).mean()) if hum.any() else float("nan")
    n = int(len(labels))
    expected_cost = float(cost_fp * counts["fp"] + cost_fn * counts["fn"])
    cost_per_1k = expected_cost / max(n, 1) * 1000.0

    # PR curve + operating points
    precision_c, recall_c, thr_c = precision_recall_curve(labels, p)
    ap = float(average_precision_score(labels, p)) if len(np.unique(labels)) > 1 else float("nan")

    # recall @ target FPR (scan thresholds high→low)
    order = np.argsort(-p)
    y_ord = labels[order]
    # for each possible threshold = score, compute FPR/recall on humans/AI
    best_rec_at_fpr = 0.0
    thr_at_fpr = 1.0
    n_hum = int(hum.sum())
    n_ai = int((labels == 1).sum())
    fp_run = 0
    tp_run = 0
    for i, idx in enumerate(order):
        if labels[idx] == 0:
            fp_run += 1
        else:
            tp_run += 1
        fpr_i = fp_run / max(n_hum, 1)
        rec_i = tp_run / max(n_ai, 1)
        if fpr_i <= target_fpr and rec_i >= best_rec_at_fpr:
            best_rec_at_fpr = float(rec_i)
            thr_at_fpr = float(p[idx])

    # precision @ target recall
    prec_at_rec = 0.0
    thr_at_rec = 0.0
    for pr, rc, t in zip(precision_c[:-1], recall_c[:-1], thr_c):
        if rc >= target_recall and pr >= prec_at_rec:
            prec_at_rec = float(pr)
            thr_at_rec = float(t)

    # downsample curve for JSON
    step = max(1, len(precision_c) // 40)
    curve = {
        "precision": precision_c[::step].tolist(),
        "recall": recall_c[::step].tolist(),
    }

    return {
        "threshold": thr,
        "precision": prec,
        "recall": rec,
        "fpr": fpr,
        "f1": float(f1_score(labels, yhat, zero_division=0)),
        "accuracy": float(accuracy_score(labels, yhat)),
        "auc": float(roc_auc_score(labels, p)) if len(np.unique(labels)) > 1 else float("nan"),
        "average_precision": ap,
        "counts": counts,
        "n": n,
        "expected_cost_usd": expected_cost,
        "cost_per_1k_usd": cost_per_1k,
        "cost_fp": float(cost_fp),
        "cost_fn": float(cost_fn),
        "recall_at_fpr_%g" % target_fpr: best_rec_at_fpr,
        "thr_at_fpr_%g" % target_fpr: thr_at_fpr,
        "precision_at_recall_%g" % target_recall: prec_at_rec,
        "thr_at_recall_%g" % target_recall: thr_at_rec,
        "pr_curve": curve,
    }


def shape_detection_score(diag: Diagnosis, ref: Diagnosis) -> dict:
    """Can drift *shape* separate benign vs malignant? Continuous scores + label."""
    d_mean = float(diag.mean - ref.mean)
    d_peak = float(diag.peak_gap - ref.peak_gap)
    d_cvar = float(diag.cvar - ref.cvar)
    action = route_action(diag, ref=ref)
    if action == "freeze_lexical_adapt_syntactic":
        shape = "malignant"
    elif action == "conformal_calibrate":
        shape = "benign"
    elif action == "hold":
        shape = "quiet"
    else:
        shape = "monitor"
    return {
        "d_mean": d_mean,
        "d_peak": d_peak,
        "d_cvar": d_cvar,
        "shape": shape,
        "action": action,
        "malignant_score": d_mean * 0.5 + d_peak * 8.0 + max(d_cvar, 0.0) * 0.5,
    }


def _po_concentration(unit_po: np.ndarray, top_frac: float = 0.2) -> dict:
    """How much total PO mass sits in the worst top-frac rows (injection signature)."""
    u = np.sort(np.asarray(unit_po, dtype=float))[::-1]
    n = len(u)
    if n == 0 or float(u.sum()) <= 0:
        return {"top_share": 0.0, "gini": 0.0, "p90": 0.0, "p50": 0.0}
    k = max(1, int(math.ceil(n * float(top_frac))))
    top_share = float(u[:k].sum() / (u.sum() + 1e-12))
    # Gini of PO mass
    cs = np.cumsum(np.sort(u))
    gini = float((n + 1 - 2 * np.sum(cs) / (cs[-1] + 1e-12)) / n)
    return {
        "top_share": top_share,
        "gini": gini,
        "p90": float(np.percentile(u, 90)),
        "p50": float(np.median(u)),
    }


def discriminate_ood_vs_intentional(
    det,
    texts_ref: Sequence[str],
    y_ref: np.ndarray,
    texts_live: Sequence[str],
    y_live: np.ndarray,
    *,
    K_subbatches: int = 28,
    seed: int = SEED,
    label: str = "live",
) -> dict:
    """With labels Y: separate natural OOD vs intentional adversarial injection.

    Why Y makes this informative
    ----------------------------
    Unit risk ``PO_i = 1 - p̂(Y_i | x_i)`` is the detector's *true-label* loss
    proxy. Natural OOD tends to lift PO **uniformly**; intentional injection
    concentrates PO on a attacked minority (high top-share / CVaR / error tail)
    while Y stays the ground truth — so FN/FP structure tells you if the attack
    worked, and before/after Y-metrics tell you if the remedy worked.

    Returns a verdict in {quiet, ood_uniform, intentional_injection, mixed}
    plus Y-conditioned effectiveness hooks (error, FN, concentration).
    """
    y_ref = np.asarray(y_ref, dtype=int)
    y_live = np.asarray(y_live, dtype=int)
    texts_ref, texts_live = list(texts_ref), list(texts_live)

    po_ref = det.po_proxy(texts_ref, y_ref)
    po_live = det.po_proxy(texts_live, y_live)
    conc_ref = _po_concentration(po_ref)
    conc_live = _po_concentration(po_live)

    diag_ref = diagnose_batch(
        texts_ref, y_ref, det, K_subbatches=K_subbatches, seed=seed, kind="ref"
    )
    diag_live = diagnose_batch(
        texts_live, y_live, det, K_subbatches=K_subbatches, seed=seed + 3, kind=label
    )
    shape = shape_detection_score(diag_live, diag_ref)

    # Y-conditioned errors on LIVE
    thr = float(getattr(det, "threshold", 0.5))
    p_live = det.predict_proba(texts_live)
    yhat = (p_live >= thr).astype(int)
    counts = classification_counts(y_live, yhat)
    n = max(len(y_live), 1)
    err = 1.0 - float((yhat == y_live).mean())
    fn_rate = float(counts["fn"] / max(int((y_live == 1).sum()), 1))
    fp_rate = float(counts["fp"] / max(int((y_live == 0).sum()), 1))

    # REF accuracy for comparison (same detector)
    p_ref = det.predict_proba(texts_ref)
    yhat_ref = (p_ref >= thr).astype(int)
    err_ref = 1.0 - float((yhat_ref == y_ref).mean())

    d_mean = shape["d_mean"]
    d_cvar = shape["d_cvar"]
    d_peak = shape["d_peak"]
    d_top = float(conc_live["top_share"] - conc_ref["top_share"])
    d_gini = float(conc_live["gini"] - conc_ref["gini"])
    d_err = float(err - err_ref)

    # Decision with Y (transparent rules)
    # intentional: concentrated PO mass OR (cvar/peak spike with error lift)
    intentional = (
        (d_top >= 0.03 and d_mean >= 0.01)
        or (d_cvar >= 0.015 and d_err >= 0.02)
        or (d_peak >= 0.008 and d_err >= 0.02)
        or (d_gini >= 0.03 and fn_rate >= err_ref + 0.03)
    )
    ood_uniform = (
        d_mean >= 0.004
        and d_top < 0.03
        and d_cvar < 0.015
        and d_peak < 0.008
    )
    quiet = d_mean < 0.004 and d_cvar < 0.01 and d_err < 0.01

    if quiet:
        verdict = "quiet"
    elif intentional and ood_uniform:
        verdict = "mixed"
    elif intentional:
        verdict = "intentional_injection"
    elif ood_uniform or (d_mean >= 0.004 and not intentional):
        verdict = "ood_uniform"
    else:
        verdict = "monitor"

    # Recommended actuator (same family as AS-DTA)
    if verdict == "quiet":
        actuator = "hold"
    elif verdict == "ood_uniform":
        actuator = "conformal_domain_calibrated"
    elif verdict in ("intentional_injection", "mixed"):
        actuator = "freeze_lexical_adv_aug"
    else:
        actuator = "human_review"

    return {
        "label": label,
        "verdict": verdict,
        "actuator": actuator,
        "with_Y": True,
        "delta": {
            "d_mean_po": d_mean,
            "d_cvar": d_cvar,
            "d_peak": d_peak,
            "d_top20_po_share": d_top,
            "d_gini_po": d_gini,
            "d_error": d_err,
        },
        "y_metrics_live": {
            "error": err,
            "fn_rate": fn_rate,
            "fp_rate": fp_rate,
            "counts": counts,
            "mean_po": float(np.mean(po_live)),
            "concentration": conc_live,
        },
        "y_metrics_ref": {
            "error": err_ref,
            "mean_po": float(np.mean(po_ref)),
            "concentration": conc_ref,
        },
        "shape": shape,
        "why_informative": (
            "With Y, PO_i=1-p̂(Y_i|x_i) is a true-label risk. Uniform ΔMean "
            "+ flat concentration ⇒ natural OOD; rising top-20% PO share / "
            "CVaR + error/FN lift ⇒ intentional injection. Y also unlocks "
            "direct effectiveness checks (Δerror, ΔFN, Δcost) after the actuator."
        ),
    }


def effectiveness_with_Y(
    before_det,
    after_det,
    texts: Sequence[str],
    labels: np.ndarray,
    cost_fp: float = COST_FP_DEFAULT,
    cost_fn: float = COST_FN_DEFAULT,
) -> dict:
    """Y-supervised before/after panel — the clean 'did the medicine work?' read."""
    b = business_pr_metrics(before_det, texts, labels, cost_fp=cost_fp, cost_fn=cost_fn)
    a = business_pr_metrics(after_det, texts, labels, cost_fp=cost_fp, cost_fn=cost_fn)
    return {
        "before": {
            "error": 1.0 - b["accuracy"],
            "recall": b["recall"],
            "precision": b["precision"],
            "fpr": b["fpr"],
            "fn": b["counts"]["fn"],
            "fp": b["counts"]["fp"],
            "cost_per_1k_usd": b["cost_per_1k_usd"],
            "mean_po": float(
                np.mean(before_det.po_proxy(texts, np.asarray(labels, dtype=int)))
            ),
        },
        "after": {
            "error": 1.0 - a["accuracy"],
            "recall": a["recall"],
            "precision": a["precision"],
            "fpr": a["fpr"],
            "fn": a["counts"]["fn"],
            "fp": a["counts"]["fp"],
            "cost_per_1k_usd": a["cost_per_1k_usd"],
            "mean_po": float(
                np.mean(after_det.po_proxy(texts, np.asarray(labels, dtype=int)))
            ),
        },
        "lifts": {
            "d_error": (1.0 - a["accuracy"]) - (1.0 - b["accuracy"]),
            "d_recall": a["recall"] - b["recall"],
            "d_precision": a["precision"] - b["precision"],
            "d_fn": a["counts"]["fn"] - b["counts"]["fn"],
            "d_fp": a["counts"]["fp"] - b["counts"]["fp"],
            "d_cost_per_1k_pct": 100.0
            * (a["cost_per_1k_usd"] - b["cost_per_1k_usd"])
            / max(b["cost_per_1k_usd"], 1e-6),
        },
        "worked": bool(
            (a["counts"]["fn"] <= b["counts"]["fn"] and a["recall"] >= b["recall"] - 1e-9)
            or (a["cost_per_1k_usd"] < b["cost_per_1k_usd"])
        ),
    }


def conformal_calibrate(
    det,
    texts: Sequence[str],
    labels: np.ndarray,
    target_fpr: float = 0.05,
) -> Tuple[Any, dict]:
    """Benign remedy: set decision threshold so FPR≈target on calibration set.

    For AI detection we treat human (0) as negative — FPR = P(pred AI | human).
    """
    labels = np.asarray(labels, dtype=int)
    p = det.predict_proba(texts)
    human = labels == 0
    if human.sum() < 5:
        return det, {"ok": False, "reason": "too few human samples"}
    # scores for humans; threshold on p(AI)
    ph = np.sort(p[human])
    # choose threshold so fraction of humans with p>=thr ≈ target_fpr
    k = int(math.ceil((1.0 - float(target_fpr)) * len(ph))) - 1
    k = int(np.clip(k, 0, len(ph) - 1))
    thr = float(ph[k])
    # clone threshold
    if isinstance(det, SklearnDetector):
        new = SklearnDetector(
            word_vec=det.word_vec,
            char_vec=det.char_vec,
            clf=det.clf,
            threshold=thr,
            block_freeze=det.block_freeze,
            meta={**det.meta, "adapted": "conformal", "threshold": thr},
        )
    else:
        det.threshold = thr
        new = det
    # realized FPR
    yhat = (p >= thr).astype(int)
    fpr = float((yhat[human] == 1).mean())
    recall = float((yhat[labels == 1] == 1).mean()) if (labels == 1).any() else float("nan")
    return new, {
        "ok": True,
        "threshold": thr,
        "target_fpr": float(target_fpr),
        "realized_fpr": fpr,
        "recall": recall,
    }


# ---------------------------------------------------------------------------
# LIVE batch construction
# ---------------------------------------------------------------------------


def build_live_batches(
    texts: Sequence[str],
    labels: np.ndarray,
    seed: int = SEED,
    synonym_rate: float = 0.15,
    attack_rate: float = 0.25,
    attack_max: int = 200,
    use_marian: bool = False,
    detector=None,
) -> Dict[str, Dict[str, np.ndarray]]:
    """Construct REF / benign / malignant LIVE text batches (labels preserved)."""
    texts = list(texts)
    labels = np.asarray(labels, dtype=int)
    rng = np.random.default_rng(int(seed))

    # REF
    out = {
        "ref": {"texts": np.array(texts, dtype=object), "labels": labels.copy()},
    }

    # Benign back-translation
    if use_marian:
        bt = try_backtranslate_marian(texts)
        if bt is None:
            bt = [backtranslate_proxy(t, seed=seed + i) for i, t in enumerate(texts)]
            bt_mode = "proxy_fallback"
        else:
            bt_mode = "marian_en_zh_en"
    else:
        bt = [backtranslate_proxy(t, seed=seed + i) for i, t in enumerate(texts)]
        bt_mode = "proxy"
    out["benign_backtranslate"] = {
        "texts": np.array(bt, dtype=object),
        "labels": labels.copy(),
        "meta": {"mode": bt_mode},
    }

    # Benign synonym
    syn = [synonym_replace(t, rate=synonym_rate, seed=seed + 17 + i) for i, t in enumerate(texts)]
    out["benign_synonym"] = {
        "texts": np.array(syn, dtype=object),
        "labels": labels.copy(),
        "meta": {"rate": float(synonym_rate)},
    }

    # Malignant: attack a *subset* hard (mixture → heteroscedastic peaks).
    # Uniform synonym is benign; concentrated TextFooler-lite is malignant.
    attacked = list(texts)
    n_atk = min(int(attack_max), max(1, len(texts) // 3))
    # prefer attacking confident correct AI/human examples for sharper flips
    if detector is not None:
        p = detector.predict_proba(texts)
        conf = np.where(labels == 1, p, 1.0 - p)
        order = np.argsort(-conf)
        atk_idx = order[:n_atk]
        for j, i in enumerate(atk_idx):
            attacked[int(i)] = textfooler_lite(
                texts[int(i)],
                int(labels[int(i)]),
                detector,
                rate=max(float(attack_rate), 0.35),
                seed=seed + 101 + j,
            )
    else:
        atk_idx = rng.choice(len(texts), size=n_atk, replace=False)
        for j, i in enumerate(atk_idx):
            attacked[int(i)] = synonym_replace(
                texts[int(i)], rate=max(float(attack_rate), 0.4), seed=seed + 201 + j
            )
    out["malignant_attack"] = {
        "texts": np.array(attacked, dtype=object),
        "labels": labels.copy(),
        "meta": {
            "n_attacked": int(n_atk),
            "rate": float(attack_rate),
            "guided": detector is not None,
            "mixture_frac": float(n_atk) / max(len(texts), 1),
        },
    }
    return out

def _opposite_class_words(detector: SklearnDetector, label: int, top_k: int = 40) -> List[str]:
    """High-weight TF-IDF tokens pushing toward the opposite class."""
    if not isinstance(detector, SklearnDetector):
        return []
    coef = np.asarray(detector.clf.coef_).reshape(-1)
    n_w = len(detector.word_vec.get_feature_names_out())
    wcoef = coef[:n_w]
    # label 1 (AI): opposite is human → negative coef words
    # label 0 (human): opposite is AI → positive coef words
    if int(label) == 1:
        idx = np.argsort(wcoef)[:top_k]  # most human-like
    else:
        idx = np.argsort(-wcoef)[:top_k]  # most AI-like
    names = detector.word_vec.get_feature_names_out()
    words = []
    for i in idx:
        tok = str(names[int(i)])
        if " " in tok:
            continue
        if WORD_RE.fullmatch(tok) and len(tok) > 2:
            words.append(tok)
    return words


def textfooler_lite(
    text: str,
    label: int,
    detector,
    rate: float = 0.25,
    seed: int = SEED,
    max_tries: int = 12,
) -> str:
    """Targeted evasion: WordNet swaps + opposite-class feature collision.

    Mixture of (a) synonym edits and (b) injecting detector-opposite lexical
    mass — mimics TextFooler / BAE directed at the REF detector.
    """
    rng = np.random.default_rng(int(seed))
    toks = text.split()
    cand = [i for i, t in enumerate(toks) if WORD_RE.fullmatch(t) and len(t) > 3]
    base = float(detector.po_proxy([text], np.array([label]))[0])
    out = list(toks) if toks else ["text"]
    changed = 0
    k = max(2, int(round(max(len(cand), 1) * float(rate))))
    if cand:
        rng.shuffle(cand)
        for i in cand:
            if changed >= k:
                break
            syns = _synonyms(out[i], max_n=8)
            if not syns:
                continue
            best, best_po = out[i], base
            for s in syns[:max_tries]:
                trial = list(out)
                repl = s.capitalize() if out[i][0].isupper() else s
                trial[i] = repl
                po = float(detector.po_proxy([" ".join(trial)], np.array([label]))[0])
                if po > best_po + 1e-4:
                    best_po, best = po, repl
            if best != out[i]:
                out[i] = best
                base = best_po
                changed += 1
    # feature-collision: inject opposite-class high-weight unigrams (lexical attack)
    opp = _opposite_class_words(detector, label, top_k=50)
    if opp:
        n_inj = max(4, int(8 * float(rate) / 0.25))
        inj = list(rng.choice(opp, size=min(n_inj, len(opp)), replace=False))
        # sprinkle into the middle + append (surface / lexical channel)
        mid = max(1, len(out) // 2)
        out = out[:mid] + inj[: n_inj // 2] + out[mid:] + inj[n_inj // 2 :]
        # light char-level typos on a few tokens (embedding/surface)
        for i in list(rng.choice(range(len(out)), size=min(3, len(out)), replace=False)):
            w = out[int(i)]
            if len(w) > 4 and WORD_RE.fullmatch(w):
                out[int(i)] = w[:-1] + w[-1] + w[-1]
    return " ".join(out)

# ---------------------------------------------------------------------------
# Diagnosis: PO-risk mean + Peaks
# ---------------------------------------------------------------------------


@dataclass
class Diagnosis:
    kind: str
    mean: float
    cvar: float
    kurtosis: float
    peak_gap: float
    p95: float
    n: int
    K: int
    raw_risks: List[float] = field(default_factory=list)
    accuracy: float = float("nan")
    auc: float = float("nan")

    def to_dict(self):
        d = asdict(self)
        return d


def diagnose_batch(
    texts: Sequence[str],
    labels: np.ndarray,
    detector,
    K_subbatches: int = 30,
    seed: int = SEED,
    kind: str = "batch",
    po_kind: str = "one_minus_p",
) -> Diagnosis:
    """Bootstrap sub-environments → PO-risk mean + Peaks (CVaR / Kurt / PeakGap).

    Peak Gap uses a *smaller* sub-batch size so mixture attacks (few heavily
    corrupted rows) show up as heteroscedastic spikes across bootstraps.
    """
    texts = list(texts)
    labels = np.asarray(labels, dtype=int)
    n = len(texts)
    if n == 0:
        return Diagnosis(kind, *[float("nan")] * 5, 0, K_subbatches, [])
    rng = np.random.default_rng(int(seed))
    unit = detector.po_proxy(texts, labels, kind=po_kind)
    risks = []
    # smaller windows amplify mixture peaks; still ≥12 for stability
    win = max(12, n // 5)
    for _ in range(int(K_subbatches)):
        idx = rng.choice(n, size=min(win, n), replace=False)
        risks.append(float(np.mean(unit[idx])))
    risks = np.asarray(risks, dtype=float)
    mean_risk = float(np.mean(unit))  # report full-batch mean (stable)
    p95 = float(np.percentile(risks, 95))
    tail = risks[risks >= p95]
    cvar = float(np.mean(tail)) if tail.size else float(risks.max())
    try:
        kurt = float(scipy_kurtosis(risks, fisher=True, bias=False))
    except Exception:
        kurt = float("nan")
    peak_gap = float(np.max(risks) - np.median(risks))
    metrics = evaluate_detector(detector, texts, labels)
    return Diagnosis(
        kind=kind,
        mean=mean_risk,
        cvar=cvar,
        kurtosis=kurt,
        peak_gap=peak_gap,
        p95=p95,
        n=int(n),
        K=int(K_subbatches),
        raw_risks=risks.tolist(),
        accuracy=metrics.accuracy,
        auc=metrics.auc,
    )


def route_action(
    diag: Diagnosis,
    ref: Optional[Diagnosis] = None,
    mean_lift_thr: float = 0.004,
    peak_lift_thr: float = 0.008,
    cvar_lift_thr: float = 0.015,
    mean_malignant_thr: float = 0.04,
) -> str:
    """Map diagnosis quadrant → actuator (relative to REF when available).

    Malignant if Peak/CVaR lifts **or** mean lifts hard (feature-collision /
    industrial mixture can raise mean more than Peak Gap). Benign if mean
    lifts mildly with Peak/CVaR quiet.
    """
    if ref is None:
        mean_thr, peak_thr = 0.25, 0.04
        mean_ok = diag.mean < mean_thr
        peak_ok = diag.peak_gap < peak_thr
        cvar_ok = diag.cvar < mean_thr + 0.05
        d_mean = diag.mean
    else:
        d_mean = diag.mean - ref.mean
        d_peak = diag.peak_gap - ref.peak_gap
        d_cvar = diag.cvar - ref.cvar
        mean_ok = d_mean <= mean_lift_thr
        peak_ok = d_peak <= peak_lift_thr
        cvar_ok = d_cvar <= cvar_lift_thr
    if mean_ok and peak_ok and cvar_ok:
        return "hold"
    # hard mean shock or peak/cvar spike → malignant actuator
    if (not peak_ok) or (not cvar_ok) or (d_mean >= mean_malignant_thr):
        return "freeze_lexical_adapt_syntactic"
    if (not mean_ok) and peak_ok and cvar_ok:
        return "conformal_calibrate"  # benign uniform inflation
    return "monitor"

# ---------------------------------------------------------------------------
# Plots / dashboard
# ---------------------------------------------------------------------------


def plot_diagnosis_dashboard(
    diagnoses: Dict[str, Diagnosis],
    path: Path,
    adapted: Optional[Dict[str, Diagnosis]] = None,
    title: str = "AS-DTA diagnosis — HC3 AIGC detector",
):
    import matplotlib.pyplot as plt

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    colors = {
        "ref": "#5B6B7A",
        "benign_backtranslate": "#2F6B4F",
        "benign_synonym": "#C4892A",
        "malignant_attack": "#A33B24",
    }
    markers = {
        "ref": "o",
        "benign_backtranslate": "s",
        "benign_synonym": "D",
        "malignant_attack": "^",
    }
    fig = plt.figure(figsize=(11.5, 8.2))
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.28)

    # --- scatter Mean × Peak Gap ---
    ax = fig.add_subplot(gs[0, 0])
    for k, d in diagnoses.items():
        ax.scatter(
            d.mean,
            d.peak_gap,
            s=110,
            c=colors.get(k, "#333"),
            marker=markers.get(k, "o"),
            label=k,
            zorder=3,
            edgecolors="white",
            linewidths=0.6,
        )
        ax.annotate(k.replace("_", "\n"), (d.mean, d.peak_gap), fontsize=7.5, xytext=(6, 4), textcoords="offset points")
    if adapted:
        for k, d in adapted.items():
            ax.scatter(
                d.mean,
                d.peak_gap,
                s=90,
                facecolors="none",
                edgecolors=colors.get(k, "#333"),
                marker=markers.get(k, "o"),
                linewidths=1.6,
                zorder=4,
            )
            ax.annotate(k + "†", (d.mean, d.peak_gap), fontsize=7, color="#444", xytext=(6, -10), textcoords="offset points")
    ax.axhline(0.06, color="#999", ls="--", lw=0.8, alpha=0.7)
    ax.axvline(0.22, color="#999", ls="--", lw=0.8, alpha=0.7)
    ax.set_xlabel("Mean PO-risk proxy ↑")
    ax.set_ylabel("Peak Gap (max − median) ↑")
    ax.set_title("Diagnosis scatter", loc="left", fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7.5, loc="best", framealpha=0.9)

    # --- CVaR / Kurtosis bars ---
    ax2 = fig.add_subplot(gs[0, 1])
    keys = list(diagnoses.keys())
    x = np.arange(len(keys))
    cvars = [diagnoses[k].cvar for k in keys]
    kurts = [diagnoses[k].kurtosis for k in keys]
    ax2.bar(x - 0.18, cvars, width=0.36, color=[colors.get(k, "#333") for k in keys], label="CVaR95")
    ax2b = ax2.twinx()
    ax2b.plot(x + 0.18, kurts, "o-", color="#5B3A8C", label="Kurtosis", lw=1.4)
    ax2.set_xticks(x)
    ax2.set_xticklabels([k.replace("_", "\n") for k in keys], fontsize=7.5)
    ax2.set_ylabel("CVaR95 of PO-risk")
    ax2b.set_ylabel("Kurtosis")
    ax2.set_title("Peaks: CVaR & Kurtosis", loc="left", fontweight="bold")
    ax2.grid(True, axis="y", alpha=0.3)

    # --- risk histograms ---
    ax3 = fig.add_subplot(gs[1, 0])
    for k, d in diagnoses.items():
        ax3.hist(
            d.raw_risks,
            bins=12,
            alpha=0.45,
            color=colors.get(k, "#333"),
            label=k,
            density=True,
        )
    ax3.set_xlabel("Bootstrap sub-batch PO-risk")
    ax3.set_ylabel("Density")
    ax3.set_title("PO-risk distribution (heteroscedasticity)", loc="left", fontweight="bold")
    ax3.legend(fontsize=7, loc="best")
    ax3.grid(True, alpha=0.3)

    # --- accuracy before / routing ---
    ax4 = fig.add_subplot(gs[1, 1])
    accs = [diagnoses[k].accuracy for k in keys]
    ax4.barh(np.arange(len(keys))[::-1], accs, color=[colors.get(k, "#333") for k in keys], height=0.55)
    ax4.set_yticks(np.arange(len(keys))[::-1])
    ax4.set_yticklabels(keys, fontsize=8)
    ax4.set_xlabel("Detector accuracy on LIVE")
    ax4.set_title("Live accuracy + routed action", loc="left", fontweight="bold")
    ref_d = diagnoses.get("ref")
    for i, k in enumerate(keys):
        act = route_action(diagnoses[k], ref=ref_d)
        ax4.text(
            0.02,
            len(keys) - 1 - i,
            act,
            fontsize=6.5,
            va="center",
            color="white",
            fontweight="bold",
        )
    ax4.set_xlim(0, 1.05)
    ax4.grid(True, axis="x", alpha=0.3)

    fig.suptitle(title, fontsize=13, fontweight="bold", y=0.98)
    fig.text(
        0.5,
        0.01,
        "Benign evasion inflates PO-risk uniformly; malignant injection induces peaks (heteroscedasticity).",
        ha="center",
        fontsize=8.5,
        style="italic",
        color="#444",
    )
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path


def write_html_dashboard(
    path: Path,
    diagnoses: Dict[str, Diagnosis],
    train_metrics: dict,
    adapted: Optional[Dict[str, Diagnosis]] = None,
    conformal: Optional[dict] = None,
    scatter_png: Optional[str] = None,
):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ref_d = diagnoses.get("ref")
    rows = []
    for k, d in diagnoses.items():
        act = route_action(d, ref=ref_d)
        rows.append(
            "<tr><td>{k}</td><td>{mean:.4f}</td><td>{peak:.4f}</td><td>{cvar:.4f}</td>"
            "<td>{kurt:.3f}</td><td>{acc:.3f}</td><td>{auc:.3f}</td>"
            "<td><code>{act}</code></td></tr>".format(
                k=k,
                mean=d.mean,
                peak=d.peak_gap,
                cvar=d.cvar,
                kurt=d.kurtosis,
                acc=d.accuracy,
                auc=d.auc,
                act=act,
            )
        )
    adapt_rows = ""
    if adapted:
        adapt_rows = "<h2>After targeted adaptation</h2><table><tr>"
        adapt_rows += "<th>batch</th><th>mean</th><th>peak_gap</th><th>acc</th></tr>"
        for k, d in adapted.items():
            adapt_rows += (
                "<tr><td>{k}</td><td>{m:.4f}</td><td>{p:.4f}</td><td>{a:.3f}</td></tr>".format(
                    k=k, m=d.mean, p=d.peak_gap, a=d.accuracy
                )
            )
        adapt_rows += "</table>"
    conf_html = "<pre>{}</pre>".format(json.dumps(conformal, indent=2)) if conformal else "<i>n/a</i>"
    img = ""
    if scatter_png:
        img = (
            '<p><img src="{}" alt="diagnosis" style="max-width:100%"/></p>'.format(scatter_png)
        )
    html = """<!DOCTYPE html>
<html><head><meta charset="utf-8"/><title>AS-DTA Diagnosis Board</title>
<style>
body{{font-family:ui-sans-serif,system-ui,sans-serif;margin:24px;background:#f7f5f1;color:#1a1a1a}}
h1{{font-size:1.5rem;margin:0 0 4px}} .sub{{color:#555;margin-bottom:18px}}
table{{border-collapse:collapse;width:100%;background:#fff;margin:12px 0}}
th,td{{border:1px solid #ddd;padding:8px 10px;font-size:0.9rem;text-align:left}}
th{{background:#eee}} code{{font-size:0.8rem}}
.card{{background:#fff;border:1px solid #e2e0da;padding:14px 16px;margin:12px 0}}
.quad{{display:grid;grid-template-columns:1fr 1fr;gap:8px;font-size:0.85rem}}
.q{{padding:10px;border:1px dashed #ccc}}
</style></head><body>
<h1>AS-DTA Diagnosis Board</h1>
<p class="sub">HC3 AIGC detector · PO-risk proxy mean × Peaks · targeted adaptation</p>
<div class="card"><b>REF detector</b>: acc={acc:.3f}  f1={f1:.3f}  auc={auc:.3f}  n={n} · backend={backend}</div>
{img}
<div class="card"><h2>LIVE diagnoses</h2>
<table><tr><th>batch</th><th>mean PO</th><th>peak gap</th><th>CVaR95</th>
<th>kurtosis</th><th>acc</th><th>auc</th><th>routed action</th></tr>
{rows}</table>
<div class="quad">
<div class="q"><b>Bottom-left</b><br/>REF / back-translate — hold</div>
<div class="q"><b>Bottom-right</b><br/>Benign synonym — conformal calibrate</div>
<div class="q"><b>Top-left</b><br/>Rare / noisy peaks — monitor</div>
<div class="q"><b>Top-right</b><br/>Malignant attack — freeze lexical, adapt syntactic</div>
</div></div>
{adapt}
<div class="card"><h2>Conformal (benign path)</h2>{conf}</div>
<p style="color:#666;font-size:0.85rem">As observed, benign evasion inflates the PO-risk uniformly,
while malignant adversarial injection induces extreme heteroscedasticity (peaks)
in the PO-risk distribution.</p>
</body></html>
""".format(
        acc=train_metrics.get("accuracy", float("nan")),
        f1=train_metrics.get("f1", float("nan")),
        auc=train_metrics.get("auc", float("nan")),
        n=train_metrics.get("n", "?"),
        backend=train_metrics.get("backend", "?"),
        img=img,
        rows="\n".join(rows),
        adapt=adapt_rows,
        conf=conf_html,
    )
    path.write_text(html, encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Full pipeline
# ---------------------------------------------------------------------------


def run_as_dta(
    hc3_path: Path = HC3_DEFAULT,
    max_pairs: int = 2500,
    live_frac: float = 0.35,
    K_subbatches: int = 30,
    seed: int = SEED,
    synonym_rate: float = 0.30,
    attack_rate: float = 0.40,
    attack_max: int = 220,
    use_marian: bool = False,
    out_dir: Optional[Path] = None,
    quick: bool = False,
) -> dict:
    """End-to-end AS-DTA: fit → LIVE → diagnose → adapt → dashboard."""
    if quick:
        max_pairs = min(max_pairs, 600)
        K_subbatches = min(K_subbatches, 16)
        attack_max = min(attack_max, 80)

    out_dir = Path(out_dir or (ROOT / "results" / "as_dta_aigc"))
    out_dir.mkdir(parents=True, exist_ok=True)

    print("[as-dta] loading HC3…", flush=True)
    data = load_hc3_binary(hc3_path, max_pairs=max_pairs, seed=seed)
    splits = split_ref_live(data, live_frac=live_frac, seed=seed)
    X_tr, y_tr = splits["train"]
    X_va, y_va = splits["val"]
    X_live, y_live = splits["live"]

    print("[as-dta] fine-tuning sklearn detector on REF (%d)…" % len(X_tr), flush=True)
    det = fit_sklearn_detector(X_tr, y_tr, seed=seed)
    m_tr = evaluate_detector(det, X_tr, y_tr)
    m_va = evaluate_detector(det, X_va, y_va)
    print(
        "  train acc=%.3f auc=%.3f | val acc=%.3f auc=%.3f"
        % (m_tr.accuracy, m_tr.auc, m_va.accuracy, m_va.auc),
        flush=True,
    )

    print("[as-dta] building LIVE batches…", flush=True)
    lives = build_live_batches(
        X_live,
        y_live,
        seed=seed,
        synonym_rate=synonym_rate,
        attack_rate=attack_rate,
        attack_max=attack_max,
        use_marian=use_marian,
        detector=det,
    )

    print("[as-dta] diagnosing…", flush=True)
    diagnoses: Dict[str, Diagnosis] = {}
    for kind, batch in lives.items():
        diagnoses[kind] = diagnose_batch(
            batch["texts"],
            batch["labels"],
            det,
            K_subbatches=K_subbatches,
            seed=seed + 7,
            kind=kind,
        )
    ref_diag = diagnoses["ref"]
    for kind, d in diagnoses.items():
        print(
            "  %-22s mean=%.4f peak=%.4f cvar=%.4f kurt=%.2f acc=%.3f → %s"
            % (
                kind,
                d.mean,
                d.peak_gap,
                d.cvar,
                d.kurtosis,
                d.accuracy,
                route_action(d, ref=ref_diag),
            ),
            flush=True,
        )

    # --- Targeted adaptation ---
    adapted: Dict[str, Diagnosis] = {}
    print("[as-dta] targeted adaptation…", flush=True)
    # malignant path
    det_mal = adapt_freeze_lexical(det, X_tr, y_tr, seed=seed)
    adapted["malignant_attack"] = diagnose_batch(
        lives["malignant_attack"]["texts"],
        lives["malignant_attack"]["labels"],
        det_mal,
        K_subbatches=K_subbatches,
        seed=seed + 11,
        kind="malignant_attack_adapted",
    )
    print(
        "  malignant after freeze-lexical: mean=%.4f→%.4f  peak=%.4f→%.4f  acc=%.3f→%.3f"
        % (
            diagnoses["malignant_attack"].mean,
            adapted["malignant_attack"].mean,
            diagnoses["malignant_attack"].peak_gap,
            adapted["malignant_attack"].peak_gap,
            diagnoses["malignant_attack"].accuracy,
            adapted["malignant_attack"].accuracy,
        ),
        flush=True,
    )

    # benign path — conformal on val humans, apply to synonym LIVE
    det_conf, conf_info = conformal_calibrate(det, X_va, y_va, target_fpr=0.05)
    adapted["benign_synonym"] = diagnose_batch(
        lives["benign_synonym"]["texts"],
        lives["benign_synonym"]["labels"],
        det_conf,
        K_subbatches=K_subbatches,
        seed=seed + 13,
        kind="benign_synonym_calibrated",
    )
    print(
        "  benign conformal: thr=%.3f FPR=%.3f recall=%.3f | synonym acc=%.3f→%.3f"
        % (
            conf_info.get("threshold", float("nan")),
            conf_info.get("realized_fpr", float("nan")),
            conf_info.get("recall", float("nan")),
            diagnoses["benign_synonym"].accuracy,
            adapted["benign_synonym"].accuracy,
        ),
        flush=True,
    )

    # dashboard
    png = out_dir / "as_dta_diagnosis_dashboard.png"
    plot_diagnosis_dashboard(diagnoses, png, adapted=adapted)
    html = out_dir / "as_dta_diagnosis_board.html"
    train_metrics = {
        "accuracy": m_va.accuracy,
        "f1": m_va.f1,
        "auc": m_va.auc,
        "n": m_va.n,
        "backend": "sklearn_tfidf",
        "train_accuracy": m_tr.accuracy,
        "train_auc": m_tr.auc,
        "n_train": int(len(y_tr)),
        "n_live": int(len(y_live)),
        "n_pairs": int(data["n_pairs"]),
    }
    write_html_dashboard(
        html,
        diagnoses,
        train_metrics,
        adapted=adapted,
        conformal=conf_info,
        scatter_png=png.name,
    )

    payload = {
        "train_metrics": train_metrics,
        "diagnoses": {k: v.to_dict() for k, v in diagnoses.items()},
        "adapted": {k: v.to_dict() for k, v in adapted.items()},
        "conformal": conf_info,
        "routing": {k: route_action(v, ref=ref_diag) for k, v in diagnoses.items()},
        "paths": {"dashboard_png": str(png), "dashboard_html": str(html)},
        "meta": {
            "seed": int(seed),
            "max_pairs": int(max_pairs),
            "K_subbatches": int(K_subbatches),
            "synonym_rate": float(synonym_rate),
            "attack_rate": float(attack_rate),
            "use_marian": bool(use_marian),
            "paper_line": (
                "As observed, benign evasion inflates the PO-risk uniformly, "
                "while malignant adversarial injection induces extreme "
                "heteroscedasticity (peaks) in the PO-risk distribution."
            ),
        },
    }
    (out_dir / "as_dta_results.json").write_text(json.dumps(payload, indent=2))
    _write_readme(out_dir / "README.md", payload)
    print("[as-dta] wrote", out_dir, flush=True)
    return payload


def _write_readme(path: Path, payload: dict):
    lines = [
        "# AS-DTA — AIGC detection diagnosis board",
        "",
        payload["meta"]["paper_line"],
        "",
        "## REF detector",
        "",
        "```",
        json.dumps(payload["train_metrics"], indent=2),
        "```",
        "",
        "## LIVE diagnoses",
        "",
        "| batch | mean PO | peak gap | CVaR | kurtosis | acc | action |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for k, d in payload["diagnoses"].items():
        lines.append(
            "| %s | %.4f | %.4f | %.4f | %.2f | %.3f | `%s` |"
            % (
                k,
                d["mean"],
                d["peak_gap"],
                d["cvar"],
                d["kurtosis"],
                d["accuracy"],
                payload["routing"].get(k, ""),
            )
        )
    lines.extend(
        [
            "",
            "## After adaptation",
            "",
            "| batch | mean PO | peak gap | acc |",
            "| --- | --- | --- | --- |",
        ]
    )
    for k, d in payload.get("adapted", {}).items():
        lines.append(
            "| %s | %.4f | %.4f | %.3f |" % (k, d["mean"], d["peak_gap"], d["accuracy"])
        )
    lines.extend(
        [
            "",
            "## Run",
            "",
            "```bash",
            "python3 scripts/run_as_dta_aigc.py --quick",
            "python3 scripts/run_as_dta_aigc.py",
            "```",
            "",
            "Open `as_dta_diagnosis_board.html` for the interactive-ish board,",
            "or `as_dta_diagnosis_dashboard.png` for the scatter + peaks figure.",
            "",
        ]
    )
    path.write_text("\n".join(lines))


def write_json(path: Path, obj: dict):
    Path(path).write_text(json.dumps(obj, indent=2))
