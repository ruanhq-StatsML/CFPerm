"""Hard-sample ranking quality for PO-risk scores.

Design (intentionally blunt)
----------------------------
On an OnlineRFPerm-reject batch O:

  truth_i  = |Y_i − μ_oracle(X_i)|
             μ_oracle fit on recent ∪ current (uniform). Diagnostic only —
             it *sees* current labels, so it is a proxy for which rows are
             actually hard, not a deployable score.

  score_i  = PO-risk from f_ref / probe / re-fit μ0
             (no access to μ_oracle; may use Y only through its own residual)

We ask: does ``score`` rank the same rows that ``truth`` says are hard?

Primary metrics
---------------
- Spearman ρ     : full-order concordance (rank–rank Pearson)
- Precision@k    : |Top_k(score) ∩ Top_k(truth)| / k
- Recall@k       : same when |Top_score|=|Top_truth|=k (≡ Precision@k here)
- NDCG@k         : graded relevance = truth hardness
- AUROC (top-q%) : treat truth top-q% as positive class; score ranks them
- Lift@k         : Precision@k / (k/n)  — vs random

k defaults to ⌈n/5⌉ (top 20%); q same for AUROC positives.
"""
from __future__ import annotations

from typing import Dict, Mapping

import numpy as np


def _safe_corr(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 3 or np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def _ndcg_at_k(relev: np.ndarray, score: np.ndarray, k: int) -> float:
    """NDCG@k with graded relevance = ``relev`` (oracle hardness)."""
    relev = np.asarray(relev, float)
    score = np.asarray(score, float)
    k = int(min(k, len(relev)))
    if k <= 0:
        return float("nan")
    order = np.argsort(score)[::-1][:k]
    gains = relev[order]
    discounts = 1.0 / np.log2(np.arange(2, k + 2))
    dcg = float(np.sum(gains * discounts))
    ideal = np.sort(relev)[::-1][:k]
    idcg = float(np.sum(ideal * discounts))
    if idcg <= 1e-12:
        return float("nan")
    return dcg / idcg


def _auroc_binary(y_pos: np.ndarray, score: np.ndarray) -> float:
    """Mann–Whitney AUROC; y_pos ∈ {0,1}."""
    y = np.asarray(y_pos, int).ravel()
    s = np.asarray(score, float).ravel()
    pos = s[y == 1]
    neg = s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    # P(score_pos > score_neg) + 0.5 P(tie)
    # vectorized via ranks
    order = np.argsort(s)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(s) + 1, dtype=float)
    # average ranks for ties
    # (simple loop-free tie break via unique)
    _, inv, counts = np.unique(s, return_inverse=True, return_counts=True)
    sum_ranks = np.bincount(inv, weights=ranks)
    avg = sum_ranks / counts
    ranks = avg[inv]
    sum_pos = float(ranks[y == 1].sum())
    n_pos, n_neg = len(pos), len(neg)
    return (sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def hard_rank_metrics(
    score: np.ndarray,
    truth_hard: np.ndarray,
    *,
    frac: float = 0.2,
    ks: tuple[float, ...] = (0.1, 0.2, 0.3),
) -> Dict[str, float]:
    """Compare a PO score against oracle hardness on one reject batch.

    Parameters
    ----------
    score :
        Candidate hard-sample score (higher = harder), e.g. PO-risk.
    truth_hard :
        Oracle hardness (higher = harder), e.g. |Y − μ_oracle(X)|.
    frac :
        Default top fraction for Precision@k / AUROC positives (0.2 = top 20%).
    ks :
        Extra Precision@k fractions to report (as ``precision_at_10pct`` …).
    """
    score = np.asarray(score, float).ravel()
    truth = np.asarray(truth_hard, float).ravel()
    n = len(score)
    out: Dict[str, float] = {
        "n": float(n),
        "score_mean": float(np.mean(score)) if n else float("nan"),
        "truth_mean": float(np.mean(truth)) if n else float("nan"),
    }
    if n < 8:
        for key in (
            "spearman",
            "pearson",
            "precision_at_k",
            "lift_at_k",
            "ndcg_at_k",
            "auroc_topk",
        ):
            out[key] = float("nan")
        for f in ks:
            out[f"precision_at_{int(round(100 * f))}pct"] = float("nan")
        return out

    # Full ranking concordance
    r_s = score.argsort().argsort().astype(float)
    r_t = truth.argsort().argsort().astype(float)
    out["spearman"] = _safe_corr(r_s, r_t)
    out["pearson"] = _safe_corr(score, truth)

    # Default k = top-frac
    k = max(1, int(round(n * frac)))
    top_s = set(np.argsort(score)[-k:])
    top_t = set(np.argsort(truth)[-k:])
    prec = float(len(top_s & top_t) / k)
    out["precision_at_k"] = prec
    out["recall_at_k"] = prec  # equal-size tops ⇒ P=R
    out["lift_at_k"] = float(prec / (k / n))
    out["ndcg_at_k"] = _ndcg_at_k(truth, score, k)

    # AUROC: truth top-frac as positive
    y_pos = np.zeros(n, dtype=int)
    y_pos[np.argsort(truth)[-k:]] = 1
    out["auroc_topk"] = _auroc_binary(y_pos, score)

    # Extra Precision@k grid
    for f in ks:
        kk = max(1, int(round(n * f)))
        ps = set(np.argsort(score)[-kk:])
        pt = set(np.argsort(truth)[-kk:])
        out[f"precision_at_{int(round(100 * f))}pct"] = float(len(ps & pt) / kk)

    # backward-compatible alias used by earlier reports
    out["topk_overlap"] = out["precision_at_k"]
    return out


def aggregate_hard_rank(
    per_batch: list[Mapping[str, float]],
) -> Dict[str, float]:
    """Mean over reject batches (skips NaNs)."""
    if not per_batch:
        return {}
    keys = set()
    for d in per_batch:
        keys.update(d.keys())
    out: Dict[str, float] = {"n_batches": float(len(per_batch))}
    for key in sorted(keys):
        if key == "n":
            out["n_rows_mean"] = float(np.mean([d.get("n", np.nan) for d in per_batch]))
            continue
        vals = [float(d[key]) for d in per_batch if key in d and d[key] == d[key]]
        out[key] = float(np.mean(vals)) if vals else float("nan")
    return out


# Alias kept for older call sites
def po_quality_vs_truth(po: np.ndarray, truth_hard: np.ndarray) -> Dict[str, float]:
    return hard_rank_metrics(po, truth_hard)
