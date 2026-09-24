"""Transfer-probe null baselines and efficiency indices.

Statistical justification
-------------------------
Adjacent-chunk HGB/LogReg AUC answers: \"selected associations on chunk t
still rank Y on chunk t+1?\"  A raw AUC near 1.0 is *not* identifiable as
skill when the label is rare, features encode volume, or the scorer is
overfit to a trivial margin.

Under H0 (no predictive association that transfers), label-permutation
on the *test* chunk destroys Y|X dependence while keeping the X geometry
and the *trained* scorer fixed.  The resulting AUC distribution is a
cheap null for the same probe pipeline.

Definitions
-----------
- ``null_auc_mean``: mean AUC after ``n_perm`` shuffles of ``y_test``
- ``excess_auc``: ``auc_obs - null_auc_mean``  (skill beyond chance)
- ``probe_eff``: ``excess_auc / relative_flops`` where relative_flops
  proxies adaptation cost of this pair (selected_k * n_train * iters)

Use across packs (DiffusionDB tokens, metro sensors, ad edges) without
locking to one business scenario — same null, different story.

Calibration (Iter2)
-------------------
AUC / excess measure *ranking* skill.  Brier is a proper scoring rule but
does not localize mis-calibration.  Equal-width ECE on the *transfer*
probabilities (fit on chunk t, scored on t+1) answers: when the probe
says p, does frequency ≈ p on the next chunk?

High excess + high ECE → transferable ranking but untrustworthy probs
(common under volume features / rare labels).  Cross-pack comparable.
"""
from __future__ import annotations

from typing import Any, Dict, Optional, Sequence

import numpy as np
from sklearn.metrics import roc_auc_score


def permute_auc(
    y_true: np.ndarray,
    proba: np.ndarray,
    *,
    n_perm: int = 5,
    seed: int = 0,
) -> Dict[str, float]:
    """Label-permutation null for a fixed score vector on one test chunk."""
    y = np.asarray(y_true).astype(int)
    p = np.asarray(proba, dtype=float)
    if len(y) < 2 or len(np.unique(y)) < 2 or n_perm <= 0:
        return {
            "null_auc_mean": float("nan"),
            "null_auc_std": float("nan"),
            "n_perm": 0.0,
        }
    rng = np.random.default_rng(seed)
    vals = []
    for i in range(n_perm):
        y_s = rng.permutation(y)
        if len(np.unique(y_s)) < 2:
            continue
        try:
            vals.append(float(roc_auc_score(y_s, p)))
        except Exception:
            continue
    if not vals:
        return {
            "null_auc_mean": float("nan"),
            "null_auc_std": float("nan"),
            "n_perm": 0.0,
        }
    return {
        "null_auc_mean": float(np.mean(vals)),
        "null_auc_std": float(np.std(vals)),
        "n_perm": float(len(vals)),
    }


def excess_auc(auc_obs: float, null_auc_mean: float) -> float:
    if not np.isfinite(auc_obs) or not np.isfinite(null_auc_mean):
        return float("nan")
    return float(auc_obs - null_auc_mean)


def relative_probe_flops(
    n_train: int,
    selected_k: int,
    *,
    hgb_iters: int = 60,
    logreg_iters: int = 200,
) -> float:
    """Order-of-magnitude adaptation FLOPs proxy (not wall-clock).

    HGB leaf updates ~ n * k * iters; LogReg ~ n * k * iters.
    We sum both because the board always fits both probes on the same X.
    """
    n = max(int(n_train), 1)
    k = max(int(selected_k), 1)
    return float(n * k * (hgb_iters + logreg_iters))


def probe_efficiency(
    excess: float,
    flops: float,
    *,
    scale: float = 1e6,
) -> float:
    """Excess AUC per million relative FLOPs — higher = more skill per compute."""
    if not np.isfinite(excess) or not np.isfinite(flops) or flops <= 0:
        return float("nan")
    return float(excess / (flops / scale))


def expected_calibration_error(
    y_true: np.ndarray,
    proba: np.ndarray,
    *,
    n_bins: int = 10,
) -> float:
    """Equal-width ECE on [0, 1] for binary labels (transfer calibration)."""
    y = np.asarray(y_true).astype(float)
    p = np.asarray(proba, dtype=float)
    if len(y) < 2 or n_bins < 2:
        return float("nan")
    # clip for binning stability
    p = np.clip(p, 0.0, 1.0)
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    n = len(y)
    for i in range(n_bins):
        lo, hi = edges[i], edges[i + 1]
        if i == n_bins - 1:
            m = (p >= lo) & (p <= hi)
        else:
            m = (p >= lo) & (p < hi)
        if not np.any(m):
            continue
        conf = float(np.mean(p[m]))
        acc = float(np.mean(y[m]))
        ece += (float(np.sum(m)) / n) * abs(acc - conf)
    return float(ece)


def enrich_row_with_null(
    row: Dict[str, Any],
    y_test: np.ndarray,
    proba: np.ndarray,
    *,
    auc_key: str = "hgb_auc",
    n_perm: int = 5,
    seed: int = 0,
    selected_k: int = 1,
    proba_lr: Optional[np.ndarray] = None,
    n_ece_bins: int = 10,
) -> Dict[str, Any]:
    """Attach null / excess / probe_eff / ECE fields to one adjacent-pair row."""
    out = dict(row)
    null = permute_auc(y_test, proba, n_perm=n_perm, seed=seed)
    out.update(null)
    auc = float(out.get(auc_key, float("nan")))
    ex = excess_auc(auc, null["null_auc_mean"])
    flops = relative_probe_flops(int(out.get("n0", 1)), selected_k)
    out["excess_auc"] = ex
    out["relative_flops"] = flops
    out["probe_eff"] = probe_efficiency(ex, flops)
    out["hgb_ece"] = expected_calibration_error(y_test, proba, n_bins=n_ece_bins)
    if proba_lr is not None:
        out["logreg_ece"] = expected_calibration_error(
            y_test, proba_lr, n_bins=n_ece_bins
        )
    else:
        out["logreg_ece"] = float("nan")
    return out


def summarize_null_pack(rows: Sequence[Dict[str, Any]]) -> Dict[str, Optional[float]]:
    def _mean(key: str) -> Optional[float]:
        vals = [
            float(r[key])
            for r in rows
            if r.get(key) is not None and np.isfinite(float(r[key]))
        ]
        return float(np.mean(vals)) if vals else None

    return {
        "mean_null_auc": _mean("null_auc_mean"),
        "mean_excess_auc": _mean("excess_auc"),
        "mean_probe_eff": _mean("probe_eff"),
        "mean_relative_flops": _mean("relative_flops"),
        "mean_hgb_ece": _mean("hgb_ece"),
        "mean_logreg_ece": _mean("logreg_ece"),
    }
