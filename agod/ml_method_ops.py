"""ML method operators that *compose* with the existing null / excess stack.

Does **not** rewrite ``transfer_null`` / ``po_eff`` / ``soft_burn``.
Adds three orthogonal ML diagnostics used by the 20-min RSI loop:

1. **Partial excess** — residualize scores vs a confounder (e.g. volume),
   then re-measure excess.  Answers: is transfer skill just the confounder?
2. **Excess learning curve** — excess vs n_train subsample.  Sample efficiency.
3. **Page–Hinkley on excess stream** — detect when *skill itself* drifts.

All quantities sit next to ``excess_auc`` / ECE / probe_eff — same H0 language.
"""
from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from sklearn.metrics import roc_auc_score

from agod.transfer_null import excess_auc, permute_auc


def residualize_vs_confounder(
    score: np.ndarray,
    confounder: np.ndarray,
) -> Dict[str, Any]:
    """Linear residualize ``score ~ 1 + confounder``; return residual scores.

    Classic partialling-out: if ranking skill collapses after removing the
    linear footprint of a volume / intensity feature, the transfer story is
    confounder-driven (still may be useful ops signal — but not \"new skill\").
    """
    s = np.asarray(score, dtype=float).ravel()
    c = np.asarray(confounder, dtype=float).ravel()
    n = min(s.size, c.size)
    s, c = s[:n], c[:n]
    if n < 3 or float(np.std(c)) < 1e-12:
        return {
            "residual": s.copy(),
            "beta": float("nan"),
            "r2": float("nan"),
            "n": float(n),
        }
    # OLS via normal equations on [1, c]
    X = np.column_stack([np.ones(n), c])
    beta, *_ = np.linalg.lstsq(X, s, rcond=None)
    fitted = X @ beta
    resid = s - fitted
    ss_tot = float(np.sum((s - s.mean()) ** 2))
    ss_res = float(np.sum(resid**2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-18 else float("nan")
    return {
        "residual": resid,
        "beta": float(beta[1]),
        "intercept": float(beta[0]),
        "r2": float(r2),
        "n": float(n),
    }


def partial_excess_auc(
    y_true: np.ndarray,
    score: np.ndarray,
    confounder: np.ndarray,
    *,
    n_perm: int = 8,
    seed: int = 0,
) -> Dict[str, float]:
    """excess_auc before vs after residualizing ``score`` vs confounder.

    Uses the same label-permutation null as ``transfer_null.permute_auc``
    — only the *score vector* changes.
    """
    y = np.asarray(y_true).astype(int).ravel()
    s = np.asarray(score, dtype=float).ravel()
    c = np.asarray(confounder, dtype=float).ravel()
    n = min(y.size, s.size, c.size)
    y, s, c = y[:n], s[:n], c[:n]
    out: Dict[str, float] = {
        "auc_raw": float("nan"),
        "excess_raw": float("nan"),
        "auc_partial": float("nan"),
        "excess_partial": float("nan"),
        "delta_excess": float("nan"),
        "confounder_r2": float("nan"),
        "n": float(n),
    }
    if n < 4 or len(np.unique(y)) < 2:
        return out

    def _pack(sc: np.ndarray) -> Tuple[float, float]:
        try:
            auc = float(roc_auc_score(y, sc))
        except Exception:
            return float("nan"), float("nan")
        null = permute_auc(y, sc, n_perm=n_perm, seed=seed)
        ex = excess_auc(auc, null["null_auc_mean"])
        return auc, ex

    auc_r, ex_r = _pack(s)
    part = residualize_vs_confounder(s, c)
    auc_p, ex_p = _pack(part["residual"])
    out.update(
        {
            "auc_raw": auc_r,
            "excess_raw": ex_r,
            "auc_partial": auc_p,
            "excess_partial": ex_p,
            "delta_excess": (
                float(ex_r - ex_p)
                if np.isfinite(ex_r) and np.isfinite(ex_p)
                else float("nan")
            ),
            "confounder_r2": float(part.get("r2", float("nan"))),
            "beta": float(part.get("beta", float("nan"))),
        }
    )
    return out


def excess_learning_curve(
    y_true: np.ndarray,
    score: np.ndarray,
    *,
    fractions: Sequence[float] = (0.2, 0.4, 0.6, 0.8, 1.0),
    n_perm: int = 5,
    seed: int = 0,
) -> Dict[str, Any]:
    """excess_auc on nested prefixes of the test chunk (proxy learning curve).

    Not a full re-fit curve (that needs the probe trainer); this answers:
    does *observed* excess stabilize with more evaluated rows — i.e. is the
    skill estimate sample-hungry?
    """
    y = np.asarray(y_true).astype(int).ravel()
    s = np.asarray(score, dtype=float).ravel()
    n = min(y.size, s.size)
    y, s = y[:n], s[:n]
    rows = []
    for f in fractions:
        m = max(4, int(np.floor(float(f) * n)))
        m = min(m, n)
        yy, ss = y[:m], s[:m]
        if len(np.unique(yy)) < 2:
            rows.append({"frac": float(f), "n": m, "excess": float("nan"), "auc": float("nan")})
            continue
        try:
            auc = float(roc_auc_score(yy, ss))
        except Exception:
            rows.append({"frac": float(f), "n": m, "excess": float("nan"), "auc": float("nan")})
            continue
        null = permute_auc(yy, ss, n_perm=n_perm, seed=seed)
        rows.append(
            {
                "frac": float(f),
                "n": m,
                "auc": auc,
                "excess": excess_auc(auc, null["null_auc_mean"]),
            }
        )
    finite = [r["excess"] for r in rows if np.isfinite(r.get("excess", np.nan))]
    return {
        "points": rows,
        "excess_at_full": rows[-1]["excess"] if rows else float("nan"),
        "excess_at_20pct": rows[0]["excess"] if rows else float("nan"),
        "curve_gain_20_to_full": (
            float(finite[-1] - finite[0]) if len(finite) >= 2 else float("nan")
        ),
        "reading": _curve_reading(rows),
    }


def _curve_reading(rows: Sequence[Mapping]) -> str:
    ex = [r.get("excess") for r in rows if np.isfinite(r.get("excess", np.nan))]
    if len(ex) < 2:
        return "insufficient points"
    gain = float(ex[-1] - ex[0])
    if abs(gain) < 0.02:
        return "excess stable early — estimate not sample-hungry"
    if gain > 0.05:
        return "excess rises with n — need more eval mass before claiming skill"
    return "excess softens / noisy with n — check rare labels / confounders"


def page_hinkley_skill(
    excess_series: Sequence[float],
    *,
    delta: float = 0.005,
    lambda_thresh: float = 0.05,
) -> Dict[str, Any]:
    """Page–Hinkley change detector on a stream of ``excess_auc`` values.

    Detects sustained *drops* in transfer skill by running classic PH on
    ``z = -excess`` (increase in skill-loss).  Meta-monitors the series —
    does not retune the probe.
    """
    xs = [float(v) for v in excess_series if np.isfinite(v)]
    if len(xs) < 3:
        return {
            "alarm": False,
            "alarm_index": None,
            "ph_stat": [],
            "n": len(xs),
            "reading": "series too short",
        }
    # PH on skill-loss z = -excess (detect upward mean shift in z)
    zs = [-x for x in xs]
    cum = 0.0
    min_cum = 0.0
    ph_stat: List[float] = []
    alarm_i = None
    running: List[float] = []
    for t, z in enumerate(zs):
        running.append(z)
        mean = float(np.mean(running))
        cum = cum + (z - mean - delta)
        min_cum = min(min_cum, cum)
        stat = float(cum - min_cum)
        ph_stat.append(stat)
        if alarm_i is None and stat > lambda_thresh:
            alarm_i = t
    return {
        "alarm": alarm_i is not None,
        "alarm_index": alarm_i,
        "ph_stat": ph_stat,
        "lambda": float(lambda_thresh),
        "delta": float(delta),
        "n": len(xs),
        "mean_excess": float(np.mean(xs)),
        "reading": (
            f"skill-drop alarm at pair index {alarm_i}"
            if alarm_i is not None
            else "no sustained excess drop detected"
        ),
    }


def ml_method_suite(
    y_true: np.ndarray,
    score: np.ndarray,
    *,
    confounder: Optional[np.ndarray] = None,
    excess_history: Optional[Sequence[float]] = None,
    n_perm: int = 6,
    seed: int = 0,
) -> Dict[str, Any]:
    """One-shot ML diagnostics beside a transfer probe (raw score fixed)."""
    y = np.asarray(y_true).astype(int).ravel()
    s = np.asarray(score, dtype=float).ravel()
    n = min(y.size, s.size)
    y, s = y[:n], s[:n]
    try:
        auc = float(roc_auc_score(y, s)) if len(np.unique(y)) > 1 else float("nan")
    except Exception:
        auc = float("nan")
    null = permute_auc(y, s, n_perm=n_perm, seed=seed)
    ex = excess_auc(auc, null["null_auc_mean"])
    curve = excess_learning_curve(y, s, n_perm=max(3, n_perm // 2), seed=seed)
    partial = None
    if confounder is not None:
        partial = partial_excess_auc(
            y, s, confounder, n_perm=n_perm, seed=seed
        )
    hist = list(excess_history or [])
    if np.isfinite(ex):
        hist = hist + [ex]
    ph = page_hinkley_skill(hist) if len(hist) >= 3 else {
        "alarm": False,
        "reading": "need ≥3 excess points for PH",
        "n": len(hist),
    }
    return {
        "auc_obs": auc,
        "excess_auc": ex,
        "null_auc_mean": null["null_auc_mean"],
        "learning_curve": curve,
        "partial_excess": partial,
        "page_hinkley": ph,
        "method_note": (
            "composes with transfer_null; does not modify null / PO / soft_burn cores"
        ),
    }
