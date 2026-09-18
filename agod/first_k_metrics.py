"""Manuscript OnlinePermOOB first1/first2/first3 + alarm-time quantiles.

Matches ``onlinePermOOB`` counting (see ``docs/biz/FIRST_K_LOGIC.py`` / OOB algo):

  first_k = first index where k consecutive rejects end  (1-based end index)
  mean / median / q25 / q75 = distribution of *all* alarm times on the trail
  DetRate (day streams) = fraction of units with ≥1 alarm
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Union

import numpy as np

ArrayLike = Union[Sequence[int], Sequence[bool], np.ndarray]


def first_k_consecutive_rej(rej: ArrayLike, k: int) -> float:
    """1-based end index of the first run of ``k`` consecutive True; else NaN."""
    det = np.asarray(rej, dtype=bool).ravel()
    n = int(det.size)
    if k <= 0 or n < k:
        return float("nan")
    for i in range(n - k + 1):
        if bool(det[i : i + k].all()):
            return float(i + k)
    return float("nan")


def alarm_time_quantiles(rej: ArrayLike) -> Dict[str, float]:
    """Distribution of all alarm times (0-based trail indices where reject=1)."""
    det = np.asarray(rej, dtype=bool).ravel()
    times = np.flatnonzero(det).astype(float)
    if times.size == 0:
        return {
            "n_alarm": 0,
            "mean": float("nan"),
            "std": float("nan"),
            "p25": float("nan"),
            "median": float("nan"),
            "p75": float("nan"),
            "min": float("nan"),
            "max": float("nan"),
        }
    qs = np.quantile(times, [0.25, 0.50, 0.75])
    return {
        "n_alarm": int(times.size),
        "mean": float(times.mean()),
        "std": float(times.std(ddof=0)),
        "p25": float(qs[0]),
        "median": float(qs[1]),
        "p75": float(qs[2]),
        "min": float(times.min()),
        "max": float(times.max()),
    }


def first_k_form(rej: ArrayLike) -> Dict[str, float]:
    """One-stream manuscript row: first1/2/3 + alarm-time quantiles + SUM."""
    det = np.asarray(rej, dtype=bool).ravel()
    q = alarm_time_quantiles(det)
    return {
        "SUM": int(det.sum()),
        "first1": first_k_consecutive_rej(det, 1),
        "first2": first_k_consecutive_rej(det, 2),
        "first3": first_k_consecutive_rej(det, 3),
        "alarm_mean": q["mean"],
        "alarm_std": q["std"],
        "alarm_p25": q["p25"],
        "alarm_median": q["median"],
        "alarm_p75": q["p75"],
        "n_alarm": q["n_alarm"],
        "n_trail": int(det.size),
    }


def aggregate_day_forms(forms: List[Dict[str, float]]) -> Dict[str, float]:
    """Aggregate per-day first_k forms → DetRate + first1 quantiles on alarmed days."""
    n_days = len(forms)
    if n_days == 0:
        return {
            "n_days": 0,
            "n_alarm_days": 0,
            "det_rate": float("nan"),
            "first1_median": float("nan"),
            "first2_median": float("nan"),
            "first3_median": float("nan"),
            "first1_mean": float("nan"),
            "first1_std": float("nan"),
            "first1_p25": float("nan"),
            "first1_p75": float("nan"),
            "mean_sum": float("nan"),
        }

    def _col(key: str) -> np.ndarray:
        return np.asarray([f.get(key, float("nan")) for f in forms], dtype=float)

    f1 = _col("first1")
    f2 = _col("first2")
    f3 = _col("first3")
    alarmed = np.isfinite(f1)
    n_alarm_days = int(alarmed.sum())
    f1_a = f1[alarmed]

    def _med(a: np.ndarray) -> float:
        a = a[np.isfinite(a)]
        return float(np.median(a)) if a.size else float("nan")

    out = {
        "n_days": n_days,
        "n_alarm_days": n_alarm_days,
        "det_rate": float(n_alarm_days / n_days) if n_days else float("nan"),
        "first1_median": _med(f1),
        "first2_median": _med(f2),
        "first3_median": _med(f3),
        "mean_sum": float(np.nanmean(_col("SUM"))),
    }
    if f1_a.size:
        qs = np.quantile(f1_a, [0.25, 0.75])
        out.update(
            {
                "first1_mean": float(f1_a.mean()),
                "first1_std": float(f1_a.std(ddof=0)),
                "first1_p25": float(qs[0]),
                "first1_p75": float(qs[1]),
            }
        )
    else:
        out.update(
            {
                "first1_mean": float("nan"),
                "first1_std": float("nan"),
                "first1_p25": float("nan"),
                "first1_p75": float("nan"),
            }
        )
    return out


def fmt_num(x: Optional[float], nd: int = 1) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "---"
    if abs(float(x) - round(float(x))) < 1e-9:
        return str(int(round(float(x))))
    return f"{float(x):.{nd}f}"
