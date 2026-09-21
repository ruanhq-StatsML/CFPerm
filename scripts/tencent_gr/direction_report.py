"""Signed direction payload for localize / FSDS JSON summaries.

Attaches outcome Δȳ and tip-feature sign(δ_j) without changing Drill gates.
"""
from __future__ import annotations

from typing import Any, Dict, Mapping, Optional, Sequence

import numpy as np
import pandas as pd


def _sign_label(x: float, *, eps: float = 1e-12) -> str:
    if not np.isfinite(x) or abs(x) <= eps:
        return "flat"
    return "pos" if x > 0 else "neg"


def _pm(x: float, *, eps: float = 1e-12) -> str:
    if not np.isfinite(x) or abs(x) <= eps:
        return "0"
    return "+" if x > 0 else "-"


def outcome_direction(
    g_ref: pd.DataFrame,
    g_cur: pd.DataFrame,
    *,
    y_col: str = "y_convert",
) -> Dict[str, Any]:
    """Δȳ = mean(y|cur) − mean(y|ref) on the localized support."""
    if y_col not in g_ref.columns or y_col not in g_cur.columns:
        return {
            "y_col": y_col,
            "y_ref": None,
            "y_cur": None,
            "Dy": None,
            "sign_Dy": "flat",
            "n_ref": int(len(g_ref)),
            "n_cur": int(len(g_cur)),
        }
    y0 = pd.to_numeric(g_ref[y_col], errors="coerce").dropna().to_numpy(dtype=float)
    y1 = pd.to_numeric(g_cur[y_col], errors="coerce").dropna().to_numpy(dtype=float)
    m0 = float(y0.mean()) if len(y0) else float("nan")
    m1 = float(y1.mean()) if len(y1) else float("nan")
    dy = float(m1 - m0) if np.isfinite(m0) and np.isfinite(m1) else float("nan")
    return {
        "y_col": y_col,
        "y_ref": m0 if np.isfinite(m0) else None,
        "y_cur": m1 if np.isfinite(m1) else None,
        "Dy": dy if np.isfinite(dy) else None,
        "sign_Dy": _sign_label(dy) if np.isfinite(dy) else "flat",
        "n_ref": int(len(y0)),
        "n_cur": int(len(y1)),
    }


def tip_direction(
    feat_diag: pd.DataFrame,
    tip_features: Sequence[str],
    *,
    mean_ref_col: str = "mean_W1",
    mean_cur_col: str = "mean_W2",
) -> Dict[str, Any]:
    """Per-tip sign(μ_cur − μ_ref) from feature shift diagnostics."""
    tip_signs: Dict[str, str] = {}
    tip_delta: Dict[str, Optional[float]] = {}
    if feat_diag is None or len(feat_diag) == 0 or "feature" not in feat_diag.columns:
        for j in tip_features:
            tip_signs[str(j)] = "0"
            tip_delta[str(j)] = None
        return {"tip_signs": tip_signs, "tip_delta": tip_delta}

    idx = feat_diag.set_index("feature")
    for j in tip_features:
        key = str(j)
        if key not in idx.index or mean_ref_col not in idx.columns or mean_cur_col not in idx.columns:
            tip_signs[key] = "0"
            tip_delta[key] = None
            continue
        d = float(idx.loc[key, mean_cur_col] - idx.loc[key, mean_ref_col])
        tip_delta[key] = d
        tip_signs[key] = _pm(d)
    return {"tip_signs": tip_signs, "tip_delta": tip_delta}


def build_direction_dict(
    g_ref: pd.DataFrame,
    g_cur: pd.DataFrame,
    tip_features: Sequence[str],
    *,
    feat_diag: Optional[pd.DataFrame] = None,
    y_col: str = "y_convert",
    extra: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    """JSON-ready direction block for summary.json."""
    out_dir = outcome_direction(g_ref, g_cur, y_col=y_col)
    tips = tip_direction(feat_diag if feat_diag is not None else pd.DataFrame(), tip_features)
    report = (
        f"[Direction] Dy={out_dir.get('Dy')} ({out_dir.get('sign_Dy')}); "
        f"tip_signs={tips['tip_signs']}"
    )
    blob: Dict[str, Any] = {
        **out_dir,
        **tips,
        "report": report,
    }
    if extra:
        blob["extra"] = dict(extra)
    return blob
