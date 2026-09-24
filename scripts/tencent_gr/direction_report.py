"""Signed direction payload for localize / FSDS JSON summaries.

Attaches outcome Δȳ and tip-feature sign(δ_j) without changing Drill gates.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

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


def _tips_from_blob(blob: Mapping[str, Any]) -> List[str]:
    for key in ("fsds_W1_holdout", "fsds_W1", "fsds_W2_temporal", "fsds_W2"):
        block = blob.get(key) or {}
        tops = block.get("top_features")
        if isinstance(tops, list) and tops:
            return [str(x) for x in tops]
    tip_signs = (blob.get("direction") or {}).get("tip_signs") or {}
    if tip_signs:
        return [str(x) for x in tip_signs.keys()]
    return []


def ensure_direction(
    blob: Mapping[str, Any],
    *,
    feat_diag: Optional[pd.DataFrame] = None,
    feat_diag_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """Attach / repair direction on a summary blob without re-running Drill.

    - tip_signs: from feature_shift_diagnostics when present
    - Dy: keep existing edge-level Dy if already set; else leave missing→flat
      (do **not** invent Dy from window-wide pos_rate — that overclaims S1)
    """
    existing = dict(blob.get("direction") or {})
    tips = _tips_from_blob(blob)
    diag = feat_diag
    if diag is None and feat_diag_path is not None and Path(feat_diag_path).exists():
        diag = pd.read_csv(feat_diag_path)

    tip_part = tip_direction(
        diag if diag is not None else pd.DataFrame(),
        tips,
    )

    has_dy = existing.get("Dy") is not None and existing.get("sign_Dy") in (
        "pos",
        "neg",
        "flat",
    )
    tip_signs_existing = existing.get("tip_signs") or {}
    tips_ok = bool(tip_signs_existing) and any(
        str(v) in ("+", "-") for v in tip_signs_existing.values()
    )

    if has_dy and tips_ok:
        out = {**existing}
        out.setdefault("report", f"[Direction] Dy={out.get('Dy')} ({out.get('sign_Dy')})")
        return out

    if has_dy:
        out = {
            **existing,
            **tip_part,
        }
        out["report"] = (
            f"[Direction] Dy={out.get('Dy')} ({out.get('sign_Dy')}); "
            f"tip_signs={out.get('tip_signs')}"
        )
        return out

    # No edge-level Dy available: honest flat + dy_missing; still fill tip signs.
    out = {
        "y_col": existing.get("y_col") or "y_convert",
        "y_ref": None,
        "y_cur": None,
        "Dy": None,
        "sign_Dy": "flat",
        "dy_missing": True,
        "dy_source": "unavailable_no_edge_y",
        "n_ref": existing.get("n_ref"),
        "n_cur": existing.get("n_cur"),
        **tip_part,
    }
    out["report"] = (
        f"[Direction] Dy=None (flat, dy_missing); tip_signs={out.get('tip_signs')}"
    )
    if existing.get("extra"):
        out["extra"] = existing["extra"]
    return out


def scenario_from_direction(direction: Mapping[str, Any]) -> Dict[str, Any]:
    """Map sign_Dy + tip_signs → S1/S2/S3 family (display clue, not auto-ban)."""
    sign = str(direction.get("sign_Dy") or "flat")
    tips = dict(direction.get("tip_signs") or {})
    dy_missing = bool(direction.get("dy_missing")) or direction.get("Dy") is None

    modifiers: List[str] = []
    if any(k.startswith("i_n_covisit") and tips.get(k) in ("+", "-") for k in tips):
        modifiers.append("S4_covisit")
    if tips.get("ui_pop_mismatch") in ("+", "-"):
        modifiers.append("S5_pop_mismatch")

    if sign == "pos":
        family, code = "刷量族", "S1_brush"
        if tips.get("i_share_last") == "+" or tips.get("i_credit_last") == "+":
            sub, sub_code = "末跳操控", "last_hop"
            read = "成功率升 + 末跳 tip+：末跳/刷量队（仍 L1，辨爆款）"
        elif tips.get("u_span_sec") == "-":
            sub, sub_code = "短窗狂点", "short_span"
            read = "成功率升 + 跨度偏短：短刷候选（仍 L1）"
        else:
            sub, sub_code = "弱刷量或真爆款", "weak_or_viral"
            read = "成功率升但无末跳 tip：弱刷量或真爆款 → 必须人工辨"
    elif sign == "neg":
        family, code = "灌入族", "S2_inject"
        if tips.get("u_span_sec") == "-":
            sub, sub_code = "劣质短会话", "short_bad_session"
            read = "成功率掉 + 跨度变短：劣质短会话灌入队"
        else:
            sub, sub_code = "差流/劫持残留", "junk_or_hijack"
            read = "成功率掉：灌入/劫持残留队列（慎当商户变差）"
    else:
        family, code = "漂移族", "S3_drift"
        sub, sub_code = "结构漂移", "structure_shift"
        if dy_missing:
            read = "Dy 缺失按 flat：X 结构漂了，成功率未证实联动 → 漂移族，禁自称刷量"
        else:
            read = "成功率平 + tip 漂：供给/推荐/活动结构漂移，慎升强动作"

    return {
        "family": family,
        "family_code": code,
        "sub": sub,
        "sub_code": sub_code,
        "modifiers": modifiers,
        "read": read,
        "auto_ban": False,
        "action_level_default": "L1_watch",
        "dy_missing": dy_missing,
    }
