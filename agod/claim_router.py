"""Global claim router: check overclaims and fall back to weaker statements.

Motivation
----------
Scorecards / RSI ticks often *overclaim*: \"skill\", \"drivers\", \"win\",
\"alarm ⇒ act\", \"stable\".  A single gate downgrades text to what the
evidence actually supports, and never lets L4 (causal / ship / ROI) through
from diagnostic boards alone.

Levels (honest ladder; aligned with coverage/confidence docs)
-------------------------------------------------------------
- L0 descriptive — raw numbers only
- L1 diagnostic  — compare / residual / curve shape (no \"skill\" word)
- L2 evidential  — excess / null / burn codes with thresholds met
- L3 action hint — soft policy suggestion (SOFTEN / investigate), not ship
- L4 forbidden   — causal, production ship, business ROI — always fallback

Usage
-----
``route_claim(proposed, evidence)`` → ``{level, text, flags, fallback}``
``harden_reading(text, evidence)``  → safe text for MD / JSON headlines.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import IntEnum
from typing import Any, Dict, List, Mapping, Optional, Sequence


class ClaimLevel(IntEnum):
    L0_DESCRIPTIVE = 0
    L1_DIAGNOSTIC = 1
    L2_EVIDENTIAL = 2
    L3_ACTION_HINT = 3
    L4_FORBIDDEN = 4


# Phrases that force a downgrade unless evidence unlocks them.
_OVERCLAIM_PATTERNS: Sequence[tuple[str, ClaimLevel, str]] = (
    (r"\bpersistent drivers?\b", ClaimLevel.L4_FORBIDDEN, "causal_driver"),
    (r"\btrue drivers?\b", ClaimLevel.L4_FORBIDDEN, "causal_driver"),
    (r"\bcausal\b", ClaimLevel.L4_FORBIDDEN, "causal"),
    (r"\bproduction\b", ClaimLevel.L4_FORBIDDEN, "ship_production"),
    (r"\bpromote\b", ClaimLevel.L4_FORBIDDEN, "ship_production"),
    (r"\bROI\b", ClaimLevel.L4_FORBIDDEN, "business_roi"),
    (r"\b省人分钟\b", ClaimLevel.L4_FORBIDDEN, "business_roi"),
    (r"\b transferable skill\b", ClaimLevel.L2_EVIDENTIAL, "skill_word"),
    (r"\bskill\b", ClaimLevel.L2_EVIDENTIAL, "skill_word"),
    (r"\bbeats?\b", ClaimLevel.L2_EVIDENTIAL, "beats_word"),
    (r"\bwins?\b", ClaimLevel.L2_EVIDENTIAL, "wins_word"),
    (r"\bstable\b", ClaimLevel.L1_DIAGNOSTIC, "stable_word"),
    (r"\balarm\b", ClaimLevel.L2_EVIDENTIAL, "alarm_word"),
)


@dataclass
class Evidence:
    """Sparse evidence bag — missing fields ⇒ cannot unlock higher claims."""

    excess_auc: Optional[float] = None
    excess_partial: Optional[float] = None
    delta_excess: Optional[float] = None
    excess_ci_lo: Optional[float] = None  # if set, must be >0 to claim skill
    rel_vs_uniform: Optional[float] = None  # <1 ⇒ gated improves MSE
    burn: Optional[bool] = None
    burn_decision: Optional[str] = None
    n_points: Optional[int] = None
    ph_alarm: Optional[bool] = None
    ph_series_len: Optional[int] = None
    mean_auc: Optional[float] = None
    mean_jaccard: Optional[float] = None
    allow_action_hint: bool = False
    extra: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        d = {
            "excess_auc": self.excess_auc,
            "excess_partial": self.excess_partial,
            "delta_excess": self.delta_excess,
            "excess_ci_lo": self.excess_ci_lo,
            "rel_vs_uniform": self.rel_vs_uniform,
            "burn": self.burn,
            "burn_decision": self.burn_decision,
            "n_points": self.n_points,
            "ph_alarm": self.ph_alarm,
            "ph_series_len": self.ph_series_len,
            "mean_auc": self.mean_auc,
            "mean_jaccard": self.mean_jaccard,
            "allow_action_hint": self.allow_action_hint,
        }
        d.update(self.extra)
        return d


def _finite(x: Optional[float]) -> bool:
    return x is not None and x == x and abs(x) != float("inf")


def scan_overclaim_text(text: str) -> List[Dict[str, str]]:
    """Return matched overclaim tags in ``text`` (order preserved)."""
    flags = []
    for pat, level, tag in _OVERCLAIM_PATTERNS:
        if re.search(pat, text, flags=re.IGNORECASE):
            flags.append({"tag": tag, "need_level": level.name, "pattern": pat})
    return flags


def max_supported_level(ev: Evidence) -> ClaimLevel:
    """Highest level the evidence bag can unlock."""
    # L2 skill/evidential unlock
    skill_ok = False
    if _finite(ev.excess_auc) and float(ev.excess_auc) >= 0.05:
        if ev.excess_ci_lo is None or (
            _finite(ev.excess_ci_lo) and float(ev.excess_ci_lo) > 0.0
        ):
            # partial check if provided: partial must not collapse to ~0
            if ev.excess_partial is None:
                skill_ok = True
            elif _finite(ev.excess_partial) and float(ev.excess_partial) >= 0.03:
                skill_ok = True

    burn_ok = bool(ev.burn) and _finite(ev.rel_vs_uniform) and float(ev.rel_vs_uniform) < 1.0
    ph_ok = bool(ev.ph_alarm) and (ev.ph_series_len or 0) >= 6
    stable_ok = (ev.n_points or 0) >= 4

    # L3 only with explicit allow + a real policy code
    action_ok = bool(ev.allow_action_hint) and bool(
        ev.burn_decision
        in ("SOFTEN_ONLY", "KEEP_UNIFORM", "BURN_SQRT", "BURN_CBRT")
    )

    if action_ok and (burn_ok or ev.burn_decision in ("SOFTEN_ONLY", "KEEP_UNIFORM")):
        return ClaimLevel.L3_ACTION_HINT
    if skill_ok or burn_ok or ph_ok:
        return ClaimLevel.L2_EVIDENTIAL
    if stable_ok or _finite(ev.mean_auc) or _finite(ev.delta_excess):
        return ClaimLevel.L1_DIAGNOSTIC
    return ClaimLevel.L0_DESCRIPTIVE


def _fallback_text(level: ClaimLevel, ev: Evidence, flags: Sequence[Mapping[str, str]]) -> str:
    tags = {f["tag"] for f in flags}
    if ClaimLevel.L4_FORBIDDEN.name in {f["need_level"] for f in flags} or tags & {
        "causal_driver",
        "causal",
        "ship_production",
        "business_roi",
    }:
        return (
            "diagnostic only — board does not support causal / ship / ROI claims; "
            "fallback to feature-family investigation"
        )
    if "skill_word" in tags and level < ClaimLevel.L2_EVIDENTIAL:
        if _finite(ev.delta_excess) and float(ev.delta_excess) > 0.05:
            return (
                "excess shrinks after partialling confounder — do not claim skill; "
                "report Δexcess only"
            )
        if _finite(ev.excess_auc) and float(ev.excess_auc) < 0.05:
            return "excess near null — little transferable association beyond chance"
        return "insufficient excess / CI / partial evidence — avoid 'skill' wording"
    if "beats_word" in tags or "wins_word" in tags:
        if not (_finite(ev.rel_vs_uniform) and float(ev.rel_vs_uniform) < 1.0):
            return "no MSE win vs uniform — cannot claim beats/wins; keep SOFTEN/KEEP"
    if "alarm_word" in tags and (ev.ph_series_len or 0) < 6:
        return "series too short for PH alarm — report statistic only"
    if "stable_word" in tags and (ev.n_points or 0) < 4:
        return "too few points to call stable — report curve values only"
    if level <= ClaimLevel.L0_DESCRIPTIVE:
        return "descriptive metrics only — no interpretive claim unlocked"
    return "claim downgraded — evidence below required level"


def route_claim(
    proposed: str,
    evidence: Optional[Evidence] = None,
    *,
    requested_level: Optional[ClaimLevel] = None,
) -> Dict[str, Any]:
    """Route a proposed claim through overclaim checks + fallback.

    Returns
    -------
    dict with keys:
      ok, level, text, flags, fallback_applied, max_supported, evidence
    """
    ev = evidence or Evidence()
    text = (proposed or "").strip()
    flags = scan_overclaim_text(text)
    supported = max_supported_level(ev)

    # implied need = max of pattern needs; default L1 if plain prose
    need = ClaimLevel.L1_DIAGNOSTIC
    for f in flags:
        need = max(need, ClaimLevel[f["need_level"]])
    if requested_level is not None:
        need = max(need, requested_level)

    # L4 always blocked
    if need >= ClaimLevel.L4_FORBIDDEN or any(
        f["need_level"] == ClaimLevel.L4_FORBIDDEN.name for f in flags
    ):
        fb = _fallback_text(ClaimLevel.L0_DESCRIPTIVE, ev, flags)
        return {
            "ok": False,
            "level": ClaimLevel.L0_DESCRIPTIVE.name,
            "text": fb,
            "proposed": text,
            "flags": flags,
            "fallback_applied": True,
            "max_supported": supported.name,
            "evidence": ev.to_dict(),
        }

    if need <= supported and text:
        return {
            "ok": True,
            "level": need.name,
            "text": text,
            "proposed": text,
            "flags": flags,
            "fallback_applied": False,
            "max_supported": supported.name,
            "evidence": ev.to_dict(),
        }

    fb = _fallback_text(supported, ev, flags)
    return {
        "ok": False,
        "level": supported.name,
        "text": fb,
        "proposed": text,
        "flags": flags,
        "fallback_applied": True,
        "max_supported": supported.name,
        "evidence": ev.to_dict(),
    }


def harden_reading(
    proposed: str,
    evidence: Optional[Evidence] = None,
    *,
    requested_level: Optional[ClaimLevel] = None,
) -> str:
    """Convenience: return routed safe text only."""
    return str(route_claim(proposed, evidence, requested_level=requested_level)["text"])


def route_pack_interpretation(
    *,
    mean_auc: Optional[float],
    mean_jaccard: Optional[float] = None,
    mean_excess: Optional[float] = None,
    excess_ci_lo: Optional[float] = None,
    excess_partial: Optional[float] = None,
) -> Dict[str, Any]:
    """Global router for adjacent-board pack readings (replaces loose wording)."""
    ev = Evidence(
        mean_auc=mean_auc,
        mean_jaccard=mean_jaccard,
        excess_auc=mean_excess,
        excess_ci_lo=excess_ci_lo,
        excess_partial=excess_partial,
        n_points=1,
    )
    if mean_auc is None:
        return route_claim("no pairs", ev, requested_level=ClaimLevel.L0_DESCRIPTIVE)

    # Build a *conservative* proposed reading, then harden.
    if mean_excess is not None and mean_excess == mean_excess and mean_excess < 0.05:
        proposed = "AUC near null (little transferable association beyond chance)"
    elif mean_auc < 0.65:
        proposed = "weak transfer (association does not travel)"
    elif mean_jaccard is not None and mean_jaccard >= 0.5:
        # was "persistent drivers" — overclaim; downgrade wording up front
        proposed = (
            "strong transfer + stable top feats "
            "(persistent correlates — not identified as drivers)"
        )
    elif mean_jaccard is not None and mean_jaccard < 0.35:
        proposed = "strong transfer but shifting top feats (regime / composition change)"
    else:
        proposed = "transfer holds; feature set partially stable"

    # \"stable\" / \"skill\" checks via router
    return route_claim(proposed, ev, requested_level=ClaimLevel.L1_DIAGNOSTIC)


def route_ml_suite_claims(suite: Mapping[str, Any]) -> Dict[str, Any]:
    """Run claim router on ml_method_suite readings / headlines."""
    partial = suite.get("partial_excess") or {}
    curve = suite.get("learning_curve") or {}
    ph = suite.get("page_hinkley") or {}
    ev = Evidence(
        excess_auc=suite.get("excess_auc"),
        excess_partial=partial.get("excess_partial"),
        delta_excess=partial.get("delta_excess"),
        n_points=len(curve.get("points") or []),
        ph_alarm=ph.get("alarm"),
        ph_series_len=ph.get("n"),
    )
    out = {}
    for key, proposed in (
        ("curve", curve.get("reading") or ""),
        ("page_hinkley", ph.get("reading") or ""),
        (
            "partial",
            (
                f"delta_excess={partial.get('delta_excess')} "
                f"(raw={partial.get('excess_raw')}, "
                f"partial={partial.get('excess_partial')})"
                if partial
                else ""
            ),
        ),
    ):
        if not proposed:
            continue
        out[key] = route_claim(str(proposed), ev)
    # headline: never say \"false skill\" as proven without Δexcess gate
    delta = partial.get("delta_excess")
    if delta is not None and delta == delta and delta > 0.05:
        head = (
            "confounder partialling: large Δexcess — report association shrink; "
            "do not claim proven skill"
        )
    elif _finite(suite.get("excess_auc")) and float(suite["excess_auc"]) >= 0.05:
        head = harden_reading(
            "excess above null under current evidence — evidential only, not a driver ID",
            ev,
            requested_level=ClaimLevel.L2_EVIDENTIAL,
        )
    else:
        head = "descriptive suite only — no evidential skill claim unlocked"
    out["headline"] = route_claim(head, ev, requested_level=ClaimLevel.L1_DIAGNOSTIC)
    return out
