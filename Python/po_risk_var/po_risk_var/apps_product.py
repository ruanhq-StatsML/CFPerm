"""Five product applications of PO-risk Mag×Var + Var(u) + Y-flip.

1. Dynamic routing & cost optimization
2. Granular canary release & segment rollback
3. Data flywheel / active-learning curation
4. Explainable compliance guardrail reports
5. Agentic workflow cascade-error blocking

All share the same sensors: severity = J + Y-flip, CF-debiased Var(u),
triage(). ROI figure $22,430 is from industrial LIVE boards
(results/as_dta_business_pr).
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

from po_risk_var._path import REPO_ROOT
from po_risk_var.product import TriageDecision, triage

OUT_DEFAULT = REPO_ROOT / "results" / "as_dta_po_risk_var_product_apps"
COST_SAVED_USD = 22430.0  # industrial LIVE boards (AS-DTA business PR)


# ---------------------------------------------------------------------------
# Shared row schema from sensors
# ---------------------------------------------------------------------------


@dataclass
class SensorRow:
    segment: str
    mu: float
    s: float
    J: float
    y_flip: float
    V_debiased: float
    V_naive: float
    V_ref: float = 0.0

    def decision(self) -> TriageDecision:
        return triage(
            mu=self.mu,
            s=self.s,
            V_debiased=self.V_debiased,
            V_naive=self.V_naive,
            J=self.J,
            y_flip=self.y_flip,
            V_ref=self.V_ref,
        )

    @property
    def severity(self) -> float:
        return float(self.J + self.y_flip)


def _demo_rows() -> List[SensorRow]:
    """Canonical industrial-style segment panel for demos / justify."""
    return [
        SensorRow("general", mu=0.04, s=0.01, J=0.04, y_flip=0.0, V_debiased=0.0, V_naive=0.01, V_ref=0.0),
        SensorRow("finance", mu=0.08, s=0.02, J=0.09, y_flip=0.0, V_debiased=0.001, V_naive=0.02, V_ref=0.0),
        SensorRow("medicine", mu=0.50, s=0.12, J=0.90, y_flip=0.32, V_debiased=0.08, V_naive=0.10, V_ref=0.01),
        SensorRow("sci", mu=0.40, s=0.10, J=0.70, y_flip=0.28, V_debiased=0.05, V_naive=0.07, V_ref=0.01),
        SensorRow("talk", mu=0.06, s=0.02, J=0.05, y_flip=0.0, V_debiased=0.0, V_naive=0.015, V_ref=0.0),
    ]


# ---------------------------------------------------------------------------
# 1. Dynamic routing & cost optimization
# ---------------------------------------------------------------------------


def app_dynamic_routing(
    rows: Optional[Sequence[SensorRow]] = None,
    *,
    cost_light: float = 0.002,
    cost_heavy: float = 0.04,
    n_requests: int = 100_000,
) -> dict:
    """Route benign → light/cache/conformal; malignant → escalate/heavy/block.

    Cost model: light path $cost_light / req, heavy $cost_heavy / req.
    Documented industrial saving: $22,430 on LIVE boards (business PR).
    """
    rows = list(rows or _demo_rows())
    routes = []
    for r in rows:
        d = r.decision()
        if d.call in ("quiet", "benign"):
            route = "light_model_or_cache_or_conformal"
            unit_cost = cost_light
        else:
            route = "escalate_heavy_or_adversarial_or_block"
            unit_cost = cost_heavy
        routes.append(
            {
                "segment": r.segment,
                "call": d.call,
                "route": route,
                "severity": r.severity,
                "unit_cost_usd": unit_cost,
                "reason": d.reason,
            }
        )
    # Simulate traffic mix from segment calls
    n_ben = sum(1 for x in routes if x["call"] in ("quiet", "benign"))
    n_mal = len(routes) - n_ben
    # share traffic: benign-like 70%, malignant-like 30% if both exist
    p_ben = 0.70 if n_ben and n_mal else (1.0 if n_ben else 0.0)
    cost_routed = n_requests * (p_ben * cost_light + (1 - p_ben) * cost_heavy)
    cost_all_heavy = n_requests * cost_heavy
    saved_sim = cost_all_heavy - cost_routed
    return {
        "app": "dynamic_routing_cost_optimization",
        "routes": routes,
        "n_requests": n_requests,
        "p_benign_share": p_ben,
        "cost_all_heavy_usd": cost_all_heavy,
        "cost_routed_usd": cost_routed,
        "cost_saved_sim_usd": saved_sim,
        "cost_saved_industrial_live_usd": COST_SAVED_USD,
        "value": (
            "PO-risk variance as confidence proxy: benign → cheap path; "
            "malignant → escalate. Industrial LIVE boards saved $%.0f."
            % COST_SAVED_USD
        ),
    }


# ---------------------------------------------------------------------------
# 2. Granular canary release & segment rollback
# ---------------------------------------------------------------------------


def app_canary_rollback(
    rows: Optional[Sequence[SensorRow]] = None,
    *,
    canary_model: str = "model_v_next",
    stable_model: str = "model_v_stable",
) -> dict:
    """Per-segment severity gate: rollback only the broken segment."""
    rows = list(rows or _demo_rows())
    actions = []
    for r in rows:
        d = r.decision()
        if d.call == "malignant":
            action = "rollback_segment_to_stable"
            model = stable_model
        else:
            action = "keep_canary"
            model = canary_model
        actions.append(
            {
                "segment": r.segment,
                "severity": r.severity,
                "y_flip": r.y_flip,
                "V_debiased": r.V_debiased,
                "call": d.call,
                "action": action,
                "serving_model": model,
                "reason": d.reason,
            }
        )
    n_roll = sum(1 for a in actions if a["action"].startswith("rollback"))
    return {
        "app": "granular_canary_release_rollback",
        "canary_model": canary_model,
        "stable_model": stable_model,
        "segments": actions,
        "n_segments": len(actions),
        "n_rolled_back": n_roll,
        "n_kept_canary": len(actions) - n_roll,
        "micro_loss": n_roll < len(actions),
        "value": (
            "Segment-level Severity=J+Y-flip + Var(u): rollback only the "
            "sick segment; other segments keep canary (micro-loss iteration)."
        ),
    }


# ---------------------------------------------------------------------------
# 3. Data flywheel & active learning curation
# ---------------------------------------------------------------------------


def app_data_flywheel(
    samples: Optional[Sequence[dict]] = None,
) -> dict:
    """Gate samples into fine-tune / RLHF sets.

    Drop high Y-flip (poison / integrity break).
    Keep Y-flip=0 with mild J as hard-but-clean gold.
    """
    if samples is None:
        rng = np.random.default_rng(2026)
        samples = []
        for i in range(40):
            yf = float(rng.choice([0.0, 0.0, 0.0, 0.25, 0.4]))
            J = float(rng.uniform(0.0, 0.4 if yf == 0 else 1.2))
            samples.append({"id": f"s{i}", "J": J, "y_flip": yf, "V": float(rng.uniform(0, 0.1))})

    keep_gold, drop_poison, quarantine = [], [], []
    for s in samples:
        yf, J = float(s["y_flip"]), float(s["J"])
        if yf >= 0.02:
            drop_poison.append({**s, "curate": "drop_poison_or_integrity_break"})
        elif J <= 0.65:
            # mild/quiet cosmetic — gold for robustness
            keep_gold.append({**s, "curate": "keep_gold_benign_hard"})
        else:
            quarantine.append({**s, "curate": "quarantine_review"})

    return {
        "app": "data_flywheel_active_learning",
        "n_in": len(samples),
        "n_gold": len(keep_gold),
        "n_drop": len(drop_poison),
        "n_quarantine": len(quarantine),
        "gold_rate": len(keep_gold) / max(len(samples), 1),
        "drop_rate": len(drop_poison) / max(len(samples), 1),
        "examples_gold": keep_gold[:5],
        "examples_drop": drop_poison[:5],
        "value": (
            "Gatekeeper for fine-tune/RLHF: filter high Y-flip poison; "
            "keep Y-flip=0 + mild J as distribution-expanding gold data."
        ),
    }


# ---------------------------------------------------------------------------
# 4. Explainable compliance / audit guardrails
# ---------------------------------------------------------------------------


def app_compliance_report(
    rows: Optional[Sequence[SensorRow]] = None,
    *,
    window: str = "2026-W44",
) -> dict:
    """Segment-level risk event report for auditors (AI Act style)."""
    rows = list(rows or _demo_rows())
    events = []
    for r in rows:
        d = r.decision()
        if d.call != "malignant":
            continue
        level = "dense" if r.y_flip >= 0.25 or r.J >= 1.0 else "sparse"
        events.append(
            {
                "window": window,
                "segment": r.segment,
                "integrity_break": True,
                "y_break": r.y_flip,
                "severity": r.severity,
                "J": r.J,
                "V_debiased": r.V_debiased,
                "level": level,
                "action_taken": "local_VAI_freeze",
                "narrative": (
                    "In %s, segment=%s detected integrity break "
                    "(Y-flip=%.2f, J=%.2f, V_debiased=%.3g) → local freeze."
                    % (window, r.segment, r.y_flip, r.J, r.V_debiased)
                ),
            }
        )
    return {
        "app": "explainable_compliance_guardrails",
        "window": window,
        "n_events": len(events),
        "events": events,
        "distinguishes": (
            "PO-risk separates cosmetic synonym / watermark scrub vs "
            "true outcome integrity break — auditable, not a black-box deny."
        ),
        "value": (
            "Turn Mag×Var + Y-flip + Var(u) into segment risk-event reports "
            "readable by compliance (EU AI Act style)."
        ),
    }


# ---------------------------------------------------------------------------
# 5. Agentic workflow cascade-error blocking
# ---------------------------------------------------------------------------


def app_agent_cascade_guard(
    steps: Optional[Sequence[dict]] = None,
) -> dict:
    """Per-step PO-risk check in a multi-step agent; abort on spike."""
    if steps is None:
        steps = [
            {"step": 1, "name": "retrieve", "J": 0.05, "y_flip": 0.0, "V_debiased": 0.0, "V_naive": 0.01},
            {"step": 2, "name": "plan", "J": 0.08, "y_flip": 0.0, "V_debiased": 0.002, "V_naive": 0.01},
            {"step": 3, "name": "tool_call", "J": 0.85, "y_flip": 0.20, "V_debiased": 0.06, "V_naive": 0.08},
            {"step": 4, "name": "synthesize", "J": 0.10, "y_flip": 0.0, "V_debiased": 0.01, "V_naive": 0.02},
        ]
    trace = []
    aborted_at = None
    for s in steps:
        d = triage(
            mu=float(s.get("mu", s["J"])),
            s=float(s.get("s", 0.1 if s["J"] > 0.5 else 0.01)),
            V_debiased=float(s["V_debiased"]),
            V_naive=float(s["V_naive"]),
            J=float(s["J"]),
            y_flip=float(s["y_flip"]),
            V_ref=float(s.get("V_ref", 0.0)),
        )
        row = {
            "step": s["step"],
            "name": s["name"],
            "call": d.call,
            "actuator": d.actuator,
            "severity": float(s["J"] + s["y_flip"]),
            "reason": d.reason,
        }
        if d.call == "malignant" and aborted_at is None:
            row["action"] = "abort_workflow_self_correct_or_human"
            aborted_at = s["step"]
            trace.append(row)
            break
        row["action"] = "continue"
        trace.append(row)

    return {
        "app": "agentic_cascade_error_blocking",
        "trace": trace,
        "aborted_at_step": aborted_at,
        "cascade_blocked": aborted_at is not None,
        "value": (
            "Embed Var(u)/Mag×Var at each agent step; abort + self-correct "
            "before early hallucination / attack cascades to the final act."
        ),
    }


# ---------------------------------------------------------------------------
# Suite
# ---------------------------------------------------------------------------


def run_product_apps_suite(
    *,
    out_dir: Optional[Path] = None,
) -> dict:
    out_dir = Path(out_dir or OUT_DEFAULT)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = _demo_rows()
    payload = {
        "scope": "Five PO-risk product applications on Mag×Var + Y-flip + CF-debiased Var(u)",
        "sensors": "severity=J+Y-flip; Var(u)=CF leaf-debiased (min leaf merge n=30/50); BCF via rpy2; bootstrap fallback",
        "apps": {
            "1_dynamic_routing": app_dynamic_routing(rows),
            "2_canary_rollback": app_canary_rollback(rows),
            "3_data_flywheel": app_data_flywheel(),
            "4_compliance": app_compliance_report(rows),
            "5_agent_cascade": app_agent_cascade_guard(),
        },
        "industrial_cost_saved_usd": COST_SAVED_USD,
        "headline": {
            "routing_saves_industrial_usd": COST_SAVED_USD,
            "canary_micro_loss": True,
            "flywheel_drops_poison": True,
            "compliance_events": True,
            "agent_cascade_blocked": True,
        },
    }
    # fill headline from apps
    a = payload["apps"]
    payload["headline"]["canary_micro_loss"] = bool(a["2_canary_rollback"]["micro_loss"])
    payload["headline"]["flywheel_drop_rate"] = a["3_data_flywheel"]["drop_rate"]
    payload["headline"]["compliance_n_events"] = a["4_compliance"]["n_events"]
    payload["headline"]["agent_cascade_blocked"] = a["5_agent_cascade"]["cascade_blocked"]

    (out_dir / "po_risk_var_product_apps.json").write_text(
        json.dumps(payload, indent=2, default=str)
    )
    _readme(out_dir / "README.md", payload)
    return payload


def _readme(path: Path, payload: dict):
    a = payload["apps"]
    lines = [
        "# PO-risk Var(u) — five product applications",
        "",
        payload["scope"],
        "",
        f"Sensors: {payload['sensors']}",
        "",
        f"Industrial LIVE cost saved (business PR): **${payload['industrial_cost_saved_usd']:.0f}**",
        "",
        "## 1. Dynamic routing & cost optimization",
        "",
        a["1_dynamic_routing"]["value"],
        "",
        f"- Simulated save vs all-heavy: **${a['1_dynamic_routing']['cost_saved_sim_usd']:.0f}**",
        f"- Documented industrial save: **${a['1_dynamic_routing']['cost_saved_industrial_live_usd']:.0f}**",
        "",
        "## 2. Granular canary & segment rollback",
        "",
        a["2_canary_rollback"]["value"],
        "",
        f"- Rolled back: **{a['2_canary_rollback']['n_rolled_back']}** / kept canary: **{a['2_canary_rollback']['n_kept_canary']}**",
        "",
        "## 3. Data flywheel / active learning",
        "",
        a["3_data_flywheel"]["value"],
        "",
        f"- Gold keep rate: **{a['3_data_flywheel']['gold_rate']:.2f}**; drop (poison) rate: **{a['3_data_flywheel']['drop_rate']:.2f}**",
        "",
        "## 4. Compliance / audit guardrails",
        "",
        a["4_compliance"]["value"],
        "",
        f"- Risk events: **{a['4_compliance']['n_events']}**",
        "",
    ]
    for e in a["4_compliance"]["events"]:
        lines.append(f"- {e['narrative']}")
    lines.extend(
        [
            "",
            "## 5. Agentic cascade-error blocking",
            "",
            a["5_agent_cascade"]["value"],
            "",
            f"- Cascade blocked: **{a['5_agent_cascade']['cascade_blocked']}** at step **{a['5_agent_cascade']['aborted_at_step']}**",
            "",
            "```bash",
            "po-risk-var products",
            "```",
            "",
        ]
    )
    path.write_text("\n".join(lines))
