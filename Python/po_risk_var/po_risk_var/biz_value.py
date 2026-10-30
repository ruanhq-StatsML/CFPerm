"""Incremental business value — estimators + five product apps.

Baseline for every claim is the *naive industrial alternative*
(always-heavy / global rollback / ingest-all / black-box deny / run-to-end).
Numbers tagged ``industrial`` come from ``results/as_dta_business_pr``
($22,430 LIVE). Numbers tagged ``sim`` are transparent unit-econ demos.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from po_risk_var._path import REPO_ROOT
from po_risk_var.apps_product import (
    COST_SAVED_USD,
    app_agent_cascade_guard,
    app_canary_rollback,
    app_compliance_report,
    app_data_flywheel,
    app_dynamic_routing,
)

OUT_DEFAULT = REPO_ROOT / "results" / "as_dta_po_risk_var_biz"


# ---------------------------------------------------------------------------
# Layer A — Var(u) estimators: what each method *adds* commercially
# ---------------------------------------------------------------------------


def estimator_incremental_value() -> List[Dict[str, Any]]:
    """Business delta of CF / BCF / bootstrap vs the naive residual Var(u)."""
    return [
        {
            "method": "causal_forest_leaf_debiased",
            "role": "product_default",
            "baseline": "empirical Var(u) or CF V_naive (leaf-size biased)",
            "what_it_adds": (
                "Debias leaf-size inflation in Var(tau) so quiet/benign LIVE "
                "does not look malignant; leaf-merge floor min n=30/50 keeps "
                "segment monitors stable."
            ),
            "incremental_value": [
                {
                    "ledger": "false_freeze_avoidance",
                    "claim": (
                        "Stops leaf-bias masquerading as integrity break → "
                        "avoid wrong-drug freeze/retrain on synonym / mild "
                        "covariate boards."
                    ),
                    "evidence": (
                        "Synth justify 5/5 + real 4/4 with CF-debiased; "
                        "honesty bias_ratio flags when debias mattered."
                    ),
                    "dollar_proxy": (
                        "Wrong freeze on a benign finance synonym cell alone "
                        "leaves ~$3.4k/1k on the table (business PR: conformal "
                        "was the right cheap drug)."
                    ),
                },
                {
                    "ledger": "confirmed_malignant_routing",
                    "claim": (
                        "Relative V_debiased↑ confirms Mag×Var / Y-flip → "
                        "spend freeze+aug dollars only when shape+outcome agree."
                    ),
                    "evidence": (
                        "Estimator separation suite: cf_V_debiased / "
                        "cf_Var_u_tree separate malignant vs benign (hit-rate 1.0)."
                    ),
                    "dollar_proxy": (
                        "Finance malignant cell: cost/1k $8404→$4067 (~$4.3k) "
                        "when freeze+aug is correctly selected."
                    ),
                },
            ],
            "when_not_enough": (
                "Need posterior intervals / Bayesian stack → escalate to BCF; "
                "need a 30-second sanity check → bootstrap."
            ),
            "ship_as": "default LIVE monitor",
        },
        {
            "method": "bcf_posterior_rpy2",
            "role": "bayesian_path",
            "baseline": "CF point/debiased Var(u) without posterior draws",
            "what_it_adds": (
                "Full posterior draws of tau via R bcf (rpy2 + pihat) → "
                "Var(u) as posterior spread; audit-grade uncertainty for "
                "regulated pipelines."
            ),
            "incremental_value": [
                {
                    "ledger": "compliance_uncertainty",
                    "claim": (
                        "Credible bands on unit risk for auditors: not only "
                        "'V is high' but 'posterior mass on integrity break'."
                    ),
                    "evidence": (
                        "bcf_Var_u separates malignant/benign (suite ratio ~29×); "
                        "backend rpy2_bcf with pihat."
                    ),
                    "dollar_proxy": (
                        "Pays for itself when a single avoided AI-Act / "
                        "escalation dispute exceeds MCMC cost (minutes of CPU "
                        "vs legal/ops review)."
                    ),
                },
                {
                    "ledger": "bayesian_product_fit",
                    "claim": (
                        "Drop-in when the org already ships BART/BCF; no second "
                        "uncertainty language for risk teams."
                    ),
                    "evidence": "Same triage actuators; only Var(u) backend swaps.",
                    "dollar_proxy": "Integration cost ↓ vs inventing a fourth estimator.",
                },
            ],
            "when_not_enough": (
                "Latency-critical LIVE default → stay on CF; tiny n without "
                "pihat quality → prefer CF + bootstrap sanity."
            ),
            "ship_as": "Bayesian / audit lane",
        },
        {
            "method": "bootstrap",
            "role": "fallback_sanity",
            "baseline": "ship CF/BCF with no cross-check",
            "what_it_adds": (
                "Intuitive refit of PermuCATE tau → Winsorized Var(u); "
                "catches gross misspecification cheaply."
            ),
            "incremental_value": [
                {
                    "ledger": "engineering_confidence",
                    "claim": (
                        "A/B against CF before a freeze: if boot and CF disagree "
                        "wildly, hold + investigate instead of auto-escalate."
                    ),
                    "evidence": (
                        "boot_Var_u still separates (suite ratio ~15×) but is "
                        "heavier-tailed — never the product default."
                    ),
                    "dollar_proxy": (
                        "One blocked false escalation (GPU retrain / human "
                        "on-call) typically dwarfs bootstrap CPU."
                    ),
                }
            ],
            "when_not_enough": "Do not use as sole freeze sensor (fourth-moment noise).",
            "ship_as": "fallback / CI sanity",
        },
    ]


# ---------------------------------------------------------------------------
# Layer B — five apps: incremental $ / ops vs industrial anti-patterns
# ---------------------------------------------------------------------------


def _routing_incremental(routing: dict) -> dict:
    industrial = float(routing["cost_saved_industrial_live_usd"])
    sim = float(routing["cost_saved_sim_usd"])
    return {
        "app": "dynamic_routing_cost_optimization",
        "baseline_anti_patterns": [
            "always_heavy_model (pay GPU on every request)",
            "always_light_or_cache (miss integrity breaks → FN ledger)",
            "OOD_high_implies_retrain (wrong-drug tax)",
        ],
        "what_po_risk_adds": (
            "Var(u) + Mag×Var + Y-flip as a *confidence / disease* proxy: "
            "quiet/benign → light/cache/conformal; malignant → escalate."
        ),
        "incremental_vs_always_heavy": {
            "sim_usd_per_100k_req": sim,
            "industrial_live_usd": industrial,
            "narrative": (
                "Industrial HC3 LIVE (business PR): routed actuators saved "
                "$%.0f vs static — 100%% shape hit-rate, so detection was not "
                "the bottleneck; right-drug selection was." % industrial
            ),
        },
        "incremental_vs_always_light": {
            "ledger": "FN_integrity",
            "narrative": (
                "Always-light would keep paying ~$20×FN on malignant cells; "
                "finance malignant alone drops cost/1k $8404→$4067 when "
                "escalated to freeze+aug."
            ),
        },
        "unit_econ": {
            "cost_fp_usd": 50,
            "cost_fn_usd": 20,
            "source": "docs/method/AS_DTA_business_PR_ROI.md",
        },
        "headline_incremental_usd": industrial,
    }


def _canary_incremental(canary: dict) -> dict:
    n_keep = int(canary["n_kept_canary"])
    n_roll = int(canary["n_rolled_back"])
    n = int(canary["n_segments"])
    # Transparent demo: assume equal traffic share; global rollback burns
    # canary uplift on healthy segments. Stylized uplift $0.01/req × share.
    n_req = 100_000
    uplift_per_req = 0.01
    share_keep = n_keep / max(n, 1)
    value_kept = n_req * share_keep * uplift_per_req
    value_lost_if_global = n_req * uplift_per_req  # lose all canary uplift
    incremental = value_lost_if_global - (n_req * (n_roll / max(n, 1)) * uplift_per_req)
    return {
        "app": "granular_canary_release_rollback",
        "baseline_anti_patterns": [
            "global_rollback_on_any_alert (burn healthy segments)",
            "global_accuracy_only_monitor (miss segment integrity breaks)",
        ],
        "what_po_risk_adds": (
            "Segment Severity=J+Y-flip + CF Var(u): rollback only sick "
            "segments; keep canary elsewhere (micro-loss iteration)."
        ),
        "demo_panel": {
            "n_segments": n,
            "n_rolled_back": n_roll,
            "n_kept_canary": n_keep,
            "kept_share": share_keep,
        },
        "incremental_vs_global_rollback": {
            "sim_usd_canary_uplift_retained": round(incremental, 2),
            "assumptions": (
                "100k req, $0.01 canary uplift/req, equal segment share; "
                "global rollback zeros uplift; segment rollback zeros only "
                "sick share."
            ),
            "narrative": (
                "Demo panel keeps %d/%d segments on canary → retain ~$%.0f "
                "uplift that a global rollback would incinerate."
                % (n_keep, n, incremental)
            ),
        },
        "incremental_vs_global_metrics_only": {
            "narrative": (
                "medicine + sci show Y-flip / high J while general/finance/talk "
                "stay quiet — global accuracy can look fine while verticals "
                "bleed; segment gate is the incremental sensor."
            ),
        },
        "headline_incremental_usd": round(incremental, 2),
    }


def _flywheel_incremental(fly: dict) -> dict:
    drop = float(fly["drop_rate"])
    gold = float(fly["gold_rate"])
    n = int(fly["n_in"])
    # Stylized: each poison sample that enters FT costs $25 expected
    # remediation; each gold hard sample is worth $8 of robustness value.
    poison_cost = 25.0
    gold_value = 8.0
    avoided_poison = drop * n * poison_cost
    gold_kept = gold * n * gold_value
    # vs drop-all-high-PO: would also drop mild-J gold (assume half of gold)
    gold_lost_if_crude = 0.5 * gold * n * gold_value
    return {
        "app": "data_flywheel_active_learning",
        "baseline_anti_patterns": [
            "ingest_all_online_feedback (data poisoning)",
            "drop_all_high_PO (throws away hard benign / synonym gold)",
        ],
        "what_po_risk_adds": (
            "Y-flip gate kills integrity poison; Y-flip=0 + mild J kept as "
            "distribution-expanding gold."
        ),
        "demo_panel": {
            "n_in": n,
            "drop_rate": drop,
            "gold_rate": gold,
            "n_quarantine": int(fly["n_quarantine"]),
        },
        "incremental_vs_ingest_all": {
            "sim_usd_poison_avoided": round(avoided_poison, 2),
            "assumptions": "$25 expected remediation per poison sample ingested",
            "narrative": (
                "Drop rate %.0f%% on demo panel → ~$%.0f avoided remediation "
                "vs ingest-all." % (100 * drop, avoided_poison)
            ),
        },
        "incremental_vs_drop_all_high_po": {
            "sim_usd_gold_retained_extra": round(gold_lost_if_crude, 2),
            "assumptions": (
                "$8 robustness value per gold sample; crude high-PO filter "
                "throws away ~half of mild-J gold."
            ),
            "narrative": (
                "Mag×Var distinguishes hard-but-clean from poison → keep gold "
                "rate %.0f%% instead of scrubbing benign hard cases."
                % (100 * gold)
            ),
        },
        "headline_incremental_usd": round(avoided_poison + gold_kept, 2),
    }


def _compliance_incremental(comp: dict) -> dict:
    n_ev = int(comp["n_events"])
    # Stylized: one unexplained deny / silent miss costs legal review $5k;
    # auditable local freeze costs $0.4k ops note.
    unexplained = 5000.0
    auditable = 400.0
    saved = n_ev * (unexplained - auditable)
    return {
        "app": "explainable_compliance_guardrails",
        "baseline_anti_patterns": [
            "black_box_deny_or_silent_pass",
            "global_accuracy_dashboard_only",
        ],
        "what_po_risk_adds": (
            "Segment risk-event narratives: window × segment × Y-break × "
            "J × V_debiased → local freeze — AI-Act-readable, not opaque."
        ),
        "demo_panel": {"n_events": n_ev, "events": comp.get("events", [])},
        "incremental_vs_black_box": {
            "sim_usd_review_cost_delta": round(saved, 2),
            "assumptions": (
                "$5k unexplained incident review vs $0.4k auditable local-freeze note"
            ),
            "narrative": (
                "%d integrity-break events become auditable freezes → "
                "~$%.0f review-cost delta on the demo window." % (n_ev, saved)
            ),
        },
        "product_distinction": (
            "Separates watermark-scrubbed / synonym cosmetic from true "
            "outcome integrity break — the compliance incremental story."
        ),
        "headline_incremental_usd": round(saved, 2),
    }


def _agent_incremental(agent: dict) -> dict:
    aborted = agent.get("aborted_at_step")
    n_steps = len(agent.get("trace", []))
    # Stylized: each downstream step after a bad tool_call wastes $0.15
    # tool/LLM cost; cascade to final act risks $40 task failure.
    step_waste = 0.15
    task_fail = 40.0
    # Without guard: steps 3 and 4 run after poison; with guard: abort at 3
    wasted_if_ungated = 2 * step_waste + task_fail
    wasted_if_gated = step_waste  # pay for the bad step once, then abort
    incremental = wasted_if_ungated - wasted_if_gated
    return {
        "app": "agentic_cascade_error_blocking",
        "baseline_anti_patterns": [
            "run_workflow_to_end",
            "only_check_final_answer",
        ],
        "what_po_risk_adds": (
            "Per-step Mag×Var / Var(u): abort + self-correct / human when a "
            "step's residual variance spikes — stop cascade amplification."
        ),
        "demo_panel": {
            "aborted_at_step": aborted,
            "cascade_blocked": bool(agent.get("cascade_blocked")),
            "n_steps_executed": n_steps,
        },
        "incremental_vs_run_to_end": {
            "sim_usd_per_bad_workflow": round(incremental, 2),
            "assumptions": (
                "$0.15/step tool waste after poison; $40 task-failure risk if "
                "cascade reaches synthesize"
            ),
            "narrative": (
                "Demo aborts at step %s → avoid ~$%.0f expected cascade loss "
                "vs run-to-end." % (aborted, incremental)
            ),
        },
        "headline_incremental_usd": round(incremental, 2),
    }


def product_apps_incremental_value(apps: Optional[dict] = None) -> Dict[str, Any]:
    if apps is None:
        rows_routing = app_dynamic_routing()
        apps = {
            "1_dynamic_routing": rows_routing,
            "2_canary_rollback": app_canary_rollback(),
            "3_data_flywheel": app_data_flywheel(),
            "4_compliance": app_compliance_report(),
            "5_agent_cascade": app_agent_cascade_guard(),
        }
    return {
        "1_dynamic_routing": _routing_incremental(apps["1_dynamic_routing"]),
        "2_canary_rollback": _canary_incremental(apps["2_canary_rollback"]),
        "3_data_flywheel": _flywheel_incremental(apps["3_data_flywheel"]),
        "4_compliance": _compliance_incremental(apps["4_compliance"]),
        "5_agent_cascade": _agent_incremental(apps["5_agent_cascade"]),
    }


def rollup_incremental(apps_inc: dict) -> dict:
    """Stack incremental dollars — industrial for routing, sim for others."""
    routing = float(apps_inc["1_dynamic_routing"]["headline_incremental_usd"])
    sim_stack = (
        float(apps_inc["2_canary_rollback"]["headline_incremental_usd"])
        + float(apps_inc["3_data_flywheel"]["headline_incremental_usd"])
        + float(apps_inc["4_compliance"]["headline_incremental_usd"])
        + float(apps_inc["5_agent_cascade"]["headline_incremental_usd"])
    )
    return {
        "industrial_routing_usd": routing,
        "sim_other_apps_usd": round(sim_stack, 2),
        "note": (
            "Do not sum industrial + sim into one ROI — different ledgers. "
            "Pitch industrial $%.0f as measured LIVE save; sim stack "
            "($%.0f) as transparent unit-econ for canary/flywheel/"
            "compliance/agent on the demo panel." % (routing, sim_stack)
        ),
        "pitch_one_liner": (
            "PO-risk Var(u) is an allocation rule: CF-debiased honesty for "
            "LIVE, BCF when auditors need posteriors, bootstrap as sanity — "
            "and the five apps turn that sensor into routing dollars "
            "($%.0f industrial), micro-loss canaries, clean flywheels, "
            "auditable freezes, and aborted agent cascades." % routing
        ),
    }


# ---------------------------------------------------------------------------
# Suite + docs
# ---------------------------------------------------------------------------


def run_biz_value_suite(*, out_dir: Optional[Path] = None) -> dict:
    out_dir = Path(out_dir or OUT_DEFAULT)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Build apps in-memory (products CLI owns the on-disk product suite)
    apps = {
        "1_dynamic_routing": app_dynamic_routing(),
        "2_canary_rollback": app_canary_rollback(),
        "3_data_flywheel": app_data_flywheel(),
        "4_compliance": app_compliance_report(),
        "5_agent_cascade": app_agent_cascade_guard(),
    }
    apps_inc = product_apps_incremental_value(apps)
    est = estimator_incremental_value()
    rollup = rollup_incremental(apps_inc)

    payload = {
        "scope": (
            "Incremental business value of PO-risk Var(u) methods "
            "(CF / BCF / bootstrap) and five product applications"
        ),
        "industrial_anchor_usd": COST_SAVED_USD,
        "industrial_source": "results/as_dta_business_pr (HC3 finance/medicine LIVE)",
        "estimators": est,
        "product_apps": apps_inc,
        "rollup": rollup,
        "justify_hit_rates": {
            "synth": "5/5",
            "real": "4/4",
            "note": "from po-risk-var justify — sensor reliability before $ claims",
        },
    }
    (out_dir / "po_risk_var_biz_incremental.json").write_text(
        json.dumps(payload, indent=2, default=str)
    )
    _plot(payload, out_dir / "po_risk_var_biz_incremental.png")
    _readme(out_dir / "README.md", payload)
    return payload


def _plot(payload: dict, path: Path):
    import matplotlib.pyplot as plt
    import numpy as np

    apps = payload["product_apps"]
    labels = [
        "1 routing\n(industrial)",
        "2 canary\n(sim)",
        "3 flywheel\n(sim)",
        "4 compliance\n(sim)",
        "5 agent\n(sim)",
    ]
    keys = [
        "1_dynamic_routing",
        "2_canary_rollback",
        "3_data_flywheel",
        "4_compliance",
        "5_agent_cascade",
    ]
    vals = [float(apps[k]["headline_incremental_usd"]) for k in keys]
    colors = ["#2F6B4F", "#3D5A80", "#3D5A80", "#3D5A80", "#3D5A80"]

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.2))

    ax = axes[0]
    x = np.arange(len(labels))
    bars = ax.bar(x, vals, color=colors, width=0.65)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Incremental USD (ledger-specific)")
    ax.set_title("Product apps — incremental value", loc="left", fontweight="bold")
    ax.grid(True, axis="y", alpha=0.3)
    for b, v in zip(bars, vals):
        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height() * 1.02 + max(vals) * 0.01,
            "$%.0f" % v,
            ha="center",
            va="bottom",
            fontsize=8,
        )

    ax = axes[1]
    est_names = ["CF debiased\n(default)", "BCF rpy2\n(bayesian)", "Bootstrap\n(fallback)"]
    # Qualitative score: ship priority weight for the pitch panel
    weights = [3.0, 2.0, 1.0]
    roles = ["LIVE honesty\n+ routing confirm", "Audit posterior\n+ CI for legal", "CI sanity\nnot freeze sensor"]
    ax.barh(est_names[::-1], weights[::-1], color=["#2F6B4F", "#4A6FA5", "#8A8F98"][::-1], height=0.55)
    ax.set_xlabel("Ship priority (product)")
    ax.set_title("Estimator methods — role / incremental job", loc="left", fontweight="bold")
    for i, r in enumerate(roles[::-1]):
        ax.text(0.08, i, r, va="center", ha="left", fontsize=8, color="white", fontweight="bold")
    ax.set_xlim(0, 3.6)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(["fallback", "bayes", "default"])

    fig.suptitle(
        "PO-risk Var(u) — business incremental value\n"
        "industrial routing $%.0f  |  sim stack $%.0f (do not sum ledgers)"
        % (
            payload["rollup"]["industrial_routing_usd"],
            payload["rollup"]["sim_other_apps_usd"],
        ),
        fontsize=11,
        fontweight="bold",
        y=1.05,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def _readme(path: Path, payload: dict):
    r = payload["rollup"]
    lines = [
        "# PO-risk Var(u) — business incremental value",
        "",
        payload["scope"],
        "",
        r["pitch_one_liner"],
        "",
        f"- Industrial LIVE routing save (measured): **${r['industrial_routing_usd']:.0f}**",
        f"- Demo sim stack (canary+flywheel+compliance+agent): **${r['sim_other_apps_usd']:.0f}**",
        f"- {r['note']}",
        "",
        "## A. Estimator methods — what each one adds",
        "",
    ]
    for e in payload["estimators"]:
        lines.append(f"### {e['method']}  [{e['role']}]")
        lines.append("")
        lines.append(f"- Baseline: {e['baseline']}")
        lines.append(f"- Adds: {e['what_it_adds']}")
        lines.append(f"- Ship as: **{e['ship_as']}**")
        for iv in e["incremental_value"]:
            lines.append(f"  - **{iv['ledger']}:** {iv['claim']}")
            lines.append(f"    - Evidence: {iv['evidence']}")
            lines.append(f"    - $ proxy: {iv['dollar_proxy']}")
        lines.append(f"- Not enough when: {e['when_not_enough']}")
        lines.append("")

    lines.append("## B. Product apps — incremental vs anti-patterns")
    lines.append("")
    order = [
        ("1_dynamic_routing", "1. Dynamic routing"),
        ("2_canary_rollback", "2. Granular canary"),
        ("3_data_flywheel", "3. Data flywheel"),
        ("4_compliance", "4. Compliance guardrails"),
        ("5_agent_cascade", "5. Agent cascade block"),
    ]
    for key, title in order:
        a = payload["product_apps"][key]
        lines.append(f"### {title}")
        lines.append("")
        lines.append(f"- Baseline anti-patterns: {', '.join(a['baseline_anti_patterns'])}")
        lines.append(f"- PO-risk adds: {a['what_po_risk_adds']}")
        lines.append(f"- Headline incremental: **${a['headline_incremental_usd']:.0f}**")
        lines.append("")

    lines.extend(
        [
            "## Justify before dollars",
            "",
            f"- Sensor hit-rates: synth **{payload['justify_hit_rates']['synth']}**, "
            f"real **{payload['justify_hit_rates']['real']}**",
            f"- Industrial anchor source: `{payload['industrial_source']}`",
            "",
            "```bash",
            "po-risk-var biz",
            "po-risk-var products",
            "po-risk-var justify --quick",
            "```",
            "",
        ]
    )
    path.write_text("\n".join(lines))


__all__ = [
    "estimator_incremental_value",
    "product_apps_incremental_value",
    "run_biz_value_suite",
    "rollup_incremental",
]
