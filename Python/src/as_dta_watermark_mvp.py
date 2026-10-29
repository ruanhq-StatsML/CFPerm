"""AS-DTA Watermark evidence MVP — real-data forms, strength, routing.

Six evidence *kinds* you actually see in production tables (not a single
boolean ``is_watermarked``):

  api_log      — generator/API provenance (strongest "watermark equivalent")
  hard_wm      — explicit bit / green-list decoder watermark
  soft_wm      — statistical soft watermark score
  fingerprint  — passive stylometric / likelihood fingerprint
  user_declare — user/plugin "AI-assisted" checkbox (policy, not crypto)
  missing      — detector not run / too short / out-of-family

Pipeline: kind → (score, robust, coverage) → strength → AS-DTA action.
Vague = gray watermark ∪ low detector margin (hook only in this MVP).
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

SEED = 2026
ROOT = Path(__file__).resolve().parents[2]

# ---------------------------------------------------------------------------
# Taxonomy
# ---------------------------------------------------------------------------

KINDS = (
    "api_log",
    "hard_wm",
    "soft_wm",
    "fingerprint",
    "user_declare",
    "missing",
)

STRENGTHS = ("strong_pos", "gray", "strong_neg", "missing")

# Capability matrix: what each kind can / cannot do (1=yes, 0.5=partial, 0=no)
CAPABILITY = {
    #                provenance  crypto  survives_paraphrase  works_open_model  sparse  false_pos_risk
    "api_log":       [1.0,       0.0,    1.0,                 0.0,              1.0,    0.1],
    "hard_wm":       [0.5,       1.0,    0.4,                 0.0,              0.7,    0.2],
    "soft_wm":       [0.3,       0.5,    0.35,                0.0,              0.6,    0.35],
    "fingerprint":   [0.1,       0.0,    0.25,                0.6,              0.2,    0.7],
    "user_declare":  [0.4,       0.0,    1.0,                 1.0,              0.8,    0.5],
    "missing":       [0.0,       0.0,    0.0,                 0.0,              0.0,    0.0],
}
CAP_LABELS = [
    "provenance",
    "cryptographic",
    "survives_paraphrase",
    "covers_open_models",
    "when_present_is_sparse",
    "false_positive_risk",
]

KIND_BLURB = {
    "api_log": "Request/model_id in logs — strongest real-data form",
    "hard_wm": "Bit/green-list watermark recovered with p-value",
    "soft_wm": "Soft statistical watermark score + robustness probe",
    "fingerprint": "Passive style/likelihood prior — weak, drifts with OOD",
    "user_declare": "Checkbox / plugin declare — policy label, not crypto",
    "missing": "Not scored (short / wrong family / detector off)",
}


@dataclass
class WmEvidence:
    doc_id: str
    kind: str
    score: float
    pvalue: float
    robust_score: float
    coverage: float
    n_tokens: int
    channel: str
    strength: str = ""
    action: str = ""
    vague: bool = False
    meta: dict = field(default_factory=dict)


def map_strength(
    kind: str,
    score: float,
    pvalue: float,
    robust_score: float,
    coverage: float,
    n_tokens: int,
) -> str:
    """Map raw detector fields → strength tier (production-shaped rules)."""
    if kind == "missing" or n_tokens < 40 or coverage < 0.15:
        return "missing"
    if kind == "api_log":
        # score encodes log confidence; pvalue unused
        if score >= 0.85 and coverage >= 0.5:
            return "strong_pos"
        if score >= 0.55:
            return "gray"
        return "strong_neg"
    if kind == "user_declare":
        return "strong_pos" if score >= 0.9 else ("gray" if score >= 0.5 else "strong_neg")
    if kind == "fingerprint":
        # never strong_pos alone — always gray at best
        if score >= 0.8 and robust_score >= 0.55:
            return "gray"
        if score < 0.4:
            return "strong_neg"
        return "gray"
    # hard_wm / soft_wm
    if pvalue <= 0.01 and score >= 0.7 and robust_score >= 0.45 and coverage >= 0.35:
        return "strong_pos"
    if pvalue <= 0.05 and score >= 0.55:
        return "gray"
    if pvalue > 0.2 and score < 0.4:
        return "strong_neg"
    return "gray"


def route_action(strength: str, kind: str, detector_margin: float = 0.25) -> dict:
    """Strength → AS-DTA hook. Vague if gray or tiny detector margin."""
    vague = strength == "gray" or detector_margin < 0.08
    if strength == "strong_pos":
        action = "provenance_short_circuit"  # dispose / audit; optional weak ai label
    elif strength == "missing" or strength == "strong_neg":
        action = "po_fn_chain"  # Magnitude → PO+FN (needs strong Y) / OOD proxy
    elif strength == "gray":
        action = "vague_human_review"
    else:
        action = "po_fn_chain"
    if vague and action == "provenance_short_circuit":
        # shouldn't happen often; keep short-circuit
        pass
    elif vague and action != "vague_human_review":
        # low margin on top of neg/missing → still prefer review if margin tiny
        if detector_margin < 0.05:
            action = "vague_human_review"
            vague = True
    return {"action": action, "vague": bool(vague or action == "vague_human_review")}


def attach_routing(ev: WmEvidence, detector_margin: float) -> WmEvidence:
    ev.strength = map_strength(
        ev.kind, ev.score, ev.pvalue, ev.robust_score, ev.coverage, ev.n_tokens
    )
    r = route_action(ev.strength, ev.kind, detector_margin=detector_margin)
    ev.action = r["action"]
    ev.vague = r["vague"]
    return ev


# ---------------------------------------------------------------------------
# Synthetic real-data-shaped LIVE (MVP — no external WM libs required)
# ---------------------------------------------------------------------------


def simulate_live_batch(n_per_kind: int = 80, seed: int = SEED) -> List[WmEvidence]:
    """Draw a LIVE table with the six kinds + realistic score geometries."""
    rng = np.random.default_rng(int(seed))
    channels = ["web_form", "api_gateway", "email_ingest", "docs_plugin"]
    rows: List[WmEvidence] = []
    i = 0
    for kind in KINDS:
        for _ in range(int(n_per_kind)):
            i += 1
            n_tok = int(rng.integers(20, 400))
            channel = str(rng.choice(channels))
            if kind == "missing":
                score = p = robust = coverage = float("nan")
                n_tok = int(rng.choice([15, 25, 80]))  # often short
                coverage = 0.0
                score = p = robust = 0.0
            elif kind == "api_log":
                # bimodal: most calls confident
                score = float(np.clip(rng.beta(8, 2) if rng.random() < 0.75 else rng.beta(2, 5), 0, 1))
                p = 1.0
                robust = score  # logs don't "wash"
                coverage = float(np.clip(rng.uniform(0.6, 1.0), 0, 1))
            elif kind == "hard_wm":
                # mixture: clear positives + washed + negatives
                u = rng.random()
                if u < 0.35:
                    score, p, robust = float(rng.uniform(0.75, 0.98)), float(rng.uniform(1e-5, 0.008)), float(rng.uniform(0.5, 0.9))
                elif u < 0.55:
                    # washed paraphrase
                    score, p, robust = float(rng.uniform(0.45, 0.7)), float(rng.uniform(0.02, 0.08)), float(rng.uniform(0.15, 0.4))
                else:
                    score, p, robust = float(rng.uniform(0.05, 0.35)), float(rng.uniform(0.25, 0.9)), float(rng.uniform(0.05, 0.35))
                coverage = float(np.clip(n_tok / 250.0, 0.1, 1.0))
            elif kind == "soft_wm":
                u = rng.random()
                if u < 0.3:
                    score, p, robust = float(rng.uniform(0.7, 0.95)), float(rng.uniform(0.001, 0.02)), float(rng.uniform(0.4, 0.75))
                elif u < 0.6:
                    score, p, robust = float(rng.uniform(0.5, 0.72)), float(rng.uniform(0.03, 0.12)), float(rng.uniform(0.25, 0.5))
                else:
                    score, p, robust = float(rng.uniform(0.1, 0.45)), float(rng.uniform(0.15, 0.8)), float(rng.uniform(0.1, 0.4))
                coverage = float(np.clip(n_tok / 280.0, 0.08, 1.0))
            elif kind == "fingerprint":
                score = float(np.clip(rng.beta(2.2, 2.2), 0, 1))
                p = float(1.0 - score)
                robust = float(np.clip(score * rng.uniform(0.3, 0.7), 0, 1))  # washes easily
                coverage = float(np.clip(n_tok / 300.0, 0.05, 1.0))
            else:  # user_declare
                # sparse declares
                score = 1.0 if rng.random() < 0.25 else (0.6 if rng.random() < 0.1 else 0.0)
                p = 0.0 if score >= 0.9 else 0.5
                robust = score
                coverage = 1.0 if score > 0 else 0.0

            margin = float(np.clip(abs(rng.normal(0.22, 0.12)), 0.01, 0.49))
            ev = WmEvidence(
                doc_id="d%05d" % i,
                kind=kind,
                score=float(score) if score == score else 0.0,
                pvalue=float(p) if p == p else 1.0,
                robust_score=float(robust) if robust == robust else 0.0,
                coverage=float(coverage) if coverage == coverage else 0.0,
                n_tokens=int(n_tok),
                channel=channel,
                meta={"detector_margin": margin},
            )
            rows.append(attach_routing(ev, detector_margin=margin))
    return rows


def summarize(rows: Sequence[WmEvidence]) -> dict:
    rows = list(rows)
    by_kind: Dict[str, dict] = {}
    for k in KINDS:
        sub = [r for r in rows if r.kind == k]
        if not sub:
            continue
        streng = {s: 0 for s in STRENGTHS}
        acts: Dict[str, int] = {}
        for r in sub:
            streng[r.strength] = streng.get(r.strength, 0) + 1
            acts[r.action] = acts.get(r.action, 0) + 1
        by_kind[k] = {
            "n": len(sub),
            "strength": streng,
            "actions": acts,
            "vague_rate": float(np.mean([r.vague for r in sub])),
            "mean_score": float(np.nanmean([r.score for r in sub])),
            "mean_robust": float(np.nanmean([r.robust_score for r in sub])),
            "blurb": KIND_BLURB[k],
        }
    # transitions kind → strength → action
    flow = []
    for r in rows:
        flow.append({"kind": r.kind, "strength": r.strength, "action": r.action})
    return {
        "n": len(rows),
        "by_kind": by_kind,
        "capability": CAPABILITY,
        "capability_labels": CAP_LABELS,
        "kind_blurbs": KIND_BLURB,
        "flow": flow,
        "same_diff": {
            "same": [
                "All kinds emit sparse evidence into the same strength enum",
                "Gray / low margin → vague → human (never hard strong-Y)",
                "strong_neg / missing → never write Y=human; hand off to PO+FN",
                "Wash/missing is observationally shared by evade vs privacy vs OOD",
            ],
            "diff": [
                "api_log: strongest provenance; robust≈score (no wash)",
                "hard/soft_wm: p-value + robust probe; wash → gray",
                "fingerprint: never strong_pos alone; high FP under OOD",
                "user_declare: policy bit, not cryptographic recovery",
                "missing: structural absence, not a negative test result",
            ],
            "negative_effects": {
                "evade_intent": (
                    "Adversary washes WM / strips api_log to spoof human — "
                    "if treated as human, FN + poisoned uplift Y"
                ),
                "privacy_intent": (
                    "User/open-model redacts provenance / strips WM for privacy — "
                    "same geometry as evade; sanctioning wash = trust/legal FP"
                ),
                "uplift_misuse": (
                    "Using WM as hard filter or pseudo-Y induces selection bias "
                    "and label noise in R_PO / ATE; WM must stay gate+stratum"
                ),
            },
            "intent_routing": {
                "evade": "attack-shaped → PO+FN / freeze-lexical / audit",
                "privacy": "stratify; hold; never force Y=AI from absence",
                "benign_ood": "conformal_calibrate — not attack actuator",
            },
        },
    }


def plot_mvp(rows: Sequence[WmEvidence], summary: dict, path: Path):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    kind_colors = {
        "api_log": "#1F6F8B",
        "hard_wm": "#2F6B4F",
        "soft_wm": "#C4892A",
        "fingerprint": "#8B5A2B",
        "user_declare": "#5B3A8C",
        "missing": "#9AA3AE",
    }
    fig = plt.figure(figsize=(12.5, 9.2))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.28)

    # --- capability heatmap ---
    ax = fig.add_subplot(gs[0, 0])
    mat = np.array([CAPABILITY[k] for k in KINDS], dtype=float)
    im = ax.imshow(mat, cmap="YlGnBu", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(CAP_LABELS)))
    ax.set_xticklabels(CAP_LABELS, rotation=30, ha="right", fontsize=7.5)
    ax.set_yticks(range(len(KINDS)))
    ax.set_yticklabels(KINDS, fontsize=8)
    ax.set_title("Kind capability matrix (same schema, different powers)", loc="left", fontweight="bold")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, "%.1f" % mat[i, j], ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.046)

    # --- score vs robust scatter ---
    ax = fig.add_subplot(gs[0, 1])
    for k in KINDS:
        if k == "missing":
            continue
        sub = [r for r in rows if r.kind == k]
        ax.scatter(
            [r.score for r in sub],
            [r.robust_score for r in sub],
            s=22,
            alpha=0.55,
            c=kind_colors[k],
            label=k,
            edgecolors="none",
        )
    ax.plot([0, 1], [0, 1], ls="--", color="#999", lw=0.8)
    ax.set_xlabel("wm_score")
    ax.set_ylabel("robust_score (after paraphrase probe)")
    ax.set_title("Score vs robustness — wash shows as drop below diagonal", loc="left", fontweight="bold")
    ax.legend(fontsize=7, loc="upper left")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.3)

    # --- strength stacked by kind ---
    ax = fig.add_subplot(gs[1, 0])
    x = np.arange(len(KINDS))
    bottom = np.zeros(len(KINDS))
    strength_colors = {
        "strong_pos": "#2F6B4F",
        "gray": "#C4892A",
        "strong_neg": "#A33B24",
        "missing": "#9AA3AE",
    }
    for s in STRENGTHS:
        vals = np.array([summary["by_kind"][k]["strength"].get(s, 0) for k in KINDS], dtype=float)
        ax.bar(x, vals, bottom=bottom, color=strength_colors[s], label=s, width=0.72)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(KINDS, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("count")
    ax.set_title("Strength mix by evidence kind", loc="left", fontweight="bold")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, axis="y", alpha=0.3)

    # --- action mix ---
    ax = fig.add_subplot(gs[1, 1])
    actions = sorted({r.action for r in rows})
    act_colors = {
        "provenance_short_circuit": "#2F6B4F",
        "vague_human_review": "#C4892A",
        "po_fn_chain": "#1F6F8B",
    }
    bottom = np.zeros(len(KINDS))
    for a in actions:
        vals = np.array(
            [summary["by_kind"][k]["actions"].get(a, 0) for k in KINDS], dtype=float
        )
        ax.bar(x, vals, bottom=bottom, color=act_colors.get(a, "#555"), label=a, width=0.72)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(KINDS, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("count")
    ax.set_title("Routed action by kind (MVP policy)", loc="left", fontweight="bold")
    ax.legend(fontsize=6.5, loc="upper right")
    ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle(
        "Watermark evidence MVP — six real-data kinds → strength → AS-DTA action",
        fontsize=12.5,
        fontweight="bold",
        y=0.995,
    )
    fig.text(
        0.5,
        0.005,
        "Same: sparse → strength enum → vague/PO+FN hooks.  Diff: provenance vs crypto vs passive vs policy vs absence.",
        ha="center",
        fontsize=8,
        style="italic",
        color="#444",
    )
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path


def run_mvp(n_per_kind: int = 80, seed: int = SEED, out_dir: Optional[Path] = None) -> dict:
    out_dir = Path(out_dir or (ROOT / "results" / "as_dta_watermark_mvp"))
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = simulate_live_batch(n_per_kind=n_per_kind, seed=seed)
    summary = summarize(rows)
    png = out_dir / "watermark_mvp_kinds.png"
    plot_mvp(rows, summary, png)
    payload = {
        "summary": summary,
        "rows": [asdict(r) for r in rows],
        "mvp_note": (
            "Six real-data watermark *forms* share one strength schema but differ in "
            "provenance power, wash robustness, and whether they may emit strong_pos. "
            "Gray→vague/human; neg/missing→PO+FN (never Y=human)."
        ),
    }
    (out_dir / "watermark_mvp.json").write_text(json.dumps(payload, indent=2))
    _write_readme(out_dir / "README.md", summary, png.name)
    return payload


def _write_readme(path: Path, summary: dict, png_name: str):
    lines = [
        "# Watermark evidence MVP — six real-data kinds",
        "",
        "## Same vs Diff",
        "",
        "### Same",
        "",
    ]
    for s in summary["same_diff"]["same"]:
        lines.append("- " + s)
    lines.extend(["", "### Diff", ""])
    for s in summary["same_diff"]["diff"]:
        lines.append("- " + s)
    neg = summary["same_diff"].get("negative_effects") or {}
    intent = summary["same_diff"].get("intent_routing") or {}
    if neg:
        lines.extend(["", "## Negative effects (uplift / policy)", ""])
        for k, v in neg.items():
            lines.append("- **%s**: %s" % (k, v))
    if intent:
        lines.extend(["", "## Same wash, three intents", ""])
        for k, v in intent.items():
            lines.append("- **%s** → %s" % (k, v))
        lines.append("")
        lines.append(
            "See also `docs/method/AS_DTA_watermark_uplift_side_effects.md`."
        )
    lines.extend(
        [
            "",
            "## Counts by kind",
            "",
            "| kind | n | strong_pos | gray | strong_neg | missing | vague_rate |",
            "| --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for k in KINDS:
        b = summary["by_kind"][k]
        st = b["strength"]
        lines.append(
            "| %s | %d | %d | %d | %d | %d | %.2f |"
            % (
                k,
                b["n"],
                st.get("strong_pos", 0),
                st.get("gray", 0),
                st.get("strong_neg", 0),
                st.get("missing", 0),
                b["vague_rate"],
            )
        )
    lines.extend(
        [
            "",
            f"![]({png_name})",
            "",
            "```bash",
            "python3 scripts/run_as_dta_watermark_mvp.py",
            "```",
            "",
        ]
    )
    path.write_text("\n".join(lines))
