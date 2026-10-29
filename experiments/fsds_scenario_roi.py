"""FSDS feature selection mapped to three business scenes.

The measured objects are X, Y, the concept-drift SSE drop, and the action
taken on the next user. Unit prices are assumptions and are labeled as such.
AUC is not computed and is not an input.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import tot_agent_probe as probe

OUT = Path(__file__).resolve().parent / "results_fsds_scenarios.json"
SEED = 2026
N = 40
DRIFT_AT = 16
REF = 12
MIN_POST = 5

OFFERS = (
    ("habit", 0.50, 0.15, 0.90),
    ("growth", 0.67, 0.12, 0.55),
    ("click", 0.20, 0.98, 0.10),
    ("save", 0.50, 0.18, 0.88),
)

# Stated unit prices. Not estimated from these episodes.
PRICES = {
    "incentive": {"value": 1.0, "unit_cost": 0.20, "value_name": "V_user", "cost_name": "incentive"},
    "churn": {"value": 1.0, "unit_cost": 0.50, "value_name": "LTV", "cost_name": "retention"},
    "payment": {"value": 1.0, "unit_cost": 0.40, "value_name": "AOV", "cost_name": "subsidy"},
}


def activates(scene: str, name: str, t: int) -> bool:
    if scene == "incentive":
        return {"habit": t % 2 == 0, "growth": t % 3 != 0, "click": t % 5 == 0}[name]
    if scene == "churn":
        # Save restores the habit retention rate. It does not clear that rate.
        return {"habit": t % 2 == 0, "growth": t % 3 != 0, "click": t % 5 == 0, "save": t % 2 == 0}[name]
    # Paid order. Click still takes the subsidy and rarely orders.
    return {"habit": t % 2 == 0, "growth": t % 3 != 0, "click": t % 5 == 0, "save": t % 2 == 0}[name]


def pick(scene: str, t: int, use_alt: bool) -> int:
    probe.CURRENT_MIX = 0.80
    probe.MIX_JITTER = 0.0
    names = [o[0] for o in OFFERS if scene != "incentive" or o[0] != "save"]
    if scene == "incentive":
        pool = [o for o in OFFERS if o[0] != "save"]
    elif scene == "churn":
        pool = [o for o in OFFERS if o[0] != "growth"]
    else:
        pool = [o for o in OFFERS if o[0] != "save"]
    scored = []
    for i, (name, stable, spurious, progress) in enumerate(pool):
        calibrated = 0.8 * progress + 0.2 * stable
        judge = probe.drifted_judge(calibrated, spurious, t, DRIFT_AT)
        if use_alt and scene == "incentive":
            key = stable
        elif use_alt and scene == "churn":
            key = 1.0 if name == "save" else 0.0
        elif use_alt and scene == "payment":
            key = stable
        else:
            key = judge
        scored.append((key, i, name, judge, stable, spurious, progress))
    scored.sort(key=lambda z: z[0], reverse=True)
    return scored[0]


def run_scene(scene: str, round_id: int) -> dict:
    """round 1 switches on the largest SSE drop. round 2 also asks for a mean increase in that feature."""
    feats = []
    rows = []
    switched_at = None
    for t in range(N):
        use_alt = False
        concept, covariate = [], []
        if len(feats) >= REF + MIN_POST:
            concept = probe.localize_fsds(feats, REF, n_perm=60, seed=SEED)
            covariate = probe.localize(feats, REF, n_perm=60, seed=SEED)
            top = concept[0]
            shift = {r["feature"]: r for r in covariate}
            moved = shift.get(top["feature"], {})
            mean_up = moved.get("abs_mean_shift", 0.0) > 0.05
            if round_id == 1:
                use_alt = top["feature"] == "spurious" and top["sse_drop"] > 0
            else:
                use_alt = top["feature"] == "spurious" and top["sse_drop"] > 0 and mean_up and top["p"] <= 0.10
            if use_alt and switched_at is None:
                switched_at = t
        key, idx, name, judge, stable, spurious, progress = pick(scene, t, use_alt)
        y = 1.0 if activates(scene, name, t) else 0.0
        paid = use_alt or scene == "payment"
        if scene == "payment":
            paid = True
        cost = PRICES[scene]["unit_cost"] if paid else 0.0
        row = {
            "t": t,
            "offer": name,
            "y": y,
            "judge": float(judge),
            "stable": float(stable),
            "spurious": float(spurious),
            "progress": float(progress),
            "alt": bool(use_alt),
            "cost": float(cost),
            "children": 3,
        }
        rows.append(row)
        feats.append({k: row[k] for k in ("judge", "stable", "spurious", "progress", "y")})
    concept = probe.localize_fsds(feats, REF, n_perm=80, seed=SEED) if len(feats) >= REF + MIN_POST else []
    covariate = probe.localize(feats, REF, n_perm=80, seed=SEED) if len(feats) >= REF + MIN_POST else []
    # Attribution on the judge-only prefix is the interpretation of the drift.
    # Recompute it from a frozen judge replay so the action does not rewrite X.
    return {
        "scene": scene,
        "round": round_id,
        "switched_at": switched_at,
        "rows": rows,
        "concept_on_acted_stream": concept,
        "covariate_on_acted_stream": covariate,
    }


def judge_replay(scene: str) -> list[dict]:
    feats = []
    for t in range(N):
        _, _, name, judge, stable, spurious, progress = pick(scene, t, False)
        y = 1.0 if activates(scene, name, t) else 0.0
        feats.append({
            "t": t, "offer": name, "y": y,
            "judge": float(judge), "stable": float(stable),
            "spurious": float(spurious), "progress": float(progress),
        })
    return feats


def summarize(scene: str, acted: dict, judge_rows: list[dict]) -> dict:
    price = PRICES[scene]
    pre_j = [r for r in judge_rows if r["t"] < DRIFT_AT]
    post_j = [r for r in judge_rows if r["t"] >= DRIFT_AT]
    post_a = [r for r in acted["rows"] if r["t"] >= DRIFT_AT]
    yj = float(np.mean([r["y"] for r in post_j]))
    ya = float(np.mean([r["y"] for r in post_a]))
    ypre = float(np.mean([r["y"] for r in pre_j]))
    extra = float(np.sum([r["y"] for r in post_a]) - np.sum([r["y"] for r in post_j]))
    cost = float(np.sum([r["cost"] for r in post_a]))
    cost_j = float(np.sum([price["unit_cost"] if scene == "payment" else 0.0 for r in post_j]))
    revenue = extra * price["value"]
    waste_j = float(np.sum([price["unit_cost"] for r in post_j if r["y"] < 0.5 and r["offer"] == "click"]))
    waste_a = float(np.sum([r["cost"] for r in post_a if r["y"] < 0.5 and r["offer"] == "click"]))
    saved = max(waste_j - waste_a, 0.0) if scene == "payment" else 0.0
    denom = cost if cost > 1e-9 else None
    roi = None if denom is None else (revenue + saved) / denom
    kind = "increment" if ya > ypre + 1e-9 else "shortfall"
    concept = probe.localize_fsds(judge_rows, REF, n_perm=80, seed=SEED)
    covariate = probe.localize(judge_rows, REF, n_perm=80, seed=SEED)
    top = concept[0] if concept else None
    x_move = []
    for key in ("judge", "stable", "spurious", "progress"):
        a = float(np.mean([r[key] for r in pre_j]))
        b = float(np.mean([r[key] for r in post_j]))
        x_move.append({"feature": key, "pre": a, "post": b, "delta": b - a})
    if top and top["feature"] == "spurious":
        if scene == "incentive":
            action = "Send the growth offer on the next user. Leave the click offer unranked."
        elif scene == "churn":
            action = "Send the save offer. Do not buy the growth offer. The target is the pre-drift retention rate."
        else:
            action = "Stop subsidizing the click offer. Subsidize the growth offer instead."
    else:
        action = "The concept-drift column is not spurious. Do not spend the scenario budget."
    return {
        "scene": scene,
        "round": acted["round"],
        "Y": "1 if the selected offer produces the scene outcome on that user, else 0.",
        "X": ["judge", "stable", "spurious", "progress"],
        "pre_Y": ypre,
        "judge_post_Y": yj,
        "acted_post_Y": ya,
        "kind": kind,
        "extra_outcomes": extra,
        "post_cost": cost,
        "judge_post_cost": cost_j,
        "fraud_cost_avoided": saved,
        "revenue_at_stated_price": revenue,
        "roi_at_stated_price": roi,
        "prices_are_assumptions": price,
        "children_per_user": 3,
        "switched_at": acted["switched_at"],
        "top_concept_on_judge_stream": top,
        "concept_on_judge_stream": concept,
        "covariate_on_judge_stream": covariate,
        "x_movement_on_judge_stream": x_move,
        "action": action,
        "n_alt": int(sum(r["alt"] for r in acted["rows"])),
    }


def main():
    report = {"rounds": []}
    for round_id in (1, 2):
        block = {"round": round_id, "scenes": {}}
        for scene in ("incentive", "churn", "payment"):
            judge_rows = judge_replay(scene)
            acted = run_scene(scene, round_id)
            block["scenes"][scene] = summarize(scene, acted, judge_rows)
        report["rounds"].append(block)
        print("round", round_id)
        for scene, s in block["scenes"].items():
            top = s["top_concept_on_judge_stream"]
            print(
                scene,
                "top", None if top is None else (top["feature"], round(top["sse_drop"], 3), round(top["p"], 3)),
                "switch", s["switched_at"],
                "Y", round(s["pre_Y"], 3), round(s["judge_post_Y"], 3), round(s["acted_post_Y"], 3),
                s["kind"],
                "extra", s["extra_outcomes"],
                "cost", round(s["post_cost"], 2),
                "roi", None if s["roi_at_stated_price"] is None else round(s["roi_at_stated_price"], 3),
            )
            print(" ", s["action"])
    OUT.write_text(json.dumps(report, indent=2))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
