"""CoT next-step scores from FSDS feature weights.

One chain, one committed step at a time. No beam.
After each episode the outcome Y is known. FSDS asks which step feature
carries the success gap between the frozen reference (first 12 episodes)
and the episodes since: y ~ 1 + batch versus y ~ 1 + batch + x + batch:x.
The SSE drop is that feature's concept-drift share.

Those shares become the next episode's weights. Mass on judge and spurious
moves the chain onto the stationary heuristic. Mass on stable and progress
leaves the chain on the judge. Spurious and progress are attribution only;
the chain score uses judge and stable, so solvability is not a steering cue.

A second modality, the covariate mean shift, is reported beside the concept
share and does not set the weight.
"""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import tot_agent_probe as probe

OUT = Path(__file__).resolve().parent / "results_cot_fsds.json"
REF = 12
MIN_POST = 5
SEED = 2026


def concept_weights(feats: list[dict]) -> tuple[dict, list[dict], list[dict]]:
    """Weights for the next episode from past episodes only."""
    concept = probe.localize_fsds(feats, REF, n_perm=80, seed=SEED)
    covariate = probe.localize(feats, REF, n_perm=80, seed=SEED)
    drops = {r["feature"]: max(float(r["sse_drop"]), 0.0) for r in concept}
    drift_mass = drops.get("judge", 0.0) + drops.get("spurious", 0.0)
    keep_mass = drops.get("stable", 0.0) + drops.get("progress", 0.0)
    z = drift_mass + keep_mass
    if z < 1e-8:
        w_stable = 0.0
    else:
        w_stable = drift_mass / z
    weights = {
        "judge": float(1.0 - w_stable),
        "stable": float(w_stable),
        "spurious": 0.0,
        "progress": 0.0,
    }
    return weights, concept, covariate


def chain_score(view, weights: dict) -> float:
    return weights["judge"] * float(view.judge) + weights["stable"] * float(view.stable)


def run_game24(n: int, drift_at: int, policy: str):
    probe.CURRENT_MIX = 0.28
    probe.MIX_JITTER = 0.0
    puzzles = probe.make_game24(random.Random(SEED), n)
    feats, successes = [], []
    chosen_feat, weight_trace = [], []

    def view_of(t, state):
        solv = 1.0 if probe._g24_solvable(state) else 0.0
        close = max(probe.math.exp(-abs(x - 24) / 8) for x in state)
        stable = close if len(state) == 1 else 0.45 * close + 0.15
        spurious = sum(state) / (13 * len(state))
        calibrated = 0.85 * solv + 0.15 * stable
        judge = probe.drifted_judge(calibrated, spurious, t, drift_at)
        return probe.StepView(judge, judge, stable, spurious, solv)

    for t, puzzle in enumerate(puzzles):
        if policy == "fsds" and len(feats) >= REF + MIN_POST:
            weights, concept, covariate = concept_weights(feats)
        else:
            weights = {"judge": 1.0, "stable": 0.0, "spurious": 0.0, "progress": 0.0}
            concept, covariate = [], []
        state = puzzle
        last = view_of(t, state)
        solved = 0.0
        for _ in range(3):
            kids = probe._g24_children(state)
            if not kids:
                break
            scored = []
            for nxt, done in kids:
                view = view_of(t, nxt)
                scored.append((chain_score(view, weights), nxt, done, view))
            scored.sort(key=lambda z: z[0], reverse=True)
            _, state, done, last = scored[0]
            if done:
                solved = 1.0
                break
        successes.append(solved)
        feats.append({
            "judge": last.judge,
            "stable": last.stable,
            "spurious": last.spurious,
            "progress": last.progress,
            "y": solved,
        })
        top_c = concept[0]["feature"] if concept else None
        top_v = covariate[0]["feature"] if covariate else None
        chosen_feat.append(top_c)
        weight_trace.append({
            "t": t,
            "w_judge": weights["judge"],
            "w_stable": weights["stable"],
            "concept_top": top_c,
            "covariate_top": top_v,
            "y": solved,
        })
    return _pack("game24", successes, weight_trace, drift_at)


def run_hotpot(n: int, drift_at: int, policy: str):
    probe.CURRENT_MIX = 0.80
    probe.MIX_JITTER = 0.0
    bank = probe.hop_questions()
    feats, successes = [], []
    weight_trace = []

    for t in range(n):
        question, gold = bank[t % len(bank)]
        need = 2
        if policy == "fsds" and len(feats) >= REF + MIN_POST:
            weights, concept, covariate = concept_weights(feats)
        else:
            weights = {"judge": 1.0, "stable": 0.0, "spurious": 0.0, "progress": 0.0}
            concept, covariate = [], []
        picked = frozenset()
        last = None
        for _ in range(2):
            cands = []
            for pid in range(len(probe.PASSAGES)):
                if pid in picked:
                    continue
                nxt = frozenset(set(picked) | {pid})
                progress = len(set(nxt) & gold) / need
                stable = float(np.mean([probe.hop_overlap(question, i) for i in nxt]))
                spurious = float(np.mean([len(probe.PASSAGES[i].split()) / 12 for i in nxt]))
                calibrated = 0.75 * progress + 0.25 * stable
                judge = probe.drifted_judge(calibrated, spurious, t, drift_at)
                view = probe.StepView(judge, judge, stable, spurious, progress)
                cands.append((chain_score(view, weights), nxt, view))
            cands.sort(key=lambda z: z[0], reverse=True)
            _, picked, last = cands[0]
        y = 1.0 if len(picked & gold) >= need else 0.0
        successes.append(y)
        feats.append({
            "judge": last.judge,
            "stable": last.stable,
            "spurious": last.spurious,
            "progress": last.progress,
            "y": y,
        })
        weight_trace.append({
            "t": t,
            "w_judge": weights["judge"],
            "w_stable": weights["stable"],
            "concept_top": concept[0]["feature"] if concept else None,
            "covariate_top": covariate[0]["feature"] if covariate else None,
            "concept": concept,
            "covariate": covariate,
            "y": y,
        })
    return _pack("hotpot_twohop", successes, weight_trace, drift_at)


def _pack(task, successes, trace, drift_at):
    post = [row for row in trace if row["t"] >= drift_at]
    guided = [row for row in post if row["t"] >= REF + MIN_POST]
    def mean_w(key):
        if not guided:
            return None
        return float(np.mean([row[key] for row in guided]))
    tops = {}
    for row in guided:
        name = row["concept_top"] or "none"
        tops[name] = tops.get(name, 0) + 1
    return {
        "task": task,
        "pre_success": float(np.mean(successes[:drift_at])),
        "post_success": float(np.mean(successes[drift_at:])),
        "post_w_judge": mean_w("w_judge"),
        "post_w_stable": mean_w("w_stable"),
        "concept_top_counts": tops,
        "last": {k: trace[-1][k] for k in ("t", "w_judge", "w_stable", "concept_top", "covariate_top", "y")},
        "trace_tail": [
            {k: row[k] for k in ("t", "w_judge", "w_stable", "concept_top", "covariate_top", "y")}
            for row in trace[-6:]
        ],
    }


def main():
    n, drift_at = 40, 16
    report = {"ref_episodes": REF, "min_post": MIN_POST, "tasks": {}}
    for task, runner in (("game24", run_game24), ("hotpot_twohop", run_hotpot)):
        report["tasks"][task] = {}
        for policy in ("judge", "fsds"):
            print(f"running {task} / {policy}")
            out = runner(n, drift_at, policy)
            report["tasks"][task][policy] = out
            print(
                f"  pre={out['pre_success']:.2f} post={out['post_success']:.2f}"
                f" w_stable={out['post_w_stable']} concept={out['concept_top_counts']}"
            )
    OUT.write_text(json.dumps(report, indent=2))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
