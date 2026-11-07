"""Selection-time clever covariate on a CoT chain.

Prior e(W) is the softmax of the stationary scores over the children.
The judge's pick is the action A = 1. The TMLE adjustment on that action is

    H = (1 - e) / (e * (1 - e)) = 1 / e.

The first 12 episodes freeze the 95% quantile of H. A later step whose H
exceeds that quantile is selected by the stationary score instead of the judge.
Both policies score every child, so the node count does not change.
"""
from __future__ import annotations

import json
import math
import random
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))
import tot_agent_probe as probe

OUT = Path(__file__).resolve().parent / "results_clever_covariate.json"
N = 40
DRIFT_AT = 16
REF = 12
SEED = 2026


def softmax(xs: list[float]) -> np.ndarray:
    a = np.asarray(xs, dtype=float)
    a = a - np.max(a)
    z = np.exp(a)
    p = z / np.sum(z)
    return np.clip(p, 1e-3, 1 - 1e-3)


def clever_h(e: float) -> float:
    e = float(np.clip(e, 1e-3, 1 - 1e-3))
    return (1.0 - e) / (e * (1.0 - e))


def pick(scored, policy: str, href: list[float]):
    """scored items are (judge, stable, payload). Returns index, H, nodes."""
    judges = [s[0] for s in scored]
    stables = [s[1] for s in scored]
    e = softmax(stables)
    j = int(np.argmax(judges))
    h = clever_h(float(e[j]))
    nodes = len(scored)
    if policy == "stable":
        return int(np.argmax(stables)), h, nodes
    if policy == "clever" and len(href) >= REF:
        thr = float(np.quantile(href[:REF], 0.95))
        if h > thr + 1e-9:
            return int(np.argmax(stables)), h, nodes
    return j, h, nodes


def run_hotpot(policy: str):
    probe.CURRENT_MIX = 0.80
    probe.MIX_JITTER = 0.0
    bank = probe.hop_questions()
    href, succ, nodes, steps, switches = [], [], [], [], []
    for t in range(N):
        question, gold = bank[t % len(bank)]
        need = 2
        picked = frozenset()
        n_nodes = 0
        n_steps = 0
        switched = 0
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
                judge = probe.drifted_judge(calibrated, spurious, t, DRIFT_AT)
                cands.append((judge, stable, nxt))
            idx, h, n_nodes_step = pick(cands, policy, href)
            j_idx = int(np.argmax([c[0] for c in cands]))
            if idx != j_idx:
                switched += 1
            if len(href) < REF:
                href.append(h)
            picked = cands[idx][2]
            n_nodes += n_nodes_step
            n_steps += 1
        y = 1.0 if len(picked & gold) >= need else 0.0
        succ.append(y)
        nodes.append(n_nodes)
        steps.append(n_steps)
        switches.append(switched)
    return _pack("hotpot_twohop", succ, nodes, steps, switches, href)


def run_game24(policy: str):
    probe.CURRENT_MIX = 0.28
    probe.MIX_JITTER = 0.0
    puzzles = probe.make_game24(random.Random(SEED), N)
    href, succ, nodes, steps, switches = [], [], [], [], []
    for t, puzzle in enumerate(puzzles):
        state = puzzle
        n_nodes = 0
        n_steps = 0
        switched = 0
        y = 0.0
        for _ in range(3):
            kids = probe._g24_children(state)
            if not kids:
                break
            cands = []
            for nxt, done in kids:
                solv = 1.0 if probe._g24_solvable(nxt) else 0.0
                close = max(math.exp(-abs(x - 24) / 8) for x in nxt)
                stable = close if len(nxt) == 1 else 0.45 * close + 0.15
                spurious = sum(nxt) / (13 * len(nxt))
                calibrated = 0.85 * solv + 0.15 * stable
                judge = probe.drifted_judge(calibrated, spurious, t, DRIFT_AT)
                cands.append((judge, stable, nxt, done))
            idx, h, n_nodes_step = pick([(c[0], c[1], c[2]) for c in cands], policy, href)
            j_idx = int(np.argmax([c[0] for c in cands]))
            if idx != j_idx:
                switched += 1
            if len(href) < REF:
                href.append(h)
            state = cands[idx][2]
            n_nodes += n_nodes_step
            n_steps += 1
            if cands[idx][3]:
                y = 1.0
                break
        succ.append(y)
        nodes.append(n_nodes)
        steps.append(n_steps)
        switches.append(switched)
    return _pack("game24", succ, nodes, steps, switches, href)


def _pack(task, succ, nodes, steps, switches, href):
    pre = succ[:DRIFT_AT]
    post = succ[DRIFT_AT:]
    return {
        "task": task,
        "pre_success": float(np.mean(pre)),
        "post_success": float(np.mean(post)),
        "nodes_per_task": float(np.mean(nodes)),
        "steps_per_task": float(np.mean(steps)),
        "post_switch_rate": float(np.mean([s > 0 for s in switches[DRIFT_AT:]])),
        "h_ref_q95": float(np.quantile(href[:REF], 0.95)) if len(href) >= REF else None,
        "success_per_node": float(np.mean(post) / max(np.mean(nodes), 1e-8)),
    }


def value_row(judge, clever, v_task: float):
    d_succ = clever["post_success"] - judge["post_success"]
    n_post = N - DRIFT_AT
    revenue = d_succ * v_task * n_post
    # Same children are scored. Node delta is numerical noise, priced at 0.
    d_nodes = clever["nodes_per_task"] - judge["nodes_per_task"]
    kind = "挽损" if clever["post_success"] <= judge["pre_success"] + 1e-9 else "用增"
    gap_left = judge["pre_success"] - clever["post_success"]
    return {
        "delta_success": d_succ,
        "n_post": n_post,
        "v_task": v_task,
        "delta_revenue": revenue,
        "delta_nodes": d_nodes,
        "delta_cost": 0.0,
        "net_value": revenue,
        "kind": kind,
        "unrecovered_success": gap_left,
    }


def main():
    report = {"policies": {}, "value": {}, "kind_rule": "post_success at or below the pre-drift judge rate is 挽损. Above that rate is 用增."}
    for name, fn in (("hotpot", run_hotpot), ("game24", run_game24)):
        report["policies"][name] = {}
        for policy in ("judge", "stable", "clever"):
            print(name, policy)
            report["policies"][name][policy] = fn(policy)
            b = report["policies"][name][policy]
            print(f"  pre={b['pre_success']:.2f} post={b['post_success']:.2f} nodes={b['nodes_per_task']:.1f} switch={b['post_switch_rate']:.2f}")
    report["value"]["hotpot"] = value_row(
        report["policies"]["hotpot"]["judge"], report["policies"]["hotpot"]["clever"], v_task=0.001
    )
    report["value"]["game24"] = value_row(
        report["policies"]["game24"]["judge"], report["policies"]["game24"]["clever"], v_task=0.10
    )
    for k, v in report["value"].items():
        print(k, v["kind"], "dS", round(v["delta_success"], 3), "net", round(v["delta_revenue"], 4))
    OUT.write_text(json.dumps(report, indent=2))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
