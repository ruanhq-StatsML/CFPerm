"""Walk the same chains until the path stops repeating.

Eight questions (or eight puzzles) are reused in order. The first visit of each
item is the reference path. A later visit agrees when the committed path is the
same. The chain is unstable at the first episode where the last eight visits
agree at most half the time.

Graph features on the committed step: depth, branch_n, sibling_gap
(judge of the chosen child minus judge of the runner-up), judge, stable,
spurious, progress. FSDS is fit once on the finished walks, split at the
detected break, and only reports which of those features carries Y.
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

OUT = Path(__file__).resolve().parent / "results_path_stability.json"
N = 100
DRIFT_AT = 16
WINDOW = 8
SEED = 2026


def _unstable_at(agree: list[float], window: int = WINDOW) -> int | None:
    for t in range(window - 1, len(agree)):
        if float(np.mean(agree[t - window + 1 : t + 1])) <= 0.5:
            return t
    return None


def walk_hotpot(t: int, question: str, gold: set[int]):
    need = 2
    picked = frozenset()
    path = []
    last = None
    graph = None
    for depth in range(2):
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
            cands.append((judge, stable, spurious, progress, pid, nxt))
        cands.sort(key=lambda z: z[0], reverse=True)
        judge, stable, spurious, progress, pid, nxt = cands[0]
        runner = cands[1][0] if len(cands) > 1 else judge
        picked = nxt
        path.append(pid)
        last = (judge, stable, spurious, progress)
        graph = {
            "depth": depth + 1,
            "branch_n": len(cands),
            "sibling_gap": float(judge - runner),
            "judge": float(judge),
            "stable": float(stable),
            "spurious": float(spurious),
            "progress": float(progress),
        }
    y = 1.0 if len(picked & gold) >= need else 0.0
    return tuple(path), y, graph


def walk_game24(t: int, puzzle):
    state = puzzle
    path = []
    graph = None
    for depth in range(3):
        kids = probe._g24_children(state)
        if not kids:
            break
        scored = []
        for nxt, _done in kids:
            solv = 1.0 if probe._g24_solvable(nxt) else 0.0
            close = max(probe.math.exp(-abs(x - 24) / 8) for x in nxt)
            stable = close if len(nxt) == 1 else 0.45 * close + 0.15
            spurious = sum(nxt) / (13 * len(nxt))
            calibrated = 0.85 * solv + 0.15 * stable
            judge = probe.drifted_judge(calibrated, spurious, t, DRIFT_AT)
            scored.append((judge, stable, spurious, solv, nxt))
        scored.sort(key=lambda z: z[0], reverse=True)
        judge, stable, spurious, solv, nxt = scored[0]
        runner = scored[1][0] if len(scored) > 1 else judge
        state = nxt
        path.append(tuple(round(x, 5) for x in nxt))
        graph = {
            "depth": depth + 1,
            "branch_n": len(scored),
            "sibling_gap": float(judge - runner),
            "judge": float(judge),
            "stable": float(stable),
            "spurious": float(spurious),
            "progress": float(solv),
        }
        if len(nxt) == 1 and abs(nxt[0] - 24) < 1e-4:
            return tuple(path), 1.0, graph
    return tuple(path), 0.0, graph


def run(name: str, items, walker):
    probe.CURRENT_MIX = 0.80 if name == "hotpot" else 0.28
    probe.MIX_JITTER = 0.0
    ref_path = {}
    agree, ys, graphs = [], [], []
    for t in range(N):
        item = items[t % len(items)]
        path, y, graph = walker(t, *item) if name == "hotpot" else walker(t, item)
        key = t % len(items)
        if key not in ref_path:
            ref_path[key] = path
            agree.append(1.0)
        else:
            agree.append(1.0 if path == ref_path[key] else 0.0)
        ys.append(y)
        graphs.append(graph)
    break_at = _unstable_at(agree)
    feats = []
    for y, g in zip(ys, graphs):
        row = dict(g)
        row["y"] = y
        feats.append(row)
    split = break_at if break_at is not None and 5 <= break_at <= N - 5 else DRIFT_AT
    concept = probe.localize_fsds(feats, split, n_perm=60, seed=SEED)
    roll = []
    for t in range(WINDOW - 1, N):
        roll.append({"t": t, "agree": float(np.mean(agree[t - WINDOW + 1 : t + 1]))})
    return {
        "task": name,
        "walks": N,
        "items": len(items),
        "unstable_at": break_at,
        "agree_before_break": None if break_at is None else float(np.mean(agree[:break_at])),
        "agree_from_break": None if break_at is None else float(np.mean(agree[break_at:])),
        "success_before_16": float(np.mean(ys[:DRIFT_AT])),
        "success_from_16": float(np.mean(ys[DRIFT_AT:])),
        "concept_at_break": concept[:4],
        "rolling_agree": roll[::4],
    }


def main():
    bank = [(q, g) for q, g in probe.hop_questions()]
    puzzles = probe.make_game24(random.Random(SEED), 8)
    report = {
        "rule": "unstable when the mean path-agreement of the last 8 walks is at most 0.5. The reference path is the first walk of that same item.",
        "hotpot": run("hotpot", bank, walk_hotpot),
        "game24": run("game24", puzzles, walk_game24),
    }
    OUT.write_text(json.dumps(report, indent=2))
    for key in ("hotpot", "game24"):
        block = report[key]
        print(key, "unstable_at", block["unstable_at"], "agree_from", block["agree_from_break"], "post16", round(block["success_from_16"], 2))
        if block["concept_at_break"]:
            top = block["concept_at_break"][0]
            print(" ", top["feature"], "sse", round(top["sse_drop"], 4), "p", round(top["p"], 3))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
