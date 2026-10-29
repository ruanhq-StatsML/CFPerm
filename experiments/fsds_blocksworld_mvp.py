"""MVP on Blocksworld: same six layers, a different search.

Four blocks, beam 4, depth 4. The state encoding is four numbers in [0, 1]:
stack count, tallest stack, goal-prefix length, and whether D is clear.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import tot_agent_probe as probe
from fsds_tot_layers import DRIFT_AT, N, SEED, episode_row, layer_report

OUT = Path(__file__).resolve().parent / "results_fsds_blocksworld.json"


def embed(state) -> tuple[float, float, float, float]:
    goal = ("A", "B", "C", "D")
    pref = 0
    d_clear = 0.0
    for st in state:
        k = 0
        for a, b in zip(st, goal):
            if a != b:
                break
            k += 1
        pref = max(pref, k)
        if st and st[-1] == "D":
            d_clear = 1.0
    tall = max(len(st) for st in state)
    return (len(state) / 4, tall / 4, pref / 4, d_clear)


def rows():
    probe.CURRENT_MIX = 0.30
    probe.MIX_JITTER = 0.0
    dist = probe.bw_distances()
    rng = probe.random.Random(SEED)
    starts = [probe.bw_random_start(rng, dist) for _ in range(N)]

    def score_at(state, t):
        d = dist.get(state, 8)
        progress = 1.0 - d / 8
        stable = probe.bw_stable(state)
        spurious = probe.bw_spurious(state)
        calibrated = 0.8 * progress + 0.2 * stable
        judge = probe.drifted_judge(calibrated, spurious, t, DRIFT_AT)
        return probe.StepView(judge, judge, stable, spurious, progress)

    out = []
    prev_y, prev_score = 0.0, 0.0
    for t, start in enumerate(starts):
        from fsds_tot_layers import search_log

        final, trace, _ = search_log(
            start, lambda s: [(n, probe.bw_success(n)) for n in probe.bw_neighbors(s)],
            lambda s, t=t: score_at(s, t),
            beam=4, depth=4, rng=probe.random.Random(SEED + 4000 + t),
        )
        y = 1.0 if final is not None else 0.0
        raw = trace[-1]["state"] if trace else start
        if trace:
            trace[-1]["state"] = (0.0, 0.0, 0.0, 0.0)
        row = episode_row(trace, y, prev_y, prev_score)
        e0, e1, e2, e3 = embed(raw)
        row["e0"], row["e1"], row["e2"], row["e3"] = e0, e1, e2, e3
        row["t"] = t
        out.append(row)
        prev_y, prev_score = y, row["roll_mean"]
    return out


def main():
    block = layer_report(rows())
    block["dataset"] = "blocksworld"
    block["mix"] = 0.30
    block["beam"] = 4
    block["depth"] = 4
    OUT.write_text(json.dumps(block, indent=2))
    print("Y", round(block["pre_Y"], 3), "->", round(block["post_Y"], 3))
    for layer in block["layers"]:
        print(
            f"{layer['layer']:12} share={layer['sse_share']:.3f} best={layer['best_feature']} "
            f"sse={layer['sse_best']:.3f} p={layer['best_p']:.3f} shift={layer['mean_shift_max']:.3f} via {layer['shift_feature']}"
        )
    print("wrote", OUT)


if __name__ == "__main__":
    main()
