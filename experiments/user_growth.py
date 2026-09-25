"""User-growth episodes on the same continuous-time prompt.

One episode is one arriving user. Three offers are scored: habit, growth, click.
Y is 1 when the selected offer activates the user.
The habit offer is what the pre-drift judge selects. Its activation rate sits near one half,
so there is room above the reference.
After drift the judge follows clicks. The growth policy leaves the judge
only after finished episodes show a higher loss than the first 12, and then
selects the growth offer on the next user.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import tot_agent_probe as probe
import tot_continuous_prompt as prompt

OUT = Path(__file__).resolve().parent / "results_user_growth.json"
OFFERS = (
    ("habit", 0.50, 0.15, 0.90),
    ("growth", 0.67, 0.12, 0.55),
    ("click", 0.20, 0.98, 0.10),
)


def activates(name: str, t: int) -> bool:
    if name == "habit":
        return t % 2 == 0
    if name == "growth":
        return t % 3 != 0
    return t % 5 == 0


def collect(n: int, drift_at: int, policy: str):
    probe.CURRENT_MIX = 0.80
    probe.MIX_JITTER = 0.0
    depth = 1
    episodes = []
    losses = []
    ref = []
    lever_at = None

    def score_at(state, t, use_growth):
        if state is None:
            return probe.StepView(0.0, 0.0, 0.0, 0.0, 0.0)
        name, stable, spurious, progress = OFFERS[state]
        calibrated = 0.8 * progress + 0.2 * stable
        judge = probe.drifted_judge(calibrated, spurious, t, drift_at)
        steer = stable if use_growth else judge
        return probe.StepView(steer, judge, stable, spurious, progress)

    for t in range(n):
        use_growth = policy == "growth" and lever_at is not None and t >= lever_at

        def expand(_state, t=t):
            return [(i, activates(OFFERS[i][0], t)) for i in range(len(OFFERS))]

        trace = []
        final, _, n_children = probe.beam_search(
            None, expand, lambda state, t=t, use_growth=use_growth: score_at(state, t, use_growth),
            beam=1, depth=depth, rng=probe.random.Random(prompt.SEED + 7000 + t), trace=trace,
        )
        xs = [prompt.step_x(view, step, depth) for step, view in trace]
        y = 1.0 if final is not None else 0.0
        episodes.append((np.vstack(xs) if xs else np.zeros((0, len(prompt.NAMES))), y, int(n_children)))
        losses.append(1.0 - y)
        if t < prompt.REF_EPISODES:
            ref.append(losses[-1])
        elif policy == "growth" and lever_at is None and len(losses) >= 4:
            if float(np.mean(losses[-4:])) > float(np.mean(ref)) + 0.15:
                lever_at = t + 1
    return episodes, lever_at


def main():
    n, drift_at = 40, 16
    report = {}
    for policy in ("judge", "growth"):
        episodes, lever_at = collect(n, drift_at, policy)
        block = prompt.run_stream(f"user_growth_{policy}", episodes, drift_at)
        ys = [item[1] for item in episodes]
        pre = ys[:drift_at]
        post = ys[drift_at:]
        block["policy"] = policy
        block["lever_at"] = lever_at
        block["pre_activation"] = float(np.mean(pre))
        block["post_activation"] = float(np.mean(post))
        block["kind"] = "用增" if block["post_activation"] > block["pre_activation"] + 1e-9 else "挽损"
        block["extra_activations_vs_pre"] = float(np.sum(post) - block["pre_activation"] * len(post))
        report[policy] = block
    judge_post = report["judge"]["post_activation"]
    growth_post = report["growth"]["post_activation"]
    report["comparison"] = {
        "pre_activation": report["judge"]["pre_activation"],
        "judge_post": judge_post,
        "growth_post": growth_post,
        "extra_vs_drifted_judge": float(
            (growth_post - judge_post) * (n - drift_at)
        ),
        "growth_kind": report["growth"]["kind"],
        "judge_kind": report["judge"]["kind"],
    }
    slim = {}
    for key, block in report.items():
        if key == "comparison":
            slim[key] = block
            continue
        slim[key] = {
            "policy": block["policy"],
            "lever_at": block["lever_at"],
            "pre_activation": block["pre_activation"],
            "post_activation": block["post_activation"],
            "kind": block["kind"],
            "extra_activations_vs_pre": block["extra_activations_vs_pre"],
            "channels": block["roi"]["channels"],
            "yongzeng_channels": block["roi"]["yongzeng_channels"],
            "w_stack": block["w_stack"],
            "trail_mse": block["trail_mse"],
            "post_anomaly": block["online_anomaly"]["post_drift"],
        }
    OUT.write_text(json.dumps(slim, indent=2))
    print(json.dumps(slim["comparison"], indent=2))
    for key in ("judge", "growth"):
        b = slim[key]
        print(key, "lever", b["lever_at"], "pre", round(b["pre_activation"], 3), "post", round(b["post_activation"], 3), b["kind"], "extra", round(b["extra_activations_vs_pre"], 2))
        print(" ", [(c["name"], round(c["delta"], 4), c["kind"]) for c in b["channels"]])
    print("wrote", OUT)


if __name__ == "__main__":
    main()
