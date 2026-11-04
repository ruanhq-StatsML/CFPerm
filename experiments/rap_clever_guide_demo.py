"""Show the next-step prompt built from a finished outcome.

The example H=(1-0.3)/(0.3*0.7) is the call in the skill note.
The stream is the three-offer judge. The prompt at user t+1 uses user t's Y.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

from rap_clever_covariate_guide import RAPCleverCovariateSkill
from user_growth import collect

OUT = Path(__file__).resolve().parent / "results_rap_clever_guide.json"


def main():
    skill = RAPCleverCovariateSkill(threshold=0.5)
    demo = skill.invoke(
        state="Current reasoning state.",
        prior=0.3,
        outcome=1,
        history=["step1", "step2"],
    )
    assert demo["hint"] == "aggressive"
    assert abs(demo["H"] - (1 - 0.3) / (0.3 * 0.7)) < 1e-9

    episodes, _ = collect(40, 16, "judge")
    history = []
    rows = []
    shown = None
    for t, (X, y, _) in enumerate(episodes):
        if t == 0:
            continue
        prev_x, prev_y, _ = episodes[t - 1]
        res = skill.invoke(
            state=f"User {t}. Three offers: habit, growth, click. Y for this user is not known yet.",
            prior=float(prev_x[0, 1]),
            outcome=float(prev_y),
            history=history[-3:],
        )
        line = f"t={t-1} Y={prev_y:.0f} e={res['prior']:.3f} H={res['H']:.3f} {res['hint']}"
        history.append(line)
        rows.append({"for_user": t, "from_user": t - 1, "hint": res["hint"], "H": res["H"], "post": t - 1 >= 16})
        if shown is None and t - 1 >= 16 and res["hint"] == "conservative":
            shown = res["prompt"]
    post = [r for r in rows if r["post"]]
    report = {
        "demo_prior_0.3_outcome_1": {"H": demo["H"], "hint": demo["hint"]},
        "post_drift_hint_counts": {
            h: int(sum(r["hint"] == h for r in post)) for h in ("aggressive", "conservative", "keep")
        },
        "prompt_after_first_post_drift_failure": shown,
    }
    OUT.write_text(json.dumps(report, indent=2))
    print("demo", round(demo["H"], 3), demo["hint"])
    print("post hints", report["post_drift_hint_counts"])
    print("--- prompt ---")
    print(shown)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
