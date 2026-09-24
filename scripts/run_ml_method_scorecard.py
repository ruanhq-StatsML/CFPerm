#!/usr/bin/env python3
"""ML method ops scorecard (partial excess / learning curve / Page–Hinkley).

Composes with transfer_null; does not rewrite cores.

  PYTHONPATH=. python3 scripts/run_ml_method_scorecard.py \\
    --out results/agod_ml_method_ops
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from agod.ml_method_ops import (
    excess_learning_curve,
    ml_method_suite,
    page_hinkley_skill,
    partial_excess_auc,
)

ROOT = Path(__file__).resolve().parents[1]


def _clean(o):
    if isinstance(o, float) and (o != o or o in (float("inf"), float("-inf"))):
        return None
    if isinstance(o, (np.floating,)):
        v = float(o)
        return None if v != v else v
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, list):
        return [_clean(v) for v in o]
    return o


def _synth(seed: int = 0) -> dict:
    rng = np.random.default_rng(seed)
    n = 800
    # Case A: volume-driven labels + volume score (false skill)
    c = rng.integers(1, 40, size=n).astype(float)
    p = 1.0 / (1.0 + np.exp(-(0.12 * c - 2.5)))
    y_vol = (rng.random(n) < p).astype(int)
    y_vol[0], y_vol[1] = 0, 1
    score_vol = c + rng.normal(scale=0.8, size=n)
    partial_vol = partial_excess_auc(y_vol, score_vol, c, n_perm=8, seed=seed)

    # Case B: real signal orthogonal to confounder
    signal = rng.normal(size=n)
    p2 = 1.0 / (1.0 + np.exp(-(1.2 * signal)))
    y_real = (rng.random(n) < p2).astype(int)
    y_real[0], y_real[1] = 0, 1
    score_real = signal + 0.05 * c + rng.normal(scale=0.4, size=n)
    partial_real = partial_excess_auc(y_real, score_real, c, n_perm=8, seed=seed + 1)

    curve = excess_learning_curve(y_real, score_real, n_perm=5, seed=seed)
    # skill stream: healthy then drop
    hist = [0.18, 0.17, 0.19, 0.18, 0.17, 0.16] + [0.03, 0.02, 0.04, 0.02, 0.03, 0.02]
    ph = page_hinkley_skill(hist, delta=0.002, lambda_thresh=0.04)
    suite = ml_method_suite(
        y_real,
        score_real,
        confounder=c,
        excess_history=hist[:-1],
        n_perm=5,
        seed=seed,
    )
    return {
        "partial_volume_driven": partial_vol,
        "partial_real_signal": partial_real,
        "learning_curve": curve,
        "page_hinkley": ph,
        "suite": suite,
        "effect_reading": {
            "volume_delta_excess": partial_vol.get("delta_excess"),
            "real_delta_excess": partial_real.get("delta_excess"),
            "curve_reading": curve.get("reading"),
            "ph_reading": ph.get("reading"),
            "headline": (
                "volume case: large Δexcess (skill collapses after partialling); "
                "real case: Δexcess small; PH alarms on injected skill drop; "
                "cores untouched."
            ),
        },
    }


def render_md(rep: dict) -> str:
    pv, pr = rep["partial_volume_driven"], rep["partial_real_signal"]
    curve, ph = rep["learning_curve"], rep["page_hinkley"]
    eff = rep["effect_reading"]

    def f(x, nd=3):
        try:
            v = float(x)
            if v != v:
                return "—"
            return f"{v:.{nd}g}"
        except Exception:
            return "—"

    lines = [
        "# ML method ops scorecard (20-min RSI)",
        "",
        f"**Headline:** {eff.get('headline')}",
        "",
        "Composes with `transfer_null` excess / null — **no core rewrite**.",
        "",
        "## Partial excess (confounder partialling)",
        "",
        "| case | excess_raw | excess_partial | Δexcess | confounder R² |",
        "|---|---:|---:|---:|---:|",
        f"| volume-driven score | {f(pv.get('excess_raw'))} | {f(pv.get('excess_partial'))} | "
        f"{f(pv.get('delta_excess'))} | {f(pv.get('confounder_r2'))} |",
        f"| real signal ⊥ volume | {f(pr.get('excess_raw'))} | {f(pr.get('excess_partial'))} | "
        f"{f(pr.get('delta_excess'))} | {f(pr.get('confounder_r2'))} |",
        "",
        "## Excess learning curve (real-signal case)",
        "",
        f"- reading: {curve.get('reading')}",
        f"- excess@20% → full: {f(curve.get('excess_at_20pct'))} → {f(curve.get('excess_at_full'))}",
        f"- gain: {f(curve.get('curve_gain_20_to_full'))}",
        "",
        "## Page–Hinkley on excess stream",
        "",
        f"- alarm: **{ph.get('alarm')}** (index={ph.get('alarm_index')})",
        f"- reading: {ph.get('reading')}",
        "",
        "## Effect (what improved)",
        "",
        "1. False skill from volume is **quantified** (Δexcess), not hand-waved.",
        "2. Sample hunger of excess estimate is visible (curve).",
        "3. Skill-drop is detectable online (PH) — dual to data-drift gates.",
        "",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=ROOT / "results/agod_ml_method_ops")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rep = _synth(args.seed)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "ml_method_scorecard.json").write_text(
        json.dumps(_clean(rep), indent=2) + "\n"
    )
    md = render_md(rep)
    (args.out / "ML_METHOD_SCORECARD.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
