"""Multi-suite B1–B6 comparison plus LLM-boost table export."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Sequence

from .chronoberg import SyntheticChronoBergConfig
from .experiment import ExperimentResult, run_baseline_comparison, save_result
from .llm_boost import results_to_latex, run_llm_boost
from .online import AGODConfig
from .suites import SUITE_NAMES, make_suite

COMPUTE_BASELINES: Sequence[str] = ("B1", "B2", "B3", "B4", "B5", "B6")


def run_multi_suite(
    *,
    suites: Sequence[str] = SUITE_NAMES,
    n_per_window: int = 56,
    seed: int = 2026,
    dim: int = 12,
    steps: int = 5,
    pretrain: int = 16,
    baselines: Sequence[str] = COMPUTE_BASELINES,
) -> Dict[str, ExperimentResult]:
    tcfg = AGODConfig(
        seed=seed,
        steps_per_window=steps,
        pretrain_steps=pretrain,
        teacher_dim=8,
        tau=0.32,
        gamma=1.6,
        momentum=0.15,
    )
    out: Dict[str, ExperimentResult] = {}
    for name in suites:
        def factory(suite=name):
            return make_suite(suite, n_per_window=n_per_window, seed=seed, dim=dim)

        # ChronoBerg uses its own dim mapping via SyntheticChronoBergConfig
        if name == "chronoberg":
            cfg = SyntheticChronoBergConfig(
                n_per_window=n_per_window,
                dims={"audio": dim, "image": dim, "text": dim},
                seed=seed,
            )
            out[name] = run_baseline_comparison(
                config=cfg, trainer_config=tcfg, baselines=baselines
            )
        else:
            out[name] = run_baseline_comparison(
                trainer_config=tcfg, baselines=baselines, stream_factory=factory
            )
    return out


def multi_to_latex(results: Dict[str, ExperimentResult]) -> str:
    header = "Suite & Method & Drift Rec@5 & Align & Rel. FLOPs & Skip \\\\"
    lines = [header, "\\midrule"]
    for suite, exp in results.items():
        first = True
        for method, m in exp.metrics.items():
            suite_cell = suite if first else ""
            first = False
            lines.append(
                f"{suite_cell} & {method} & {m['drift_subgroup_recall']:.3f} & "
                f"{m['audio_alignment_on_drift']:.3f} & {m.get('rel_flops', 1.0):.2f} & "
                f"{m.get('skip_rate', 0.0):.2f} \\\\"
            )
        lines.append("\\midrule")
    if lines[-1] == "\\midrule":
        lines.pop()
    body = "\n".join(lines)
    return f"""\\begin{{table}}[t]
\\centering
\\small
\\caption{{AGOD vs static (B1), covariate-only (B2), sparse (B4), budgeted (B5) and cascade (B6) distillation on three shift suites. Rel.\\ FLOPs is distill MAC-count vs.\\ B1; Skip is the fraction of windows where cascade idles the student.}}
\\label{{tab:agod-suites}}
\\begin{{tabular}}{{llcccc}}
\\toprule
{body}
\\bottomrule
\\end{{tabular}}
\\end{{table}}
"""


def write_reports(out_dir: Path, *, quick: bool = False) -> Dict[str, Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    n = 40 if quick else 56
    steps = 3 if quick else 5
    pretrain = 8 if quick else 16
    multi = run_multi_suite(n_per_window=n, steps=steps, pretrain=pretrain)
    payload = {k: {"metrics": v.metrics, "ranking": v.ranking} for k, v in multi.items()}
    (out_dir / "multi_suite.json").write_text(json.dumps(payload, indent=2, default=float))
    table1 = multi_to_latex(multi)
    (out_dir / "tab_agod_suites.tex").write_text(table1)

    boost, diag = run_llm_boost(n_id=64 if quick else 96, n_ood=64 if quick else 96)
    table2 = results_to_latex(boost)
    (out_dir / "tab_llm_boost.tex").write_text(table2)
    (out_dir / "llm_boost.json").write_text(
        json.dumps(
            {
                k: vars(v) for k, v in boost.items()
            } | {"diagnostics": diag},
            indent=2,
        )
    )
    return {
        "multi": out_dir / "multi_suite.json",
        "tab_suites": out_dir / "tab_agod_suites.tex",
        "tab_boost": out_dir / "tab_llm_boost.tex",
        "boost": out_dir / "llm_boost.json",
    }
