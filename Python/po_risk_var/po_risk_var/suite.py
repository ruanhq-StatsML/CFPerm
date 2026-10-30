"""Run synthetic + real-world justify suite; write results + plot."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import numpy as np

from po_risk_var._path import REPO_ROOT
from po_risk_var.apps_realworld import run_realworld_apps
from po_risk_var.apps_synthetic import run_synthetic_apps
from po_risk_var.product import estimator_priority

OUT_DEFAULT = REPO_ROOT / "results" / "as_dta_po_risk_var_apps"


def run_justify_suite(
    *,
    seed: int = 2026,
    out_dir: Optional[Path] = None,
    quick: bool = True,
    method: str = "cf",
) -> dict:
    out_dir = Path(out_dir or OUT_DEFAULT)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("[suite] synthetic …", flush=True)
    synth = run_synthetic_apps(seed=seed, quick=quick, method=method)
    print("[suite] real-world …", flush=True)
    real = run_realworld_apps(seed=seed, quick=quick, method=method)

    payload = {
        "scope": (
            "Product landing justify: CF-debiased Var(u) + Mag×Var triage "
            "on synthetic DGPs and real-world text boards."
        ),
        "product_logic": {
            "steps": [
                "fit φ, τ̂",
                "u=(φ−τ̂)²",
                "Mag×Var (μ, s)",
                "Var(u) CF-debiased (default)",
                "Y-flip if text integrity",
                "actuator via triage()",
                "honesty: bias_ratio",
            ],
            "estimator_priority": estimator_priority(),
        },
        "synthetic": synth,
        "realworld": real,
        "headline": {
            "synth_hit_rate": synth["hit_rate"],
            "real_hit_rate": real["hit_rate"],
            "synth_n": f"{synth['n_correct']}/{synth['n_boards']}",
            "real_n": f"{real['n_pair_ok']}/{real['n_boards']}",
            "meets_expectation": bool(
                synth["hit_rate"] >= 0.8 and real["hit_rate"] >= 0.75
            ),
        },
        "recommend": (
            "Ship CF leaf-size-debiased + leaf merge as default Var(u); "
            "gate with Mag×Var + Y-flip; BCF optional; bootstrap fallback."
        ),
    }
    (out_dir / "po_risk_var_apps.json").write_text(
        json.dumps(payload, indent=2, default=str)
    )
    _plot(payload, out_dir / "po_risk_var_apps.png")
    _readme(out_dir / "README.md", payload)
    return payload


def _plot(payload: dict, path: Path):
    import matplotlib.pyplot as plt

    synth = payload["synthetic"]["boards"]
    real = payload["realworld"]["boards"]

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.5))

    ax = axes[0]
    names = [r["board"] for r in synth]
    x = np.arange(len(names))
    colors = ["#2F6B4F" if r["correct"] else "#A33B24" for r in synth]
    vals = [r["V_debiased"] for r in synth]
    ax.bar(x, vals, color=colors, width=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("mean V_debiased (LIVE)")
    ax.set_title(
        "Synthetic — hit %.0f%% (%s)"
        % (100 * payload["headline"]["synth_hit_rate"], payload["headline"]["synth_n"]),
        loc="left",
        fontweight="bold",
    )
    ax.grid(True, axis="y", alpha=0.3)

    ax = axes[1]
    names = [r["board"] for r in real]
    x = np.arange(len(names))
    gaps = [r["V_gap"] for r in real]
    colors = ["#2F6B4F" if r["pair_ok"] else "#A33B24" for r in real]
    ax.bar(x, gaps, color=colors, width=0.7)
    ax.axhline(0, color="#666", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("V_debiased gap (mal − ben)")
    ax.set_title(
        "Real-world — pair hit %.0f%% (%s)"
        % (100 * payload["headline"]["real_hit_rate"], payload["headline"]["real_n"]),
        loc="left",
        fontweight="bold",
    )
    ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle(
        "PO-risk Var(u) product justify — synthetic + real-world",
        fontsize=12,
        fontweight="bold",
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def _readme(path: Path, payload: dict):
    h = payload["headline"]
    lines = [
        "# PO-risk Var(u) — synthetic + real-world justify",
        "",
        payload["scope"],
        "",
        "**Recommend:** " + payload["recommend"],
        "",
        f"- Synthetic hit-rate: **{h['synth_hit_rate']:.2f}** ({h['synth_n']})",
        f"- Real-world pair hit-rate: **{h['real_hit_rate']:.2f}** ({h['real_n']})",
        f"- Meets expectation: **{h['meets_expectation']}**",
        "",
        "## Product logic (next step)",
        "",
    ]
    for i, s in enumerate(payload["product_logic"]["steps"], 1):
        lines.append(f"{i}. {s}")
    lines.extend(
        [
            "",
            "## Synthetic",
            "",
            "| board | oracle | pred | V_debiased | bias_ratio | OK? |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
    )
    for r in payload["synthetic"]["boards"]:
        lines.append(
            "| %s | %s | %s | %.4g | %.2f | %s |"
            % (
                r["board"],
                r["oracle"],
                r["pred"],
                r["V_debiased"],
                r["bias_ratio"],
                "yes" if r["correct"] else "no",
            )
        )
    lines.extend(
        [
            "",
            "## Real-world",
            "",
            "| board | ben pred | mal pred | V_gap | pair OK? |",
            "| --- | --- | --- | --- | --- |",
        ]
    )
    for r in payload["realworld"]["boards"]:
        lines.append(
            "| %s | %s | %s | %.4g | %s |"
            % (
                r["board"],
                r["benign"]["pred"],
                r["malignant"]["pred"],
                r["V_gap"],
                "yes" if r["pair_ok"] else "no",
            )
        )
    lines.extend(
        [
            "",
            "```bash",
            "pip install -e Python/po_risk_var",
            "po-risk-var justify --quick",
            "# or: python3 -m po_risk_var.cli justify --quick",
            "```",
            "",
        ]
    )
    path.write_text("\n".join(lines))
