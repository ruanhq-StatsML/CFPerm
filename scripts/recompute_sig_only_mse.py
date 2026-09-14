#!/usr/bin/env python3
"""Recompute next-MSE on OnlineRFPerm-significant batches only.

Non-significant steps follow uniform for gated modes — exclude them from the
comparison. Uses ``gate_on`` from RFPerm refit results (same reject sequence
as the cbrt compare when seed/hparams match).

  PYTHONPATH=. python3 scripts/recompute_sig_only_mse.py
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from agod.sig_batch_metrics import annotate_results_with_sig, pick_gate_mask, significant_only_stats

REFIT = Path("results/agod_rfperm_po_refit/summary.json")
CBRT = Path("results/agod_po_cbrt/summary.json")
OUT = Path("results/agod_sig_only")

REFIT_MODES = ("uniform", "sqrt", "sqrt_gated", "sqrt_gated_refit", "dre")
CBRT_MODES = ("uniform", "sqrt", "cbrt", "gated_sqrt", "gated_cbrt", "dre")
REFIT_COLORS = {
    "uniform": "#4C566A",
    "sqrt": "#A3BE8C",
    "sqrt_gated": "#88C0D0",
    "sqrt_gated_refit": "#5E81AC",
    "dre": "#BF616A",
}
CBRT_COLORS = {
    "uniform": "#4C566A",
    "sqrt": "#A3BE8C",
    "cbrt": "#EBCB8B",
    "gated_sqrt": "#88C0D0",
    "gated_cbrt": "#5E81AC",
    "dre": "#BF616A",
}


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def inject_gate_from_refit(cbrt: dict, refit: dict) -> None:
    """Copy RFPerm gate_on into cbrt results (duties already match)."""
    for ds, block in cbrt["datasets"].items():
        gate = pick_gate_mask(
            refit["datasets"][ds]["results"],
            preferred_modes=("sqrt_gated_refit", "sqrt_gated"),
        )
        if gate is None:
            continue
        # store full-length gate from refit (pre-align); annotate will re-align
        src = refit["datasets"][ds]["results"]["sqrt_gated_refit"]["gate_on"]
        for mode, r in block["results"].items():
            r["gate_on"] = list(src)


def summarize_table(
    all_sig: Dict[str, dict],
    modes: List[str],
    highlight: str,
) -> str:
    header = "| dataset | n_sig | duty | " + " | ".join(modes) + " | best |"
    sep = "|---|---:|---:|" + "|".join(["---:"] * len(modes)) + "|---|"
    lines = [header, sep]
    wins = {m: 0 for m in modes}
    for ds, blob in all_sig.items():
        meta = blob["_meta"]
        means = {m: blob[m]["mse_mean_sig"] for m in modes}
        best = min(means, key=lambda k: means[k] if means[k] == means[k] else 1e99)
        wins[best] += 1
        cells = []
        for m in modes:
            v = means[m]
            s = f"{v:.4g}" if v == v else "nan"
            if m == highlight:
                s = f"**{s}**"
            cells.append(s)
        lines.append(
            f"| `{ds}` | {meta['n_significant']}/{meta['n_eval']} | "
            f"{meta['duty']:.2f} | " + " | ".join(cells) + f" | `{best}` |"
        )
    lines.append("")
    lines.append("**Wins (sig-only):** " + ", ".join(f"`{m}`={wins[m]}" for m in modes))
    return "\n".join(lines)


def rel_to_uniform_table(all_sig: Dict[str, dict], modes: List[str]) -> str:
    header = "| dataset | " + " | ".join(m for m in modes if m != "uniform") + " |"
    sep = "|---|" + "|".join(["---:"] * (len(modes) - 1)) + "|"
    lines = [header, sep]
    for ds, blob in all_sig.items():
        cells = []
        for m in modes:
            if m == "uniform":
                continue
            v = blob[m]["rel_mse_vs_uniform_sig"]
            cells.append(f"{v:.3f}×" if v == v else "nan")
        lines.append(f"| `{ds}` | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def plot_rel(all_sig: Dict[str, dict], modes: List[str], colors: dict, title: str, path: Path) -> None:
    ds = list(all_sig.keys())
    fig, ax = plt.subplots(figsize=(max(8, 1.7 * len(ds)), 4.6))
    x = np.arange(len(ds))
    w = 0.8 / len(modes)
    for i, mode in enumerate(modes):
        vals = [all_sig[d][mode]["rel_mse_vs_uniform_sig"] for d in ds]
        ax.bar(x + (i - (len(modes) - 1) / 2) * w, vals, w, label=mode, color=colors[mode])
    ax.axhline(1.0, color="k", ls="--", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ds, rotation=15, ha="right")
    ax.set_ylabel("Next-MSE / uniform (significant batches only)")
    ax.set_title(title)
    ax.legend(ncol=3, fontsize=7)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def build_sig(datasets: dict, preferred: tuple) -> Dict[str, dict]:
    out = {}
    for ds, block in datasets.items():
        stats = significant_only_stats(
            block["results"], preferred_gate_modes=preferred
        )
        out[ds] = stats
        annotate_results_with_sig(block["results"], preferred_gate_modes=preferred)
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    refit = load(REFIT)
    cbrt = load(CBRT)
    inject_gate_from_refit(cbrt, refit)

    refit_sig = build_sig(refit["datasets"], ("sqrt_gated_refit", "sqrt_gated"))
    cbrt_sig = build_sig(cbrt["datasets"], ("gated_cbrt", "gated_sqrt"))

    payload = {
        "note": (
            "Next-batch MSE averaged only on OnlineRFPerm reject steps. "
            "Non-significant batches stay uniform and are excluded."
        ),
        "refit": refit_sig,
        "cbrt": cbrt_sig,
    }
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    # also refresh source summaries with mse_mean_sig fields
    REFIT.write_text(json.dumps(refit, indent=2), encoding="utf-8")
    CBRT.write_text(json.dumps(cbrt, indent=2), encoding="utf-8")

    plot_rel(
        refit_sig,
        list(REFIT_MODES),
        REFIT_COLORS,
        "Significant batches only — RFPerm + T0/T1 √PO re-fit",
        OUT / "mse_rel_sig_refit.png",
    )
    plot_rel(
        cbrt_sig,
        list(CBRT_MODES),
        CBRT_COLORS,
        "Significant batches only — PO^{1/2} vs PO^{1/3} gated",
        OUT / "mse_rel_sig_cbrt.png",
    )

    md = "\n".join(
        [
            "# Significant-batch-only next-MSE",
            "",
            "Non-reject steps: gated modes = **uniform** → exclude from the mean.",
            "Compare methods only where OnlineRFPerm opens the gate.",
            "",
            "## RFPerm → T0/T1 √PO re-fit",
            "",
            summarize_table(refit_sig, list(REFIT_MODES), "sqrt_gated_refit"),
            "",
            "Relative to uniform (sig-only, ↓ better):",
            "",
            rel_to_uniform_table(refit_sig, list(REFIT_MODES)),
            "",
            "## PO^{1/3} vs PO^{1/2} (same gate)",
            "",
            summarize_table(cbrt_sig, list(CBRT_MODES), "gated_cbrt"),
            "",
            "Relative to uniform (sig-only, ↓ better):",
            "",
            rel_to_uniform_table(cbrt_sig, list(CBRT_MODES)),
            "",
            "```",
            "metric = mean(mse_next[t] for t where reject_t)",
            "# non-significant t: all gated modes ≡ uniform → drop",
            "```",
            "",
        ]
    )
    (OUT / "SIG_ONLY_REPORT.md").write_text(md, encoding="utf-8")
    Path("docs/agod").mkdir(parents=True, exist_ok=True)
    Path("docs/agod/AGOD_sig_only_mse.md").write_text(md, encoding="utf-8")

    # refresh refit/cbrt report sections
    refit_md = Path("results/agod_rfperm_po_refit/RFPERM_PO_REFIT_REPORT.md")
    if refit_md.exists():
        extra = (
            "\n## Significant batches only\n\n"
            + summarize_table(refit_sig, list(REFIT_MODES), "sqrt_gated_refit")
            + "\n"
        )
        text = refit_md.read_text(encoding="utf-8")
        if "## Significant batches only" not in text:
            refit_md.write_text(text.rstrip() + "\n" + extra, encoding="utf-8")
            Path("docs/agod/AGOD_rfperm_po_refit.md").write_text(
                refit_md.read_text(encoding="utf-8"), encoding="utf-8"
            )

    cbrt_md = Path("results/agod_po_cbrt/PO_CBRT_REPORT.md")
    if cbrt_md.exists():
        extra = (
            "\n## Significant batches only\n\n"
            + summarize_table(cbrt_sig, list(CBRT_MODES), "gated_cbrt")
            + "\n"
        )
        text = cbrt_md.read_text(encoding="utf-8")
        if "## Significant batches only" not in text:
            cbrt_md.write_text(text.rstrip() + "\n" + extra, encoding="utf-8")
            Path("docs/agod/AGOD_po_cbrt.md").write_text(
                cbrt_md.read_text(encoding="utf-8"), encoding="utf-8"
            )

    print(md)


if __name__ == "__main__":
    main()
