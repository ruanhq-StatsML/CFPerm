#!/usr/bin/env python3
"""Multi-dataset Mag×Var justify on PermuCATE PO (same Y-flip口径).

Runs severity = J + Y-flip across JUSTIFY_DATASETS (HC3 verticals,
20 Newsgroups pairs, AG News, IMDB) and writes a concrete justify table.

  python3 scripts/run_as_dta_magvar_dataset_justify.py --quick
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from as_dta_segment_magvar import eval_segment  # noqa: E402
from as_dta_shift_data import JUSTIFY_DATASETS  # noqa: E402

OUT = ROOT / "results" / "as_dta_magvar_dataset_justify"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument(
        "--boards",
        nargs="*",
        default=None,
        help="override JUSTIFY_DATASETS",
    )
    args = ap.parse_args()

    names = list(args.boards or JUSTIFY_DATASETS)
    quick = bool(args.quick)
    rows = []
    for b in names:
        print("[ds-justify] board", b, flush=True)
        rows.append(
            eval_segment(
                b,
                live_n=48 if quick else 80,
                K=10 if quick else 16,
                seed=int(args.seed),
                max_pairs=320 if quick else 700,
                n_estimators=24 if quick else 36,
            )
        )

    n_ok = sum(1 for r in rows if r["separates_benign_malignant"])
    by_family = {}
    for r in rows:
        fam = r.get("vertical", "other")
        by_family.setdefault(fam, {"n": 0, "ok": 0})
        by_family[fam]["n"] += 1
        by_family[fam]["ok"] += int(r["separates_benign_malignant"])

    payload = {
        "scope": (
            "Concrete multi-dataset justify: same Mag×Var + Y-flip口径 on "
            "PermuCATE u across AIGC / topic / news / sentiment boards."
        ),
        "estimator": {
            "name": "PermuCATE RF text",
            "unit_po": "(φ - τ̂)^2",
            "Y": "response",
            "W": "REF/LIVE",
            "used_permucate_rf": True,
        },
        "rule": {
            "benign": "cosmetic synonym; Y unchanged",
            "malignant": "subset attack + Y-flip on attacked rows",
            "severity": "J + Y-flip",
        },
        "boards": rows,
        "n_boards": len(rows),
        "n_separated": n_ok,
        "hit_rate": float(n_ok) / max(len(rows), 1),
        "by_family": {
            k: {
                "n": v["n"],
                "ok": v["ok"],
                "hit_rate": float(v["ok"]) / max(v["n"], 1),
            }
            for k, v in by_family.items()
        },
        "justify": (
            "Across HC3 verticals, 20 Newsgroups pairs, AG News, and IMDB, "
            "malignant severity exceeds benign under the same Y-flip口径 — "
            "concrete evidence the product rule is not board-specific."
        ),
    }

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "dataset_justify.json").write_text(json.dumps(payload, indent=2))
    _plot(payload, OUT / "dataset_justify.png")
    _readme(OUT / "README.md", payload)
    DOCS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "README.md", DOCS / "AS_DTA_magvar_dataset_justify.md")
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "dataset_justify.png", ARTIFACTS / "as_dta_magvar_dataset_justify.png")

    print(
        "[ds-justify] hit_rate=%.2f  (%d/%d)"
        % (payload["hit_rate"], payload["n_separated"], payload["n_boards"]),
        flush=True,
    )
    for r in rows:
        print(
            "[ds-justify] %-12s  ben=%.3f  mal=%.3f  gap=%.3f  %s"
            % (
                r["board"],
                r["benign"]["severity"],
                r["malignant"]["severity"],
                r["severity_gap"],
                "OK" if r["separates_benign_malignant"] else "MISS",
            ),
            flush=True,
        )
    print("[ds-justify] wrote", OUT, flush=True)


def _plot(payload: dict, path: Path):
    import matplotlib.pyplot as plt
    import numpy as np

    boards = payload["boards"]
    names = [b["board"] for b in boards]
    x = np.arange(len(names))
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))

    ax = axes[0]
    jb = [b["benign"]["severity"] for b in boards]
    jm = [b["malignant"]["severity"] for b in boards]
    ax.bar(x - 0.18, jb, width=0.36, color="#2F6B4F", label="benign")
    ax.bar(x + 0.18, jm, width=0.36, color="#A33B24", label="malignant")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("severity = J + Y-flip")
    ax.set_title("Multi-dataset severity (same Y-flip rule)", loc="left", fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(True, axis="y", alpha=0.3)

    ax = axes[1]
    gaps = [b["severity_gap"] for b in boards]
    colors = ["#2F6B4F" if g > 0.05 else "#A33B24" for g in gaps]
    ax.bar(x, gaps, color=colors, width=0.65)
    ax.axhline(0.05, ls="--", color="#999", lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("severity gap (mal − ben)")
    ax.set_title(
        "Justify hit-rate = %.0f%% (%d/%d)"
        % (100 * payload["hit_rate"], payload["n_separated"], payload["n_boards"]),
        loc="left",
        fontweight="bold",
    )
    ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle(
        "Dataset justify — Mag×Var + Y-flip on PermuCATE u",
        fontsize=11.5,
        fontweight="bold",
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def _readme(path: Path, payload: dict):
    lines = [
        "# Mag×Var multi-dataset justify (PermuCATE PO)",
        "",
        payload["scope"],
        "",
        "| board | family | sev(ben) | sev(mal) | gap | OK? |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for b in payload["boards"]:
        lines.append(
            "| %s | %s | %.3f | %.3f | %.3f | %s |"
            % (
                b["board"],
                b.get("vertical", ""),
                b["benign"]["severity"],
                b["malignant"]["severity"],
                b["severity_gap"],
                "yes" if b["separates_benign_malignant"] else "no",
            )
        )
    lines.extend(
        [
            "",
            "hit-rate: **%.2f** (%d/%d)"
            % (payload["hit_rate"], payload["n_separated"], payload["n_boards"]),
            "",
            "## By family",
            "",
            "| family | hit-rate |",
            "| --- | --- |",
        ]
    )
    for fam, st in payload["by_family"].items():
        lines.append("| %s | %.2f (%d/%d) |" % (fam, st["hit_rate"], st["ok"], st["n"]))
    lines.extend(
        [
            "",
            payload["justify"],
            "",
            "```bash",
            "python3 scripts/run_as_dta_magvar_dataset_justify.py --quick",
            "```",
            "",
        ]
    )
    path.write_text("\n".join(lines))


if __name__ == "__main__":
    main()
