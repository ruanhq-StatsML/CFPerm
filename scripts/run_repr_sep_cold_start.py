#!/usr/bin/env python3
"""Representation-first cold-start board (4-step pipeline).

  python3 scripts/run_repr_sep_cold_start.py --quick
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from category_shift_boards import run_cgcd  # noqa: E402

OUT = ROOT / "results" / "repr_sep_cold_start"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")


def _fmt(x):
    try:
        return "%.3f" % float(x)
    except Exception:
        return "—"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--expand-min", type=int, default=10)
    ap.add_argument("--tta", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    modes = {}
    for mode in ("repr_sep", "weighted", "uniform"):
        modes[mode] = run_cgcd(
            seed=args.seed, mode=mode, quick=bool(args.quick),
            expand_min=args.expand_min, do_tta=bool(args.tta) and mode == "repr_sep",
        )

    r, w, u = modes["repr_sep"], modes["weighted"], modes["uniform"]
    payload = {
        "seed": args.seed,
        "quick": bool(args.quick),
        "expand_min": args.expand_min,
        "tta": bool(args.tta),
        "modes": modes,
        "delta_novel_repr_vs_uniform": r["mean_novel_acc"] - u["mean_novel_acc"],
        "delta_novel_repr_vs_weighted": r["mean_novel_acc"] - w["mean_novel_acc"],
        "delta_forget_repr_vs_uniform": u["mean_forgetting"] - r["mean_forgetting"],
        "takeaway": (
            "repr_sep = attribute → margin separation on S → expand head when "
            "buffer≥N* / footprint stable. Instant learn = sep; learn = expand."
        ),
    }
    (OUT / "repr_sep_cold_start.json").write_text(json.dumps(payload, indent=2))

    rows_diag = []
    for row in r.get("stages") or []:
        rows_diag.append(
            "| %s | %s | %s | %s | %s | %s |" % (
                row.get("stage"), _fmt(row.get("novel_acc")), _fmt(row.get("known_acc")),
                row.get("phase"), _fmt(row.get("gap_lift")), row.get("expand_fitted"),
            )
        )
    md = "\n".join([
        "# Representation-first cold-start",
        "",
        payload["takeaway"],
        "",
        "| mode | known | novel | forgetting | feat_hit |",
        "| --- | --- | --- | --- | --- |",
        "| repr_sep | %s | %s | %s | %s |" % (
            _fmt(r["mean_known_acc"]), _fmt(r["mean_novel_acc"]),
            _fmt(r["mean_forgetting"]), _fmt(r["mean_feat_hit"]),
        ),
        "| weighted | %s | %s | %s | %s |" % (
            _fmt(w["mean_known_acc"]), _fmt(w["mean_novel_acc"]),
            _fmt(w["mean_forgetting"]), _fmt(w["mean_feat_hit"]),
        ),
        "| uniform | %s | %s | %s | — |" % (
            _fmt(u["mean_known_acc"]), _fmt(u["mean_novel_acc"]),
            _fmt(u["mean_forgetting"]),
        ),
        "",
        "Δ novel (repr − uniform): **%s**" % _fmt(payload["delta_novel_repr_vs_uniform"]),
        "Δ novel (repr − weighted): **%s**" % _fmt(payload["delta_novel_repr_vs_weighted"]),
        "",
        "## repr_sep per-stage",
        "",
        "| stage | novel | known | phase | gap_lift | expanded |",
        "| --- | --- | --- | --- | --- | --- |",
        *rows_diag,
        "",
        "```bash",
        "python3 scripts/run_repr_sep_cold_start.py --quick",
        "```",
        "",
    ])
    (OUT / "README.md").write_text(md)
    DOCS.mkdir(parents=True, exist_ok=True)
    (DOCS / "Repr_sep_cold_start.md").write_text(md)
    try:
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        shutil.copy2(OUT / "repr_sep_cold_start.json", ARTIFACTS / "repr_sep_cold_start.json")
        shutil.copy2(OUT / "README.md", ARTIFACTS / "Repr_sep_cold_start.md")
    except Exception:
        pass
    print(md)
    print("wrote", OUT / "repr_sep_cold_start.json")


if __name__ == "__main__":
    main()
