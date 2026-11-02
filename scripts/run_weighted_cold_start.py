#!/usr/bin/env python3
"""Weighted cold-start under category shift (C-GCD + Office-Home PDA proxies).

  python3 scripts/run_weighted_cold_start.py --quick
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from category_shift_boards import run_suite  # noqa: E402

OUT = ROOT / "results" / "weighted_cold_start"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")


def _fmt(x):
    try:
        return "%.3f" % float(x)
    except Exception:
        return "—"


def _readme(payload: dict) -> str:
    w, m, u = payload["cgcd_weighted"], payload["cgcd_membership"], payload["cgcd_uniform"]
    eff, oh = payload["efficiency"], payload["officehome"]
    lines = [
        "# Weighted cold-start under category shift (v2)",
        "",
        payload["takeaway"],
        "",
        "## Protocol (four axes)",
        "",
    ]
    for a in payload["protocol"]["axes"]:
        lines.append("- " + a)
    lines.extend([
        "",
        "### Modes",
        "",
        "| mode | meaning |",
        "| --- | --- |",
    ])
    for k, v in payload["protocol"]["modes"].items():
        lines.append("| %s | %s |" % (k, v))
    lines.extend([
        "",
        "## CIFAR-100 C-GCD proxy",
        "",
        "| mode | known | novel | forgetting | feat_hit | vimp_top_mass |",
        "| --- | --- | --- | --- | --- | --- |",
        "| weighted | %s | %s | %s | %s | %s |" % (
            _fmt(w["mean_known_acc"]), _fmt(w["mean_novel_acc"]), _fmt(w["mean_forgetting"]),
            _fmt(w["mean_feat_hit"]), _fmt(w["mean_vimp_top_mass"]),
        ),
        "| membership | %s | %s | %s | %s | %s |" % (
            _fmt(m["mean_known_acc"]), _fmt(m["mean_novel_acc"]), _fmt(m["mean_forgetting"]),
            _fmt(m["mean_feat_hit"]), _fmt(m["mean_vimp_top_mass"]),
        ),
        "| uniform | %s | %s | %s | — | — |" % (
            _fmt(u["mean_known_acc"]), _fmt(u["mean_novel_acc"]), _fmt(u["mean_forgetting"]),
        ),
        "",
        "Δ novel (weighted − uniform): **%s**" % _fmt(payload["delta_novel_acc_vs_uniform"]),
        "",
        "Δ novel (membership − uniform): **%s**" % _fmt(payload["delta_novel_acc_membership"]),
        "",
        "Δ forgetting (uniform − weighted; &gt;0 ⇒ weighted forgets less): **%s**" % _fmt(
            payload["delta_forgetting_vs_uniform"]
        ),
        "",
        "## Cold-start efficiency (n* novel samples → novel_acc ≥ %.2f)" % eff["target"],
        "",
        "| mode | n* | curve (budget:acc) |",
        "| --- | --- | --- |",
    ])
    for mode, cur in eff["curves"].items():
        curve = ", ".join("%d:%.2f" % (p["budget"], p["novel_acc"]) for p in cur["points"])
        lines.append("| %s | %s | %s |" % (mode, cur["n_star"] if cur["n_star"] is not None else "—", curve))
    lines.extend([
        "",
        "efficiency winner: **%s**" % (eff.get("winner") or "—"),
        "",
        "## Office-Home PDA proxy",
        "",
        "| overlap | shared acc | incomparable rate | adapt rate |",
        "| --- | --- | --- | --- |",
    ])
    for ov, v in oh["by_overlap"].items():
        lines.append("| %s | %s | %s | %s |" % (
            ov, _fmt(v["mean_shared_acc"]), _fmt(v["incomparable_rate"]), _fmt(v["adapt_rate"]),
        ))
    lines.extend([
        "",
        "Directed adapt = detect gate + Jaccard boundary + `w=1+λ·membership` + `partial_fit`.",
        "",
        "```bash",
        "python3 scripts/run_weighted_cold_start.py --quick",
        "```",
        "",
    ])
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    payload = run_suite(seed=int(args.seed), quick=bool(args.quick))
    (OUT / "weighted_cold_start.json").write_text(json.dumps(payload, indent=2, default=str))
    readme = _readme(payload)
    (OUT / "README.md").write_text(readme)
    DOCS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "README.md", DOCS / "Weighted_cold_start_category_shift.md")
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "README.md", ARTIFACTS / "weighted_cold_start_README.md")
    shutil.copy(OUT / "weighted_cold_start.json", ARTIFACTS / "weighted_cold_start.json")

    print("[wcs] Δnovel(w-u)=%s  Δnovel(m-u)=%s  Δforget=%s  eff_winner=%s" % (
        _fmt(payload["delta_novel_acc_vs_uniform"]),
        _fmt(payload["delta_novel_acc_membership"]),
        _fmt(payload["delta_forgetting_vs_uniform"]),
        payload["efficiency"].get("winner"),
    ), flush=True)
    for mode, cur in payload["efficiency"]["curves"].items():
        print("[wcs] n*(%s)=%s" % (mode, cur["n_star"]), flush=True)
    for ov, v in payload["officehome"]["by_overlap"].items():
        print("[wcs] overlap=%s  shared=%s  incomparable=%s  adapt=%s" % (
            ov, _fmt(v["mean_shared_acc"]), _fmt(v["incomparable_rate"]), _fmt(v["adapt_rate"]),
        ), flush=True)
    print("[wcs] wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
