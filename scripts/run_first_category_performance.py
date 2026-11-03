#!/usr/bin/env python3
"""Per-dataset first-category indices → performance (weights vs no-weights).

  python3 scripts/run_first_category_performance.py --quick
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from first_category_performance import (  # noqa: E402
    _fit_scheme,
    evaluate_board_first_category,
    first_indices_per_category,
    predict_scaled,
    record_first_category_performance,
    run_all_boards_first_category,
    split_first_indices,
    union_indices,
)
from category_shift_boards import make_cgcd_stream  # noqa: E402
from weighted_cold_start import SEED  # noqa: E402

OUT = ROOT / "results" / "first_category_performance"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")


def _fmt(x):
    try:
        return "%.3f" % float(x)
    except Exception:
        return "—"


def eval_cgcd_stream_first_category(*, seed: int, quick: bool, frac: float) -> dict:
    """Streaming C-GCD: concat stages in time order; first indices per novel class."""
    stream = make_cgcd_stream(
        n_known=10 if quick else 12,
        n_novel_per_stage=2 if quick else 3,
        n_stages=3 if quick else 4,
        n_per=24 if quick else 40,
        p=24 if quick else 32,
        novel_dims=6 if quick else 8,
        seed=seed,
    )
    # build time-ordered stream: REF known, then stages in order
    Xs = [stream["X_ref"]]
    ys = [stream["y_ref"]]
    for st in stream["stages"]:
        Xs.append(st["X"])
        ys.append(st["y"])
    X_all = np.vstack(Xs)
    y_all = np.concatenate(ys)
    # novel cats across stages
    novel = []
    for st in stream["stages"]:
        novel.extend(st["novel"])
    known = stream["known"]

    first = first_indices_per_category(y_all, frac=frac, categories=novel)
    first_tr, first_ev = split_first_indices(first, seed=seed + 7, eval_frac=0.35)
    known_ctx = np.flatnonzero(np.isin(y_all, known))
    rng = np.random.default_rng(seed)
    if len(known_ctx) > 80:
        known_ctx = rng.choice(known_ctx, size=80, replace=False)
    tr = np.unique(np.concatenate([known_ctx, union_indices(first_tr)]))

    classes = np.asarray(stream["all_classes"], int)
    rows = {}
    for sch in ("uniform", "weighted", "both_fused"):
        clf = _fit_scheme(
            stream["X_ref"], stream["y_ref"], X_all[tr], y_all[tr],
            scheme=sch, seed=seed, classes=classes, top_k=8,
        )
        pred = predict_scaled(clf, X_all)
        rec = record_first_category_performance(y_all, pred, first_ev, categories=novel)
        rows[sch] = {
            "first_micro_acc": rec["micro_acc"],
            "first_n_total": rec["n_total"],
            "per_category": {
                str(c): {"n": v["n"], "acc": v["acc"]} for c, v in rec["per_category"].items()
            },
        }
    base = rows["uniform"]["first_micro_acc"]
    deltas = {
        k: (float(v["first_micro_acc"] - base) if v["first_micro_acc"] == v["first_micro_acc"] and base == base else None)
        for k, v in rows.items() if k != "uniform"
    }
    return {
        "dataset": "cgcd_stream",
        "novel_labels": novel,
        "first_indices": {str(c): v.tolist() for c, v in first.items()},
        "schemes": rows,
        "delta_first_vs_uniform": deltas,
        "n_first_total": int(sum(len(v) for v in first.values())),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--frac", type=float, default=0.25, help="first fraction per category")
    ap.add_argument("--n-first", type=int, default=None, help="override frac with fixed count")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    payload = run_all_boards_first_category(
        seed=args.seed, quick=bool(args.quick),
        frac=args.frac, n_first=args.n_first, include_legacy=True,
    )
    cg = eval_cgcd_stream_first_category(
        seed=args.seed, quick=bool(args.quick), frac=args.frac,
    )
    payload["cgcd_stream"] = cg

    (OUT / "first_category_performance.json").write_text(json.dumps(payload, indent=2))

    lines = [
        "# First-category performance (per dataset)",
        "",
        "For each novel category, take the **earliest** stream indices "
        f"(frac={args.frac}" + (f", n_first={args.n_first}" if args.n_first else "") + "), "
        "then record accuracy on that cold-start window. Compare weights vs no-weights.",
        "",
        "| dataset | uniform first_acc | weighted | both_fused | Δ weighted | Δ fused |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for r in payload["boards"]:
        sc = r["schemes"]
        d = r.get("delta_first_vs_uniform") or {}
        lines.append(
            "| %s | %s | %s | %s | %s | %s |" % (
                r["dataset"],
                _fmt(sc["uniform"]["first_micro_acc"]),
                _fmt(sc["weighted"]["first_micro_acc"]),
                _fmt(sc["both_fused"]["first_micro_acc"]),
                _fmt(d.get("weighted")),
                _fmt(d.get("both_fused")),
            )
        )
    # cgcd row
    sc = cg["schemes"]
    d = cg.get("delta_first_vs_uniform") or {}
    lines.append(
        "| %s | %s | %s | %s | %s | %s |" % (
            cg["dataset"],
            _fmt(sc["uniform"]["first_micro_acc"]),
            _fmt(sc["weighted"]["first_micro_acc"]),
            _fmt(sc["both_fused"]["first_micro_acc"]),
            _fmt(d.get("weighted")),
            _fmt(d.get("both_fused")),
        )
    )
    lines.extend([
        "",
        f"mean Δ weighted (tabular boards): **{_fmt(payload.get('mean_delta_weighted'))}**",
        f"mean Δ both_fused (tabular boards): **{_fmt(payload.get('mean_delta_both_fused'))}**",
        "",
        "## Per-category (both_fused) detail",
        "",
    ])
    for r in payload["boards"] + [cg]:
        lines.append(f"### {r['dataset']}")
        lines.append("")
        lines.append("| category | n_first_eval | acc (fused) | acc (uniform) |")
        lines.append("| --- | --- | --- | --- |")
        fu = r["schemes"]["both_fused"]["per_category"]
        un = r["schemes"]["uniform"]["per_category"]
        for c in sorted(fu.keys(), key=lambda x: int(x)):
            lines.append(
                "| %s | %s | %s | %s |" % (
                    c, fu[c]["n"], _fmt(fu[c]["acc"]), _fmt(un.get(c, {}).get("acc")),
                )
            )
        lines.append("")

    lines.extend([
        "```bash",
        "python3 scripts/run_first_category_performance.py --quick",
        "```",
        "",
    ])
    md = "\n".join(lines)
    (OUT / "README.md").write_text(md)
    DOCS.mkdir(parents=True, exist_ok=True)
    (DOCS / "First_category_performance.md").write_text(md)
    try:
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        shutil.copy2(OUT / "first_category_performance.json", ARTIFACTS / "first_category_performance.json")
        shutil.copy2(OUT / "README.md", ARTIFACTS / "First_category_performance.md")
    except Exception:
        pass
    print(md)
    print("wrote", OUT / "first_category_performance.json")


if __name__ == "__main__":
    main()
