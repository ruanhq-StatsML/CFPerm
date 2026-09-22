#!/usr/bin/env python3
"""Run post-hoc localization directly on the on-disk audit tables.

No new statistics. Load the csv, pull subset indices, look at the mean,
then call MMD() and po_risk() on every pair.

    HH two-stream          T = helpful / harmless queue
    HH multi-step          T = hop 0 / 1 / 2+
    HH helpful hop         T = batch before/after the policy cut
    real two-stream        T = BeaverTails / ToxicChat

HH chosen is not Y. Raw text is not stored.

Usage::

    PYTHONPATH=. python3 scripts/run_posthoc_on_dataset.py
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.posthoc_localization import (  # noqa: E402
    localize,
    markdown_conditional_means,
)

AUDIT = ROOT / "results" / "manuscript" / "llm_audit"
OUT = ROOT / "results" / "manuscript" / "posthoc_dataset"
LEAK = ("chosen", "rejected")
CUT_BATCH = 4
SEED = 2026
N_PERM = 25
TOP_X = "x_n_toks"


def x_columns(fieldnames) -> list[str]:
    return [c for c in fieldnames if c.startswith("x_")]


def load_dataset(path: Path):
    with path.open() as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"empty {path}")
    keys = {k.lower() for k in rows[0]}
    for bad in LEAK:
        if bad in keys:
            raise RuntimeError(f"{path} still has {bad} — that is not Y")
    cols = x_columns(rows[0].keys())
    ycol = "Y" if "Y" in rows[0] else "y"
    Y = np.asarray([int(float(r[ycol])) for r in rows], dtype=int)
    X = np.asarray([[float(r[c]) for c in cols] for r in rows], dtype=float)
    labels = None
    if "T" in rows[0]:
        labels = np.asarray([int(r["T"]) for r in rows], dtype=int)
        how = "column T"
    elif "step" in rows[0]:
        labels = np.clip(np.asarray([int(r["step"]) for r in rows], dtype=int), 0, 2)
        how = "step clipped to 0/1/2+"
    elif "batch" in rows[0]:
        batch = np.asarray([int(r["batch"]) for r in rows], dtype=int)
        labels = (batch >= CUT_BATCH).astype(int)
        how = f"batch >= {CUT_BATCH}"
    else:
        raise RuntimeError(f"{path} has no T, step, or batch to pull subsets from")
    return X, Y, labels, cols, how


def top_feature(X: np.ndarray, cols: list[str]) -> np.ndarray:
    j = cols.index(TOP_X) if TOP_X in cols else 0
    return X[:, j], cols[j]


DATASETS = [
    {
        "name": "hh_two_stream",
        "title": "HH two-stream",
        "path": AUDIT / "xy_hh_two_stream.csv",
        "note": "T = helpful vs harmless queue",
    },
    {
        "name": "hh_multistep",
        "title": "HH multi-step hops",
        "path": AUDIT / "xy_hh_multistep_consistent.csv",
        "note": "T = hop 0 / 1 / 2+",
    },
    {
        "name": "hh_helpful_hop",
        "title": "HH helpful policy hop",
        "path": AUDIT / "xy_hh_helpful_hop.csv",
        "note": "T = batch before/after the cut",
    },
    {
        "name": "real_two_stream",
        "title": "BeaverTails vs ToxicChat",
        "path": AUDIT / "xy_real_two_stream.csv",
        "note": "T = two real judge queues",
    },
]


def jsonable(obj):
    if isinstance(obj, dict):
        return {k: jsonable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return jsonable(obj.tolist())
    if isinstance(obj, (np.floating, float)):
        x = float(obj)
        return None if not np.isfinite(x) else x
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    return obj


def compact(rec: dict, spec: dict) -> dict:
    loc = rec["localization"]
    groups = loc.get("groups") or {}
    bins = loc.get("bins") or {}
    return {
        "name": spec["name"],
        "title": spec["title"],
        "note": spec["note"],
        "how": rec["how"],
        "n": rec["n"],
        "top_x": rec["top_x"],
        "path": str(spec["path"].relative_to(ROOT)),
        "loc_sig_pairs": groups.get("n_sig_pairs"),
        "loc_pairwise": groups.get("pairwise"),
        "loc_means": groups.get("conditional_mean"),
        "loc_bin_means": bins.get("conditional_mean"),
    }


def plot_means(rows: list[dict], path: Path) -> None:
    fig, axes = plt.subplots(1, len(rows), figsize=(3.4 * len(rows), 3.2), sharey=False)
    if len(rows) == 1:
        axes = [axes]
    for ax, r in zip(axes, rows):
        means = r.get("loc_means") or {}
        keys = sorted(means)
        ys = [means[k]["mean_Y"] for k in keys]
        ax.bar(keys, ys, color="#1f4e79")
        ax.set_title(r["title"], fontsize=9)
        ax.set_ylabel("mean Y")
        ax.set_ylim(0, 1)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def render_report(rows: list[dict]) -> str:
    lines = [
        "# Post-hoc localization on the audit tables",
        "",
        "Directly on the csv: pull subset indices, **look at the mean**, then call `MMD()` and `po_risk()`.",
        "HH chosen is not Y.",
        "",
        "```bash",
        "PYTHONPATH=. python3 scripts/run_posthoc_on_dataset.py",
        "```",
        "",
        "| Dataset | subsets from | n | sig pairs |",
        "|---|---|---:|---:|",
    ]
    for r in rows:
        lines.append(
            f"| {r['title']} | {r['how']} | {r['n']} | {r.get('loc_sig_pairs', '—')} |"
        )
    lines += ["", "Look at the mean first (already computed). Then pairwise MMD / PO-risk.", ""]
    lines += markdown_conditional_means(rows)
    lines += [
        "## Pairwise subset MMD / PO-risk",
        "",
        "| Dataset | pair | n | MMD | MMD p | PO-risk | PO p | mean Y |",
        "|---|---|---|---:|---:|---:|---:|---|",
    ]
    for r in rows:
        for p in r.get("loc_pairwise") or []:
            lines.append(
                "| {title} | {a} vs {b} | {na}/{nb} | {mmd:.3g} | {mp:.3g} | {po:.3g} | {pp:.3g} | {ya:.3f}/{yb:.3f} |".format(
                    title=r["title"],
                    a=p["a"],
                    b=p["b"],
                    na=p["n_a"],
                    nb=p["n_b"],
                    mmd=p["mmd"],
                    mp=p["mmd_p"],
                    po=p.get("po_risk", float("nan")),
                    pp=p.get("po_p", float("nan")),
                    ya=p.get("mean_Y_a", float("nan")),
                    yb=p.get("mean_Y_b", float("nan")),
                )
            )
    lines += ["", "Quartile splits of `x_n_toks` are in the mean table (Q0–Q3). Nothing else.", ""]
    return "\n".join(lines) + "\n"


def run_one(spec: dict) -> dict:
    X, Y, labels, cols, how = load_dataset(spec["path"])
    feat, feat_name = top_feature(X, cols)
    loc = localize(X, Y, group_labels=labels, feature=feat, n_perm=N_PERM, seed=SEED)
    return {
        "name": spec["name"],
        "title": spec["title"],
        "how": how,
        "n": int(len(Y)),
        "top_x": feat_name,
        "localization": loc,
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    runs = []
    skipped = []
    for spec in DATASETS:
        if not spec["path"].exists():
            skipped.append(str(spec["path"].relative_to(ROOT)))
            continue
        rec = run_one(spec)
        runs.append((rec, spec))
        (OUT / f"{spec['name']}.json").write_text(
            json.dumps(jsonable(rec), indent=2) + "\n", encoding="utf-8"
        )
    if not runs:
        raise RuntimeError("no audit tables on disk")
    summary = [compact(rec, spec) for rec, spec in runs]
    plot_means(summary, OUT / "mean_by_subset.png")
    report = render_report(summary)
    if skipped:
        report += "Skipped missing files:\n\n" + "\n".join(f"- `{p}`" for p in skipped) + "\n"
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (OUT / "REPORT.md").write_text(report, encoding="utf-8")
    print(report)
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
