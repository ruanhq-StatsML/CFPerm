#!/usr/bin/env python3
"""Correlation buckets + R-run frequency on Amazon / MSR-VTT.

Replays grad-cos online trajectories and reports:

  - low / mid / high correlation buckets (ρ̄ + erank)
  - R-run frequency (sticky redundant streaks, min_streak=2)
  - length-1 R events
  - leader-swap after uniqueness discount
  - LR ratio under soft_decorr geometry

  PYTHONPATH=. python3 scripts/run_agod_decorr_buckets.py
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from agod.ensemble_decorr import characterize_buckets_and_rruns

OUT = ROOT / "results" / "agod_decorr_buckets"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_decorr_buckets")
GRADCOS = ROOT / "results" / "agod_gradcos_lr" / "agod_gradcos_lr.json"

KEYS = (
    "amazon:equal",
    "amazon:soft",
    "amazon:soft_gradcos",
    "msrvtt:equal",
    "msrvtt:soft",
    "msrvtt:soft_gradcos",
)


def _f(x, nd=3):
    if x is None or (isinstance(x, float) and (np.isnan(x) or np.isinf(x))):
        return "—"
    return f"{float(x):.{nd}f}"


def load_gradcos_cells(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(path)
    payload = json.loads(path.read_text())
    traj = payload.get("trajectory") or {}
    out = {}
    for key in KEYS:
        rows = traj.get(key) or []
        if not rows:
            continue
        mods = list(rows[0]["alpha"].keys())
        norm = [
            {
                "t": r.get("t"),
                "alpha": r.get("alpha", {}),
                "mean_pair_cos": r.get("mean_pair_cos"),
                "pair_cos": r.get("pair_cos"),
                "acc_lift": r.get("acc_lift"),
            }
            for r in rows
        ]
        report = characterize_buckets_and_rruns(norm, mods, min_streak=2)
        windows = report.pop("windows")
        report["dataset"] = key.split(":")[0]
        report["source_scheduler"] = key.split(":")[1]
        report["n_source_windows"] = len(rows)
        lifts = [float(r["acc_lift"]) for r in rows if r.get("acc_lift") is not None]
        report["mean_acc_lift_source"] = float(np.mean(lifts)) if lifts else float("nan")
        out[key] = {"summary": report, "windows": windows}
    return out


def plot_board(cells: dict, path: Path):
    datasets = sorted({c["summary"]["dataset"] for c in cells.values()})
    fig, axes = plt.subplots(
        1, len(datasets), figsize=(4.4 * len(datasets), 4.0), facecolor="#f7f5f1"
    )
    if len(datasets) == 1:
        axes = [axes]
    for ax, ds in zip(axes, datasets):
        keys = [k for k in KEYS if k.startswith(ds + ":") and k in cells]
        x = np.arange(len(keys))
        low = [cells[k]["summary"]["buckets"]["low"]["frac_windows"] for k in keys]
        mid = [cells[k]["summary"]["buckets"]["mid"]["frac_windows"] for k in keys]
        high = [cells[k]["summary"]["buckets"]["high"]["frac_windows"] for k in keys]
        rhit = [cells[k]["summary"]["r_runs"]["frac_windows_with_any_R"] for k in keys]
        sticky = [
            cells[k]["summary"]["r_runs"]["runs_per_high_bucket_window"] for k in keys
        ]
        w = 0.18
        ax.bar(x - 1.5 * w, low, w, label="bucket low", color="#7a9e7e")
        ax.bar(x - 0.5 * w, mid, w, label="bucket mid", color="#c4a35a")
        ax.bar(x + 0.5 * w, high, w, label="bucket high", color="#b85c38")
        ax.plot(x, rhit, "o--", color="#333", label="frac any-R | decorr", lw=1.5)
        ax.plot(
            x, sticky, "s-", color="#1d3557", label="sticky R-runs / high-win", lw=1.5
        )
        ax.set_xticks(x)
        ax.set_xticklabels([k.split(":")[1] for k in keys], rotation=15)
        ax.set_ylim(0, 1.05)
        ax.set_title(f"{ds}: corr buckets + R-run rate")
        ax.set_ylabel("fraction / rate")
        ax.legend(frameon=False, fontsize=7, loc="upper right")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def write_docs(cells: dict, path: Path):
    rows = []
    for key in KEYS:
        if key not in cells:
            continue
        s = cells[key]["summary"]
        b = s["buckets"]
        rr = s["r_runs"]
        rows.append(
            f"| {s['dataset']} | {s['source_scheduler']} | "
            f"{_f(s['overall']['mean_rho'])} | {_f(s['overall']['mean_erank'])} | "
            f"{_f(b['low']['frac_windows'])}/{_f(b['mid']['frac_windows'])}/{_f(b['high']['frac_windows'])} | "
            f"{_f(s['overall']['frac_decorr_active'])} | "
            f"{_f(rr['frac_windows_with_any_R'])} | "
            f"{rr['n_runs']} ({_f(rr['runs_per_high_bucket_window'])}/high) | "
            f"{rr.get('n_r_events_len1', 0)} | "
            f"{_f(b['high']['mean_lr_ratio'])} | "
            f"{_f(s['mean_acc_lift_source'])} |"
        )
    high_role = []
    for key in KEYS:
        if key not in cells:
            continue
        s = cells[key]["summary"]
        h = s["buckets"]["high"]
        if h["n_windows"] == 0:
            continue
        rf = h["role_frac"]
        high_role.append(
            f"| {s['dataset']} | {s['source_scheduler']} | {h['n_windows']} | "
            f"{_f(rf.get('leader'))}/{_f(rf.get('diversifier'))}/{_f(rf.get('redundant'))} | "
            f"{_f(h['frac_leader_swap'])} | {_f(h['mean_acc_lift'])} |"
        )
    md = f"""# Correlation buckets + R-run frequency (Amazon / MSR-VTT)

## Protocol

Per online window from grad-cos trajectories:

1. rebuild pair geometry from `mean_pair_cos` (or `pair_cos` if logged)
2. `soft_decorr` role split → `ρ̄`, `erank`, roles, LR
3. **bucket** by correlation:
   - **high**: `ρ̄ ≥ 0.55` OR `erank ≤ 1.55` (decorr trigger-aligned)
   - **low**: `ρ̄ < 0.35` AND `erank > 2.20` (independent voters)
   - **mid**: borderline
4. **R-run**: contiguous `role=redundant` while decorr on, length ≥ 2 (sticky);
   also count length-1 R events

Budget actuator stays **LR** (`γ_role`); FWD always on. This script characterizes
adjust-correlation geometry on logged streams (no retrain).

## Summary

| Dataset | source | mean ρ | erank | bucket L/M/H | frac decorr | frac any-R | sticky R-runs (rate/high) | R events (len≥1) | LR ratio (high) | Acc↑ source |
|---|---|---:|---:|---|---:|---:|---|---:|---:|---:|
{chr(10).join(rows)}

## High-bucket role mix

| Dataset | source | n_high | role frac L/D/R | leader-swap | Acc↑ in high |
|---|---|---:|---|---:|---:|
{chr(10).join(high_role) if high_role else '| — | — | — | — | — | — |'}

## Readout (honest)

- **Amazon**: windows sit in **high** bucket (ρ̄≈0.67–0.71, erank≈1.5); decorr always on.
  Length-1 R events appear; sticky R-runs (k≥2) are rare on this short smoke (n=6).
  High-bucket LR ratio is large → soft_decorr reallocates step budget.
- **MSR-VTT**: mostly **low** (+ occasional mid); decorr off; R-run rate = 0
  → negative control for adjust-correlation.

```bash
PYTHONPATH=. python3 scripts/run_agod_decorr_buckets.py
```
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(md)


def write_latex(cells: dict, path: Path):
    lines = []
    for key in KEYS:
        if key not in cells:
            continue
        s = cells[key]["summary"]
        b = s["buckets"]
        rr = s["r_runs"]
        lines.append(
            f"{s['dataset']} & {s['source_scheduler'].replace('_', '\\_')} & "
            f"{_f(s['overall']['mean_rho'])} & {_f(s['overall']['mean_erank'])} & "
            f"{_f(b['high']['frac_windows'])} & "
            f"{_f(rr['frac_windows_with_any_R'])} & "
            f"{rr['n_runs']} & {_f(b['high']['mean_lr_ratio'])} \\\\"
        )
    tex = (
        "% Correlation buckets + R-run frequency\n"
        "\\begin{table}[t]\\centering\n"
        "\\caption{Amazon (high-corr) vs MSR-VTT (low-corr) buckets and R-run rates "
        "under soft\\_decorr role geometry.}\n"
        "\\label{tab:agod-decorr-buckets}\n"
        "\\begin{tabular}{llrrrrr}\\toprule\n"
        "dataset & source & mean $\\rho$ & erank & frac high & frac any-R & sticky $n$ & LR ratio \\\\\n"
        "\\midrule\n"
        + "\n".join(lines)
        + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(tex)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)
    cells = load_gradcos_cells(GRADCOS)
    payload = {
        "agod_version": "0.1.0",
        "focus": "corr buckets + R-run frequency on Amazon / MSR-VTT",
        "bucket_rule": {
            "high": "mean_rho>=0.55 OR erank<=1.55",
            "low": "mean_rho<0.35 AND erank>2.20",
            "mid": "otherwise",
        },
        "r_run": "contiguous redundant while decorr_active, min_streak=2",
        "cells": {k: v["summary"] for k, v in cells.items()},
        "windows": {k: v["windows"] for k, v in cells.items()},
    }
    (OUT / "agod_decorr_buckets.json").write_text(json.dumps(payload, indent=2))
    plot_board(cells, OUT / "AGOD_Decorr_Buckets_Board.png")
    write_docs(cells, OUT / "README.md")
    write_latex(cells, OUT / "AGOD_decorr_buckets_tables_only.tex")
    write_docs(cells, DOCS / "AGOD_decorr_buckets_rruns.md")
    write_latex(cells, DOCS / "AGOD_decorr_buckets_tables_only.tex")
    shutil.copy2(OUT / "AGOD_Decorr_Buckets_Board.png", ART / "AGOD_Decorr_Buckets_Board.png")
    shutil.copy2(OUT / "README.md", ART / "README.md")

    print("=== corr buckets + R-runs ===")
    for key, cell in cells.items():
        s = cell["summary"]
        b = s["buckets"]
        rr = s["r_runs"]
        print(
            f"  {key}: rho={s['overall']['mean_rho']:.3f} "
            f"erank={s['overall']['mean_erank']:.2f} "
            f"L/M/H={b['low']['frac_windows']:.2f}/"
            f"{b['mid']['frac_windows']:.2f}/"
            f"{b['high']['frac_windows']:.2f} "
            f"decorr={s['overall']['frac_decorr_active']:.2f} "
            f"anyR={rr['frac_windows_with_any_R']:.2f} "
            f"sticky={rr['n_runs']} "
            f"(rate/high={rr['runs_per_high_bucket_window']:.2f}) "
            f"R1={rr.get('n_r_events_len1', 0)} "
            f"LRratio_high={_f(b['high']['mean_lr_ratio'])}"
        )
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
