#!/usr/bin/env python3
"""Sweep DiffusionDB time-window schemes (+ optional hyperparam concat).

Shows that equal_count vs equal_time vs width=Nd are *not* the same cut:
n_by_T / span_hours / ΔȲ / HGB AUC / top tokens all move.

  PYTHONPATH=. python3 scripts/run_diffusiondb_window_sweep.py \\
    --n-sample 3000 --out results/diffusiondb_window_sweep
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from scripts.run_diffusiondb_temporal_fsds import load_subset, run_once

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", type=Path, default=ROOT / "data/diffusiondb/metadata.parquet")
    ap.add_argument("--n-sample", type=int, default=3000)
    ap.add_argument("--max-features", type=int, default=256)
    ap.add_argument("--select-k", type=int, default=25)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=ROOT / "results/diffusiondb_window_sweep")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    df = load_subset(args.meta, n_sample=args.n_sample, seed=args.seed)
    configs = [
        {"name": "equal_count_k2", "window_scheme": "equal_count", "n_windows": 2, "width_hours": None, "concat_hyperparams": False},
        {"name": "equal_count_k4", "window_scheme": "equal_count", "n_windows": 4, "width_hours": None, "concat_hyperparams": False},
        {"name": "equal_time_k2", "window_scheme": "equal_time", "n_windows": 2, "width_hours": None, "concat_hyperparams": False},
        {"name": "equal_time_k4", "window_scheme": "equal_time", "n_windows": 4, "width_hours": None, "concat_hyperparams": False},
        {"name": "width_48h", "window_scheme": "width", "n_windows": 2, "width_hours": 48.0, "concat_hyperparams": False},
        {"name": "equal_count_k2_hp", "window_scheme": "equal_count", "n_windows": 2, "width_hours": None, "concat_hyperparams": True},
    ]
    rows = []
    for cfg in configs:
        print("run", cfg["name"], flush=True)
        out = run_once(
            df,
            embed="tfidf",
            max_features=args.max_features,
            select_k=args.select_k,
            y_quantile=0.7,
            seed=args.seed,
            concat_hyperparams_flag=cfg["concat_hyperparams"],
            n_windows=cfg["n_windows"],
            window_scheme=cfg["window_scheme"],
            width_hours=cfg["width_hours"],
        )
        if not out.get("ok"):
            rows.append({"name": cfg["name"], "ok": False, "reason": out.get("reason")})
            continue
        out.pop("_tables", None)
        hgb = ((out.get("fsds") or {}).get("models") or {}).get("hgb") or {}
        rows.append(
            {
                "name": cfg["name"],
                "ok": True,
                "scheme": cfg["window_scheme"],
                "n_windows": cfg["n_windows"],
                "width_hours": cfg["width_hours"],
                "concat_hyperparams": cfg["concat_hyperparams"],
                "n_by_T": out["n_by_T"],
                "span_hours_by_T": out["span_hours_by_T"],
                "delta_Y": out["delta_Y_image_nsfw"],
                "direction": out["direction"],
                "hp_vimp_share": out["hp_vimp_share"],
                "hgb_auc_late": hgb.get("auc"),
                "top_blend": [r["feature"] for r in (out.get("top_blend") or [])[:8]],
                "top_cov": [r["feature"] for r in (out.get("top_covariate") or [])[:5]],
            }
        )

    blob = {
        "n_sample": args.n_sample,
        "sec": float(time.time() - t0),
        "claim": (
            "Window scheme/size is part of the estimand: equal_count vs equal_time "
            "vs fixed width change n_by_T, calendar span, ΔȲ and selected tokens. "
            "Concatenating cfg/step/sampler into X surfaces hyperparam drift in VIMP."
        ),
        "rows": rows,
    }
    (args.out / "summary.json").write_text(json.dumps(blob, indent=2, default=str) + "\n")
    lines = [
        "# DiffusionDB window-scheme sweep",
        "",
        blob["claim"],
        "",
        "| config | ΔȲ | HGB AUC | hp_vimp | n_by_T | span_h | top_blend |",
        "|---|---:|---:|---:|---|---|---|",
    ]
    for r in rows:
        if not r.get("ok"):
            lines.append(f"| {r['name']} | — | — | — | fail | — | {r.get('reason')} |")
            continue
        lines.append(
            f"| `{r['name']}` | {r['delta_Y']:.4f} | {r['hgb_auc_late']:.3f} | "
            f"{r['hp_vimp_share']:.3f} | `{r['n_by_T']}` | `{r['span_hours_by_T']}` | "
            f"{', '.join(r['top_blend'][:5])} |"
        )
    lines += [
        "",
        "## Read",
        "- `equal_count_*`: balanced n, unequal calendar width",
        "- `equal_time_*`: equal calendar span, unequal n",
        "- `width_*`: fixed hours; early/late = first/last occupied bin",
        "- `*_hp`: TF-IDF ⊕ cfg/step/sampler — check `hp_vimp_share` and top_cov",
        "",
    ]
    (args.out / "WINDOW_SWEEP_REPORT.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({"n_ok": sum(1 for r in rows if r.get("ok")), "out": str(args.out)}, indent=2))


if __name__ == "__main__":
    main()
