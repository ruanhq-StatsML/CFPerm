#!/usr/bin/env python3
"""Predict adapt FLOPs ratio from modality gradient correlation (green).

Uses existing grad-cos trajectories when present; also emits a synthetic
correlation sweep to show the closed-form predictor.

  PYTHONPATH=. python3 scripts/run_agod_compute_ratio_pred.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from agod.compute_ratio import (
    characterize_compute_ratio,
    predict_flops_from_correlation,
)

OUT = ROOT / "results" / "agod_compute_ratio"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_compute_ratio")
GRADCOS = ROOT / "results" / "agod_gradcos_lr" / "agod_gradcos_lr.json"


def synthetic_sweep(mods=("video", "text", "audio"), n: int = 21):
    """Sweep mean pairwise cos → predicted flops_rel (equal α)."""
    alpha = {m: 1.0 / len(mods) for m in mods}
    rows = []
    for c in np.linspace(-0.2, 0.95, n):
        pair = {
            f"{a}|{b}": float(c)
            for i, a in enumerate(mods)
            for b in mods[i + 1 :]
        }
        pred = predict_flops_from_correlation(alpha, pair, mods, rho_skip=0.70)
        rows.append(
            {
                "pair_cos": float(c),
                "mean_redundancy": pred["mean_redundancy"],
                "flops_rel_pred": pred["flops_rel_pred"],
                "flops_rel_proxy": pred["flops_rel_proxy"],
                "n_keep": pred["n_keep"],
                "keep": pred["keep"],
            }
        )
    return rows


def from_gradcos_json(path: Path) -> dict:
    if not path.exists():
        return {}
    payload = json.loads(path.read_text())
    out = {}
    traj = payload.get("trajectory") or {}
    for key, rows in traj.items():
        # key like "msrvtt:soft_gradcos"
        if ":" not in key:
            continue
        ds, sched = key.split(":", 1)
        if not rows:
            continue
        mods = list(rows[0]["alpha"].keys())
        # map field names from gradcos runner
        norm = []
        for r in rows:
            norm.append(
                {
                    "alpha": r.get("alpha", {}),
                    "mean_pair_cos": r.get("mean_pair_cos", 0.0),
                    "pair_cos": r.get("pair_cos"),
                    "flops_rel": r.get("flops_rel", 1.0),  # equal/soft often 1.0
                    "active": r.get("active"),
                }
            )
        out[key] = characterize_compute_ratio(norm, mods)
        out[key]["dataset"] = ds
        out[key]["scheduler"] = sched
    return out


def plot_sweep(rows, path: Path):
    xs = [r["pair_cos"] for r in rows]
    fig, ax = plt.subplots(figsize=(7.2, 4.2), facecolor="#f7f5f1")
    ax.plot(xs, [r["flops_rel_proxy"] for r in rows], label="flops_rel proxy (closed form)", lw=2)
    ax.plot(xs, [r["flops_rel_pred"] for r in rows], label="flops_rel pred (keep-set)", lw=2, ls="--")
    ax.axhline(1.0, color="#999", ls=":", lw=0.8)
    ax.set_xlabel("pairwise grad cos (synthetic, equal α)")
    ax.set_ylabel("predicted flops_rel")
    ax.set_title("Green compute: high modality correlation → lower adapt FLOPs")
    ax.legend(frameon=False, fontsize=8)
    ax.set_ylim(0.4, 1.05)
    fig.text(
        0.5,
        0.02,
        "FWD always on · skip adapt-BWD on collinear (high-cos) modalities · savings ∝ ρ̄·c_pb",
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def write_docs(sweep, from_real, path: Path):
    # pick a few sweep anchors
    anchors = []
    for target in (-0.2, 0.0, 0.5, 0.7, 0.9):
        r = min(sweep, key=lambda x: abs(x["pair_cos"] - target))
        anchors.append(
            f"| {r['pair_cos']:+.2f} | {r['mean_redundancy']:.2f} | "
            f"{r['flops_rel_proxy']:.3f} | {r['flops_rel_pred']:.3f} | {r['n_keep']} |"
        )
    real_rows = []
    for k, s in sorted(from_real.items()):
        real_rows.append(
            f"| {s.get('dataset','')} | {s.get('scheduler','')} | "
            f"{s.get('mean_redundancy', float('nan')):.3f} | "
            f"{s.get('mean_flops_rel_proxy', float('nan')):.3f} | "
            f"{s.get('mean_flops_rel_pred', float('nan')):.3f} | "
            f"{s.get('mean_keep_frac', float('nan')):.2f} |"
        )
    md = f"""# Computational ratio from modality correlation (green)

## Justification

Adapt efficiency claim is **adapt FLOPs**, not inference latency:

```
C_fwd = |M|·c_pf + c_sf          # always paid (FWD on)
C_bwd(A) = |A|·c_pb + c_sb       # A = adapt-active set
flops_rel(A) = (C_fwd + C_bwd(A)) / (C_fwd + C_bwd(M))
```

When modality gradients are **highly correlated** (`cos(g_m,g_{{m'}})↑`):

1. updates are nearly **collinear** → second tower's proj-BWD buys little new direction
2. define redundancy `ρ_m = mean_{{m'≠m}} max(cos_{{mm'}}, 0)` (conflict cos<0 does *not* justify skip)
3. unique mass `u_m = α_m·(1−η·ρ_m)_+` then keep top unique / drop collinear losers
4. **predict** `flops_rel` from the keep set *before* paying BWD

Closed-form savings proxy (no keep set needed):

```
savings_proxy ≈ ρ̄ · (|M|·c_pb) / C_full
flops_rel_proxy ≈ max(1 − savings_proxy, (C_fwd+c_sb)/C_full)
```

This is green by construction: FWD unchanged; only redundant adapt-BWD is skipped.

## Synthetic sweep (equal α)

| pair_cos | mean ρ | flops_proxy | flops_pred | n_keep |
|---|---:|---:|---:|---:|
{chr(10).join(anchors)}

## From grad-cos smoke trajectories

| Dataset | scheduler | mean ρ | flops_proxy | flops_pred | keep_frac |
|---|---|---:|---:|---:|---:|
{chr(10).join(real_rows) if real_rows else '| — | — | — | — | — | — |'}

```bash
PYTHONPATH=. python3 scripts/run_agod_compute_ratio_pred.py
```
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(md)


def write_latex(sweep, from_real, path: Path):
    lines = []
    for target in (0.0, 0.5, 0.7, 0.9):
        r = min(sweep, key=lambda x: abs(x["pair_cos"] - target))
        lines.append(
            f"{r['pair_cos']:+.2f} & {r['mean_redundancy']:.2f} & "
            f"{r['flops_rel_proxy']:.3f} & {r['flops_rel_pred']:.3f} & {r['n_keep']} \\\\"
        )
    tex = (
        "% Computational ratio vs modality gradient correlation\n"
        "\\begin{table}[t]\\centering\n"
        "\\caption{Predicted adapt FLOPs ratio from pairwise gradient cosine. "
        "FWD always on; high correlation justifies skipping redundant proj-BWD.}\n"
        "\\label{tab:agod-compute-ratio-corr}\n"
        "\\begin{tabular}{rrrrr}\\toprule\n"
        "pair\\_cos & mean $\\rho$ & flops\\_proxy & flops\\_pred & $n_{\\mathrm{keep}}$ \\\\\n"
        "\\midrule\n"
        + "\n".join(lines)
        + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(tex)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)
    sweep = synthetic_sweep()
    real = from_gradcos_json(GRADCOS)
    payload = {
        "agod_version": "0.1.0",
        "focus": "predict adapt flops_rel from modality gradient correlation",
        "accounting": {
            "C_fwd": "|M|*c_pf + c_sf (always on)",
            "C_bwd": "|A|*c_pb + c_sb",
            "flops_rel": "(C_fwd+C_bwd(A))/(C_fwd+C_bwd(M))",
            "redundancy": "rho_m = mean max(cos,0) over other mods",
            "unique": "u_m = alpha_m * (1 - eta * rho_m)_+",
            "savings_proxy": "rho_bar * (|M|*c_pb) / C_full",
        },
        "synthetic_sweep": sweep,
        "from_gradcos": real,
    }
    (OUT / "agod_compute_ratio.json").write_text(json.dumps(payload, indent=2))
    plot_sweep(sweep, OUT / "AGOD_Compute_Ratio_Corr_Board.png")
    write_docs(sweep, real, OUT / "README.md")
    write_latex(sweep, real, OUT / "AGOD_compute_ratio_tables_only.tex")
    write_docs(sweep, real, DOCS / "AGOD_compute_ratio_correlation.md")
    write_latex(sweep, real, DOCS / "AGOD_compute_ratio_tables_only.tex")
    import shutil

    shutil.copy2(
        OUT / "AGOD_Compute_Ratio_Corr_Board.png",
        ART / "AGOD_Compute_Ratio_Corr_Board.png",
    )
    shutil.copy2(OUT / "README.md", ART / "README.md")
    print("=== synthetic anchors ===")
    for target in (0.0, 0.5, 0.7, 0.9):
        r = min(sweep, key=lambda x: abs(x["pair_cos"] - target))
        print(
            f"  cos={r['pair_cos']:+.2f} rho={r['mean_redundancy']:.2f} "
            f"proxy={r['flops_rel_proxy']:.3f} pred={r['flops_rel_pred']:.3f} "
            f"keep={r['n_keep']}"
        )
    if real:
        print("=== from grad-cos traj ===")
        for k, s in sorted(real.items()):
            print(
                f"  {k}: rho={s['mean_redundancy']:.3f} "
                f"proxy={s['mean_flops_rel_proxy']:.3f} "
                f"pred={s['mean_flops_rel_pred']:.3f} "
                f"keep_frac={s['mean_keep_frac']:.2f}"
            )
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
