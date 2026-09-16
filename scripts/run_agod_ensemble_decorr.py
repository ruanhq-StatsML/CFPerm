#!/usr/bin/env python3
"""High correlation → ensemble decorrelation → per-modality LR.

When modality grads are collinear, the ensemble is near-rank-1. Instead of
only *skipping* adapt-BWD, decompose roles and reshape LR:

  leader      — owns shared update direction → raise LR
  diversifier — residual uniqueness after leader → keep / mild boost
  redundant   — collinear with leader → damp toward β

  PYTHONPATH=. python3 scripts/run_agod_ensemble_decorr.py
  PYTHONPATH=. python3 scripts/run_agod_ensemble_decorr.py --online msrvtt
  PYTHONPATH=. python3 scripts/run_agod_ensemble_decorr.py --no-online
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from agod.ensemble_decorr import (
    assign_ensemble_roles,
    characterize_ensemble_traj,
    soft_decorr_lr,
)
from agod.lr_controller import alpha_to_lr, schedule_modality_lr

OUT = ROOT / "results" / "agod_ensemble_decorr"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_ensemble_decorr")
GRADCOS = ROOT / "results" / "agod_gradcos_lr" / "agod_gradcos_lr.json"

MODS3 = ("video", "text", "audio")


def _iso_pair(c: float, mods=MODS3) -> dict[str, float]:
    return {
        f"{a}|{b}": float(c)
        for i, a in enumerate(mods)
        for b in mods[i + 1 :]
    }


def _asym_pair(leader_cos: float, other_cos: float, mods=MODS3) -> dict[str, float]:
    """Leader (mods[0]) highly aligned with others; others less aligned."""
    a, b, c = mods
    return {
        f"{a}|{b}": float(leader_cos),
        f"{a}|{c}": float(leader_cos),
        f"{b}|{c}": float(other_cos),
    }


def synthetic_sweep(mods=MODS3, n: int = 21):
    """Isotropic pair_cos sweep → roles + soft_decorr LR (equal α)."""
    alpha = {m: 1.0 / len(mods) for m in mods}
    rows = []
    for c in np.linspace(-0.1, 0.95, n):
        pair = _iso_pair(c, mods)
        packed = soft_decorr_lr(alpha, pair, mods)
        soft = alpha_to_lr(alpha, mods, beta=0.10)
        d = packed["decomp"]
        rows.append(
            {
                "pair_cos": float(c),
                "geometry": "isotropic",
                "mean_redundancy": d["mean_redundancy"],
                "effective_rank": d["spectral"]["effective_rank"],
                "decorr_active": d["decorr_active"],
                "leader": d["leader"],
                "roles": d["roles"],
                "role_gain": packed["role_gain"],
                "lr_decorr": {m: packed["lr"][m] for m in mods},
                "lr_soft": {m: soft[m] for m in mods},
                "lr_ratio_decorr": float(
                    max(packed["lr"][m] for m in mods)
                    / max(min(packed["lr"][m] for m in mods), 1e-9)
                ),
            }
        )
    return rows


def synthetic_asym_cases(mods=MODS3):
    """Asymmetric geometries: one shared leader + a residual diversifier pair."""
    alpha = {m: float(v) for m, v in zip(mods, (0.45, 0.35, 0.20))}
    cases = [
        ("low_corr", _iso_pair(0.15, mods)),
        ("high_iso", _iso_pair(0.80, mods)),
        ("leader_plus_residual", _asym_pair(0.85, 0.20, mods)),
        ("full_collinear", _iso_pair(0.92, mods)),
    ]
    out = []
    for name, pair in cases:
        packed = soft_decorr_lr(alpha, pair, mods)
        d = packed["decomp"]
        out.append(
            {
                "case": name,
                "pair_cos": pair,
                "alpha": alpha,
                "mean_redundancy": d["mean_redundancy"],
                "effective_rank": d["spectral"]["effective_rank"],
                "decorr_active": d["decorr_active"],
                "leader": d["leader"],
                "roles": d["roles"],
                "residual": d["residual"],
                "role_gain": packed["role_gain"],
                "lr": {m: packed["lr"][m] for m in mods},
            }
        )
    return out


def from_gradcos_json(path: Path) -> dict:
    if not path.exists():
        return {}
    payload = json.loads(path.read_text())
    traj = payload.get("trajectory") or {}
    out = {}
    for key, rows in traj.items():
        if ":" not in key or not rows:
            continue
        ds, sched = key.split(":", 1)
        mods = list(rows[0]["alpha"].keys())
        norm = []
        decorr_lrs = []
        for r in rows:
            alpha = r.get("alpha", {})
            pair = r.get("pair_cos") or {}
            if not pair:
                mpc = float(r.get("mean_pair_cos", 0.0))
                pair = _iso_pair(mpc, mods)
            packed = soft_decorr_lr(alpha, pair, mods)
            norm.append(
                {
                    "alpha": alpha,
                    "pair_cos": pair,
                    "mean_pair_cos": r.get("mean_pair_cos"),
                    "acc_lift": r.get("acc_lift"),
                }
            )
            decorr_lrs.append(
                {
                    "t": r.get("t"),
                    "roles": packed["decomp"]["roles"],
                    "leader": packed["decomp"]["leader"],
                    "decorr_active": packed["decomp"]["decorr_active"],
                    "effective_rank": packed["decomp"]["spectral"]["effective_rank"],
                    "role_gain": packed["role_gain"],
                    "lr_decorr": {m: packed["lr"][m] for m in mods},
                    "lr_logged": r.get("lr_mult"),
                    "acc_lift_logged": r.get("acc_lift"),
                }
            )
        summary = characterize_ensemble_traj(norm, mods)
        summary["dataset"] = ds
        summary["source_scheduler"] = sched
        summary["windows"] = decorr_lrs
        # Acc context from logged source scheduler (not soft_decorr Acc)
        lifts = [float(r["acc_lift"]) for r in rows if r.get("acc_lift") is not None]
        if lifts:
            summary["mean_acc_lift_source"] = float(np.mean(lifts))
        out[key] = summary
    return out


def plot_board(sweep, asym, path: Path):
    xs = [r["pair_cos"] for r in sweep]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0), facecolor="#f7f5f1")
    ax = axes[0]
    ax.plot(xs, [r["effective_rank"] for r in sweep], lw=2, label="effective rank")
    ax.axhline(1.55, color="#999", ls=":", lw=0.9, label="erank trigger")
    ax.set_xlabel("pairwise grad cos (isotropic)")
    ax.set_ylabel("effective rank")
    ax.set_title("High corr → near-rank-1 ensemble")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    active = [1.0 if r["decorr_active"] else 0.0 for r in sweep]
    ax.step(xs, active, where="mid", lw=2, label="decorr active")
    ax.plot(
        xs,
        [r["lr_ratio_decorr"] for r in sweep],
        lw=2,
        ls="--",
        label="LR max/min (decorr)",
    )
    ax.set_xlabel("pairwise grad cos")
    ax.set_ylabel("flag / LR ratio")
    ax.set_title("Decorrelate: damp redundant, boost leader")
    ax.legend(frameon=False, fontsize=8)

    # annotate asymmetric cases under figure
    bits = []
    for c in asym:
        roles = ",".join(f"{m[0]}:{c['roles'][m][:3]}" for m in MODS3)
        bits.append(f"{c['case']}: erank={c['effective_rank']:.2f} [{roles}]")
    fig.text(0.5, 0.02, " · ".join(bits), ha="center", fontsize=7.5)
    fig.tight_layout(rect=[0, 0.08, 1, 1])
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def write_docs(sweep, asym, from_real, online, path: Path):
    anchors = []
    for target in (0.0, 0.4, 0.55, 0.75, 0.90):
        r = min(sweep, key=lambda x: abs(x["pair_cos"] - target))
        roles = "/".join(r["roles"][m][0].upper() for m in MODS3)
        anchors.append(
            f"| {r['pair_cos']:+.2f} | {r['mean_redundancy']:.2f} | "
            f"{r['effective_rank']:.2f} | {str(r['decorr_active'])} | "
            f"{r['leader']} | {roles} | {r['lr_ratio_decorr']:.2f} |"
        )
    asym_rows = []
    for c in asym:
        roles = ", ".join(f"{m}={c['roles'][m]}" for m in MODS3)
        lrs = ", ".join(f"{m}={c['lr'][m]:.2f}" for m in MODS3)
        asym_rows.append(
            f"| `{c['case']}` | {c['mean_redundancy']:.2f} | "
            f"{c['effective_rank']:.2f} | {c['leader']} | {roles} | {lrs} |"
        )
    real_rows = []
    for k, s in sorted(from_real.items()):
        rf = s.get("role_frac") or {}
        real_rows.append(
            f"| {s.get('dataset','')} | {s.get('source_scheduler','')} | "
            f"{s.get('mean_redundancy', float('nan')):.3f} | "
            f"{s.get('mean_effective_rank', float('nan')):.2f} | "
            f"{s.get('frac_decorr_active', float('nan')):.2f} | "
            f"{rf.get('leader', 0):.2f}/{rf.get('diversifier', 0):.2f}/"
            f"{rf.get('redundant', 0):.2f} |"
        )
    online_rows = []
    for k, s in sorted((online or {}).items()):
        online_rows.append(
            f"| {s.get('dataset','')} | {s.get('scheduler','')} | "
            f"{s.get('mean_pair_cos', float('nan')):.3f} | "
            f"{s.get('frac_decorr_active', float('nan')):.2f} | "
            f"{s.get('mean_acc_lift', float('nan')):+.3f} | "
            f"{s.get('mean_lr_ratio', float('nan')):.2f} |"
        )
    md = f"""# Ensemble decorrelation under high modality correlation

## Justification

High pairwise grad cosine ⇒ modalities are **not** independent voters; the
cosine Gram is near-rank-1. Skip alone is blunt. Decorrelate the *ensemble*:

1. `effective_rank(G)` from pairwise `cos(g_m,g_{{m'}})` — near 1 ⇒ collinear
2. **leader** = argmax unique mass `α_m·(1−η·ρ_m)_+` (owns shared direction)
3. **diversifier** = residual uniqueness `r_m=(1−max(cos,0))_+` after leader
4. **redundant** = collinear with leader → damp adapt LR (FWD still on)

LR map (`soft_decorr`):

```
LR_m = α_to_lr(α_m · role_gain_m)
role_gain: leader↑ · diversifier·(1+boost·r_m) · redundant→floor
```

When `mean ρ` low / effective rank high → decorr inactive → reduces to soft LR.

## Synthetic isotropic sweep

| pair_cos | mean ρ | erank | decorr | leader | roles(v/t/a) | LR ratio |
|---|---:|---:|---|---|---|---:|
{chr(10).join(anchors)}

## Asymmetric cases (unequal α)

| case | mean ρ | erank | leader | roles | LR |
|---|---:|---:|---|---|---|
{chr(10).join(asym_rows)}

## From grad-cos trajectories (role replay)

| Dataset | source sched | mean ρ | erank | frac decorr | role frac L/D/R |
|---|---|---:|---:|---:|---|
{chr(10).join(real_rows) if real_rows else '| — | — | — | — | — | — |'}

## Online smoke (optional)

| Dataset | scheduler | mean pair_cos | frac decorr | mean Acc↑ | mean LR ratio |
|---|---|---:|---:|---:|---:|
{chr(10).join(online_rows) if online_rows else '| — | — | — | — | — | — |'}

```bash
PYTHONPATH=. python3 scripts/run_agod_ensemble_decorr.py
PYTHONPATH=. python3 scripts/run_agod_ensemble_decorr.py --online msrvtt
```
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(md)


def write_latex(sweep, asym, path: Path):
    lines = []
    for target in (0.0, 0.55, 0.75, 0.90):
        r = min(sweep, key=lambda x: abs(x["pair_cos"] - target))
        lines.append(
            f"{r['pair_cos']:+.2f} & {r['mean_redundancy']:.2f} & "
            f"{r['effective_rank']:.2f} & "
            f"{'Y' if r['decorr_active'] else 'N'} & "
            f"{r['leader']} & {r['lr_ratio_decorr']:.2f} \\\\"
        )
    asym_lines = []
    for c in asym:
        case = c["case"].replace("_", r"\_")
        asym_lines.append(
            f"{case} & {c['mean_redundancy']:.2f} & "
            f"{c['effective_rank']:.2f} & {c['leader']} & "
            f"{c['roles']['video'][0]}/{c['roles']['text'][0]}/{c['roles']['audio'][0]} \\\\"
        )
    tex = (
        "% Ensemble decorrelation vs modality gradient correlation\n"
        "\\begin{table}[t]\\centering\n"
        "\\caption{High correlation triggers ensemble role split "
        "(leader / diversifier / redundant) and soft\\_decorr LR.}\n"
        "\\label{tab:agod-ensemble-decorr}\n"
        "\\begin{tabular}{rrrlrl}\\toprule\n"
        "pair\\_cos & mean $\\rho$ & erank & decorr & leader & LR ratio \\\\\n"
        "\\midrule\n"
        + "\n".join(lines)
        + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n\n"
        "\\begin{table}[t]\\centering\n"
        "\\caption{Asymmetric geometries: residual diversifier vs full collinear.}\n"
        "\\label{tab:agod-ensemble-decorr-asym}\n"
        "\\begin{tabular}{lrrll}\\toprule\n"
        "case & mean $\\rho$ & erank & leader & roles (v/t/a) \\\\\n"
        "\\midrule\n"
        + "\n".join(asym_lines)
        + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(tex)


def run_online_msrvtt(device: str = "cpu"):
    """Compact MSR-VTT smoke: equal / soft / soft_decorr with pair_cos wired."""
    # reuse building blocks from gradcos runner
    sys.path.insert(0, str(ROOT / "scripts"))
    import run_agod_gradcos_lr as g

    schedulers = ("equal", "soft", "soft_decorr")
    cells = {}
    traj_all = {}
    for sched in schedulers:
        feats, y = g.load_msrvtt()
        mods = list(g.MSRVTT_MODS)
        stream = g.make_stream_y(y, g.N_CUR_MSRVTT)
        import torch

        torch.manual_seed(g.SEED)
        np.random.seed(g.SEED)
        model = g.Fusion({m: feats[m].shape[1] for m in mods}, mods).to(device)
        opt = g.build_optim(model, mods)
        router = g.EMARouter(mods, ema=g.EMA)
        ref_idx = stream["ref_idx"].copy()
        ref_b, y_ref = g.blocks(feats, y, ref_idx, mods)
        warm, _ = g.split_hold(ref_idx, g.SEED)
        g.train_window(
            model,
            opt,
            feats,
            y,
            warm,
            mods,
            device,
            lr_mult={**{m: 1.0 for m in mods}, "shared": 1.0},
            steps=g.STEPS_M,
        )
        traj = []
        decorr_flags = []
        for w in stream["windows"]:
            adapt, hold = g.split_hold(w["idx"], g.SEED + 13 * w["t"])
            cur_b, y_cur = g.blocks(feats, y, adapt, mods)
            msg = g.domain_msg(
                ref_b, cur_b, y_ref, y_cur, mods, seed=g.SEED + 10 * w["t"]
            )
            raw, _ = g.select_b5_raw(
                msg, ref_b, cur_b, y_ref, y_cur, mods, seed=g.SEED + 10 * w["t"]
            )
            alpha = router.update(raw)
            grads, shared, _ = g.probe_grads(model, feats, y, adapt, mods, device)
            ginfo = g.modality_grad_cosine(grads, mods, shared=shared)
            pair = ginfo.get("pair_cos") or _iso_pair(ginfo["mean_pair_cos"], mods)
            decomp = assign_ensemble_roles(alpha, pair, mods)
            lr_mult = schedule_modality_lr(
                sched,
                alpha,
                mods,
                t=w["t"],
                t_max=g.T_WIN,
                beta=g.BETA,
                align_gain=ginfo.get("align_gain"),
                pair_cos=pair,
            )
            pre = g.eval_acc(model, feats, y, hold, mods, device)
            loss, wall = g.train_window(
                model,
                opt,
                feats,
                y,
                adapt,
                mods,
                device,
                lr_mult=lr_mult,
                steps=g.STEPS_M,
            )
            post = g.eval_acc(model, feats, y, hold, mods, device)
            dacc = post["acc"] - pre["acc"]
            disp = g.soft_lr_dispersion(lr_mult, mods)
            keep = g.N_REF // 2
            ref_idx = np.concatenate(
                [ref_idx[-keep:], adapt[: min(keep, len(adapt))]]
            )
            ref_b, y_ref = g.blocks(feats, y, ref_idx, mods)
            decorr_flags.append(1.0 if decomp["decorr_active"] else 0.0)
            row = {
                "t": w["t"],
                "alpha": {m: float(alpha[m]) for m in mods},
                "lr_mult": {k: float(v) for k, v in lr_mult.items()},
                "acc_pre": pre["acc"],
                "acc_post": post["acc"],
                "acc_lift": dacc,
                "train_loss": loss,
                "wall_ms": wall,
                "mean_pair_cos": ginfo["mean_pair_cos"],
                "pair_cos": {k: float(v) for k, v in pair.items()},
                "roles": decomp["roles"],
                "leader": decomp["leader"],
                "decorr_active": decomp["decorr_active"],
                "effective_rank": decomp["spectral"]["effective_rank"],
                **disp,
            }
            traj.append(row)
            print(
                f"[msrvtt/{sched}] t={w['t']} cos={ginfo['mean_pair_cos']:+.2f} "
                f"decorr={decomp['decorr_active']} leader={decomp['leader']} "
                f"acc {pre['acc']:.3f}->{post['acc']:.3f} (d={dacc:+.3f})",
                flush=True,
            )
        lifts = [r["acc_lift"] for r in traj]
        cells[f"msrvtt:{sched}"] = {
            "dataset": "msrvtt",
            "scheduler": sched,
            "mean_pair_cos": float(np.mean([r["mean_pair_cos"] for r in traj])),
            "frac_decorr_active": float(np.mean(decorr_flags)),
            "mean_acc_lift": float(np.mean(lifts)),
            "mean_lr_ratio": float(np.mean([r["lr_ratio"] for r in traj])),
            "mean_effective_rank": float(
                np.mean([r["effective_rank"] for r in traj])
            ),
        }
        traj_all[f"msrvtt:{sched}"] = traj
    return cells, traj_all


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--online",
        nargs="?",
        const="msrvtt",
        default="msrvtt",
        help="online smoke dataset (default: msrvtt). Use --no-online to skip.",
    )
    ap.add_argument("--no-online", action="store_true")
    ap.add_argument("--device", default="cpu")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)

    sweep = synthetic_sweep()
    asym = synthetic_asym_cases()
    real = from_gradcos_json(GRADCOS)

    online_cells, online_traj = {}, {}
    if not args.no_online and args.online:
        if args.online == "msrvtt":
            online_cells, online_traj = run_online_msrvtt(device=args.device)
        else:
            print(f"online dataset {args.online!r} not implemented; skipping")

    payload = {
        "agod_version": "0.1.0",
        "focus": "high correlation → ensemble decorrelation → soft_decorr LR",
        "roles": ["leader", "diversifier", "redundant"],
        "design": {
            "trigger": "mean_rho >= rho_high OR effective_rank <= erank_collinear",
            "leader": "argmax unique mass after redundancy discount",
            "diversifier": "residual uniqueness after removing leader direction",
            "redundant": "collinear with leader → damp LR (FWD on)",
            "scheduler": "soft_decorr",
        },
        "synthetic_isotropic": sweep,
        "synthetic_asymmetric": asym,
        "from_gradcos": {
            k: {kk: vv for kk, vv in s.items() if kk != "windows"}
            for k, s in real.items()
        },
        "from_gradcos_windows": {k: s.get("windows") for k, s in real.items()},
        "online_cells": online_cells,
        "online_trajectory": online_traj,
    }
    (OUT / "agod_ensemble_decorr.json").write_text(json.dumps(payload, indent=2))
    plot_board(sweep, asym, OUT / "AGOD_Ensemble_Decorr_Board.png")
    write_docs(sweep, asym, real, online_cells, OUT / "README.md")
    write_latex(sweep, asym, OUT / "AGOD_ensemble_decorr_tables_only.tex")
    write_docs(sweep, asym, real, online_cells, DOCS / "AGOD_ensemble_decorr.md")
    write_latex(sweep, asym, DOCS / "AGOD_ensemble_decorr_tables_only.tex")
    shutil.copy2(
        OUT / "AGOD_Ensemble_Decorr_Board.png", ART / "AGOD_Ensemble_Decorr_Board.png"
    )
    shutil.copy2(OUT / "README.md", ART / "README.md")

    print("=== synthetic anchors ===")
    for target in (0.0, 0.55, 0.75, 0.90):
        r = min(sweep, key=lambda x: abs(x["pair_cos"] - target))
        print(
            f"  cos={r['pair_cos']:+.2f} rho={r['mean_redundancy']:.2f} "
            f"erank={r['effective_rank']:.2f} decorr={r['decorr_active']} "
            f"leader={r['leader']} lr_ratio={r['lr_ratio_decorr']:.2f}"
        )
    print("=== asymmetric ===")
    for c in asym:
        print(
            f"  {c['case']}: erank={c['effective_rank']:.2f} "
            f"roles={c['roles']} lr={[round(c['lr'][m],2) for m in MODS3]}"
        )
    if real:
        print("=== from grad-cos traj ===")
        for k, s in sorted(real.items()):
            print(
                f"  {k}: rho={s['mean_redundancy']:.3f} "
                f"erank={s['mean_effective_rank']:.2f} "
                f"decorr_frac={s['frac_decorr_active']:.2f}"
            )
    if online_cells:
        print("=== online smoke ===")
        for k, s in sorted(online_cells.items()):
            print(
                f"  {k}: cos={s['mean_pair_cos']:.3f} "
                f"decorr={s['frac_decorr_active']:.2f} "
                f"acc_lift={s['mean_acc_lift']:+.3f} "
                f"lr_ratio={s['mean_lr_ratio']:.2f}"
            )
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
