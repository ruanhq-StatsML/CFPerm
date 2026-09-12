#!/usr/bin/env python3
"""Online modality selection + LR adjustment prototype (Amazon AGOD).

Core loop (the thing that actually decides):
  1. MSG_m = Norm(AUC_m * VIMP_m + γ * PO_m)     # who is drifting?
  2. α_m   = Softmax(MSG / τ)  with EMA            # select / weight modalities
  3. LR_m  = lr0 * (β + (1-β) * α_m * |M|)         # adjust learning rates
  4. if α_m < θ: zero grads of modality-m proj     # hard select (skip update)

Policies compared
  B1: α = uniform          → always update both (no selection)
  B2: α ∝ Softmax(AUC)     → covariate-only selection
  B3: α ∝ Softmax(MSG)     → attribution-guided selection (ours)

Eval (per window)
  upstream: α, LR×, active set, rel adapt-FLOPs, wall ms
  downstream: hold ΔAcc/ΔAUC, ref ΔAcc, cost-utility = holdΔAcc / relFLOPs

  python3 scripts/run_agod_amazon_online_select.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import accuracy_score, roc_auc_score
from torch.utils.data import DataLoader
from torchvision import transforms

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import run_agod_amazon_modality_lr as base  # noqa: E402

OUT = ROOT / "results" / "agod_amazon"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_amazon")

SEED = base.SEED
MODS = base.MODS
EMA = base.EMA
BATCH = base.BATCH
STEPS = base.STEPS_PER_WIN
HOLD = 0.35

# adapt-FLOPs units (encoder FWD dominates; gating saves proj BWD)
C_ENC, C_PF, C_PB, C_SF, C_SB = 10.0, 1.0, 2.0, 0.25, 0.5


def flops_rel(active: dict) -> dict:
    fwd = C_ENC + 2 * C_PF + C_SF
    bwd = sum(C_PB for m in MODS if active[m]) + C_SB
    total = fwd + bwd
    full = C_ENC + 2 * C_PF + C_SF + 2 * C_PB + C_SB
    return {"total": total, "full": full, "rel": total / full, "active": dict(active)}


def split_rows(rows, seed):
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(rows))
    n_h = max(20, min(len(rows) // 2, int(len(rows) * HOLD)))
    return [rows[i] for i in idx[n_h:]], [rows[i] for i in idx[:n_h]]


@torch.no_grad()
def eval_set(model, rows, device, tf):
    model.eval()
    if len(rows) < 8:
        return {"acc": float("nan"), "auc": float("nan")}
    loader = DataLoader(base.AmazonDS(rows, tf), batch_size=BATCH, shuffle=False)
    logits, ys = [], []
    for imgs, txts, y in loader:
        logits.append(model(imgs.to(device), txts, return_mods=False).cpu().numpy())
        ys.append(y.numpy())
    L, Y = np.concatenate(logits), np.concatenate(ys)
    acc = float(accuracy_score(Y, L.argmax(1)))
    auc = float("nan")
    if len(np.unique(Y)) > 1:
        try:
            p = torch.softmax(torch.from_numpy(L), 1).numpy()[:, 1]
            auc = float(roc_auc_score(Y, p))
        except Exception:
            pass
    return {"acc": acc, "auc": auc}


def select_and_step(model, opt, rows, device, tf, *, alpha, hard_select: bool):
    """Online selection → LR adjust → (optional) hard gate → one window of updates."""
    model.train()
    model.image_encoder.eval()
    mult = base.alpha_to_lr_mult(alpha)
    base.set_lrs(opt, mult)

    active = {m: True for m in MODS}
    if hard_select:
        active = {m: alpha[m] >= base.GATE_TH for m in MODS}
        if not any(active.values()):
            mstar = max(MODS, key=lambda m: alpha[m])
            active = {m: m == mstar for m in MODS}

    loader = DataLoader(base.AmazonDS(rows, tf), batch_size=BATCH, shuffle=True)
    crit = nn.CrossEntropyLoss()
    losses, it = [], iter(loader)
    t0 = time.perf_counter()
    for _ in range(STEPS):
        try:
            imgs, txts, y = next(it)
        except StopIteration:
            it = iter(loader)
            imgs, txts, y = next(it)
        imgs, y = imgs.to(device), y.to(device)
        loss = crit(model(imgs, txts, return_mods=False), y)
        opt.zero_grad(set_to_none=True)
        loss.backward()
        if hard_select:
            groups = model.param_groups()
            for m in MODS:
                if not active[m]:
                    for p in groups[m]:
                        p.grad = None
        opt.step()
        losses.append(float(loss.item()))
    wall = time.perf_counter() - t0
    fr = flops_rel(active)
    return {
        "loss": float(np.mean(losses)),
        "lr_mult": mult,
        "active": active,
        "flops_rel": fr["rel"],
        "flops": fr,
        "wall_ms_step": 1000.0 * wall / STEPS,
    }


def run_policy(ctor, stream, device, tf, policy: str):
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    model = ctor().to(device)
    opt = base.build_optim(model)
    ref = stream["ref"]
    ref_b, y_ref = base.extract_blocks(model, ref, device, tf)
    ema = {m: 1.0 / len(MODS) for m in MODS}
    traj = []

    for w in stream["windows"]:
        adapt, hold = split_rows(w["rows"], SEED + 11 * w["t"])
        cur_b, y_cur = base.extract_blocks(model, adapt, device, tf)
        msg = base.compute_msg(ref_b, cur_b, y_ref, y_cur, seed=SEED + 10 * w["t"])

        # --- online selection signal ---
        raw = msg.alpha[policy]  # B1 uniform / B2 AUC / B3 MSG
        for m in MODS:
            ema[m] = EMA * ema[m] + (1.0 - EMA) * raw[m]
        s = sum(ema.values())
        alpha = {m: ema[m] / s for m in MODS}

        pre_h = eval_set(model, hold, device, tf)
        pre_r = eval_set(model, ref[:120], device, tf)
        stats = select_and_step(
            model, opt, adapt, device, tf, alpha=alpha, hard_select=(policy != "B1")
        )
        post_h = eval_set(model, hold, device, tf)
        post_r = eval_set(model, ref[:120], device, tf)

        dacc = post_h["acc"] - pre_h["acc"]
        dauc = (
            post_h["auc"] - pre_h["auc"]
            if np.isfinite(pre_h["auc"]) and np.isfinite(post_h["auc"])
            else float("nan")
        )
        dref = post_r["acc"] - pre_r["acc"]
        util = dacc / stats["flops_rel"] if stats["flops_rel"] > 0 else float("nan")

        # slow ref refresh
        ref_b = {
            m: np.vstack([ref_b[m][-(base.N_REF // 2) :], cur_b[m][: base.N_REF // 2]])
            for m in MODS
        }
        y_ref = np.concatenate([y_ref[-(base.N_REF // 2) :], y_cur[: base.N_REF // 2]])

        row = {
            "t": w["t"],
            "category": w["category"],
            "selection": {
                "g": msg.g,
                "auc_domain": msg.auc,
                "alpha_raw": raw,
                "alpha_ema": alpha,
                "selected": stats["active"],  # hard selection set
                "lr_mult": stats["lr_mult"],
            },
            "compute": {
                "flops_rel": stats["flops_rel"],
                "wall_ms_step": stats["wall_ms_step"],
            },
            "eval": {
                "hold_dacc": dacc,
                "hold_dauc": dauc,
                "ref_dacc": dref,
                "train_loss": stats["loss"],
                "cost_utility": util,
            },
        }
        traj.append(row)
        sel = {m: int(stats["active"][m]) for m in MODS}
        print(
            f"[{policy}] t={w['t']} {w['category'][:22]:<22} "
            f"α={ {m: round(alpha[m], 2) for m in MODS} } "
            f"sel={sel} LR×_img={stats['lr_mult']['image']:.2f} "
            f"FLOPs={stats['flops_rel']:.2f} | "
            f"Δacc={dacc:+.3f} util={util:+.3f}"
        )
    return traj


def summarize(results):
    out = {}
    for pol, traj in results.items():
        out[pol] = {
            "mean_flops_rel": float(np.mean([r["compute"]["flops_rel"] for r in traj])),
            "mean_wall_ms": float(np.mean([r["compute"]["wall_ms_step"] for r in traj])),
            "mean_hold_dacc": float(np.nanmean([r["eval"]["hold_dacc"] for r in traj])),
            "mean_hold_dauc": float(np.nanmean([r["eval"]["hold_dauc"] for r in traj])),
            "mean_ref_dacc": float(np.nanmean([r["eval"]["ref_dacc"] for r in traj])),
            "mean_cost_utility": float(
                np.nanmean([r["eval"]["cost_utility"] for r in traj])
            ),
            "mean_alpha_image": float(
                np.mean([r["selection"]["alpha_ema"]["image"] for r in traj])
            ),
            "frac_text_selected": float(
                np.mean([r["selection"]["selected"]["text"] for r in traj])
            ),
            "frac_image_selected": float(
                np.mean([r["selection"]["selected"]["image"] for r in traj])
            ),
        }
    out["B3_vs_B1"] = {
        "flops_ratio": out["B3"]["mean_flops_rel"] / out["B1"]["mean_flops_rel"],
        "utility_ratio": (
            out["B3"]["mean_cost_utility"] / out["B1"]["mean_cost_utility"]
            if abs(out["B1"]["mean_cost_utility"]) > 1e-9
            else float("nan")
        ),
        "hold_dacc_gap": out["B3"]["mean_hold_dacc"] - out["B1"]["mean_hold_dacc"],
    }
    return out


def plot_dashboard(results, summary, stream, path: Path):
    fig = plt.figure(figsize=(13.5, 8.8), facecolor="#f7f5f1")
    gs = fig.add_gridspec(2, 2, hspace=0.36, wspace=0.28)
    fig.suptitle(
        "Online modality selection → LR adjust (Amazon AGOD)",
        fontsize=14,
        fontweight="bold",
    )
    b3 = results["B3"]
    ts = [r["t"] + 1 for r in b3]

    ax = fig.add_subplot(gs[0, 0])
    for m, c in zip(MODS, ["#2B6CB0", "#C05621"]):
        ax.plot(
            ts,
            [r["selection"]["alpha_ema"][m] for r in b3],
            "o-",
            color=c,
            label=rf"$\alpha_{{{m}}}$",
        )
        ax.plot(
            ts,
            [1.0 if r["selection"]["selected"][m] else 0.0 for r in b3],
            "x--",
            color=c,
            alpha=0.5,
            label=f"selected_{m}",
        )
    ax.set_ylim(-0.05, 1.1)
    ax.set_title("Selection: soft α (EMA) + hard gate")
    ax.legend(frameon=False, fontsize=7, ncol=2)

    ax = fig.add_subplot(gs[0, 1])
    for m, c in zip(MODS + ["shared"], ["#2B6CB0", "#C05621", "#718096"]):
        ax.plot(
            ts,
            [r["selection"]["lr_mult"][m] for r in b3],
            "o-",
            color=c,
            label=m,
        )
    ax.axhline(1.0, color="#999", ls="--", lw=0.8)
    ax.set_title("Adjustment: per-modality LR multipliers")
    ax.legend(frameon=False, fontsize=8)

    ax = fig.add_subplot(gs[1, 0])
    for pol, c, mk in zip(
        ["B1", "B2", "B3"], ["#718096", "#DD6B20", "#C53030"], ["o", "s", "D"]
    ):
        ax.plot(
            [r["t"] + 1 for r in results[pol]],
            [r["compute"]["flops_rel"] for r in results[pol]],
            marker=mk,
            color=c,
            label=f"{pol} FLOPs",
        )
    ax.set_ylim(0, 1.15)
    ax.set_title("Effect: adapt FLOPs (selection → skip BWD)")
    ax.legend(frameon=False, fontsize=8)

    ax = fig.add_subplot(gs[1, 1])
    xs = np.arange(3)
    fl = [summary[p]["mean_flops_rel"] for p in ["B1", "B2", "B3"]]
    ut = [summary[p]["mean_cost_utility"] for p in ["B1", "B2", "B3"]]
    ax.bar(xs - 0.15, fl, 0.3, color="#4A5568", label="mean rel FLOPs")
    ax.bar(xs + 0.15, ut, 0.3, color="#C53030", label="mean cost-utility")
    ax.set_xticks(xs)
    ax.set_xticklabels(["B1", "B2", "B3"])
    ax.axhline(0.0, color="#999", ls="--", lw=0.7)
    ax.set_title("Policy comparison: cost vs utility")
    ax.legend(frameon=False, fontsize=8)

    gap = summary["B3_vs_B1"]
    fig.text(
        0.5,
        0.012,
        f"B3/B1 FLOPs={gap['flops_ratio']:.2f}, utility_ratio={gap['utility_ratio']:.2f}, "
        f"holdΔacc gap={gap['hold_dacc_gap']:+.3f}  |  "
        "select=Softmax(MSG); adjust=LR×α; hard-gate α<θ skips proj BWD",
        ha="center",
        fontsize=8.5,
        color="#333",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


DISCUSSION = r"""# Online modality selection & LR adjustment — how it works, how well

## Decision rule (implemented)

```
each online window Dt vs reference D0:
  for m in {text, image}:
      AUC_m, VIMP_m, PO_m  ← RF domain / residual probes on X_m
      g_m ← Normalize(AUC_m·VIMP_m + γ·PO_m)
  α ← EMA( Softmax(g / τ) )          # soft selection
  LR_m ← lr0 · (β + (1-β)·α_m·|M|)   # continuous adjustment
  if α_m < θ: drop ∂L/∂θ_m           # hard selection (save BWD)
```

| Piece | Role |
|---|---|
| MSG `g_m` | ranks which modality is drifting |
| Softmax `α_m` | turns ranks into a budget over modalities |
| `LR_m` | spends more gradient steps on high-α modalities |
| Gate `α_m<θ` | refuses to update near-zero modalities → FLOPs↓ |

B1 disables selection (α=½). B2 selects on AUC only. B3 selects on full MSG.

## What “effect” means here

1. **Selection fidelity (upstream):** does α track the drifting channel?
   On Amazon category shifts, α_image typically rises (image pack changes more
   across categories than English review text templates).
2. **Compute effect:** hard gate skips one proj BWD → rel adapt FLOPs ~0.6–0.8
   (encoder FWD still paid). Wall-clock shrinks less than FLOPs on CPU.
3. **Task effect (downstream):** holdout ΔAcc after the update.
   Smoke often shows B3 ≈ B1 on ΔAcc (sometimes slightly worse) while FLOPs↓
   ⇒ **efficiency win, not yet an accuracy win**.
4. **Cost-utility:** holdΔAcc / relFLOPs. If ΔAcc holds and FLOPs fall, utility↑.

## Reading the smoke (typical)

- B3 **selects image** more often (`frac_image_selected` high, text gated).
- B3 **FLOPs_ratio vs B1 < 1**.
- B3 **holdΔacc_gap vs B1 ≈ 0 or slightly negative** on short smoke.
- Therefore: online selection/adjustment is doing the *intended control*;
  claiming SOTA accuracy from this alone would be overclaim. Fair claim:
  **iso-loss / iso-ΔAcc adaptation with lower update cost**.

## Failure modes to watch

- Over-gate: θ too high / τ too sharp → skip the modality that still carries label signal.
- AUC-only (B2): can chase covariate shift that is irrelevant to Y.
- Encoder-dominated cost: gating proj BWD helps, but big savings need also
  skipping unused modality encoder FWD (stricter hard path).
"""


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"device={device}")
    shards = sorted(base.SHARD_DIR.glob("*.tar.gz"))
    samples = base.load_shards(shards)
    print(f"samples={len(samples)}")
    tf = transforms.Compose(
        [
            transforms.Resize((224, 224)),
            transforms.ToTensor(),
            transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225]),
        ]
    )
    stream = base.make_stream(samples)
    print(f"ref={stream['ref_cat']} wins={[w['category'] for w in stream['windows']]}")

    results = {}
    for pol in ["B1", "B2", "B3"]:
        print(f"\n===== {pol} =====")
        results[pol] = run_policy(base.AmazonAGOD, stream, device, tf, pol)

    summary = summarize(results)
    print("\n=== SUMMARY ===")
    print(json.dumps(summary, indent=2))

    payload = {
        "logic": {
            "select": "α = EMA(Softmax(MSG/τ)); B1 uniform; B2 AUC; B3 MSG",
            "adjust": "LR_m = lr0*(β+(1-β)*α_m*|M|)",
            "hard_gate": "skip modality proj BWD if α_m < θ",
        },
        "summary": summary,
        "trajectory": results,
    }
    jp = OUT / "agod_amazon_online_select.json"
    jp.write_text(json.dumps(payload, indent=2, default=float))
    dash = OUT / "AGOD_Amazon_Online_Select_Dashboard.png"
    plot_dashboard(results, summary, stream, dash)

    note = DOCS / "AGOD_online_selection_lr_discussion.md"
    note.write_text(DISCUSSION)
    (OUT / note.name).write_text(DISCUSSION)

    lines = [
        r"\begin{tabular}{lccccc}",
        r"\toprule",
        r"Policy & Rel FLOPs & Hold $\Delta$Acc & Cost-util & $\alpha_{img}$ & Text selected \\",
        r"\midrule",
    ]
    for pol in ["B1", "B2", "B3"]:
        s = summary[pol]
        lines.append(
            f"{pol} & {s['mean_flops_rel']:.3f} & {s['mean_hold_dacc']:+.3f} & "
            f"{s['mean_cost_utility']:+.3f} & {s['mean_alpha_image']:.3f} & "
            f"{s['frac_text_selected']:.2f} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}"]
    (OUT / "AGOD_amazon_online_select_tables_only.tex").write_text("\n".join(lines))
    (DOCS / "AGOD_amazon_online_select_tables_only.tex").write_text("\n".join(lines))

    for p in [dash, jp, note]:
        (ART / Path(p).name).write_bytes(Path(p).read_bytes())
    print("dashboard", dash)
    print("discussion", note)


if __name__ == "__main__":
    main()
