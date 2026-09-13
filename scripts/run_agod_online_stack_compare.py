#!/usr/bin/env python3
"""Attribution-guided online stacking vs baselines.

Datasets: Amazon, MSR-VTT, COCO ImgTxt, Affec, Fashion-IQ, Food-101 (HF CLIP pack).

Variants (same towers / steps / holdout; only the weight socket changes):
  mean_ce / stack_ce / stack_alpha / stack_uniform / stack_fixed / stack_temp / stack_erank
  (see agod.online_stack.WEIGHT_MODES)

alpha = EMA of attribution (MSG-B3 / B5 hybrid / ImgTxt hybrid / Affec VIMP).
Primary metrics: holdout Brier/MSE drop and Acc lift.

  PYTHONPATH=. python3 scripts/run_agod_online_stack_compare.py
  PYTHONPATH=. python3 scripts/run_agod_online_stack_compare.py \\
      --datasets food101 fashion_iq --merge-existing
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import accuracy_score
from torchvision import transforms

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import run_agod_amazon_modality_lr as amazon
import run_agod_gradcos_lr as msrvtt
import run_agod_imgtxt_mmd_lr as imgtxt
from agod.lr_controller import EMARouter
from agod.online_stack import (
    MeanFusion,
    StackFusion,
    WEIGHT_MODES,
    freezes_stack_psi,
    probs_mse,
    stack_temperature_for,
    stack_weight_aux,
    uses_stack_fusion,
)

OUT = ROOT / "results" / "agod_online_stack_compare"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_online_stack_compare")
AFFEC_CACHE = ROOT / "results" / "affec_fsds" / "affec_fsds_xyw_cache.npz"

SEED = 2026
VARIANTS = WEIGHT_MODES
LAMBDA_KL = 0.50
FUSE = 128
HOLD = 0.35
DATASETS_ALL = ("msrvtt", "amazon", "coco", "affec", "fashion_iq", "food101")
IMGTXT_PACKS = {
    "coco": "coco_outdoor_indoor",
    "fashion_iq": "fashion_iq",
    "food101": "food101",
}


def split_hold(idx, seed):
    rng = np.random.default_rng(seed)
    idx = np.asarray(idx)
    perm = rng.permutation(len(idx))
    n_h = max(16, min(len(idx) // 2, int(len(idx) * HOLD)))
    return idx[perm[n_h:]], idx[perm[:n_h]]


def to_batch(feats, y, idx, mods, device):
    batch = {
        m: torch.tensor(feats[m][idx], dtype=torch.float32, device=device) for m in mods
    }
    yy = torch.tensor(np.asarray(y)[idx], dtype=torch.long, device=device)
    return batch, yy


@torch.no_grad()
def eval_hold(model, feats, y, idx, mods, device, batch_size):
    model.eval()
    if len(idx) < 8:
        return {"acc": float("nan"), "mse": float("nan"), "n": len(idx)}
    logits, ys = [], []
    for s in range(0, len(idx), batch_size):
        sl = idx[s : s + batch_size]
        batch, yy = to_batch(feats, y, sl, mods, device)
        logits.append(model(batch).cpu().numpy())
        ys.append(yy.cpu().numpy())
    L, Y = np.concatenate(logits), np.concatenate(ys)
    return {
        "acc": float(accuracy_score(Y, L.argmax(1))),
        "mse": float(probs_mse(L, Y)),
        "n": int(len(Y)),
    }


def train_window(
    model,
    opt,
    feats,
    y,
    idx,
    mods,
    device,
    *,
    kind,
    steps,
    batch_size,
    alpha,
    seed,
):
    model.train()
    crit = nn.CrossEntropyLoss()
    losses, kls = [], []
    rng = np.random.default_rng(seed)
    t0 = time.perf_counter()
    for _ in range(steps):
        sel = rng.choice(
            idx, size=min(batch_size, len(idx)), replace=len(idx) < batch_size
        )
        batch, yy = to_batch(feats, y, sel, mods, device)
        if kind == "mean_ce":
            logits = model(batch)
            loss = crit(logits, yy)
            kl_v = float("nan")
        else:
            fixed = alpha if kind == "stack_fixed" and alpha is not None else None
            logits, hiddens, _lm, stack_w = model(
                batch, return_parts=True, fixed_w=fixed
            )
            loss = crit(logits, yy)
            pack = stack_weight_aux(
                kind,
                stack_w=stack_w,
                alpha=alpha,
                mods=mods,
                hiddens=hiddens,
                lambda_kl=LAMBDA_KL,
            )
            loss = loss + pack["loss"]
            kl_raw = pack.get("kl", float("nan"))
            kl_v = (
                float(kl_raw.detach().cpu())
                if hasattr(kl_raw, "detach")
                else float(kl_raw)
            )
        opt.zero_grad(set_to_none=True)
        loss.backward()
        opt.step()
        losses.append(float(loss.item()))
        kls.append(kl_v)
    arr = np.asarray(kls, float)
    return {
        "train_loss": float(np.mean(losses)),
        "mean_kl": float(np.nanmean(arr)) if np.isfinite(arr).any() else float("nan"),
        "wall_ms": (time.perf_counter() - t0) * 1000.0,
    }


def alpha_msrvtt(feats, y, ref_idx, adapt_idx, mods, *, seed):
    b0 = {m: feats[m][ref_idx] for m in mods}
    b1 = {m: feats[m][adapt_idx] for m in mods}
    y0 = y[ref_idx].astype(float)
    y1 = y[adapt_idx].astype(float)
    try:
        msg = msrvtt.domain_msg(b0, b1, y0, y1, mods, seed=seed)
        raw, _ = msrvtt.select_b5_raw(msg, b0, b1, y0, y1, mods, seed=seed)
        return {m: float(raw[m]) for m in mods}
    except Exception:
        return {m: 1.0 / len(mods) for m in mods}


def alpha_amazon(feats, y, ref_idx, adapt_idx, mods, *, seed):
    b0 = {m: feats[m][ref_idx] for m in mods}
    b1 = {m: feats[m][adapt_idx] for m in mods}
    y0 = y[ref_idx].astype(float)
    y1 = y[adapt_idx].astype(float)
    try:
        msg = amazon.compute_msg(b0, b1, y0, y1, seed=seed)
        return {m: float(msg.alpha["B3"][m]) for m in mods}
    except Exception:
        return {m: 1.0 / len(mods) for m in mods}


def alpha_imgtxt(feats, y, ref_idx, adapt_idx, mods, *, seed):
    b0 = {m: feats[m][ref_idx] for m in mods}
    b1 = {m: feats[m][adapt_idx] for m in mods}
    y0 = y[ref_idx].astype(float)
    y1 = y[adapt_idx].astype(float)
    try:
        msg = imgtxt.domain_msg(b0, b1, y0, y1, mods, seed=seed)
        raw, _ = imgtxt.select_raw("B5", msg, b0, b1, y0, y1, mods, seed=seed)
        return {m: float(raw[m]) for m in mods}
    except Exception:
        return {m: 1.0 / len(mods) for m in mods}


def alpha_affec(feats, y, ref_idx, adapt_idx, mods, *, seed):
    """Per-modality RF domain VIMP → softmax alpha (same spirit as ImgTxt B5)."""
    return alpha_imgtxt(feats, y, ref_idx, adapt_idx, mods, seed=seed)


def load_msrvtt_pack():
    feats, y = msrvtt.load_msrvtt()
    mods = list(msrvtt.MSRVTT_MODS)
    stream = msrvtt.make_stream_y(y, msrvtt.N_CUR_MSRVTT)
    windows = [
        {"t": int(w["t"]), "idx": np.asarray(w["idx"])} for w in stream["windows"]
    ]
    return {
        "name": "msrvtt",
        "mods": mods,
        "feats": feats,
        "y": y.astype(int),
        "ref_idx": np.asarray(stream["ref_idx"]),
        "windows": windows,
        "batch": int(msrvtt.BATCH),
        "steps": int(msrvtt.STEPS_M),
        "lr": float(msrvtt.LR0),
        "ema": float(msrvtt.EMA),
        "alpha_fn": alpha_msrvtt,
        "n_ref_keep": int(msrvtt.N_REF // 2),
    }


@torch.no_grad()
def extract_amazon_features(samples, device):
    image_tf = transforms.Compose(
        [
            transforms.Resize((224, 224)),
            transforms.ToTensor(),
            transforms.Normalize(
                mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]
            ),
        ]
    )
    model = amazon.AmazonAGOD().to(device)
    model.eval()
    feats = {"text": [], "image": []}
    ys = []
    chunk = 64
    for s in range(0, len(samples), chunk):
        rows = samples[s : s + chunk]
        blocks, y = amazon.extract_blocks(model, rows, device, image_tf)
        feats["text"].append(blocks["text"])
        feats["image"].append(blocks["image"])
        ys.append(y)
    return (
        {
            "text": np.concatenate(feats["text"]).astype(np.float32),
            "image": np.concatenate(feats["image"]).astype(np.float32),
        },
        np.concatenate(ys).astype(int),
    )


def load_amazon_pack(device):
    shards = sorted(amazon.SHARD_DIR.glob("*.tar.gz"))[:3]
    samples = amazon.load_shards(shards)
    stream = amazon.make_stream(samples)
    all_rows = list(stream["ref"])
    for w in stream["windows"]:
        all_rows.extend(w["rows"])
    uniq, seen = [], set()
    for r in all_rows:
        key = id(r)
        if key in seen:
            continue
        seen.add(key)
        uniq.append(r)
    print(f"amazon: extracting features for {len(uniq)} rows...", flush=True)
    feats, y = extract_amazon_features(uniq, device)
    id_to_i = {id(r): i for i, r in enumerate(uniq)}

    def rows_to_idx(rows):
        return np.asarray([id_to_i[id(r)] for r in rows], dtype=int)

    windows = [
        {"t": int(w["t"]), "idx": rows_to_idx(w["rows"])} for w in stream["windows"]
    ]
    return {
        "name": "amazon",
        "mods": list(amazon.MODS),
        "feats": feats,
        "y": y,
        "ref_idx": rows_to_idx(stream["ref"]),
        "windows": windows,
        "batch": 32,
        "steps": 24,
        "lr": 3e-3,
        "ema": 0.40,
        "alpha_fn": alpha_amazon,
        "n_ref_keep": int(amazon.N_REF // 2),
    }


def load_imgtxt_pack(alias: str):
    """Generic image+text pack under data/img_txt/<folder>/."""
    folder = IMGTXT_PACKS[alias]
    feats, y, meta = imgtxt.load_dataset(folder)
    # Food-101 has ~101 balanced classes: mode-vs-rest is tiny-pos.
    # Use coarse even/odd class split so the CE socket is well-posed.
    if alias == "food101":
        y_raw = np.load(imgtxt.DATA_ROOT / folder / "labels.npy").astype(np.int64)
        y = (y_raw % 2 == 0).astype(np.int64)
        meta = {
            **meta,
            "pos_rate": float(y.mean()),
            "mode_label": "even_class_id",
            "label_rule": "y = 1{class_id % 2 == 0}",
        }
    # image+text only (bbox optional; keep two-mod stack comparable to Amazon)
    mods = ["image", "text"]
    feats = {m: feats[m] for m in mods}
    stream = imgtxt.make_stream(y)
    windows = [
        {"t": int(w["t"]), "idx": np.asarray(w["idx"])} for w in stream["windows"]
    ]
    print(
        f"{alias}: folder={folder} n={meta['n']} pos_rate={meta['pos_rate']:.3f} "
        f"mode_label={meta['mode_label']}",
        flush=True,
    )
    return {
        "name": alias,
        "mods": mods,
        "feats": feats,
        "y": y.astype(int),
        "ref_idx": np.asarray(stream["ref_idx"]),
        "windows": windows,
        "batch": int(imgtxt.BATCH),
        "steps": int(imgtxt.STEPS),
        "lr": float(imgtxt.LR0),
        "ema": float(imgtxt.EMA),
        "alpha_fn": alpha_imgtxt,
        "n_ref_keep": int(imgtxt.N_REF // 2),
    }


def load_coco_pack():
    return load_imgtxt_pack("coco")


def load_affec_pack():
    if not AFFEC_CACHE.exists():
        raise FileNotFoundError(f"missing Affec cache: {AFFEC_CACHE}")
    z = np.load(AFFEC_CACHE, allow_pickle=True)
    X = z["X"].astype(np.float32)
    y_cont = z["Y"].astype(np.float64)
    block_slices = z["block_slices"].item()
    mods = [str(m) for m in z["mods"].tolist()]
    # classification socket: median split of continuous affect score
    thr = float(np.median(y_cont))
    y = (y_cont >= thr).astype(np.int64)
    feats = {}
    for m in mods:
        a, b = block_slices[m]
        feats[m] = X[:, int(a) : int(b)]
    stream = imgtxt.make_stream(y)
    windows = [
        {"t": int(w["t"]), "idx": np.asarray(w["idx"])} for w in stream["windows"]
    ]
    print(
        f"affec: n={len(y)} mods={mods} thr={thr:.3f} pos_rate={y.mean():.3f}",
        flush=True,
    )
    return {
        "name": "affec",
        "mods": mods,
        "feats": feats,
        "y": y,
        "ref_idx": np.asarray(stream["ref_idx"]),
        "windows": windows,
        "batch": 64,
        "steps": 28,
        "lr": 3e-3,
        "ema": 0.40,
        "alpha_fn": alpha_affec,
        "n_ref_keep": int(imgtxt.N_REF // 2),
    }


def run_variant(pack, device, kind):
    mods = list(pack["mods"])
    feats, y = pack["feats"], pack["y"]
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    dims = {m: int(feats[m].shape[1]) for m in mods}
    if not uses_stack_fusion(kind):
        model = MeanFusion(dims, mods, fuse=FUSE).to(device)
        opt = torch.optim.Adam(model.parameters(), lr=pack["lr"])
    else:
        tau = stack_temperature_for(kind, temp=0.5)
        model = StackFusion(dims, mods, fuse=FUSE, temperature=tau).to(device)
        if freezes_stack_psi(kind):
            # towers only — fusion weights locked to attribution alpha
            params = [p for n, p in model.named_parameters() if "stack_logits" not in n]
            opt = torch.optim.Adam(params, lr=pack["lr"])
        else:
            opt = torch.optim.Adam(model.parameters(), lr=pack["lr"])
    router = EMARouter(mods, ema=pack["ema"])
    ref_idx = pack["ref_idx"].copy()
    warm, _ = split_hold(ref_idx, SEED)
    train_window(
        model,
        opt,
        feats,
        y,
        warm,
        mods,
        device,
        kind=kind,
        steps=pack["steps"],
        batch_size=pack["batch"],
        alpha=None,
        seed=SEED,
    )

    traj = []
    for win in pack["windows"]:
        adapt, hold = split_hold(win["idx"], SEED + 13 * win["t"])
        raw = pack["alpha_fn"](
            feats, y, ref_idx, adapt, mods, seed=SEED + 10 * win["t"]
        )
        if kind == "mean_ce":
            alpha = {m: 1.0 / len(mods) for m in mods}
        else:
            alpha = router.update(raw)

        with torch.no_grad():
            model.eval()
            sel = adapt[: min(pack["batch"], len(adapt))]
            batch, _ = to_batch(feats, y, sel, mods, device)
            if kind == "mean_ce":
                w_np = {m: 1.0 / len(mods) for m in mods}
            else:
                _logits, _h, _lm, stack_w = model(batch, return_parts=True)
                w_np = {m: float(stack_w[i].cpu()) for i, m in enumerate(mods)}

        pre = eval_hold(model, feats, y, hold, mods, device, pack["batch"])
        tr = train_window(
            model,
            opt,
            feats,
            y,
            adapt,
            mods,
            device,
            kind=kind,
            steps=pack["steps"],
            batch_size=pack["batch"],
            alpha=alpha if kind in ("stack_alpha", "stack_uniform", "stack_fixed", "stack_erank") else None,
            seed=SEED + 7 * win["t"],
        )
        post = eval_hold(model, feats, y, hold, mods, device, pack["batch"])

        keep = pack["n_ref_keep"]
        ref_idx = np.concatenate(
            [ref_idx[-keep:], adapt[: min(keep, len(adapt))]]
        )
        prop_mae = float(np.mean([abs(w_np[m] - alpha[m]) for m in mods]))
        row = {
            "t": int(win["t"]),
            "alpha": {m: float(alpha[m]) for m in mods},
            "stack_w": w_np,
            "prop_mae": prop_mae,
            "acc_pre": pre["acc"],
            "acc_post": post["acc"],
            "acc_lift": post["acc"] - pre["acc"],
            "mse_pre": pre["mse"],
            "mse_post": post["mse"],
            "mse_drop": pre["mse"] - post["mse"],
            "train_loss": tr["train_loss"],
            "mean_kl": tr["mean_kl"],
            "wall_ms": tr["wall_ms"],
        }
        traj.append(row)
        print(
            f"[{pack['name']}/{kind}] t={win['t']} "
            f"MSE {pre['mse']:.4f}->{post['mse']:.4f} (drop={row['mse_drop']:+.4f}) "
            f"Acc {pre['acc']:.3f}->{post['acc']:.3f} (d={row['acc_lift']:+.3f}) "
            f"prop_mae={prop_mae:.3f}",
            flush=True,
        )
    return traj


def summarize(traj, *, dataset, kind, mods):
    return {
        "dataset": dataset,
        "variant": kind,
        "n_windows": len(traj),
        "mean_mse_drop": float(np.mean([r["mse_drop"] for r in traj])),
        "mean_mse_post": float(np.mean([r["mse_post"] for r in traj])),
        "mean_acc_lift": float(np.mean([r["acc_lift"] for r in traj])),
        "mean_acc_post": float(np.mean([r["acc_post"] for r in traj])),
        "mean_prop_mae": float(np.mean([r["prop_mae"] for r in traj])),
        "mean_stack_w": {
            m: float(np.mean([r["stack_w"][m] for r in traj])) for m in mods
        },
        "mean_alpha": {
            m: float(np.mean([r["alpha"][m] for r in traj])) for m in mods
        },
    }


def plot_board(cells_by_ds, path: Path):
    datasets = list(cells_by_ds.keys())
    n = len(datasets)
    if n <= 2:
        nrows, ncols = 1, max(n, 1)
        figsize = (5.2 * ncols, 4.0)
    else:
        ncols = 2
        nrows = int(np.ceil(n / 2))
        figsize = (10.5, 3.6 * nrows)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=figsize, facecolor="#f7f5f1", squeeze=False
    )
    flat = axes.ravel()
    for ax, ds in zip(flat, datasets):
        cells = cells_by_ds[ds]
        names = [c["variant"] for c in cells]
        x = np.arange(len(names))
        ax.bar(
            x - 0.15,
            [c["mean_mse_drop"] for c in cells],
            0.3,
            label="MSE drop",
            color="#b85c38",
        )
        ax.bar(
            x + 0.15,
            [c["mean_acc_lift"] for c in cells],
            0.3,
            label="Acc lift",
            color="#2f5d50",
        )
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=15, ha="right")
        ax.axhline(0, color="#999", ls=":", lw=0.9)
        ax.set_title(f"{ds}: stacking vs baselines")
        ax.legend(frameon=False, fontsize=8)
    for ax in flat[n:]:
        ax.axis("off")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def write_docs(cells, path: Path):
    rows = []
    for c in cells:
        sw = " / ".join(f"{m}:{c['mean_stack_w'][m]:.2f}" for m in c["mean_stack_w"])
        rows.append(
            f"| `{c['dataset']}` | `{c['variant']}` | {c['mean_mse_drop']:+.4f} | "
            f"{c['mean_acc_lift']:+.3f} | {c['mean_acc_post']:.3f} | "
            f"{c['mean_prop_mae']:.3f} | {sw} |"
        )
    best = max(cells, key=lambda c: (c["mean_mse_drop"], c["mean_acc_lift"]))
    md = f"""# Attribution-guided online stacking

## Method

Online stacking learns fusion weights `w = softmax(psi)` over modality logits.

Weight socket: attribution alpha (MSG-B3 / B5 hybrid / ImgTxt hybrid / Affec VIMP, EMA) plugs in as

```
L = CE(stack_w · logits, y) + lambda * KL(stack_w || alpha)
```

Classical online stacking + an FSDS/MSG attribution prior — not a new fusion architecture.

## Variants

| variant | fusion | alpha used? |
|---|---|---|
| `mean_ce` | mean-pool | no |
| `stack_ce` | learned stack_w | logged only |
| `stack_alpha` | stack_w + KL(w||alpha) | yes |

## Smoke

| dataset | variant | MSE drop | Acc lift | Acc post | |w-alpha| | mean stack_w |
|---|---|---:|---:|---:|---:|---|
{chr(10).join(rows)}

Best (MSE drop, Acc lift): **`{best['dataset']}/{best['variant']}`**
(MSE drop={best['mean_mse_drop']:+.4f}, Acc lift={best['mean_acc_lift']:+.3f}).

```bash
PYTHONPATH=. python3 scripts/run_agod_online_stack_compare.py
PYTHONPATH=. python3 scripts/run_agod_online_stack_compare.py --datasets coco affec
```
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(md)


def _merge_cells(old_cells, new_cells):
    """Keep previous dataset rows; overwrite matching dataset/variant."""
    key = lambda c: (c["dataset"], c["variant"])
    by = {key(c): c for c in old_cells}
    for c in new_cells:
        by[key(c)] = c
    order = []
    for ds in DATASETS_ALL:
        for v in VARIANTS:
            k = (ds, v)
            if k in by:
                order.append(by.pop(k))
    order.extend(by.values())
    return order


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cpu")
    ap.add_argument(
        "--datasets",
        nargs="+",
        default=["msrvtt", "amazon"],
        choices=list(DATASETS_ALL),
    )
    ap.add_argument(
        "--variants",
        nargs="+",
        default=list(VARIANTS),
        choices=list(VARIANTS),
    )
    ap.add_argument(
        "--merge-existing",
        action="store_true",
        help="merge new cells into existing JSON (keep other datasets)",
    )
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)
    device = args.device

    packs = {}
    if "msrvtt" in args.datasets:
        print("loading MSR-VTT...", flush=True)
        packs["msrvtt"] = load_msrvtt_pack()
    if "amazon" in args.datasets:
        print("loading Amazon...", flush=True)
        packs["amazon"] = load_amazon_pack(device)
    if "coco" in args.datasets:
        print("loading COCO ImgTxt...", flush=True)
        packs["coco"] = load_coco_pack()
    if "fashion_iq" in args.datasets:
        print("loading Fashion-IQ ImgTxt...", flush=True)
        packs["fashion_iq"] = load_imgtxt_pack("fashion_iq")
    if "food101" in args.datasets:
        print("loading Food-101 ImgTxt (HF CLIP pack)...", flush=True)
        packs["food101"] = load_imgtxt_pack("food101")
    if "affec" in args.datasets:
        print("loading Affec...", flush=True)
        packs["affec"] = load_affec_pack()

    cells, trajs, cells_by_ds = [], {}, {}
    for ds, pack in packs.items():
        cells_by_ds[ds] = []
        for kind in args.variants:
            traj = run_variant(pack, device, kind)
            key = f"{ds}/{kind}"
            trajs[key] = traj
            cell = summarize(
                traj, dataset=ds, kind=kind, mods=pack["mods"]
            )
            cells.append(cell)
            cells_by_ds[ds].append(cell)

    out_json = OUT / "agod_online_stack_compare.json"
    if args.merge_existing and out_json.exists():
        prev = json.loads(out_json.read_text())
        cells = _merge_cells(prev.get("cells", []), cells)
        old_traj = prev.get("trajectory", {})
        old_traj.update(trajs)
        trajs = old_traj
        cells_by_ds = {}
        for c in cells:
            cells_by_ds.setdefault(c["dataset"], []).append(c)

    payload = {
        "agod_version": "0.1.0",
        "method": "attribution-guided online stacking",
        "rule": "L = CE(stack(w), y) + lambda KL(w || alpha); alpha = EMA(FSDS/MSG)",
        "lambda_kl": LAMBDA_KL,
        "variants": list(VARIANTS),
        "datasets_run": list(args.datasets),
        "cells": cells,
        "trajectory": trajs,
    }
    out_json.write_text(json.dumps(payload, indent=2))
    plot_board(cells_by_ds, OUT / "AGOD_Online_Stack_Compare_Board.png")
    write_docs(cells, OUT / "README.md")
    write_docs(cells, DOCS / "AGOD_online_stack_compare.md")
    shutil.copy2(
        OUT / "AGOD_Online_Stack_Compare_Board.png",
        ART / "AGOD_Online_Stack_Compare_Board.png",
    )
    shutil.copy2(OUT / "README.md", ART / "README.md")

    print("=== summary ===")
    for c in cells:
        print(
            f"  {c['dataset']}/{c['variant']}: MSE drop={c['mean_mse_drop']:+.4f} "
            f"Acc+={c['mean_acc_lift']:+.3f} Acc_post={c['mean_acc_post']:.3f} "
            f"|w-a|={c['mean_prop_mae']:.3f}"
        )
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
