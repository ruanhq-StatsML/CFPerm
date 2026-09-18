#!/usr/bin/env python3
"""Image-OOD bench: PO-risk is NOT an obs-level detector.

Obs-level PO-risk AUROC is too noisy — if you want binary OOD detection,
an RF binary classifier (oracle ID vs OOD) is the right tool.

This bench therefore:
  1. Primary: **batch-mean** AUROC (OnlineRFPerm unit) + **hard-rank**
  2. Secondary: obs-level AUROC only as a noisy diagnostic
  3. Ceiling: oracle RF-binary (trained with OOD labels) — not a fair PO rival

  PYTHONPATH=. python3 scripts/run_agod_image_ood_bench.py \\
    --packs food101_vit coco_outdoor_indoor coco_time_order indiana_cxr fashion_iq
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from agod.hard_rank_metrics import hard_rank_metrics
from agod.image_ood import (
    batch_aggregate_scores,
    binary_ood_metrics,
    eval_scores_on_split,
    fit_mahalanobis,
    fit_po_classifier,
    fit_po_regressor,
    fit_rf_binary,
    list_available_packs,
    load_image_ood_pack,
    score_centroid,
    score_knn,
    score_mahalanobis,
    score_po_energy,
    score_po_msp,
    score_po_nll,
    score_po_resid,
    score_rf_binary,
)

# Methods we compare as *deployable* (no OOD labels at train)
DEPLOY = ("po_msp", "po_energy", "maha", "knn", "centroid")
# Oracle ceiling (needs OOD labels)
ORACLE = ("rf_binary",)
SCORE_COLORS = {
    "po_msp": "#5E81AC",
    "po_energy": "#EBCB8B",
    "maha": "#B48EAD",
    "knn": "#A3BE8C",
    "centroid": "#D08770",
    "rf_binary": "#BF616A",
}


def _batch_metrics(scores: Dict[str, np.ndarray], y: np.ndarray, *, batch_size: int, seed: int):
    out = {}
    for name, sc in scores.items():
        bs, by = batch_aggregate_scores(sc, y, batch_size=batch_size, seed=seed)
        out[name] = {
            **binary_ood_metrics(by, bs),
            "n_batches": int(len(by)),
            "batch_size": int(batch_size),
        }
    return out


def run_pack(
    root: Path,
    pack: str,
    *,
    pca_d: int,
    seed: int,
    extract_max_n: int,
    batch_size: int,
) -> dict:
    kwargs = {"pca_d": pca_d, "seed": seed}
    if pack == "food101_vit":
        kwargs["extract_max_n"] = extract_max_n
    split = load_image_ood_pack(root, pack, **kwargs)

    clf = fit_po_classifier(split.X_id_train, split.y_id_train, seed=seed)
    reg = fit_po_regressor(split.X_id_train, split.y_id_train, seed=seed)
    maha = fit_mahalanobis(split.X_id_train, split.y_id_train)

    X_eval = np.vstack([split.X_id_test, split.X_ood])
    y_bin = np.concatenate(
        [np.zeros(len(split.X_id_test), int), np.ones(len(split.X_ood), int)]
    )
    y_eval = np.concatenate([split.y_id_test, split.y_ood])

    scores = {
        "po_msp": np.concatenate(
            [score_po_msp(clf, split.X_id_test), score_po_msp(clf, split.X_ood)]
        ),
        "po_energy": np.concatenate(
            [score_po_energy(clf, split.X_id_test), score_po_energy(clf, split.X_ood)]
        ),
        "po_nll": np.concatenate(
            [
                score_po_nll(clf, split.X_id_test, split.y_id_test),
                score_po_nll(clf, split.X_ood, split.y_ood),
            ]
        ),
        "maha": np.concatenate(
            [
                score_mahalanobis(maha, split.X_id_test),
                score_mahalanobis(maha, split.X_ood),
            ]
        ),
        "knn": score_knn(split.X_id_train, X_eval, k=5),
        "centroid": score_centroid(split.X_id_train, X_eval, split.y_id_train),
    }

    # Oracle RF-binary ceiling: train on ID-train + half OOD, eval on ID-test + held-out OOD
    rf_bin, ood_tr_idx, ood_te_idx = fit_rf_binary(
        split.X_id_train, split.X_ood, seed=seed, ood_train_frac=0.5
    )
    # Build eval slice that excludes OOD rows used for RF training
    X_ood_te = split.X_ood[ood_te_idx]
    y_ood_te = split.y_ood[ood_te_idx]
    X_rf_eval = np.vstack([split.X_id_test, X_ood_te])
    y_rf_bin = np.concatenate(
        [np.zeros(len(split.X_id_test), int), np.ones(len(X_ood_te), int)]
    )
    rf_score = score_rf_binary(rf_bin, X_rf_eval)

    # Obs-level (noisy diagnostic) — deployable methods on full eval
    obs_metrics = eval_scores_on_split(split, {k: scores[k] for k in DEPLOY})
    # RF-binary obs metrics on its held-out eval
    obs_metrics["rf_binary"] = binary_ood_metrics(y_rf_bin, rf_score)
    obs_metrics["rf_binary"]["note"] = "oracle: trained with OOD labels"

    # Batch-level (primary for PO / OnlineRFPerm)
    batch_metrics = _batch_metrics(
        {k: scores[k] for k in DEPLOY}, y_bin, batch_size=batch_size, seed=seed
    )
    # RF-binary batches on its own eval slice
    bs, by = batch_aggregate_scores(rf_score, y_rf_bin, batch_size=batch_size, seed=seed)
    batch_metrics["rf_binary"] = {
        **binary_ood_metrics(by, bs),
        "n_batches": int(len(by)),
        "batch_size": int(batch_size),
        "note": "oracle: trained with OOD labels",
    }

    class_holdout = split.meta.get("protocol") == "class_holdout"
    if class_holdout:
        truth = y_bin.astype(float)
        rank_names = DEPLOY
    else:
        truth = score_po_resid(reg, X_eval, y_eval)
        rank_names = DEPLOY
    hard = {name: hard_rank_metrics(scores[name], truth) for name in rank_names}
    hard["rf_binary"] = hard_rank_metrics(
        rf_score,
        np.concatenate(
            [np.zeros(len(split.X_id_test), float), np.ones(len(X_ood_te), float)]
        )
        if class_holdout
        else score_po_resid(
            reg,
            X_rf_eval,
            np.concatenate([split.y_id_test, y_ood_te]),
        ),
    )

    return {
        "pack": pack,
        "id_domain": split.id_domain,
        "ood_domain": split.ood_domain,
        "meta": {
            **split.meta,
            "n_ood_rf_train": int(len(ood_tr_idx)),
            "n_ood_rf_eval": int(len(ood_te_idx)),
            "batch_size": int(batch_size),
        },
        "obs_metrics": obs_metrics,
        "batch_metrics": batch_metrics,
        "hard_rank": {
            k: {
                "spearman": v.get("spearman"),
                "precision_at_k": v.get("precision_at_k"),
                "lift_at_k": v.get("lift_at_k"),
                "ndcg_at_k": v.get("ndcg_at_k"),
                "auroc_topk": v.get("auroc_topk"),
            }
            for k, v in hard.items()
        },
    }


def _bar(
    all_res: dict,
    out_dir: Path,
    section: str,
    key: str,
    ylabel: str,
    title: str,
    fname: str,
    methods: Tuple[str, ...],
    ylim=None,
):
    packs = list(all_res.keys())
    fig, ax = plt.subplots(figsize=(max(9, 1.8 * len(packs)), 4.8))
    x = np.arange(len(packs))
    w = 0.12
    mid = (len(methods) - 1) / 2.0
    for i, m in enumerate(methods):
        vals = [all_res[p][section][m][key] for p in packs]
        ax.bar(x + (i - mid) * w, vals, w, label=m, color=SCORE_COLORS.get(m, "#888"))
    if key == "auroc":
        ax.axhline(0.5, color="k", ls="--", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(packs, rotation=15, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if ylim:
        ax.set_ylim(*ylim)
    ax.legend(ncol=3, fontsize=7)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    p = out_dir / fname
    fig.savefig(p, dpi=140)
    plt.close(fig)
    return p


def report(all_res: dict) -> str:
    def _fmt(v):
        if v is None or (isinstance(v, float) and (v != v)):
            return "—"
        return f"{v:.3f}"

    methods = DEPLOY + ORACLE
    lines = [
        "# Image-OOD: PO-risk is batch / hard-rank — not obs-level detection",
        "",
        "**Point:** observation-level PO-risk AUROC is too noisy for detection.",
        "If you want a binary OOD detector, use an **RF binary classifier** (oracle ceiling).",
        "PO-risk belongs to AGOD as a **batch / hard-sample** score after OnlineRFPerm.",
        "",
        "- `food101_vit`: frozen ViT + class-holdout (far-OOD) — primary.",
        "- CLIP packs: cached embeddings + domain shift.",
        "- `rf_binary`: oracle ID vs OOD RF (**uses OOD labels at train**) — ceiling, not a fair PO rival.",
        "",
        "## Primary: batch-mean AUROC (batch_size=32)",
        "",
        "| pack | po_msp | po_energy | maha | knn | centroid | rf_binary† | best deploy |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    wins = {m: 0 for m in DEPLOY}
    for pack, blob in all_res.items():
        m = blob["batch_metrics"]
        best = max(DEPLOY, key=lambda k: m[k]["auroc"] if m[k]["auroc"] == m[k]["auroc"] else -1)
        if m[best]["auroc"] == m[best]["auroc"]:
            wins[best] += 1
        lines.append(
            f"| `{pack}` | {_fmt(m['po_msp']['auroc'])} | {_fmt(m['po_energy']['auroc'])} | "
            f"{_fmt(m['maha']['auroc'])} | {_fmt(m['knn']['auroc'])} | {_fmt(m['centroid']['auroc'])} | "
            f"{_fmt(m['rf_binary']['auroc'])} | `{best}` |"
        )
    lines += [
        "",
        f"**Deployable batch-AUROC wins:** " + ", ".join(f"`{k}`={v}" for k, v in wins.items()),
        "",
        "† `rf_binary` = oracle ceiling (OOD labels at train).",
        "",
        "## Diagnostic: obs-level AUROC (noisy — do not prefer)",
        "",
        "| pack | po_msp | po_energy | maha | knn | centroid | rf_binary† |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for pack, blob in all_res.items():
        m = blob["obs_metrics"]
        lines.append(
            f"| `{pack}` | {_fmt(m['po_msp']['auroc'])} | {_fmt(m['po_energy']['auroc'])} | "
            f"{_fmt(m['maha']['auroc'])} | {_fmt(m['knn']['auroc'])} | {_fmt(m['centroid']['auroc'])} | "
            f"{_fmt(m['rf_binary']['auroc'])} |"
        )
    lines += [
        "",
        "## Hard-rank (PO job): spearman / P@20%",
        "",
        "| pack | po_msp | po_energy | maha | knn | centroid | rf_binary† |",
        "|---|---|---|---|---|---|---|",
    ]
    for pack, blob in all_res.items():
        h = blob["hard_rank"]
        cells = []
        for name in methods:
            if name not in h:
                cells.append("—")
                continue
            cells.append(f"{_fmt(h[name]['spearman'])}/{_fmt(h[name]['precision_at_k'])}")
        lines.append(f"| `{pack}` | " + " | ".join(cells) + " |")
    lines += [
        "",
        "### Takeaway",
        "",
        "1. Obs-level PO-risk AUROC is a bad primary — too noisy; RF-binary dominates if OOD labels exist.",
        "2. Batch-mean AUROC is closer to the OnlineRFPerm reject unit.",
        "3. Hard-rank is the AGOD-native PO question: does μ-risk surface the hard rows?",
        "",
        "See `docs/agod/AGOD_image_ood_bench.md`.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--packs", nargs="+", default=None)
    ap.add_argument("--pca-d", type=int, default=64)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--extract-max-n", type=int, default=8000)
    ap.add_argument("--batch-size", type=int, default=32)
    ap.add_argument("--out", type=Path, default=Path("results/agod_image_ood"))
    args = ap.parse_args()

    available = list_available_packs(args.root)
    packs = args.packs or available
    packs = [p for p in packs if p in available or p == "food101_vit"]
    if not packs:
        raise SystemExit(f"no packs under {args.root}; found {available}")

    args.out.mkdir(parents=True, exist_ok=True)
    all_res: Dict[str, dict] = {}
    for pack in packs:
        print(f"=== {pack} ===", flush=True)
        try:
            blob = run_pack(
                args.root,
                pack,
                pca_d=args.pca_d,
                seed=args.seed,
                extract_max_n=args.extract_max_n,
                batch_size=args.batch_size,
            )
        except Exception as e:
            print(f"[skip] {pack}: {e}", flush=True)
            continue
        all_res[pack] = blob
        b, o = blob["batch_metrics"], blob["obs_metrics"]
        print(
            f"  batch AUROC po_msp={b['po_msp']['auroc']:.3f} knn={b['knn']['auroc']:.3f} "
            f"rf_binary†={b['rf_binary']['auroc']:.3f} | "
            f"obs po_msp={o['po_msp']['auroc']:.3f} rf_binary†={o['rf_binary']['auroc']:.3f}",
            flush=True,
        )

    payload = {
        "pca_d": args.pca_d,
        "seed": args.seed,
        "batch_size": args.batch_size,
        "deploy_methods": list(DEPLOY),
        "oracle_methods": list(ORACLE),
        "packs": list(all_res.keys()),
        "results": all_res,
        "note": "Primary=batch AUROC + hard-rank; obs AUROC diagnostic; rf_binary=oracle ceiling",
    }
    (args.out / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md = report(all_res)
    (args.out / "IMAGE_OOD_REPORT.md").write_text(md, encoding="utf-8")
    Path("docs/agod").mkdir(parents=True, exist_ok=True)
    Path("docs/agod/AGOD_image_ood_bench.md").write_text(md, encoding="utf-8")
    plots = []
    if all_res:
        methods_plot = DEPLOY + ORACLE
        plots.append(
            _bar(
                all_res,
                args.out,
                "batch_metrics",
                "auroc",
                "batch AUROC",
                "Primary: batch-mean AUROC (PO vs classical vs oracle RF-binary)",
                "image_ood_batch_auroc.png",
                methods_plot,
                ylim=(0.35, 1.02),
            )
        )
        plots.append(
            _bar(
                all_res,
                args.out,
                "obs_metrics",
                "auroc",
                "obs AUROC (noisy)",
                "Diagnostic: obs-level AUROC (prefer batch / hard-rank)",
                "image_ood_obs_auroc.png",
                methods_plot,
                ylim=(0.35, 1.02),
            )
        )
    print(md)
    print("plots:", [str(p) for p in plots])


if __name__ == "__main__":
    main()
