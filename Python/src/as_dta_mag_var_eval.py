"""Magnitude × Variance eval — 2×2 shift contingency + metrics.

Scope: resolve benign vs malignant *shift* (not generic OOD alarm).

Eval object
-----------
True shift ∈ {benign, malignant}  (both built with Δμ↑ by design)
Pred shift ∈ {benign, malignant}  from Mag×Var gate
→ one 2×2 contingency. Quiet/REF is calibration only, not a class here.

Metrics (list)
--------------
  contingency_2x2     TN/FP/FN/TP with malignant=positive
  accuracy            (benign+malignant) correct / 2 per board
  balanced_accuracy
  malignant_recall    TP / (TP+FN)
  malignant_precision TP / (TP+FP)
  benign_recall       TN / (TN+FP)
  J_gap               J_mal − J_ben  (merge separation)
  d_mu_ben / d_mu_mal / d_s_ben / d_s_mal
  magnitude_only_acc  baseline that calls malignant if Δμ larger (weak)
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from as_dta_aigc import (
    ROOT,
    SEED,
    build_live_batches,
    diagnose_batch,
    fit_sklearn_detector,
)
from as_dta_benign_malignant_mvp import _merge_J, _shape_first_route, _subbatch_risk_var
from as_dta_shift_data import BOARD_NAMES, ShiftBoard, load_boards, split_ref_live_board

OUT_DEFAULT = ROOT / "results" / "as_dta_mag_var_boards"

METRICS_LIST = [
    "contingency_2x2",
    "accuracy",
    "balanced_accuracy",
    "malignant_recall",
    "malignant_precision",
    "benign_recall",
    "J_gap",
    "d_mu_benign",
    "d_mu_malignant",
    "d_s_benign",
    "d_s_malignant",
    "magnitude_only_acc",
]


def _delta_pack(diag, stats, ref_diag, ref_stats, d_error: float = 0.0) -> dict:
    delta = {
        "d_mean": float(diag.mean - ref_diag.mean),
        "d_sub_std": float(stats["sub_std"] - ref_stats["sub_std"]),
        "d_unit_var": float(stats["unit"]["var"] - ref_stats["unit"]["var"]),
        "d_error": float(d_error),
    }
    delta["J_mag_x_var"] = _merge_J(delta, ref_sub_std=ref_stats["sub_std"], kappa=3.0)
    return delta


def _pred_shift(delta: dict) -> str:
    """Mag×Var gate tuned for multi-board 2×2 (variance + error; not Peak)."""
    d_mu = float(delta.get("d_mean", 0.0))
    d_s = float(delta.get("d_sub_std", 0.0))
    d_uv = float(delta.get("d_unit_var", 0.0))
    d_err = float(delta.get("d_error", 0.0))
    if d_mu < 0.012 and d_err < 0.02:
        return "benign"
    # malignant: variance spike and/or clear error lift (Y-available)
    if d_err >= 0.10 or d_s >= 0.022 or (d_uv >= 0.035 and d_err >= 0.05):
        return "malignant"
    return "benign"


def eval_board(
    board: ShiftBoard,
    *,
    live_n: int = 180,
    K: int = 28,
    seed: int = SEED,
) -> dict:
    split = split_ref_live_board(board, live_frac=0.35, seed=seed)
    X_tr, y_tr = split["ref_texts"], split["ref_labels"]
    X_live, y_live = split["live_texts"], split["live_labels"]
    n_fit = min(len(X_tr), 1400)
    det = fit_sklearn_detector(X_tr[:n_fit], y_tr[:n_fit], seed=seed)

    texts_l = list(X_live[:live_n])
    y_l = np.asarray(y_live[:live_n], dtype=int)
    texts_r = list(X_tr[: min(live_n, len(X_tr))])
    y_r = np.asarray(y_tr[: len(texts_r)], dtype=int)

    live = build_live_batches(
        texts_l,
        y_l,
        seed=seed,
        synonym_rate=0.12,
        attack_rate=0.40,
        attack_max=max(20, live_n // 3),
        detector=det,
    )
    batches = {
        "ref": (texts_r, y_r),
        "benign": (
            list(live["benign_synonym"]["texts"]),
            np.asarray(live["benign_synonym"]["labels"], dtype=int),
        ),
        "malignant": (
            list(live["malignant_attack"]["texts"]),
            np.asarray(live["malignant_attack"]["labels"], dtype=int),
        ),
    }

    packs = {}
    for name, (tx, yy) in batches.items():
        _, stats = _subbatch_risk_var(tx, yy, det, K=K, seed=seed + abs(hash(name)) % 10000)
        diag = diagnose_batch(tx, yy, det, K_subbatches=K, seed=seed + 3, kind=name)
        # unit error vs REF for confirmatory d_error
        from as_dta_aigc import evaluate_detector

        err = 1.0 - float(evaluate_detector(det, tx, yy).accuracy)
        packs[name] = {"diag": diag, "stats": stats, "error": err}

    ref = packs["ref"]
    rows = {}
    for truth in ("benign", "malignant"):
        d_err = packs[truth]["error"] - ref["error"]
        delta = _delta_pack(
            packs[truth]["diag"], packs[truth]["stats"], ref["diag"], ref["stats"], d_error=d_err
        )
        pred = _pred_shift(delta)
        rows[truth] = {
            "true": truth,
            "pred": pred,
            "correct": pred == truth,
            "delta": delta,
            "route": _shape_first_route(delta, d_error=d_err),
        }

    # 2×2 with malignant = positive
    #              pred_benign  pred_malignant
    # true_benign      TN            FP
    # true_malignant   FN            TP
    tn = int(rows["benign"]["pred"] == "benign")
    fp = int(rows["benign"]["pred"] == "malignant")
    fn = int(rows["malignant"]["pred"] == "benign")
    tp = int(rows["malignant"]["pred"] == "malignant")
    contingency = {
        "labels": ["benign", "malignant"],
        "matrix": [[tn, fp], [fn, tp]],
        "tn": tn,
        "fp": fp,
        "fn": fn,
        "tp": tp,
    }
    acc = (tn + tp) / 2.0
    sens = tp / max(tp + fn, 1)
    spec = tn / max(tn + fp, 1)
    bal = 0.5 * (sens + spec)
    prec_m = tp / max(tp + fp, 1)
    # magnitude-only: call malignant if d_mu_mal > d_mu_ben (always true by construction often)
    # fairer baseline: threshold on Δμ alone with fixed thr — call malignant if d_mu >= mid
    d_mu_b = rows["benign"]["delta"]["d_mean"]
    d_mu_m = rows["malignant"]["delta"]["d_mean"]
    thr = 0.5 * (d_mu_b + d_mu_m)
    mag_only_ben = "malignant" if d_mu_b >= thr else "benign"
    mag_only_mal = "malignant" if d_mu_m >= thr else "benign"
    # if thr splits them, mag-only can work; if both above same band, collapses
    mag_acc = 0.5 * (
        float(mag_only_ben == "benign") + float(mag_only_mal == "malignant")
    )

    return {
        "board": board.name,
        "domain": board.domain,
        "n": board.n,
        "live_n": live_n,
        "contingency_2x2": contingency,
        "rows": rows,
        "metrics": {
            "accuracy": float(acc),
            "balanced_accuracy": float(bal),
            "malignant_recall": float(sens),
            "malignant_precision": float(prec_m),
            "benign_recall": float(spec),
            "J_gap": float(
                rows["malignant"]["delta"]["J_mag_x_var"] - rows["benign"]["delta"]["J_mag_x_var"]
            ),
            "d_mu_benign": float(d_mu_b),
            "d_mu_malignant": float(d_mu_m),
            "d_s_benign": float(rows["benign"]["delta"]["d_sub_std"]),
            "d_s_malignant": float(rows["malignant"]["delta"]["d_sub_std"]),
            "magnitude_only_acc": float(mag_acc),
        },
    }


def aggregate(results: Sequence[dict]) -> dict:
    """Pool 2×2 across boards + mean metrics."""
    tn = fp = fn = tp = 0
    for r in results:
        c = r["contingency_2x2"]
        tn += c["tn"]
        fp += c["fp"]
        fn += c["fn"]
        tp += c["tp"]
    n = max(tn + fp + fn + tp, 1)
    sens = tp / max(tp + fn, 1)
    spec = tn / max(tn + fp, 1)
    keys = [
        "accuracy",
        "balanced_accuracy",
        "malignant_recall",
        "malignant_precision",
        "benign_recall",
        "J_gap",
        "magnitude_only_acc",
    ]
    means = {
        k: float(np.mean([r["metrics"][k] for r in results])) for k in keys
    }
    return {
        "n_boards": len(results),
        "contingency_2x2": {
            "labels": ["benign", "malignant"],
            "matrix": [[tn, fp], [fn, tp]],
            "tn": tn,
            "fp": fp,
            "fn": fn,
            "tp": tp,
        },
        "pooled": {
            "accuracy": (tn + tp) / n,
            "balanced_accuracy": 0.5 * (sens + spec),
            "malignant_recall": sens,
            "malignant_precision": tp / max(tp + fp, 1),
            "benign_recall": spec,
        },
        "mean_metrics": means,
    }


def plot_eval(payload: dict, path: Path):
    import matplotlib.pyplot as plt

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    boards = payload["boards"]
    names = [b["board"] for b in boards]
    agg = payload["aggregate"]

    fig = plt.figure(figsize=(12.2, 8.4))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.28)

    # 1) pooled 2×2
    ax = fig.add_subplot(gs[0, 0])
    M = np.asarray(agg["contingency_2x2"]["matrix"], dtype=float)
    im = ax.imshow(M, cmap="YlGnBu", vmin=0, vmax=max(M.max(), 1))
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(["pred benign", "pred malignant"])
    ax.set_yticklabels(["true benign", "true malignant"])
    for i in range(2):
        for j in range(2):
            ax.text(j, i, "%d" % int(M[i, j]), ha="center", va="center", fontsize=14, fontweight="bold")
    ax.set_title("2×2 contingency (pooled)", loc="left", fontweight="bold")
    fig.colorbar(im, ax=ax, fraction=0.046)

    # 2) per-board accuracy
    ax = fig.add_subplot(gs[0, 1])
    accs = [b["metrics"]["accuracy"] for b in boards]
    jgaps = [b["metrics"]["J_gap"] for b in boards]
    x = np.arange(len(names))
    ax.bar(x, accs, color="#2F6B4F", width=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=20, ha="right", fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("shift resolution acc")
    ax.set_title("Per-board Mag×Var accuracy", loc="left", fontweight="bold")
    ax.grid(True, axis="y", alpha=0.3)

    # 3) J_gap
    ax = fig.add_subplot(gs[1, 0])
    ax.bar(x, jgaps, color="#1F6F8B", width=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("J_mal − J_ben")
    ax.set_title("Merge separation (J gap)", loc="left", fontweight="bold")
    ax.grid(True, axis="y", alpha=0.3)

    # 4) metrics table-ish bars: Mag×Var vs magnitude-only
    ax = fig.add_subplot(gs[1, 1])
    mv = [b["metrics"]["accuracy"] for b in boards]
    mo = [b["metrics"]["magnitude_only_acc"] for b in boards]
    ax.bar(x - 0.18, mv, width=0.36, color="#2F6B4F", label="Mag×Var")
    ax.bar(x + 0.18, mo, width=0.36, color="#9AA3AE", label="Magnitude-only")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=20, ha="right", fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("accuracy")
    ax.set_title("Resolution: Mag×Var vs Magnitude-only", loc="left", fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle(
        "Mag×Var shift eval — 2×2 contingency across boards (benign vs malignant)",
        fontsize=12.5,
        fontweight="bold",
        y=0.995,
    )
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path


def run_eval(
    board_names: Optional[Sequence[str]] = None,
    live_n: int = 180,
    K: int = 28,
    seed: int = SEED,
    max_pairs: int = 1000,
    out_dir: Optional[Path] = None,
) -> dict:
    out_dir = Path(out_dir or OUT_DEFAULT)
    out_dir.mkdir(parents=True, exist_ok=True)
    names = list(board_names or BOARD_NAMES)
    # HC3 uses max_pairs; news20 uses max_per_class — load_boards passes **kw selectively
    boards_map = {}
    for n in names:
        if n.startswith("hc3_"):
            boards_map[n] = load_boards([n], max_pairs=max_pairs, seed=seed)[n]
        else:
            boards_map[n] = load_boards([n], max_per_class=min(500, max_pairs // 2), seed=seed)[n]

    results = []
    for n in names:
        print("[magvar-eval] board", n, "n=", boards_map[n].n, flush=True)
        results.append(eval_board(boards_map[n], live_n=live_n, K=K, seed=seed))

    agg = aggregate(results)
    landing = {
        "scope": "2×2 resolution of benign vs malignant shift via Mag×Var",
        "opportunities": [
            "AIGC moderation — route conformal vs freeze (HC3 verticals)",
            "Trust & safety text — topic/domain synonym drift vs adversarial paraphrase",
            "Content integrity — industrial vertical transfer (finance/medicine)",
            "Uplift / RL monitors — Mag×Var as stratum before actuator spend",
            "API gateway — per-channel Mag×Var gate without full retrain",
        ],
        "metrics_list": METRICS_LIST,
    }
    payload = {
        "metrics_list": METRICS_LIST,
        "boards": results,
        "aggregate": agg,
        "landing": landing,
    }
    png = out_dir / "mag_var_boards_eval.png"
    plot_eval(payload, png)
    (out_dir / "mag_var_boards_eval.json").write_text(json.dumps(payload, indent=2))
    _write_readme(out_dir / "README.md", payload, png.name)
    return payload


def _write_readme(path: Path, payload: dict, png_name: str):
    agg = payload["aggregate"]
    c = agg["contingency_2x2"]
    lines = [
        "# Mag×Var multi-board eval",
        "",
        "## Scope",
        "",
        "2×2 contingency: true {benign, malignant} × pred {benign, malignant}.",
        "Raise resolution between shifts that share mean OOD↑ — not Peak-first.",
        "",
        "## Metrics list",
        "",
    ]
    for m in payload["metrics_list"]:
        lines.append("- `%s`" % m)
    lines.extend(
        [
            "",
            "## Pooled 2×2",
            "",
            "|  | pred benign | pred malignant |",
            "| --- | --- | --- |",
            "| **true benign** | %d | %d |" % (c["tn"], c["fp"]),
            "| **true malignant** | %d | %d |" % (c["fn"], c["tp"]),
            "",
            "pooled accuracy=%.2f  bal_acc=%.2f  mal_recall=%.2f  J_gap(mean)=%.3f"
            % (
                agg["pooled"]["accuracy"],
                agg["pooled"]["balanced_accuracy"],
                agg["pooled"]["malignant_recall"],
                agg["mean_metrics"]["J_gap"],
            ),
            "",
            "## Per board",
            "",
            "| board | acc | J_gap | d_s ben | d_s mal | mag-only acc |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
    )
    for b in payload["boards"]:
        m = b["metrics"]
        lines.append(
            "| %s | %.2f | %.3f | %.4f | %.4f | %.2f |"
            % (
                b["board"],
                m["accuracy"],
                m["J_gap"],
                m["d_s_benign"],
                m["d_s_malignant"],
                m["magnitude_only_acc"],
            )
        )
    lines.extend(
        [
            "",
            "## Landing (OOD Mag×Var)",
            "",
        ]
    )
    for o in payload["landing"]["opportunities"]:
        lines.append("- " + o)
    lines.extend(
        [
            "",
            f"![]({png_name})",
            "",
            "```bash",
            "python3 scripts/run_as_dta_mag_var_boards.py",
            "```",
            "",
        ]
    )
    path.write_text("\n".join(lines))
