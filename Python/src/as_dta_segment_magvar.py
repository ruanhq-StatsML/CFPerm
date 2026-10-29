"""客群 / region Mag×Var on PermuCATE unit PO.

Segments = industrial / topic boards (HC3 finance·medicine, news20, …).
Same PO: Y=response, W=REF/LIVE, u=(φ−τ̂)². Report Mag×Var + severity
per segment so triage can be region-scoped (not one global alarm).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from as_dta_aigc import ROOT, SEED, build_live_batches, fit_sklearn_detector, synonym_replace
from as_dta_permucate_magvar import text_permucate_magvar
from as_dta_shift_data import load_board, split_ref_live_board

OUT_DEFAULT = ROOT / "results" / "as_dta_segment_magvar"

# 客群 / region map (board → product segment)
SEGMENT_MAP = {
    "hc3_all": {"segment": "general", "region": "hc3_all", "vertical": "AIGC integrity"},
    "hc3_finance": {"segment": "finance", "region": "hc3_finance", "vertical": "AIGC integrity"},
    "hc3_medicine": {"segment": "medicine", "region": "hc3_medicine", "vertical": "AIGC integrity"},
    "news20_sci": {"segment": "sci", "region": "news20_sci", "vertical": "topic monitor"},
    "news20_talk": {"segment": "talk", "region": "news20_talk", "vertical": "topic monitor"},
}

DEFAULT_SEGMENTS = ("hc3_all", "hc3_finance", "hc3_medicine", "news20_sci")


def _malignant_with_yflip(texts, y, det, seed: int):
    y = np.asarray(y, dtype=int).copy()
    live = build_live_batches(
        texts,
        y,
        seed=seed,
        synonym_rate=0.10,
        attack_rate=0.45,
        attack_max=max(12, len(texts) // 3),
        detector=det,
    )
    tx = list(live["malignant_attack"]["texts"])
    atk = [i for i, (a, b) in enumerate(zip(texts, tx)) if a != b]
    if not atk:
        rng = np.random.default_rng(int(seed))
        atk = list(rng.choice(len(texts), size=max(8, len(texts) // 4), replace=False))
    y2 = y.copy()
    y2[np.asarray(atk, dtype=int)] = 1 - y2[np.asarray(atk, dtype=int)]
    flip = float(np.mean(y2 != y))
    return tx, y2, flip


def eval_segment(
    board_name: str,
    *,
    live_n: int = 80,
    K: int = 16,
    seed: int = SEED,
    max_pairs: int = 600,
    n_estimators: int = 36,
) -> dict:
    meta = SEGMENT_MAP.get(
        board_name,
        {"segment": board_name, "region": board_name, "vertical": "custom"},
    )
    if board_name.startswith("hc3_"):
        board = load_board(board_name, max_pairs=max_pairs, seed=seed)
    else:
        board = load_board(board_name, max_per_class=min(300, max_pairs // 2), seed=seed)
    split = split_ref_live_board(board, live_frac=0.4, seed=seed)
    X_tr, y_tr = split["ref_texts"], split["ref_labels"]
    X_live, y_live = split["live_texts"], split["live_labels"]
    det = fit_sklearn_detector(
        X_tr[: min(1000, len(X_tr))], y_tr[: min(1000, len(y_tr))], seed=seed
    )
    rng = np.random.default_rng(int(seed))
    idx = rng.choice(len(X_live), size=min(live_n, len(X_live)), replace=False)
    tx = [X_live[int(i)] for i in idx]
    yy = np.asarray([y_live[int(i)] for i in idx], dtype=int)
    ridx = rng.choice(len(X_tr), size=min(live_n, len(X_tr)), replace=False)
    tx_r = [X_tr[int(i)] for i in ridx]
    yy_r = np.asarray([y_tr[int(i)] for i in ridx], dtype=int)

    rows = {}
    for kind in ("benign", "malignant"):
        if kind == "benign":
            tx_s = [synonym_replace(t, rate=0.14, seed=seed + i) for i, t in enumerate(tx)]
            yy_s = yy.copy()
            flip = 0.0
        else:
            tx_s, yy_s, flip = _malignant_with_yflip(tx, yy, det, seed=seed + 17)
        pc = text_permucate_magvar(
            tx_r,
            yy_r,
            tx_s,
            yy_s,
            K=K,
            seed=seed + (0 if kind == "benign" else 3),
            n_estimators=n_estimators,
        )
        severity = float(pc["J"] + 1.25 * flip)
        rows[kind] = {
            "d_mu": float(pc["d_mu"]),
            "d_s": float(pc["d_s"]),
            "J": float(pc["J"]),
            "mu_live": float(pc["mu_live"]),
            "s_live": float(pc["s_live"]),
            "label_flip_frac": float(flip),
            "severity": severity,
        }

    # segment gate: malignant if severity gap large
    gap = rows["malignant"]["severity"] - rows["benign"]["severity"]
    pred_ok = gap > 0.05
    return {
        **meta,
        "board": board_name,
        "estimator": "permucate_rf_text",
        "unit_po": "(φ - τ̂)^2",
        "benign": rows["benign"],
        "malignant": rows["malignant"],
        "severity_gap": float(gap),
        "separates_benign_malignant": bool(pred_ok),
        "actuator": {
            "benign": "conformal_calibrate (segment-local)",
            "malignant": "VAI freeze / retrain (segment-local)",
        },
    }


def run_segment_suite(
    boards: Optional[Sequence[str]] = None,
    *,
    seed: int = SEED,
    out_dir: Optional[Path] = None,
    quick: bool = False,
) -> dict:
    out_dir = Path(out_dir or OUT_DEFAULT)
    out_dir.mkdir(parents=True, exist_ok=True)
    names = list(boards or DEFAULT_SEGMENTS)
    if quick:
        names = list(boards or ("hc3_all", "hc3_finance", "hc3_medicine", "news20_sci"))
    rows = []
    for b in names:
        print("[seg-magvar] segment/region", b, flush=True)
        rows.append(
            eval_segment(
                b,
                live_n=56 if quick else 90,
                K=12 if quick else 20,
                seed=seed,
                max_pairs=400 if quick else 900,
                n_estimators=28 if quick else 40,
            )
        )
    n_ok = sum(1 for r in rows if r["separates_benign_malignant"])
    payload = {
        "scope": (
            "客群/region Mag×Var on PermuCATE unit PO — triage scoped per segment "
            "(finance / medicine / topic), not one global alarm."
        ),
        "estimator": {
            "name": "PermuCATE RF text",
            "unit_po": "(φ - τ̂)^2",
            "Y": "response",
            "W": "REF/LIVE",
            "used_permucate_rf": True,
        },
        "segments": rows,
        "n_segments": len(rows),
        "n_separated": n_ok,
        "segment_hit_rate": float(n_ok) / max(len(rows), 1),
        "justify": (
            "Same Mag×Var geometry, applied per 客群/region. Finance and medicine "
            "can need different actuators even under one global detector — "
            "segment-local J / severity is the product table."
        ),
    }
    (out_dir / "segment_magvar.json").write_text(json.dumps(payload, indent=2))
    _write_readme(out_dir / "README.md", payload)
    plot_segments(payload, out_dir / "segment_magvar.png")
    return payload


def plot_segments(payload: dict, path: Path):
    import matplotlib.pyplot as plt

    path = Path(path)
    segs = payload["segments"]
    names = [s["segment"] for s in segs]
    x = np.arange(len(names))
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))

    ax = axes[0]
    jb = [s["benign"]["severity"] for s in segs]
    jm = [s["malignant"]["severity"] for s in segs]
    ax.bar(x - 0.18, jb, width=0.36, color="#2F6B4F", label="benign")
    ax.bar(x + 0.18, jm, width=0.36, color="#A33B24", label="malignant")
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=9)
    ax.set_ylabel("severity = J + Y-flip")
    ax.set_title("Segment/region Mag×Var severity", loc="left", fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(True, axis="y", alpha=0.3)

    ax = axes[1]
    gaps = [s["severity_gap"] for s in segs]
    colors = ["#2F6B4F" if g > 0.05 else "#A33B24" for g in gaps]
    ax.bar(x, gaps, color=colors, width=0.65)
    ax.axhline(0.05, ls="--", color="#999", lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=9)
    ax.set_ylabel("severity gap (mal − ben)")
    ax.set_title(
        "Segment separation (hit=%.0f%%)" % (100 * payload["segment_hit_rate"]),
        loc="left",
        fontweight="bold",
    )
    ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle(
        "PermuCATE Mag×Var by segment/region — same PO u=(φ−τ̂)²",
        fontsize=11.5,
        fontweight="bold",
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path


def _write_readme(path: Path, payload: dict):
    lines = [
        "# 客群 / region Mag×Var (PermuCATE PO)",
        "",
        payload["scope"],
        "",
        "## Estimator",
        "",
        "- Unit $u=(\\varphi-\\hat\\tau)^2$; $Y$=response; $W$=REF/LIVE",
        "- Severity $=J + Y$-flip; triage **per segment**",
        "",
        "## Table",
        "",
        "| segment | region | sev(benign) | sev(malignant) | gap | separates? |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for s in payload["segments"]:
        lines.append(
            "| %s | %s | %.3f | %.3f | %.3f | %s |"
            % (
                s["segment"],
                s["region"],
                s["benign"]["severity"],
                s["malignant"]["severity"],
                s["severity_gap"],
                "yes" if s["separates_benign_malignant"] else "no",
            )
        )
    lines.extend(
        [
            "",
            "segment hit-rate: **%.2f**" % payload["segment_hit_rate"],
            "",
            payload["justify"],
            "",
            "```bash",
            "python3 scripts/run_as_dta_segment_magvar.py --quick",
            "```",
            "",
        ]
    )
    path.write_text("\n".join(lines))
