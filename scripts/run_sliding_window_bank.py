#!/usr/bin/env python3
"""Put the sliding-window bank on a gradual concept vs covariate walk."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
for p in (SRC, ROOT):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from sliding_window_bank import SlidingWindowBank  # noqa: E402
from stream_dgps import (  # noqa: E402
    TRIMODAL_GROUPS,
    make_trimodal_gradual_concept,
    make_trimodal_stream,
)

OUT = ROOT / "results" / "sliding_window_bank"


def _fmt(x, n=3):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    return f"{v:.{n}g}"


def _md_row(cells) -> str:
    return "| " + " | ".join(str(c) for c in cells) + " |"


def jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return float(obj) if isinstance(obj, np.floating) else int(obj)
    return obj


def run_concept():
    X, Y, alpha, _, meta = make_trimodal_gradual_concept(
        n_ref=400, n_new=40, n_batches=12, onset_batch=3, seed=2026
    )
    return X, Y, meta, "concept"


def run_covariate():
    X, Y, sl, meta = make_trimodal_stream(
        n_ref=400,
        n_new=40,
        n_batches=12,
        onset_batch=3,
        kind="covariate_audio",
        seed=2026,
    )
    return X, Y, meta, "covariate"


def walk(kind: str):
    if kind == "concept":
        X, Y, meta, name = run_concept()
    else:
        X, Y, meta, name = run_covariate()
    n_ref, n_new = meta["n_ref"], meta["n_new"]
    n_batches = int(meta.get("n_batches") or (len(Y) - n_ref) // n_new)
    bank = SlidingWindowBank(
        X[:n_ref],
        Y[:n_ref],
        window=120,
        n_components=4,
        groups=TRIMODAL_GROUPS,
        with_po=False,
    )
    rows = []
    for t in range(n_batches):
        lo = n_ref + t * n_new
        feat = bank.step(X[lo : lo + n_new], Y[lo : lo + n_new])
        feat["t"] = t
        feat["onset"] = t >= int(meta["onset_batch"])
        feat["vector"] = bank.vector(feat).tolist()
        rows.append(feat)
    return {"kind": name, "meta": meta, "rows": rows}


def plot_walks(results: dict, dest: Path) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.6), sharey=False)
    for ax, kind in zip(axes, ("concept", "covariate")):
        rows = results[kind]["rows"]
        t = [r["t"] for r in rows]
        ax.plot(t, [r["rfperm_T"] for r in rows], label="T", color="#1f4e79")
        ax.plot(t, [r["brier_excess"] for r in rows], label="Brier excess", color="#555")
        ax.plot(t, [r["pca_recon_excess"] for r in rows], label="PCA recon vs ref", color="#d48b16")
        g = [r.get("pca_recon_group") or {} for r in rows]
        if g and g[-1]:
            ax.plot(t, [x.get("audio", 0) for x in g], label="PCA audio", color="#6b4ea0", ls="--")
            ax.plot(t, [x.get("video", 0) for x in g], label="PCA video", color="#b33", ls=":")
        ax.axvline(2.5, color="#999", ls=":", lw=1)
        ax.set_title(kind)
        ax.set_xlabel("batch")
        ax.grid(alpha=0.3)
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].set_ylabel("feature")
    fig.suptitle("sliding-window bank · same D_ref", fontsize=11)
    fig.tight_layout()
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    results = {k: walk(k) for k in ("concept", "covariate")}
    plot_walks(results, OUT / "bank_walk.png")
    lines = [
        "# Sliding-window candidate feature bank",
        "",
        "Cheap rolling stats vs frozen D_ref. Online PCA is the fast P(X) pointer.",
        "MMD still vs D_ref, not pairwise history.",
        "",
        _md_row(["kind", "t", "onset", "T", "Brier excess", "PCA recon", "PCA video", "PCA audio", "PCA text", "argmax"]),
        _md_row(["---"] * 10),
    ]
    for kind in ("concept", "covariate"):
        for r in results[kind]["rows"]:
            g = r.get("pca_recon_group") or {}
            lines.append(
                _md_row(
                    [
                        kind,
                        r["t"],
                        "yes" if r["onset"] else "",
                        _fmt(r["rfperm_T"]),
                        _fmt(r["brier_excess"]),
                        _fmt(r["pca_recon_excess"]),
                        _fmt(g.get("video")),
                        _fmt(g.get("audio")),
                        _fmt(g.get("text")),
                        r.get("pca_recon_group_argmax") or "",
                    ]
                )
            )
    (OUT / "TABLES.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    slim = {
        k: {
            "meta": v["meta"],
            "rows": [
                {
                    key: rec[key]
                    for key in (
                        "t",
                        "onset",
                        "n_window",
                        "mean_l2",
                        "pca_recon_excess",
                        "pca_score_l2",
                        "pca_subspace_gap",
                        "mmd_vs_ref",
                        "brier_excess",
                        "rfperm_T",
                        "pca_recon_group",
                        "pca_recon_group_argmax",
                    )
                    if key in rec
                }
                for rec in v["rows"]
            ],
        }
        for k, v in results.items()
    }
    (OUT / "summary.json").write_text(json.dumps(jsonable(slim), indent=2), encoding="utf-8")
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
