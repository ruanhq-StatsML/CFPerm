#!/usr/bin/env python3
"""Order-grain FSDS + two-layer subset attribution (MMD, PO-risk, Conditional Mean)."""
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

from graph_fsds_localize import run_pipeline  # noqa: E402
from stream_dgps import make_order_graph_stream  # noqa: E402

OUT = ROOT / "results" / "graph_fsds_localize"
KINDS = ("covariate_south", "concept_south", "both")
N_REF = 360
N_NEW = 100
N_BATCHES = 3
ONSET = 1


def _fmt(x, n=3):
    if x is None:
        return ""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if v != v:
        return ""
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


def slim_rank(rank, k=4):
    rows = []
    for r in rank[:k]:
        rows.append(
            {
                "feature": r["feature"],
                "score": r["score"],
                "mmd": r["mmd"],
                "cmean_x": r["cmean_x"],
                "cmean_y": r["cmean_y"],
                "loud": r["loud"],
            }
        )
    return rows


def run_kind(kind: str, seed: int = 2026) -> dict:
    tables = make_order_graph_stream(
        n_ref=N_REF,
        n_new=N_NEW,
        n_batches=N_BATCHES,
        onset_batch=ONSET,
        kind=kind,
        seed=seed,
    )
    rows = []
    for t in range(N_BATCHES):
        use_logo = t == N_BATCHES - 1
        rec = run_pipeline(
            tables,
            t=t,
            grain="order",
            mode="localize",
            top_k=4,
            seed=seed,
            min_n=20,
            with_logo=use_logo,
            with_po=True,
        )
        rec["onset"] = t >= ONSET
        rows.append(rec)
    return {"meta": tables["meta"], "rows": rows}


def portrait_map(portraits):
    return {p["subset"]: p for p in portraits}


def write_markdown(results: dict, dest: Path) -> None:
    lines = [
        "# Graph localization → FSDS unify → two-layer subset (order grain)",
        "",
        "Shares are localization proxies, not a unique decomposition. "
        "T is the batch label. Y is never a feature and never the subset key.",
        "",
        f"n_ref={N_REF}, n_new={N_NEW}, batches={N_BATCHES}, onset={ONSET}. "
        "Planted subset is always **south**.",
        "",
        "## Per-batch order-grain portraits (MMD + PO-risk + Conditional Mean)",
        "",
    ]
    header = [
        "kind",
        "t",
        "onset",
        "FSDS selected",
        "loud subset",
        "south π_MMD",
        "south π_PO",
        "south π_CMean",
        "south mix",
        "south MMD",
        "south PO",
        "south ΔE[Y]",
        "south ‖ΔE[X]‖",
        "gap MMD vs other",
        "fingerprint",
        "LOGO π_MMD o/m/u",
    ]
    lines.append(_md_row(header))
    lines.append(_md_row(["---"] * len(header)))
    for kind in KINDS:
        for rec in results[kind]["rows"]:
            pm = portrait_map(rec["portraits"])
            s = pm.get("south") or {}
            logo = rec.get("logo_full") or {}
            pi = logo.get("pi_mmd") or {}
            logo_trip = "/".join(_fmt(pi.get(g, 0)) for g in ("order", "merchant", "user"))
            mix = "/".join(
                _fmt(s.get(k))
                for k in ("mix_mmd", "mix_po", "mix_cmean")
            )
            lines.append(
                _md_row(
                    [
                        kind,
                        rec["t"],
                        "yes" if rec["onset"] else "",
                        ",".join(rec["fsds"]["selected_names"]),
                        rec["loud_subset"],
                        _fmt(s.get("pi_mmd")),
                        _fmt(s.get("pi_po")),
                        _fmt(s.get("pi_cmean")),
                        mix,
                        _fmt(s.get("mmd")),
                        _fmt(s.get("po")),
                        _fmt(s.get("cmean_y")),
                        _fmt(s.get("cmean_x")),
                        _fmt((s.get("gap_vs_other") or {}).get("mmd")),
                        s.get("fingerprint") or "",
                        logo_trip,
                    ]
                )
            )
    lines += [
        "",
        "## FSDS ranking (last batch, top of each grain)",
        "",
        _md_row(["kind", "grain", "feature", "score", "mmd", "cmean_x", "cmean_y", "loud"]),
        _md_row(["---"] * 8),
    ]
    for kind in KINDS:
        last = results[kind]["rows"][-1]
        for grain, rank in last["fsds"]["rank"].items():
            for r in slim_rank(rank, k=3):
                lines.append(
                    _md_row(
                        [
                            kind,
                            grain,
                            r["feature"],
                            _fmt(r["score"]),
                            _fmt(r["mmd"]),
                            _fmt(r["cmean_x"]),
                            _fmt(r["cmean_y"]),
                            "yes" if r["loud"] else "",
                        ]
                    )
                )
    dest.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot_portraits(results: dict, dest: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), sharey=True)
    metrics = [("pi_mmd", "π_MMD"), ("pi_po", "π_PO"), ("pi_cmean", "π_CMean")]
    colors = {"south": "#c0392b", "north": "#2980b9"}
    for ax, kind in zip(axes, KINDS):
        rec = results[kind]["rows"][-1]
        pm = portrait_map(rec["portraits"])
        x = np.arange(len(metrics))
        w = 0.35
        for i, lab in enumerate(("south", "north")):
            p = pm.get(lab) or {}
            vals = [float(p.get(k) or 0.0) for k, _ in metrics]
            ax.bar(x + (i - 0.5) * w, vals, w, label=lab, color=colors[lab])
        ax.set_xticks(x)
        ax.set_xticklabels([t for _, t in metrics])
        ax.set_title(kind)
        ax.set_ylim(0, 1.05)
        ax.grid(axis="y", alpha=0.3)
    axes[0].set_ylabel("subset share")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Order grain · last batch · subset shares (localization, not Shapley)")
    fig.tight_layout()
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    results = {kind: run_kind(kind) for kind in KINDS}
    plot_portraits(results, OUT / "subset_shares.png")
    write_markdown(results, OUT / "TABLES.md")
    slim = {}
    for kind, r in results.items():
        slim[kind] = {
            "meta": r["meta"],
            "rows": [
                {
                    "t": rec["t"],
                    "onset": rec["onset"],
                    "selected": rec["fsds"]["selected_names"],
                    "loud_subset": rec["loud_subset"],
                    "portraits": rec["portraits"],
                    "logo_full": rec.get("logo_full"),
                    "logo_subset": rec.get("logo_subset"),
                    "leakage": rec["leakage"],
                    "fsds_rank": {
                        g: slim_rank(rank, k=4)
                        for g, rank in rec["fsds"]["rank"].items()
                    },
                }
                for rec in r["rows"]
            ],
        }
    (OUT / "summary.json").write_text(json.dumps(jsonable(slim), indent=2), encoding="utf-8")
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
