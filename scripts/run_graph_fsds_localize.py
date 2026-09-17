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
N_REF = 480
N_NEW = 160
N_BATCHES = 3
ONSET = 1
N_MERCHANTS = 16
N_USERS = 48


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
        n_merchants=N_MERCHANTS,
        n_users=N_USERS,
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
            subset_by="level_set",
            top_k=4,
            seed=seed,
            min_n=16,
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
        "Package: **networkx** for the incidence graph. Default cut = **own-ref subset scan** "
        "(level set `{φ≥τ}` / coverage prefix). The graph lifts a layer's Ŝ onto orders; "
        "it does not run community detection. Bundled / structural Louvain are contrast only.",
        "Shares are localization proxies, not a unique decomposition. Y is never a feature.",
        "",
        f"n_ref={N_REF}, n_new={N_NEW}, merchants={N_MERCHANTS}, batches={N_BATCHES}, onset={ONSET}. "
        "Planted region is **south** (second half of merchant ids).",
        "",
        "## Subset scan by layer (lift to orders; Louvain is contrast only)",
        "",
    ]
    gheader = [
        "kind",
        "t",
        "onset",
        "cut",
        "loud south_frac",
        "J(scan mer,south)",
        "J(coverage mer,south)",
        "J(mass mer,south)",
        "J(user scan,south)",
        "J(Louvain mer,south)",
        "n_loud mer",
        "FSDS selected",
        "fingerprint",
    ]
    lines.append(_md_row(gheader))
    lines.append(_md_row(["---"] * len(gheader)))
    for kind in KINDS:
        for rec in results[kind]["rows"]:
            g = rec.get("graph") or {}
            ls = g.get("level_set") or {}
            L = g.get("layers") or {}
            layer2 = rec.get("loud_subset")
            fp = ""
            south_frac = ""
            for p in rec.get("portraits") or []:
                if str(p.get("subset")) == str(layer2):
                    fp = p.get("fingerprint") or ""
                    south_frac = p.get("south_frac")
                    break
            mer = ls.get("merchant") or {}
            lines.append(
                _md_row(
                    [
                        kind,
                        rec["t"],
                        "yes" if rec["onset"] else "",
                        g.get("cut") or "",
                        _fmt(south_frac),
                        _fmt(ls.get("jaccard_merchant_vs_south")),
                        _fmt(ls.get("jaccard_coverage_merchant_vs_south")),
                        _fmt(ls.get("jaccard_mass_merchant_vs_south")),
                        _fmt(ls.get("jaccard_user_vs_south")),
                        _fmt(L.get("jaccard_merchant_vs_south")),
                        mer.get("n_loud_nodes"),
                        ",".join(rec["fsds"]["selected_names"]),
                        fp,
                    ]
                )
            )
    lines += [
        "",
        "## Own-ref vs full-ref (node clock) and multi-layer lift-to-order",
        "",
        _md_row(
            [
                "kind",
                "t",
                "south own MMD",
                "north own MMD",
                "south own ‖ΔX‖",
                "north own ‖ΔX‖",
                "J(scan mer,south)",
                "J(coverage mer,south)",
                "J(user scan,south)",
                "MMD-slice south_frac",
                "ΔY-slice south_frac",
            ]
        ),
        _md_row(["---"] * 11),
    ]
    for kind in KINDS:
        for rec in results[kind]["rows"]:
            if not rec.get("onset"):
                continue
            g = rec.get("graph") or {}
            o = g.get("own_vs_full") or {}
            ls = g.get("level_set") or {}
            mer = ls.get("merchant") or {}
            slices = mer.get("slices") or {}
            lines.append(
                _md_row(
                    [
                        kind,
                        rec["t"],
                        _fmt(o.get("south_own_mmd")),
                        _fmt(o.get("north_own_mmd")),
                        _fmt(o.get("south_own_cmean_x")),
                        _fmt(o.get("north_own_cmean_x")),
                        _fmt(ls.get("jaccard_merchant_vs_south")),
                        _fmt(ls.get("jaccard_coverage_merchant_vs_south")),
                        _fmt(ls.get("jaccard_user_vs_south")),
                        _fmt((slices.get("mmd") or {}).get("south_frac")),
                        _fmt((slices.get("cmean_y") or {}).get("south_frac")),
                    ]
                )
            )
    lines += [
        "",
        "## Subset portraits (scan loud vs other · MMD + PO + CMean)",
        "",
    ]
    header = [
        "kind",
        "t",
        "subset",
        "n",
        "south_frac",
        "π_MMD",
        "π_PO",
        "π_CMean",
        "MMD",
        "PO",
        "ΔE[Y]",
        "‖ΔE[X]‖",
        "fingerprint",
    ]
    lines.append(_md_row(header))
    lines.append(_md_row(["---"] * len(header)))
    for kind in KINDS:
        for rec in results[kind]["rows"]:
            if not rec.get("onset"):
                continue
            for p in rec.get("portraits") or []:
                lines.append(
                    _md_row(
                        [
                            kind,
                            rec["t"],
                            p.get("subset"),
                            p.get("n"),
                            _fmt(p.get("south_frac")),
                            _fmt(p.get("pi_mmd")),
                            _fmt(p.get("pi_po")),
                            _fmt(p.get("pi_cmean")),
                            _fmt(p.get("mmd")),
                            _fmt(p.get("po")),
                            _fmt(p.get("cmean_y")),
                            _fmt(p.get("cmean_x")),
                            p.get("fingerprint") or "",
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

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6), sharey=True)
    metrics = [("pi_mmd", "π_MMD"), ("pi_po", "π_PO"), ("pi_cmean", "π_CMean")]
    for ax, kind in zip(axes, KINDS):
        rec = results[kind]["rows"][-1]
        portraits = rec.get("portraits") or []
        if not portraits:
            ax.set_title(kind)
            continue
        x = np.arange(len(metrics))
        w = 0.8 / max(len(portraits), 1)
        for i, p in enumerate(portraits[:4]):
            vals = [float(p.get(k) or 0.0) for k, _ in metrics]
            frac = float(p.get("south_frac") or 0.0)
            color = (0.75, 0.18, 0.15, 0.35 + 0.65 * frac)
            ax.bar(
                x + (i - 0.5 * (len(portraits[:4]) - 1)) * w,
                vals,
                w * 0.9,
                label=f"{p['subset']} s={frac:.2f}",
                color=color,
            )
        ax.set_xticks(x)
        ax.set_xticklabels([t for _, t in metrics])
        ax.set_title(kind)
        ax.set_ylim(0, 1.05)
        ax.grid(axis="y", alpha=0.3)
        ax.legend(frameon=False, fontsize=7)
    axes[0].set_ylabel("subset share")
    fig.suptitle("Subset-scan loud vs other · last batch · redder = higher south fraction")
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
                    "graph": rec.get("graph"),
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
