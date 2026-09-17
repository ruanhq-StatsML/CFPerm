#!/usr/bin/env python3
"""Multimodal attribution (CFPerm): permute T, report modality blocks.

Group attribution asks which *coordinates* X_j carry a group CATE.
This object uses the same permute-T loop, but the unit is a **modality block**:

    X = [text | vision-or-dense | graph | fusion]
    T = pack version / community / queue
    permute T, not X, not pixels
    I_m = sum_{j in modality m} importance_j
    reject a block iff it clears the across-block threshold

Live tables already on disk (no raw text, no images):

    hybrid dense hop     T = dense-channel pack after the cut
                         blocks = query / sparse / dense / fusion
    Graph-RAG community  T = local pack vs community pack
                         blocks = query-seed / graph topology

Y is still the product label. HH chosen is not Y.

Usage::

    PYTHONPATH=. python3 scripts/prototype_multimodal_attribution.py
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.prototype_group_attribution import (  # noqa: E402
    LEAK,
    LEVEL_ACROSS,
    LEVEL_FEATURE,
    N_PERM,
    SEED,
    TOP_K,
    interaction_vimp,
    zscore,
)
from scripts.posthoc_localization import localize  # noqa: E402

HYBRID = ROOT / "results" / "manuscript" / "hybrid_retrieval"
GRAPH = ROOT / "results" / "manuscript" / "graph_rag_batches"
OUT = ROOT / "results" / "manuscript" / "multimodal_attribution"
CUT_BATCH = 4

HYBRID_BLOCKS = {
    "query": ["x_q_toks", "x_q_chars", "x_qmark", "x_q_ents"],
    "sparse": ["x_sparse_margin", "x_ks"],
    "dense": ["x_dense_margin", "x_kd"],
    "fusion": [
        "x_n_cand",
        "x_js_overlap",
        "x_rank_corr",
        "x_rrf_top1_mass",
        "x_fuse_uniq",
        "x_mean_q_overlap",
    ],
}

GRAPH_BLOCKS = {
    "query": ["x_n_q_seeds", "x_seed_frac"],
    "graph": ["x_n_nodes", "x_n_edges", "x_mean_deg", "x_n_cc", "x_lcc_frac"],
}


def x_columns(fieldnames) -> list[str]:
    return [c for c in fieldnames if c.startswith("x_")]


def block_index(names: list[str], blocks: dict[str, list[str]]) -> dict[str, np.ndarray]:
    pos = {n: i for i, n in enumerate(names)}
    out = {}
    for mod, cols in blocks.items():
        missing = [c for c in cols if c not in pos]
        if missing:
            raise KeyError(f"{mod} missing {missing}")
        out[mod] = np.asarray([pos[c] for c in cols], dtype=int)
    return out


def cfperm_blocks(
    X,
    Y,
    T,
    blocks: dict[str, list[str]] | dict[str, np.ndarray],
    *,
    names: list[str] | None = None,
    n_perm: int = N_PERM,
    level_feature: float = LEVEL_FEATURE,
    level_across_feature: float = LEVEL_ACROSS,
    top_k: int = TOP_K,
    seed: int = SEED,
):
    """Same CFPerm decision as groups, then sum importance inside each modality."""
    X = zscore(np.asarray(X, dtype=float))
    Y = np.asarray(Y, dtype=float).ravel()
    T = np.asarray(T, dtype=int).ravel()
    p = X.shape[1]
    names = list(names) if names is not None else [f"x{j}" for j in range(p)]
    sample = next(iter(blocks.values()))
    first = sample[0] if len(sample) else None
    if isinstance(first, str):
        idx = block_index(names, blocks)  # type: ignore[arg-type]
    else:
        idx = {m: np.asarray(v, dtype=int) for m, v in blocks.items()}
    mods = list(idx)
    rng = np.random.default_rng(seed)
    imp = interaction_vimp(X, Y, T)
    perm = np.zeros((p, n_perm), dtype=float)
    for b in range(n_perm):
        perm[:, b] = interaction_vimp(X, Y, rng.permutation(T))

    def _decide(obs: np.ndarray, null: np.ndarray, labels: list[str]) -> dict:
        pvals = (1.0 + np.sum(null >= obs[:, None], axis=1)) / (1.0 + n_perm)
        q1_upper = np.quantile(null, 1.0 - level_feature, axis=1)
        threshold = float(np.quantile(q1_upper, 1.0 - level_across_feature))
        hits = obs > threshold
        order = np.argsort(-obs)
        return {
            "names": labels,
            "imp": obs,
            "pvals": pvals,
            "q1_upper": q1_upper,
            "threshold": threshold,
            "hits": hits.astype(int),
            "rejected": int(int(hits.sum()) >= int(top_k)),
            "n_hits": int(hits.sum()),
            "top": [labels[int(j)] for j in order[: min(5, len(labels))]],
            "hit_names": [labels[i] for i, h in enumerate(hits) if h],
        }

    feat = _decide(imp, perm, names)
    block_obs = np.asarray([float(imp[idx[m]].sum()) for m in mods], dtype=float)
    block_null = np.asarray(
        [[float(perm[idx[m], b].sum()) for b in range(n_perm)] for m in mods],
        dtype=float,
    )
    blk = _decide(block_obs, block_null, mods)
    return {
        "n": int(len(Y)),
        "n_groups": int(len(np.unique(T))),
        "group_counts": {str(int(g)): int(np.sum(T == g)) for g in np.unique(T)},
        "y_by_group": {str(int(g)): float(np.mean(Y[T == g])) for g in np.unique(T)},
        "n_perm": int(n_perm),
        "level_feature": float(level_feature),
        "level_across_feature": float(level_across_feature),
        "top_k": int(top_k),
        "feature": feat,
        "block": blk,
        "blocks": {m: [names[int(j)] for j in idx[m]] for m in mods},
    }


def load_xy(path: Path):
    with path.open() as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"empty {path}")
    keys = {k.lower() for k in rows[0]}
    for bad in LEAK:
        if bad in keys:
            raise RuntimeError(f"{path} still has {bad} — that is not Y")
    cols = x_columns(rows[0].keys())
    ycol = "Y" if "Y" in rows[0] else "y"
    Y = np.asarray([int(float(r[ycol])) for r in rows], dtype=int)
    batch = np.asarray([int(r["batch"]) for r in rows], dtype=int)
    X = np.asarray([[float(r[c]) for c in cols] for r in rows], dtype=float)
    return X, Y, batch, cols


def load_hybrid_dense_hop(path: Path):
    """T = dense pack version (after the labeled cut). Blocks are retrieval modalities."""
    X, Y, batch, cols = load_xy(path)
    T = (batch >= CUT_BATCH).astype(int)
    return X, Y, T, cols, HYBRID_BLOCKS


def load_graph_community(quiet: Path, hop: Path):
    """T = local pack vs community pack. Blocks = query-seed vs graph topology."""
    X0, Y0, _b0, cols = load_xy(quiet)
    X1, Y1, _b1, cols1 = load_xy(hop)
    if cols != cols1:
        raise RuntimeError("graph quiet/hop X mismatch")
    X = np.vstack([X0, X1])
    Y = np.concatenate([Y0, Y1])
    T = np.asarray([0] * len(Y0) + [1] * len(Y1), dtype=int)
    return X, Y, T, cols, GRAPH_BLOCKS


def planted_multimodal(n=600, seed=0):
    """Three modalities. Only vision interacts with T."""
    rng = np.random.default_rng(seed)
    text = rng.normal(size=(n, 4))
    vision = rng.normal(size=(n, 4))
    graph = rng.normal(size=(n, 4))
    T = rng.integers(0, 2, size=n)
    Y = 0.35 * text[:, 0] + 2.5 * T.astype(float) * vision[:, 1] + 0.25 * rng.normal(size=n)
    X = np.hstack([text, vision, graph])
    names = [f"text{j}" for j in range(4)] + [f"vision{j}" for j in range(4)] + [f"graph{j}" for j in range(4)]
    blocks = {
        "text": names[:4],
        "vision": names[4:8],
        "graph": names[8:],
    }
    return X, Y, T, names, blocks


def null_multimodal(n=600, seed=1):
    rng = np.random.default_rng(seed)
    text = rng.normal(size=(n, 4))
    vision = rng.normal(size=(n, 4))
    graph = rng.normal(size=(n, 4))
    T = rng.integers(0, 2, size=n)
    Y = 0.8 * text[:, 0] - 0.4 * graph[:, 2] + 0.25 * rng.normal(size=n)
    X = np.hstack([text, vision, graph])
    names = [f"text{j}" for j in range(4)] + [f"vision{j}" for j in range(4)] + [f"graph{j}" for j in range(4)]
    blocks = {
        "text": names[:4],
        "vision": names[4:8],
        "graph": names[8:],
    }
    return X, Y, T, names, blocks


STREAMS = [
    {
        "name": "hybrid_dense_hop",
        "title": "hybrid: dense pack hop",
        "t_meaning": "T=0 before dense cut, T=1 after",
        "y_meaning": "fused answer valid",
        "kind": "hybrid_hop",
        "path": HYBRID / "xy_hotpot_hybrid_hop.csv",
    },
    {
        "name": "graph_community",
        "title": "Graph-RAG: local vs community",
        "t_meaning": "T=0 local pack, T=1 community pack",
        "y_meaning": "supporting nodes in the pack",
        "kind": "graph_stack",
        "quiet": GRAPH / "xy_graph_query.csv",
        "hop": GRAPH / "xy_graph_query_hop.csv",
    },
]


def jsonable(obj):
    if isinstance(obj, dict):
        return {k: jsonable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return jsonable(obj.tolist())
    if isinstance(obj, (np.floating, float)):
        x = float(obj)
        return None if not np.isfinite(x) else x
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    return obj


def compact(rec: dict, spec: dict) -> dict:
    blk = rec["block"]
    loc = rec.get("localization") or {}
    groups = loc.get("groups") or {}
    return {
        "name": spec["name"],
        "title": spec["title"],
        "t_meaning": spec["t_meaning"],
        "y_meaning": spec["y_meaning"],
        "n": rec["n"],
        "n_groups": rec["n_groups"],
        "rejected": blk["rejected"],
        "n_hits": blk["n_hits"],
        "hits": blk["hit_names"],
        "top": blk["top"],
        "feature_top": rec["feature"]["top"][:3],
        "blocks": rec["blocks"],
        "loc_sig_pairs": groups.get("n_sig_pairs"),
        "loc_pairwise": groups.get("pairwise"),
    }


def plot_blocks(runs: list[tuple[dict, dict]], path: Path) -> None:
    n = len(runs)
    fig, axes = plt.subplots(n, 1, figsize=(8.4, 2.6 * n))
    if n == 1:
        axes = [axes]
    for ax, (rec, spec) in zip(axes, runs):
        blk = rec["block"]
        names = blk["names"]
        imp = blk["imp"]
        colors = ["#9b2c2c" if h else "#1f4e79" for h in blk["hits"]]
        ax.bar(range(len(names)), imp, color=colors)
        ax.axhline(blk["threshold"], color="#9b2c2c", ls="--", lw=1.0)
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names, fontsize=11)
        flag = "REJECT" if blk["rejected"] else "quiet"
        ax.set_title(f"{spec['title']}  ·  {flag}  ·  top {', '.join(blk['top'][:3])}", fontsize=10)
        ax.set_ylabel("block imp")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def render_report(rows: list[dict], extras: list[dict]) -> str:
    lines = [
        "# Multimodal attribution — permute T, report blocks",
        "",
        "Same CFPerm loop as group attribution. The unit is a **modality block**, not a single x_*.",
        "X = concat(text, vision/dense, graph, fusion). T is the pack/group. Permute T.",
        "No pixels, no raw prompts, no HH chosen.",
        "",
        "```bash",
        "PYTHONPATH=. python3 scripts/prototype_multimodal_attribution.py",
        "```",
        "",
        "| Stream | T | reject | block hits | top blocks | loc sig pairs |",
        "|---|---|---|---|---|---:|",
    ]
    for r in rows:
        hits = ", ".join(r["hits"]) if r["hits"] else "—"
        top = ", ".join(r["top"][:3])
        lines.append(
            f"| {r['title']} | {r['t_meaning']} | "
            f"{'yes' if r['rejected'] else 'no'} | {hits} | {top} | {r.get('loc_sig_pairs', '—')} |"
        )
    lines += [
        "",
        "Post-hoc: subset indices from T (and quartiles of the top coordinate), then pairwise **MMD** and **PO-risk**.",
        "Conditional means stay in the JSON; they are not the localization test.",
        "",
        "## Pairwise subset MMD / PO-risk",
        "",
        "| Stream | pair | n | MMD | MMD p | PO-risk | PO p | mean Y |",
        "|---|---|---|---:|---:|---:|---:|---|",
    ]
    for r in rows:
        for p in r.get("loc_pairwise") or []:
            lines.append(
                "| {title} | {a} vs {b} | {na}/{nb} | {mmd:.3g} | {mp:.3g} | {po:.3g} | {pp:.3g} | {ya:.3f}/{yb:.3f} |".format(
                    title=r["title"],
                    a=p["a"],
                    b=p["b"],
                    na=p["n_a"],
                    nb=p["n_b"],
                    mmd=p["mmd"],
                    mp=p["mmd_p"],
                    po=p.get("po_risk", float("nan")),
                    pp=p.get("po_p", float("nan")),
                    ya=p.get("mean_Y_a", float("nan")),
                    yb=p.get("mean_Y_b", float("nan")),
                )
            )
    lines += [
        "",
        "## Synthetic check",
        "",
        "| DGP | reject | top blocks |",
        "|---|---|---|",
    ]
    for e in extras:
        lines.append(
            f"| {e['title']} | {'yes' if e['rejected'] else 'no'} | {', '.join(e['top'][:3])} |"
        )
    lines += [
        "",
        "Planted: only the vision block interacts with T. Null: Y depends on text/graph, not on T.",
        "",
        "Not this object: time-gate hop, single-channel Recall, HH chosen, raw image/audio.",
        "",
    ]
    return "\n".join(lines) + "\n"


def run_stream(spec: dict) -> dict:
    if spec["kind"] == "hybrid_hop":
        X, Y, T, cols, blocks = load_hybrid_dense_hop(spec["path"])
    elif spec["kind"] == "graph_stack":
        X, Y, T, cols, blocks = load_graph_community(spec["quiet"], spec["hop"])
    else:
        raise ValueError(spec["kind"])
    rec = cfperm_blocks(X, Y, T, blocks, names=cols, seed=SEED)
    rec["name"] = spec["name"]
    rec["title"] = spec["title"]
    top_name = rec["feature"]["top"][0] if rec["feature"]["top"] else None
    feat = X[:, cols.index(top_name)] if top_name in cols else None
    rec["localization"] = localize(X, Y, group_labels=T, feature=feat, n_perm=25, seed=SEED)
    return rec


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    runs = []
    skipped = []
    for spec in STREAMS:
        paths = []
        if spec.get("path"):
            paths.append(spec["path"])
        for k in ("quiet", "hop"):
            if spec.get(k):
                paths.append(spec[k])
        if any(not p.exists() for p in paths):
            skipped.extend(str(p.relative_to(ROOT)) for p in paths if not p.exists())
            continue
        rec = run_stream(spec)
        runs.append((rec, spec))
        (OUT / f"{spec['name']}.json").write_text(
            json.dumps(jsonable(rec), indent=2) + "\n", encoding="utf-8"
        )
    extras = []
    for title, fn, seed in (
        ("planted CATE on vision block", planted_multimodal, 0),
        ("null: Y depends on X, not T", null_multimodal, 1),
    ):
        X, Y, T, names, blocks = fn(n=700, seed=seed)
        rec = cfperm_blocks(X, Y, T, blocks, names=names, n_perm=25, seed=seed)
        rec["localization"] = localize(X, Y, group_labels=T, n_perm=20, seed=seed)
        extras.append(
            {
                "title": title,
                "rejected": rec["block"]["rejected"],
                "top": rec["block"]["top"],
                "hits": rec["block"]["hit_names"],
            }
        )
        (OUT / ("planted.json" if "planted" in title else "null.json")).write_text(
            json.dumps(jsonable(rec), indent=2) + "\n", encoding="utf-8"
        )
    summary = [compact(rec, spec) for rec, spec in runs]
    if runs:
        plot_blocks(runs, OUT / "importance_by_block.png")
    report = render_report(summary, extras)
    if skipped:
        report += "Skipped missing files:\n\n" + "\n".join(f"- `{p}`" for p in skipped) + "\n"
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (OUT / "REPORT.md").write_text(report, encoding="utf-8")
    print(report)
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
