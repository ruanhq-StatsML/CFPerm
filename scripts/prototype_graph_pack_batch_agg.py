#!/usr/bin/env python3
"""Graph-RAG 图包 / community 换代: one serving batch = aggregated graph features.

Per query (Hotpot distractor pool, 10 titles):
  nodes = candidate titles
  edges = shared-token co-mention, or rewired after the cut
  X     = n_nodes, n_edges, mean_deg, n_cc, n_q_seeds, seed_frac, lcc_frac
  Y     = 1 iff every supporting title sits in the served pack
        (local: seed ∪ 1-hop; after community 换代: largest CC)

Per arrival window:
  X_batch = mean (and std) of those seven coordinates
  Y_batch = pack-usable rate

Hotpot file order is not wall-clock time. This script prototypes the table
shape. A live serving gate still needs a timestamped request log.

Usage::

    PYTHONPATH=. python3 scripts/prototype_graph_pack_batch_agg.py
"""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.build_hybrid_retrieval_xy import (  # noqa: E402
    GRAPH_COLS,
    GRAPH_DIMS,
    N_CAND,
    N_PER,
    N_ROWS,
    load_hotpot,
    pad_pool,
    query_seeds,
    tokenize,
    write_xy,
)
from scripts.llm_audit_online_bootstrap_prototype import freeze_and_score, last_two_hops

OUT = ROOT / "results" / "manuscript" / "graph_rag_batches"
CUT_BATCH = 4
N_REF = 4
GATE = 1.5
SEED = 2026
STD_COLS = [f"{c}_std" for c in GRAPH_COLS]
BATCH_COLS = GRAPH_COLS + STD_COLS


def co_mention_edges(titles: list[str]) -> list[tuple[int, int]]:
    tok_sets = [set(tokenize(t)) for t in titles]
    n = len(titles)
    edges: list[tuple[int, int]] = []
    for i in range(n):
        for j in range(i + 1, n):
            if tok_sets[i] & tok_sets[j]:
                edges.append((i, j))
    return edges


def rewire_edges(n: int, n_edges: int, rng: np.random.Generator) -> list[tuple[int, int]]:
    """Same edge count, random pairs: a cheap community recompute."""
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    if n <= 1 or n_edges <= 0 or not pairs:
        return []
    k = min(int(n_edges), len(pairs))
    pick = rng.choice(len(pairs), size=k, replace=False)
    return [pairs[int(i)] for i in np.atleast_1d(pick)]


def graph_features(n: int, edges: list[tuple[int, int]], seeds: np.ndarray) -> dict:
    parent = list(range(n))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    deg = np.zeros(n, dtype=float)
    for i, j in edges:
        deg[i] += 1.0
        deg[j] += 1.0
        pi, pj = find(i), find(j)
        if pi != pj:
            parent[pj] = pi
    roots = [find(i) for i in range(n)] if n else []
    n_cc = len(set(roots)) if n else 0
    lcc = max(Counter(roots).values()) if roots else 0
    return {
        "x_n_nodes": float(n) / 10.0,
        "x_n_edges": float(len(edges)) / 20.0,
        "x_mean_deg": float(deg.mean()) / 6.0 if n else 0.0,
        "x_n_cc": float(n_cc) / 10.0,
        "x_n_q_seeds": float(np.asarray(seeds).sum()) / 6.0,
        "x_seed_frac": float(np.asarray(seeds).mean()) if n else 0.0,
        "x_lcc_frac": float(lcc / n) if n else 0.0,
    }


def largest_cc(n: int, edges: list[tuple[int, int]]) -> set[int]:
    parent = list(range(n))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j in edges:
        pi, pj = find(i), find(j)
        if pi != pj:
            parent[pj] = pi
    roots = [find(i) for i in range(n)] if n else []
    if not roots:
        return set()
    top = Counter(roots).most_common(1)[0][0]
    return {i for i, r in enumerate(roots) if r == top}


def served_pack(
    n: int,
    edges: list[tuple[int, int]],
    seeds: np.ndarray,
    *,
    community: bool = False,
) -> set[int]:
    """Local pack = seed ∪ 1-hop. Community pack = largest connected component."""
    if community:
        return largest_cc(n, edges)
    adj = [set() for _ in range(n)]
    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)
    pack = {i for i, s in enumerate(np.asarray(seeds).ravel()) if s}
    extra: set[int] = set()
    for i in pack:
        extra |= adj[i]
    return pack | extra


def pack_usable(
    n: int,
    edges: list[tuple[int, int]],
    seeds: np.ndarray,
    gold_idx: set[int],
    *,
    community: bool = False,
) -> int:
    """Y = 1 iff every supporting title sits in the served pack."""
    pack = served_pack(n, edges, seeds, community=community)
    return int(bool(gold_idx) and gold_idx <= pack)


def query_graph_row(
    question, context, supporting, *, batch: int, edges=None, community: bool = False
) -> dict:
    titles, _docs = pad_pool(
        context["title"],
        [" ".join(s) for s in context["sentences"]],
        N_CAND,
    )
    gold = set(supporting["title"])
    gold_idx = {i for i, t in enumerate(titles) if t in gold and t}
    seeds = query_seeds(question, titles)
    if edges is None:
        edges = co_mention_edges(titles)
    feats = graph_features(len(titles), edges, seeds)
    return {
        "y": pack_usable(len(titles), edges, seeds, gold_idx, community=community),
        "batch": int(batch),
        "n_raw_edges": len(edges),
        **feats,
    }


def aggregate_batch_rows(query_rows: list[dict]) -> list[dict]:
    """One serving window = mean/std of graph features + pack-usable rate."""
    by_batch: dict[int, list[dict]] = defaultdict(list)
    for r in query_rows:
        by_batch[int(r["batch"])].append(r)
    out = []
    for b in sorted(by_batch):
        block = by_batch[b]
        rec = {
            "y": float(np.mean([int(r["y"]) for r in block])),
            "batch": int(b),
            "n_queries": len(block),
        }
        for c in GRAPH_COLS:
            xs = np.asarray([float(r[c]) for r in block], dtype=float)
            rec[c] = float(xs.mean())
            rec[f"{c}_std"] = float(xs.std(ddof=0))
        out.append(rec)
    return out


def assign_window(i: int, n_per: int = N_PER) -> int:
    return int(i) // int(n_per)


def build_query_tables(questions, supporting, contexts, *, hop: bool, rng: np.random.Generator):
    rows = []
    for i, (q, sf, ctx) in enumerate(zip(questions, supporting, contexts)):
        batch = assign_window(i)
        titles, _ = pad_pool(ctx["title"], [" ".join(s) for s in ctx["sentences"]], N_CAND)
        native = co_mention_edges(titles)
        if hop and batch >= CUT_BATCH:
            edges = rewire_edges(len(titles), len(native), rng)
            community = True
        else:
            edges = native
            community = False
        rows.append(
            query_graph_row(q, ctx, sf, batch=batch, edges=edges, community=community)
        )
    n_use = (len(rows) // N_PER) * N_PER
    return rows[:n_use]


def write_batch_xy(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["y", "batch", "n_queries", *BATCH_COLS]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(
                {
                    "y": f"{float(r['y']):.6g}",
                    "batch": int(r["batch"]),
                    "n_queries": int(r["n_queries"]),
                    **{c: f"{float(r[c]):.6g}" for c in BATCH_COLS},
                }
            )


def matrix_from_query(rows: list[dict]):
    X = np.asarray([[float(r[c]) for c in GRAPH_COLS] for r in rows], dtype=float)
    y = np.asarray([int(r["y"]) for r in rows], dtype=int)
    batch = np.asarray([int(r["batch"]) for r in rows], dtype=int)
    return X, y, batch


def summarize_regime(name: str, query_rows: list[dict], batch_rows: list[dict]) -> dict:
    X, y, batch = matrix_from_query(query_rows)
    frozen = freeze_and_score(X, y, batch, n_ref_batches=N_REF, seed=SEED)
    hops = last_two_hops(X, y, batch, gate=GATE, seed=SEED)
    cut_hops = [h for h in hops if h["abs_batch"] == CUT_BATCH]
    y_pre = [int(r["y"]) for r in query_rows if int(r["batch"]) < CUT_BATCH]
    y_post = [int(r["y"]) for r in query_rows if int(r["batch"]) >= CUT_BATCH]
    edges_pre = [int(r["n_raw_edges"]) for r in query_rows if int(r["batch"]) < CUT_BATCH]
    edges_post = [int(r["n_raw_edges"]) for r in query_rows if int(r["batch"]) >= CUT_BATCH]
    trail = list(frozen["trail_batches"])
    delta = [float(s - frozen["mu_ref"]) for s in frozen["scores"]]
    return {
        "name": name,
        "n_query": len(query_rows),
        "n_batch": len(batch_rows),
        "y_rate": float(y.mean()),
        "y_pre": float(np.mean(y_pre)) if y_pre else None,
        "y_post": float(np.mean(y_post)) if y_post else None,
        "mean_edges_pre": float(np.mean(edges_pre)) if edges_pre else None,
        "mean_edges_post": float(np.mean(edges_post)) if edges_post else None,
        "mu_ref": float(frozen["mu_ref"]),
        "mean_delta": float(np.mean(delta)) if delta else None,
        "delta_by_trail": [{"batch": int(b), "delta": float(d)} for b, d in zip(trail, delta)],
        "hop_fire_at_cut": bool(cut_hops and cut_hops[0]["fired"]),
        "n_hop_fires": int(sum(1 for h in hops if h["fired"])),
        "cut_hop": cut_hops[0] if cut_hops else None,
        "hops": hops,
    }


def plot_regimes(native: dict, hop: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6))
    for ax, rec, title in (
        (axes[0], native, "native local pack"),
        (axes[1], hop, "community pack after cut"),
    ):
        xs = [d["batch"] for d in rec["delta_by_trail"]]
        ys = [d["delta"] for d in rec["delta_by_trail"]]
        ax.axhline(0.0, color="#888", lw=0.8)
        ax.axvline(CUT_BATCH - 0.5, color="#9b2c2c", ls="--", lw=1.0, label="cut")
        ax.plot(xs, ys, marker="o", color="#1f4e79")
        ax.set_title(title)
        ax.set_xlabel("batch")
        ax.set_ylabel(r"$\Delta_t$")
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def write_report(native: dict, hop: dict, path: Path) -> None:
    def hop_row(h):
        prev = "—" if h["e_prev"] is None else f"{h['e_prev']:.3f}"
        fired = "yes" if h["fired"] else "no"
        return f"| {h['abs_batch']} | {fired} | {prev} | {h['e_now']:.3f} | {h['ratio']:.3f} |"

    lines = [
        "# Graph-RAG batch aggregate prototype",
        "",
        "One serving batch = mean/std of title-graph features. "
        "Native Y: supporting titles sit in seed ∪ 1-hop. "
        "Hop Y: supporting titles sit in the largest CC after an edge rewire "
        "(community 换代). Hotpot file order is not wall-clock time.",
        "",
        "| Regime | y pre | y post | mean Δ | hop@cut | n hop fires |",
        "|---|---:|---:|---:|---|---:|",
        (
            f"| native | {native['y_pre']:.3f} | {native['y_post']:.3f} | "
            f"{native['mean_delta']:.3f} | "
            f"{'yes' if native['hop_fire_at_cut'] else 'no'} | {native['n_hop_fires']} |"
        ),
        (
            f"| community 换代 | {hop['y_pre']:.3f} | {hop['y_post']:.3f} | "
            f"{hop['mean_delta']:.3f} | "
            f"{'yes' if hop['hop_fire_at_cut'] else 'no'} | {hop['n_hop_fires']} |"
        ),
        "",
        "## Last-two around the cut (community 换代)",
        "",
        "| abs batch | fire | e_prev | e_now | ratio |",
        "|---:|---|---:|---:|---:|",
    ]
    for h in hop["hops"]:
        if 2 <= h["abs_batch"] <= 6:
            lines.append(hop_row(h))
    lines += [
        "",
        "Rebuild: `PYTHONPATH=. python3 scripts/prototype_graph_pack_batch_agg.py`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    print("load HotpotQA distractor validation")
    questions, supporting, contexts = load_hotpot(N_ROWS)
    native_q = build_query_tables(
        questions, supporting, contexts, hop=False, rng=np.random.default_rng(SEED)
    )
    hop_q = build_query_tables(
        questions, supporting, contexts, hop=True, rng=np.random.default_rng(SEED)
    )
    native_b = aggregate_batch_rows(native_q)
    hop_b = aggregate_batch_rows(hop_q)
    write_xy(OUT / "xy_graph_query.csv", native_q, GRAPH_COLS)
    write_xy(OUT / "xy_graph_query_hop.csv", hop_q, GRAPH_COLS)
    write_batch_xy(OUT / "xy_graph_batch.csv", native_b)
    write_batch_xy(OUT / "xy_graph_batch_hop.csv", hop_b)

    native_sum = summarize_regime("native", native_q, native_b)
    hop_sum = summarize_regime("community_rewire", hop_q, hop_b)
    plot_regimes(native_sum, hop_sum, OUT / "delta_by_batch.png")
    write_report(native_sum, hop_sum, OUT / "REPORT.md")
    native_manifest = {k: native_sum[k] for k in native_sum if k not in {"delta_by_trail", "hops", "cut_hop"}}
    hop_manifest = {k: hop_sum[k] for k in hop_sum if k not in {"delta_by_trail", "hops", "cut_hop"}}

    manifest = {
        "note": (
            "Graph-RAG serving table: one batch is the aggregate of subgraph "
            "features. Native Y: supporting titles in seed ∪ 1-hop. "
            "Hop: rewire edges and serve the largest CC after CUT_BATCH "
            "(graph-pack / community refresh). Hotpot order is not wall-clock time."
        ),
        "cut_batch": CUT_BATCH,
        "n_per": N_PER,
        "graph_x": GRAPH_DIMS,
        "y": "1 iff supporting titles sit in the served pack (local 1-hop vs LCC)",
        "files": {
            "xy_graph_query.csv": {"n": len(native_q), "shape": "one query, graph X"},
            "xy_graph_query_hop.csv": {
                "n": len(hop_q),
                "hop": "rewire edges + serve largest CC after cut",
            },
            "xy_graph_batch.csv": {"n": len(native_b), "shape": "one window, aggregated X"},
            "xy_graph_batch_hop.csv": {"n": len(hop_b)},
        },
        "native": native_manifest,
        "community_rewire": hop_manifest,
    }
    (OUT / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (OUT / "PREPARED_XY.md").write_text(
        "\n".join(
            [
                "# Graph-RAG serving windows (batch-aggregated graph features)",
                "",
                "Hotpot distractor validation. Questions and wiki text are not stored.",
                "",
                "This is a **table-shape prototype**, not a wall-clock stream. "
                "Validation queries have no arrival order; window index follows file order.",
                "",
                "## Query table",
                "",
                "`xy_graph_query.csv`: one row per query. X is title-graph geometry. "
                "Y=1 iff every supporting title sits in the served pack "
                "(local: seed ∪ 1-hop).",
                "",
                "## Batch table",
                "",
                "`xy_graph_batch.csv`: **one row per serving window**. "
                "X is mean and std of the seven graph coordinates on that window. "
                "Y is the pack-usable rate. A graph-pack / community refresh rewires "
                "edges and serves the largest connected component, then aggregates again.",
                "",
                "Hop files apply that community pack after `batch >= 4`.",
                "",
                "Rebuild:",
                "",
                "```bash",
                "PYTHONPATH=. python3 scripts/prototype_graph_pack_batch_agg.py",
                "```",
                "",
            ]
        ),
        encoding="utf-8",
    )
    print(json.dumps(hop_manifest, indent=2))
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
