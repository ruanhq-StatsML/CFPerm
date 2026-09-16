#!/usr/bin/env python3
"""Graph continuity: user-id → list → aggregate, then GraphSAGE ⊕ hand features.

Serving log
    user_id mapsto  list of title-graphs in an arrival window
    aggregate the list (mean pool)
    X_batch = [ GraphSAGE-mean(list)  |  hand graph features(list) ]
    Y       = pack usable (supporting titles in the served pack)

OnlineRFPerm (Algorithm 1 in the Sep 15 writeup), on this table only:
    fit f_ref once on D_ref
    T_t = MSE_t − E_ref
    rank p-value against the historical T pool
    last-two hop_fires is the adjacent-window gate

GraphSAGE here is the frozen mean aggregator (no online fine-tune), matching
the reference-model protocol. Hand features are the seven title-graph
coordinates already on disk.

Hotpot has no real user_id; query index hashed into N_USERS prototypes the
log shape. File order is not wall-clock time.

Usage::

    PYTHONPATH=. python3 scripts/prototype_graph_continuity.py
"""
from __future__ import annotations

import csv
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from agod.rf_probe import (  # noqa: E402
    brier_score,
    error_floor,
    fit_online_rf,
    hop_fires,
    probe_err,
    shift_ratio,
)
from scripts.build_hybrid_retrieval_xy import (  # noqa: E402
    GRAPH_COLS,
    N_CAND,
    N_PER,
    N_ROWS,
    load_hotpot,
    pad_pool,
    query_seeds,
)
from scripts.prototype_graph_pack_batch_agg import (  # noqa: E402
    CUT_BATCH,
    co_mention_edges,
    graph_features,
    pack_usable,
    rewire_edges,
)

OUT = ROOT / "results" / "manuscript" / "graph_continuity"
N_USERS = 20
N_REF = 4
GATE = 1.5
SEED = 2026
SAGE_LAYERS = 2
NODE_DIM = 4  # seed, deg, slot, in_lcc
SAGE_DIM = NODE_DIM * (2**SAGE_LAYERS)  # concat each layer
HAND_DIM = len(GRAPH_COLS)
X_DIM = SAGE_DIM + HAND_DIM
SAGE_COLS = [f"x_sage_{i}" for i in range(SAGE_DIM)]
ALL_COLS = SAGE_COLS + GRAPH_COLS
ALPHA = 0.05


def adj_list(n: int, edges: list[tuple[int, int]]) -> list[list[int]]:
    adj = [[] for _ in range(n)]
    for i, j in edges:
        adj[i].append(j)
        adj[j].append(i)
    return adj


def in_lcc(n: int, edges: list[tuple[int, int]]) -> np.ndarray:
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
        return np.zeros(n, dtype=float)
    from collections import Counter

    top = Counter(roots).most_common(1)[0][0]
    return np.asarray([1.0 if r == top else 0.0 for r in roots], dtype=float)


def node_features(n: int, edges: list[tuple[int, int]], seeds: np.ndarray) -> np.ndarray:
    deg = np.zeros(n, dtype=float)
    for i, j in edges:
        deg[i] += 1.0
        deg[j] += 1.0
    slot = (np.arange(n, dtype=float) + 1.0) / max(n, 1)
    lcc = in_lcc(n, edges)
    seeds = np.asarray(seeds, dtype=float).ravel()
    if len(seeds) < n:
        seeds = np.pad(seeds, (0, n - len(seeds)))
    return np.column_stack(
        [
            seeds[:n],
            deg / 6.0,
            slot,
            lcc,
        ]
    )


def graphsage_mean(node_x: np.ndarray, edges: list[tuple[int, int]], layers: int = SAGE_LAYERS) -> np.ndarray:
    """Frozen GraphSAGE-mean: concat(self, mean-neighbors), then mean-pool."""
    h = np.asarray(node_x, dtype=float)
    n = len(h)
    if n == 0:
        dim = NODE_DIM * (2 ** int(layers))
        return np.zeros(dim, dtype=float)
    adj = adj_list(n, edges)
    for _ in range(int(layers)):
        neigh = np.zeros_like(h)
        for v in range(n):
            if adj[v]:
                neigh[v] = h[np.asarray(adj[v], dtype=int)].mean(axis=0)
            else:
                neigh[v] = h[v]
        h = np.concatenate([h, neigh], axis=1)
    return h.mean(axis=0)


def query_record(question, context, supporting, *, community: bool, edges) -> dict:
    titles, _ = pad_pool(context["title"], [" ".join(s) for s in context["sentences"]], N_CAND)
    gold = set(supporting["title"])
    gold_idx = {i for i, t in enumerate(titles) if t in gold and t}
    seeds = query_seeds(question, titles)
    n = len(titles)
    hand = graph_features(n, edges, seeds)
    sage = graphsage_mean(node_features(n, edges, seeds), edges)
    y = pack_usable(n, edges, seeds, gold_idx, community=community)
    rec = {"y": int(y), "n_raw_edges": len(edges)}
    for i, v in enumerate(sage):
        rec[f"x_sage_{i}"] = float(v)
    rec.update(hand)
    return rec


def user_of(i: int) -> int:
    return int(i) % N_USERS


def aggregate_user_lists(query_rows: list[dict]) -> list[dict]:
    """user_id → list of graphs in a window → mean pool SAGE ⊕ hand, Y = pack rate."""
    buckets: dict[tuple[int, int], list[dict]] = defaultdict(list)
    for r in query_rows:
        buckets[(int(r["batch"]), int(r["user_id"]))].append(r)
    out = []
    for (batch, user), block in sorted(buckets.items()):
        rec = {
            "y": float(np.mean([int(r["y"]) for r in block])),
            "batch": int(batch),
            "user_id": int(user),
            "n_list": len(block),
        }
        for c in ALL_COLS:
            rec[c] = float(np.mean([float(r[c]) for r in block]))
        out.append(rec)
    return out


def build_query_rows(questions, supporting, contexts, *, hop: bool, rng: np.random.Generator):
    rows = []
    for i, (q, sf, ctx) in enumerate(zip(questions, supporting, contexts)):
        batch = i // N_PER
        titles, _ = pad_pool(ctx["title"], [" ".join(s) for s in ctx["sentences"]], N_CAND)
        native = co_mention_edges(titles)
        if hop and batch >= CUT_BATCH:
            edges = rewire_edges(len(titles), len(native), rng)
            community = True
        else:
            edges = native
            community = False
        rec = query_record(q, ctx, sf, community=community, edges=edges)
        rec["batch"] = int(batch)
        rec["user_id"] = user_of(i)
        rec["query_index"] = int(i)
        rows.append(rec)
    n_use = (len(rows) // N_PER) * N_PER
    return rows[:n_use]


def matrix(rows: list[dict], y_int: bool = False):
    X = np.asarray([[float(r[c]) for c in ALL_COLS] for r in rows], dtype=float)
    y = np.asarray([r["y"] for r in rows], dtype=float)
    if y_int:
        y = (y >= 0.5).astype(int)
    batch = np.asarray([int(r["batch"]) for r in rows], dtype=int)
    return X, y, batch


def rank_pvalues(T: np.ndarray) -> np.ndarray:
    """OnlineRFPerm S3: p_k = #{i < k: T_k ≤ T_i} / k. Small p ⇒ current T is extreme."""
    T = np.asarray(T, dtype=float).ravel()
    p = np.ones(len(T), dtype=float)
    for k in range(1, len(T)):
        p[k] = float(np.sum(T[k] <= T[:k]) + 1) / float(k + 1)
    return p


def freeze_rf_mse(X, y, batch, *, n_ref: int, seed: int):
    """Algorithm 1 S1–S3 with a frozen RF. Y may be 0/1 or a pack-usable rate."""
    batch = np.asarray(batch, dtype=int)
    y = np.asarray(y, dtype=float).ravel()
    X = np.asarray(X, dtype=float)
    ref = batch < int(n_ref)
    y_ref = y[ref]
    # 0/1 → classifier Brier; rate → regressor MSE
    binary = set(np.unique(np.round(y_ref, 6))) <= {0.0, 1.0}
    if binary:
        probe = fit_online_rf(X[ref], (y_ref >= 0.5).astype(int), seed=seed, task="acc")

        def mse(mask):
            return brier_score(probe, X[mask], (y[mask] >= 0.5).astype(int))

    else:
        probe = fit_online_rf(X[ref], y_ref, seed=seed, task="mse")

        def mse(mask):
            pred = np.asarray(probe.predict(X[mask]), dtype=float).ravel()
            return float(np.mean((y[mask] - pred) ** 2))

    batches = sorted(int(b) for b in np.unique(batch))
    ref_b = [b for b in batches if b < int(n_ref)]
    trail_b = [b for b in batches if b >= int(n_ref)]
    ref_scores = [mse(batch == b) for b in ref_b]
    e_ref = float(np.mean(ref_scores))
    scores = np.asarray([mse(batch == b) for b in trail_b], dtype=float)
    T = scores - e_ref
    hist_T = np.asarray(ref_scores, dtype=float) - e_ref
    p_trail = []
    pool = list(hist_T)
    for tval in T:
        pool_arr = np.asarray(pool, dtype=float)
        p_trail.append(float(np.sum(tval <= pool_arr) + 1) / float(len(pool) + 1))
        pool.append(float(tval))
    return {
        "e_ref": e_ref,
        "ref_batches": ref_b,
        "trail_batches": trail_b,
        "scores": scores,
        "T": T,
        "p": np.asarray(p_trail, dtype=float),
        "binary": binary,
    }


def last_two_on_rates(X, y, batch, *, gate: float, seed: int):
    batch = np.asarray(batch, dtype=int)
    y = np.asarray(y, dtype=float).ravel()
    X = np.asarray(X, dtype=float)
    batches = sorted(int(b) for b in np.unique(batch))
    e_prev = None
    hist = []
    for t in batches[1:]:
        prev = batch == (t - 1)
        cur = batch == t
        probe = fit_online_rf(X[prev], y[prev], seed=int(seed) + t, task="mse")
        e_now = probe_err(probe, X[cur], y[cur], task="mse")
        e_fl = error_floor("mse", int(cur.sum()))
        fired = hop_fires(e_now, e_prev, gate=gate, e_floor=e_fl)
        ratio = 1.0 if e_prev is None else shift_ratio(e_now, e_prev, e_floor=e_fl)
        hist.append(
            {
                "abs_batch": int(t),
                "fired": bool(fired),
                "e_now": float(e_now),
                "e_prev": None if e_prev is None else float(e_prev),
                "ratio": float(ratio),
            }
        )
        e_prev = e_now
    return hist


def write_xy(path: Path, rows: list[dict], extra: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["y", "batch", *extra, *ALL_COLS]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(
                {
                    "y": f"{float(r['y']):.6g}",
                    "batch": int(r["batch"]),
                    **{k: int(r[k]) if k in extra else f"{float(r[k]):.6g}" for k in extra},
                    **{c: f"{float(r[c]):.6g}" for c in ALL_COLS},
                }
            )


def plot_pipeline(path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.6, 3.4))
    ax.set_xlim(0, 11)
    ax.set_ylim(0, 3.2)
    ax.axis("off")
    boxes = [
        (0.15, 1.15, 1.7, 1.1, "user_id\n→ list of graphs"),
        (2.15, 1.15, 1.9, 1.1, "aggregate\nthe list"),
        (4.35, 1.85, 1.85, 0.85, "GraphSAGE-mean\n(frozen)"),
        (4.35, 0.55, 1.85, 0.85, "hand features\n7-d graph X"),
        (6.55, 1.15, 1.6, 1.1, "concat\nX = SAGE ⊕ hand"),
        (8.45, 1.15, 2.25, 1.1, "OnlineRFPerm\nT = MSE − E_ref"),
    ]
    for x, y, w, h, txt in boxes:
        ax.add_patch(
            FancyBboxPatch(
                (x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.12",
                facecolor="#fffdf8", edgecolor="#1f4e79", linewidth=1.4,
            )
        )
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=8.5, color="#1b1f24")
    arrows = [
        ((1.85, 1.7), (2.15, 1.7)),
        ((4.05, 1.7), (4.35, 2.27)),
        ((4.05, 1.7), (4.35, 0.97)),
        ((6.2, 2.27), (6.55, 1.85)),
        ((6.2, 0.97), (6.55, 1.55)),
        ((8.15, 1.7), (8.45, 1.7)),
    ]
    for (x0, y0), (x1, y1) in arrows:
        ax.add_patch(
            FancyArrowPatch(
                (x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=10,
                color="#1f4e79", lw=1.2, shrinkA=0, shrinkB=0,
            )
        )
    ax.set_title("Graph continuity → OnlineRFPerm (this facet only)", loc="left", fontsize=11, color="#1f4e79")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_y_rate(native_rows: list[dict], hop_rows: list[dict], path: Path) -> None:
    def series(rows):
        by = defaultdict(list)
        for r in rows:
            by[int(r["batch"])].append(float(r["y"]))
        xs = sorted(by)
        return xs, [float(np.mean(by[b])) for b in xs]

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.2), sharey=True)
    for ax, rows, title, hop in (
        (axes[0], native_rows, "pack-usable rate · quiet", False),
        (axes[1], hop_rows, "pack-usable rate · community hop", True),
    ):
        xs, ys = series(rows)
        if hop:
            ax.axvline(CUT_BATCH - 0.5, color="#9b2c2c", ls="--", lw=1.0)
        ax.plot(xs, ys, marker="o", color="#1f4e79")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("batch")
        ax.set_ylabel(r"$Y$ (pack usable)")
        ax.set_ylim(0.2, 0.75)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_T_and_p(native: dict, hop: dict, path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 6.4), sharex="col")
    specs = [
        (axes[0, 0], native, "T_t", "quiet"),
        (axes[0, 1], hop, "T_t", "community hop"),
        (axes[1, 0], native, "p", "quiet"),
        (axes[1, 1], hop, "p", "community hop"),
    ]
    for ax, rec, which, title in specs:
        xs = rec["trail_batches"]
        if which == "T_t":
            ys = rec["T"]
            ax.axhline(0.0, color="#888", lw=0.8)
            ax.set_ylabel(r"$T_t=\mathrm{MSE}_t-E_{\mathrm{ref}}$")
        else:
            ys = rec["p"]
            ax.axhline(ALPHA, color="#9b2c2c", ls=":", lw=1.0, label=r"$\alpha=0.05$")
            ax.set_ylabel("rank $p_t$")
            ax.set_ylim(-0.05, 1.05)
        if rec.get("regime") == "hop":
            ax.axvline(CUT_BATCH - 0.5, color="#9b2c2c", ls="--", lw=1.0)
        ax.plot(xs, ys, marker="o", color="#1f4e79")
        ax.set_title(f"{which} · {title}", fontsize=10)
        ax.set_xlabel("batch")
        if which == "p":
            ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def summarize(name: str, user_rows: list[dict]) -> dict:
    X, y, batch = matrix(user_rows)
    frozen = freeze_rf_mse(X, y, batch, n_ref=N_REF, seed=SEED)
    hops = last_two_on_rates(X, y, batch, gate=GATE, seed=SEED)
    cut = next((h for h in hops if h["abs_batch"] == CUT_BATCH), None)
    y_pre = [r["y"] for r in user_rows if int(r["batch"]) < CUT_BATCH]
    y_post = [r["y"] for r in user_rows if int(r["batch"]) >= CUT_BATCH]
    p = frozen["p"]
    n_p_fire = int(np.sum(p < ALPHA)) if len(p) else 0
    first = None
    for b, pv in zip(frozen["trail_batches"], p):
        if pv < ALPHA:
            first = int(b)
            break
    return {
        "name": name,
        "regime": "hop" if "hop" in name or name == "community" else "quiet",
        "n": len(user_rows),
        "n_users": N_USERS,
        "x_dim": X_DIM,
        "sage_dim": SAGE_DIM,
        "hand_dim": HAND_DIM,
        "y_pre": float(np.mean(y_pre)) if y_pre else None,
        "y_post": float(np.mean(y_post)) if y_post else None,
        "e_ref": frozen["e_ref"],
        "mean_T": float(np.mean(frozen["T"])) if len(frozen["T"]) else None,
        "n_rank_p_lt_alpha": n_p_fire,
        "first_p_fire": first,
        "hop_fire_at_cut": None if cut is None else bool(cut["fired"]),
        "n_hop_fires": int(sum(1 for h in hops if h["fired"])),
        "trail_batches": frozen["trail_batches"],
        "T": [float(v) for v in frozen["T"]],
        "p": [float(v) for v in frozen["p"]],
        "hops": hops,
        "cut_hop": cut,
    }


def render_report(native: dict, hop: dict) -> str:
    def hop_row(h):
        prev = "—" if h["e_prev"] is None else f"{h['e_prev']:.3f}"
        return f"| {h['abs_batch']} | {'yes' if h['fired'] else 'no'} | {prev} | {h['e_now']:.3f} | {h['ratio']:.3f} |"

    lines = [
        "# Graph continuity → OnlineRFPerm",
        "",
        f"X = GraphSAGE-mean ({SAGE_DIM}-d, frozen) ⊕ hand graph features ({HAND_DIM}-d). "
        f"user_id → list in the window → mean-pool. Y = pack-usable rate.",
        "",
        "OnlineRFPerm Algorithm 1: frozen RF, $T_t=\\mathrm{MSE}_t-E_{\\mathrm{ref}}$, rank $p_t$. "
        "Last-two `hop_fires` is the adjacent-window gate.",
        "",
        "| Regime | y pre | y post | mean T | rank-p < 0.05 | hop@cut | n hop |",
        "|---|---:|---:|---:|---:|---|---:|",
        (
            f"| quiet local pack | {native['y_pre']:.3f} | {native['y_post']:.3f} | "
            f"{native['mean_T']:.3f} | {native['n_rank_p_lt_alpha']} | "
            f"{'yes' if native['hop_fire_at_cut'] else 'no'} | {native['n_hop_fires']} |"
        ),
        (
            f"| community hop | {hop['y_pre']:.3f} | {hop['y_post']:.3f} | "
            f"{hop['mean_T']:.3f} | {hop['n_rank_p_lt_alpha']} | "
            f"{'yes' if hop['hop_fire_at_cut'] else 'no'} | {hop['n_hop_fires']} |"
        ),
        "",
        "## Last-two around the cut (community hop)",
        "",
        "| abs batch | fire | e_prev | e_now | ratio |",
        "|---:|---|---:|---:|---:|",
    ]
    for h in hop["hops"]:
        if 2 <= h["abs_batch"] <= 6:
            lines.append(hop_row(h))
    lines += [
        "",
        "Hotpot file order is not wall-clock time. Figures: `pipeline.png`, `y_rate.png`, `T_and_rank_p.png`. "
        "Rebuild: `PYTHONPATH=. python3 scripts/prototype_graph_continuity.py`.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    print("load HotpotQA distractor validation")
    questions, supporting, contexts = load_hotpot(N_ROWS)
    native_q = build_query_rows(
        questions, supporting, contexts, hop=False, rng=np.random.default_rng(SEED)
    )
    hop_q = build_query_rows(
        questions, supporting, contexts, hop=True, rng=np.random.default_rng(SEED)
    )
    native_u = aggregate_user_lists(native_q)
    hop_u = aggregate_user_lists(hop_q)
    write_xy(OUT / "xy_user_window.csv", native_u, ["user_id", "n_list"])
    write_xy(OUT / "xy_user_window_hop.csv", hop_u, ["user_id", "n_list"])
    native_sum = summarize("quiet", native_u)
    hop_sum = summarize("community", hop_u)
    native_sum["regime"] = "quiet"
    hop_sum["regime"] = "hop"
    plot_pipeline(OUT / "pipeline.png")
    plot_T_and_p(native_sum, hop_sum, OUT / "T_and_rank_p.png")
    plot_y_rate(native_u, hop_u, OUT / "y_rate.png")
    report = render_report(native_sum, hop_sum)
    (OUT / "REPORT.md").write_text(report, encoding="utf-8")
    slim_n = {k: native_sum[k] for k in native_sum if k not in {"T", "p", "hops", "cut_hop", "trail_batches"}}
    slim_h = {k: hop_sum[k] for k in hop_sum if k not in {"T", "p", "hops", "cut_hop", "trail_batches"}}
    manifest = {
        "x": "concat(GraphSAGE-mean, hand graph features)",
        "sage_dim": SAGE_DIM,
        "hand_dim": HAND_DIM,
        "n_users": N_USERS,
        "y": "pack-usable rate on the user list",
        "onlinerfperm": "frozen RF, T=MSE-E_ref, rank p; last-two hop_fires",
        "note": "Hotpot has no real user_id; query index % N_USERS prototypes the log.",
        "quiet": slim_n,
        "community_hop": slim_h,
    }
    (OUT / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(report)
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
