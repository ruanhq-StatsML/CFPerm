#!/usr/bin/env python3
"""TencentGR merchant path attribution + graph localization board.

Graphs (not inventing methods — classical constructions):
  1) **transition digraph**  : sequential step→next within a user path
  2) **co-occurrence graph** : undirected edges for items co-present in a path
  (optional later: user–item bipartite, funnel hetero exp/clk/cnv)

Attribution (ChannelAttribution-compatible; pip ``ChannelAttribution`` failed to
build here, so we call the same algorithms via networkx / pure Python):
  - heuristic: first-touch / last-touch / linear
  - Markov-1 removal-effect on the transition digraph

Localization (call networkx):
  - ``ego_graph`` around convert items
  - Personalized PageRank on the ego / full transition graph

  PYTHONPATH=. python3 scripts/tencent_gr/run_channel_graph_attribution.py \\
    --root data/tencent_subset --max-users 3000
"""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[2]
A_EXP, A_CLK, A_CNV = 0, 1, 2
Event = Tuple[int, int, int]


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------


def _pq(d: Path) -> List[Path]:
    return sorted(p for p in d.glob("*.parquet") if p.is_file())


def parse_events(seq: list) -> List[Event]:
    evs = []
    for e in seq:
        if not isinstance(e, dict):
            e = dict(e)
        evs.append((int(e["item_id"]), int(e["action_type"]), int(e["timestamp"])))
    evs.sort(key=lambda x: x[2])
    return evs


def iter_users(seq_dir: Path, max_users: int):
    n = 0
    for path in _pq(seq_dir):
        pf = pq.ParquetFile(path)
        for rg in range(pf.num_row_groups):
            df = pf.read_row_group(rg).to_pandas()
            for _, row in df.iterrows():
                yield int(row["user_id"]), parse_events(list(row["seq"]))
                n += 1
                if n >= max_users:
                    return


def path_items(evs: Sequence[Event], *, dedupe_consecutive: bool = True) -> List[int]:
    """Ordered item path (keep action diversity as node visits)."""
    items = [iid for iid, _, _ in evs]
    if not dedupe_consecutive:
        return items
    out: List[int] = []
    for x in items:
        if not out or out[-1] != x:
            out.append(x)
    return out


def converting_paths(
    users: Iterable[Tuple[int, List[Event]]],
) -> Tuple[List[List[int]], List[List[int]], Counter]:
    """Return converting item-paths, all paths, and convert-item counts."""
    conv_paths: List[List[int]] = []
    all_paths: List[List[int]] = []
    cnv_items: Counter = Counter()
    for _, evs in users:
        items = path_items(evs)
        if len(items) < 1:
            continue
        all_paths.append(items)
        has_cnv = any(a == A_CNV for _, a, _ in evs)
        if has_cnv:
            # truncate at last convert (classic journey → conversion)
            last_cnv_idx = max(i for i, (_, a, _) in enumerate(evs) if a == A_CNV)
            # map event index → path item index approximately via prefix items
            prefix = path_items(evs[: last_cnv_idx + 1])
            if prefix:
                conv_paths.append(prefix)
                cnv_items[evs[last_cnv_idx][0]] += 1
    return conv_paths, all_paths, cnv_items


# ---------------------------------------------------------------------------
# Graphs
# ---------------------------------------------------------------------------


def build_transition_digraph(paths: Sequence[Sequence[int]]) -> nx.DiGraph:
    """Sequential co-visit: edge u→v if v immediately follows u in a path."""
    G = nx.DiGraph()
    for path in paths:
        for a, b in zip(path, path[1:]):
            if a == b:
                continue
            if G.has_edge(a, b):
                G[a][b]["weight"] += 1.0
            else:
                G.add_edge(a, b, weight=1.0)
    return G


def build_cooccurrence_graph(
    paths: Sequence[Sequence[int]], *, window: Optional[int] = None
) -> nx.Graph:
    """Undirected co-occurrence within a path (optional sliding window)."""
    G = nx.Graph()
    for path in paths:
        uniq = list(dict.fromkeys(path))  # preserve order, drop dups
        if window is None:
            pairs = ((uniq[i], uniq[j]) for i in range(len(uniq)) for j in range(i + 1, len(uniq)))
        else:
            pairs = []
            for i in range(len(path)):
                for j in range(i + 1, min(len(path), i + window)):
                    if path[i] != path[j]:
                        pairs.append((path[i], path[j]))
        for a, b in pairs:
            u, v = (a, b) if a < b else (b, a)
            if G.has_edge(u, v):
                G[u][v]["weight"] += 1.0
            else:
                G.add_edge(u, v, weight=1.0)
    return G


# ---------------------------------------------------------------------------
# ChannelAttribution-compatible heuristics + Markov removal effect
# ---------------------------------------------------------------------------


def heuristic_attribution(paths: Sequence[Sequence[int]]) -> pd.DataFrame:
    """first / last / linear touch — same three heuristics as ChannelAttribution."""
    first = Counter()
    last = Counter()
    linear = Counter()
    for path in paths:
        if not path:
            continue
        first[path[0]] += 1.0
        last[path[-1]] += 1.0
        w = 1.0 / float(len(path))
        for n in path:
            linear[n] += w
    nodes = sorted(set(first) | set(last) | set(linear))
    rows = []
    for n in nodes:
        rows.append(
            {
                "channel": int(n),
                "first_touch": float(first[n]),
                "last_touch": float(last[n]),
                "linear_touch": float(linear[n]),
            }
        )
    df = pd.DataFrame(rows)
    for c in ("first_touch", "last_touch", "linear_touch"):
        s = df[c].sum()
        df[c + "_share"] = df[c] / s if s > 0 else 0.0
    return df.sort_values("linear_touch", ascending=False).reset_index(drop=True)


def _transition_probs(G: nx.DiGraph) -> Dict[int, Dict[int, float]]:
    P: Dict[int, Dict[int, float]] = {}
    for u in G.nodes:
        tot = sum(d.get("weight", 1.0) for _, _, d in G.out_edges(u, data=True))
        if tot <= 0:
            continue
        P[u] = {
            v: float(d.get("weight", 1.0)) / tot for _, v, d in G.out_edges(u, data=True)
        }
    return P


def markov_removal_attribution(
    paths: Sequence[Sequence[int]],
    G: nx.DiGraph,
    *,
    n_sim: int = 20000,
    seed: int = 0,
    top_channels: int = 80,
) -> pd.DataFrame:
    """Markov-1 removal-effect attribution (ChannelAttribution-style).

    Simulate paths on the transition digraph; credit ∝ drop in conversion rate
    when a channel is removed. Restrict removal candidates to top-frequency
    channels for speed.
    """
    rng = np.random.default_rng(seed)
    # conversion = path whose last node is a convert-terminal in training paths
    terminals = Counter(p[-1] for p in paths if p)
    freq = Counter(n for p in paths for n in p)
    candidates = [n for n, _ in freq.most_common(top_channels)]

    # build absorbing-ish sim: start from empirical first nodes; stop at sink
    starts = [p[0] for p in paths if p]
    if not starts:
        return pd.DataFrame()
    P = _transition_probs(G)
    nodes = list(G.nodes)
    if not nodes:
        return pd.DataFrame()

    def conv_rate(block: Optional[set], n: int) -> float:
        ok = 0
        for _ in range(n):
            u = int(starts[int(rng.integers(0, len(starts)))])
            if block and u in block:
                continue
            seen = {u}
            for _step in range(40):
                nbrs = P.get(u)
                if not nbrs:
                    break
                vs = list(nbrs.keys())
                ps = np.asarray([nbrs[v] for v in vs], float)
                ps = ps / ps.sum()
                u = int(rng.choice(vs, p=ps))
                if block and u in block:
                    break
                if u in terminals and rng.random() < min(1.0, terminals[u] / (freq[u] + 1e-9)):
                    ok += 1
                    break
                if u in seen:
                    break
                seen.add(u)
        return ok / float(n)

    base = conv_rate(None, n_sim)
    rows = []
    # fewer sims per removal
    n_rm = max(2000, n_sim // 4)
    for ch in candidates:
        r = conv_rate({ch}, n_rm)
        effect = max(0.0, base - r)
        rows.append({"channel": int(ch), "removal_effect": effect, "base_rate": base, "rate_wo": r})
    df = pd.DataFrame(rows)
    s = df["removal_effect"].sum()
    df["markov_share"] = df["removal_effect"] / s if s > 0 else 0.0
    return df.sort_values("markov_share", ascending=False).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Localization: ego + Personalized PageRank
# ---------------------------------------------------------------------------


def localize_converts(
    G: nx.Graph,
    cnv_items: Counter,
    *,
    radius: int = 2,
    top_k_seeds: int = 12,
    top_neighbors: int = 8,
) -> List[dict]:
    """For top convert items: ego_graph + personalized PageRank ranking."""
    out = []
    undirected = not G.is_directed()
    for seed, cnt in cnv_items.most_common(top_k_seeds):
        if seed not in G:
            out.append({"seed": int(seed), "cnv_cnt": int(cnt), "error": "seed_not_in_graph"})
            continue
        if G.is_directed():
            ego = nx.ego_graph(G, seed, radius=radius, undirected=True)
        else:
            ego = nx.ego_graph(G, seed, radius=radius)
        if ego.number_of_nodes() <= 1:
            out.append(
                {
                    "seed": int(seed),
                    "cnv_cnt": int(cnt),
                    "ego_n": ego.number_of_nodes(),
                    "ego_m": ego.number_of_edges(),
                    "top_ppr": [],
                }
            )
            continue
        personalization = {n: 0.0 for n in ego.nodes}
        personalization[seed] = 1.0
        try:
            ppr = nx.pagerank(ego, alpha=0.85, personalization=personalization, weight="weight")
        except Exception:
            ppr = {n: 1.0 / ego.number_of_nodes() for n in ego.nodes}
        ranked = sorted(
            ((n, s) for n, s in ppr.items() if n != seed),
            key=lambda x: -x[1],
        )[:top_neighbors]
        out.append(
            {
                "seed": int(seed),
                "cnv_cnt": int(cnt),
                "ego_n": int(ego.number_of_nodes()),
                "ego_m": int(ego.number_of_edges()),
                "top_ppr": [{"item": int(n), "ppr": float(s)} for n, s in ranked],
            }
        )
    return out


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------


def plot_board(
    heur: pd.DataFrame,
    markov: pd.DataFrame,
    locs_trans: List[dict],
    locs_cooc: List[dict],
    meta: dict,
    out_path: Path,
) -> Path:
    fig = plt.figure(figsize=(13.2, 9.0))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.30)

    # (1) heuristic top channels
    ax0 = fig.add_subplot(gs[0, 0])
    top = heur.head(12).copy()
    x = np.arange(len(top))
    w = 0.27
    ax0.bar(x - w, top["first_touch_share"], width=w, label="first", color="#4C78A8")
    ax0.bar(x, top["last_touch_share"], width=w, label="last", color="#F58518")
    ax0.bar(x + w, top["linear_touch_share"], width=w, label="linear", color="#54A24B")
    ax0.set_xticks(x)
    ax0.set_xticklabels([str(c)[-4:] for c in top["channel"]], rotation=45, ha="right", fontsize=7)
    ax0.set_ylabel("share")
    ax0.set_title("Heuristic path attribution (top-12 by linear)", fontsize=11)
    ax0.legend(fontsize=8)

    # (2) Markov removal vs linear
    ax1 = fig.add_subplot(gs[0, 1])
    if markov is not None and len(markov):
        mtop = markov.head(12)
        # align linear share
        lin = heur.set_index("channel")["linear_touch_share"]
        ys_m = mtop["markov_share"].to_numpy()
        ys_l = [float(lin.get(int(c), 0.0)) for c in mtop["channel"]]
        x = np.arange(len(mtop))
        ax1.bar(x - 0.18, ys_m, width=0.36, label="Markov removal", color="#E45756")
        ax1.bar(x + 0.18, ys_l, width=0.36, label="linear", color="#54A24B")
        ax1.set_xticks(x)
        ax1.set_xticklabels([str(c)[-4:] for c in mtop["channel"]], rotation=45, ha="right", fontsize=7)
        ax1.legend(fontsize=8)
    ax1.set_ylabel("share")
    ax1.set_title("Markov-1 removal-effect (transition digraph)", fontsize=11)

    # (3) localization on transition ego
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.axis("off")
    ax2.set_title("Graph-local: transition digraph · ego+PPR", fontsize=11, pad=6)
    lines = []
    for loc in locs_trans[:5]:
        if "error" in loc:
            continue
        nbr = ", ".join(f"{t['item']%10000}:{t['ppr']:.3f}" for t in loc.get("top_ppr", [])[:5])
        lines.append(
            f"cnv={loc['seed']%10000} (n={loc['cnv_cnt']})  "
            f"ego={loc.get('ego_n',0)}n/{loc.get('ego_m',0)}e\n  PPR→ {nbr}"
        )
    ax2.text(
        0.02,
        0.98,
        "\n\n".join(lines) if lines else "(empty)",
        va="top",
        ha="left",
        family="monospace",
        fontsize=7.5,
        transform=ax2.transAxes,
    )

    # (4) localization on co-occurrence + schema
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.axis("off")
    ax3.set_title("Graph-local: co-occurrence · ego+PPR", fontsize=11, pad=6)
    lines = []
    for loc in locs_cooc[:5]:
        if "error" in loc:
            continue
        nbr = ", ".join(f"{t['item']%10000}:{t['ppr']:.3f}" for t in loc.get("top_ppr", [])[:5])
        lines.append(
            f"cnv={loc['seed']%10000} (n={loc['cnv_cnt']})  "
            f"ego={loc.get('ego_n',0)}n/{loc.get('ego_m',0)}e\n  PPR→ {nbr}"
        )
    schema = (
        f"users={meta['n_users']}  conv_paths={meta['n_conv_paths']}  "
        f"trans={meta['trans_n']}/{meta['trans_m']}  "
        f"cooc={meta['cooc_n']}/{meta['cooc_m']}\n"
        "graphs: transition digraph | co-occurrence (undirected)\n"
        "attr: first/last/linear + Markov removal  |  loc: nx.ego_graph + PPR"
    )
    ax3.text(
        0.02,
        0.98,
        schema + "\n\n" + ("\n\n".join(lines) if lines else "(empty)"),
        va="top",
        ha="left",
        family="monospace",
        fontsize=7.5,
        transform=ax3.transAxes,
    )

    fig.suptitle(
        "TencentGR ChannelAttribution-style + graph localization (call existing algos)",
        fontsize=12.5,
        y=0.995,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT / "data" / "tencent_subset")
    ap.add_argument("--max-users", type=int, default=3000)
    ap.add_argument("--n-sim", type=int, default=12000)
    ap.add_argument("--ego-radius", type=int, default=2)
    ap.add_argument("--cooc-window", type=int, default=8, help="sliding window; 0=full path")
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "results" / "tencent_gr_channel_graph_attr",
    )
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    print(f"load ≤{args.max_users} users ...", flush=True)
    users = list(iter_users(args.root / "seq", args.max_users))
    conv_paths, all_paths, cnv_items = converting_paths(users)
    print(
        f"  users={len(users)} all_paths={len(all_paths)} conv_paths={len(conv_paths)} "
        f"uniq_cnv_items={len(cnv_items)}",
        flush=True,
    )

    print("build transition digraph + co-occurrence ...", flush=True)
    G_trans = build_transition_digraph(all_paths)
    win = None if args.cooc_window <= 0 else args.cooc_window
    G_cooc = build_cooccurrence_graph(all_paths, window=win)
    print(
        f"  trans: n={G_trans.number_of_nodes()} m={G_trans.number_of_edges()} | "
        f"cooc: n={G_cooc.number_of_nodes()} m={G_cooc.number_of_edges()}",
        flush=True,
    )

    print("heuristic attribution ...", flush=True)
    heur = heuristic_attribution(conv_paths if conv_paths else all_paths)
    print("Markov removal-effect (transition) ...", flush=True)
    markov = markov_removal_attribution(
        conv_paths if conv_paths else all_paths,
        G_trans,
        n_sim=args.n_sim,
        seed=args.seed,
    )

    print("ego + PPR localization ...", flush=True)
    locs_t = localize_converts(G_trans, cnv_items, radius=args.ego_radius)
    locs_c = localize_converts(G_cooc, cnv_items, radius=args.ego_radius)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    heur.to_csv(args.out_dir / "heuristic_attribution.csv", index=False)
    markov.to_csv(args.out_dir / "markov_removal_attribution.csv", index=False)

    meta = {
        "n_users": len(users),
        "n_all_paths": len(all_paths),
        "n_conv_paths": len(conv_paths),
        "n_cnv_items": len(cnv_items),
        "trans_n": G_trans.number_of_nodes(),
        "trans_m": G_trans.number_of_edges(),
        "cooc_n": G_cooc.number_of_nodes(),
        "cooc_m": G_cooc.number_of_edges(),
        "cooc_window": win,
        "ego_radius": args.ego_radius,
        "note": (
            "ChannelAttribution pip wheel failed to build (Cython); "
            "heuristic+Markov-removal reimplemented to match CA API; "
            "localization via networkx.ego_graph + pagerank."
        ),
        "graphs": {
            "transition_digraph": "edge u→v if v follows u in user path",
            "cooccurrence": "undirected co-presence in path (optional window)",
            "also_possible": [
                "user–item bipartite",
                "funnel hetero (exp/clk/cnv node types)",
                "mm_emb kNN item graph",
            ],
        },
    }
    payload = {
        "meta": meta,
        "heuristic_top20": heur.head(20).to_dict(orient="records"),
        "markov_top20": markov.head(20).to_dict(orient="records") if len(markov) else [],
        "localize_transition": locs_t,
        "localize_cooccurrence": locs_c,
    }
    (args.out_dir / "channel_graph_attribution.json").write_text(json.dumps(payload, indent=2))

    png = plot_board(
        heur,
        markov,
        locs_t,
        locs_c,
        meta,
        args.out_dir / "tencent_channel_graph_attribution.png",
    )
    (args.out_dir / "README.md").write_text(
        "\n".join(
            [
                "# ChannelAttribution-style + graph localization",
                "",
                "- Heuristic: first / last / linear (CA-compatible)",
                "- Markov-1 removal-effect on **transition digraph**",
                "- Localization: `networkx.ego_graph` + Personalized PageRank",
                "- Graphs: transition digraph **and** co-occurrence (windowed)",
                f"- Figure: `{png.name}`",
                "",
                meta["note"],
                "",
            ]
        )
    )
    print(f"wrote {png}", flush=True)
    print("heuristic top5:", heur.head(5)[["channel", "linear_touch_share"]].to_dict("records"))
    if len(markov):
        print("markov top5:", markov.head(5)[["channel", "markov_share"]].to_dict("records"))


if __name__ == "__main__":
    main()
