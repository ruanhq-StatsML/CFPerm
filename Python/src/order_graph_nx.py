"""Order-graph localization on **networkx** (not PyG / DGL / GraphRAG / igraph).

Package: networkx (3.6.x in this env).

The changing subset is a **Fast Subset Scan** of a node potential, not a
community. Rank entities by own-ref φ (or mass-weighted φ·n) and take a
prefix:

    level set   Ŝ = { i : ψ_i ≥ τ },  τ = max(floor, α · max ψ)
    coverage    smallest prefix of the ranking with Σψ ≥ β · Σψ

φ_i is the bundled own-ref score (MMD, CMean_X, CMean_Y, PO). The multi-layer
graph is incidence only: which orders hang on which merchant / user. Scan
each layer independently, lift Ŝ to the order grain, then Jaccard. Do not
mix merchant-MMD and user-MMD in one simplex. Do not run Louvain.

Louvain is an optional contrast, not the cut:

  structural  co-order + kNN(X). Modularity Q finds dense subgraphs
              ("who shares users"). That is **not** "who shifted".
  bundled     same-direction affinity of the multi-metric bundle.
              Louvain here groups loud nodes; still not required.

Frozen GraphSAGE-mean is a readout of the structural graph, not the cut.
Y never enters node features. Y may score the bundle (monitoring outcome).
"""
from __future__ import annotations

from collections import defaultdict
from typing import Mapping, Sequence

import networkx as nx
import numpy as np

GRAPH_PACKAGE = "networkx"
GRAPH_PACKAGE_VERSION = nx.__version__
ENCODER = "graphsage-mean-frozen"
CUT_STRUCTURAL = "networkx.community.louvain_communities/structural"
CUT_BUNDLED = "networkx.community.louvain_communities/bundled-shift"
CUT_LEVEL_SET = "level_set/{phi>=tau}"
CUT_SCAN_COVERAGE = "subset_scan/{coverage}"
CUT_SCAN_MASS = "subset_scan/{mass_coverage}"
SCAN_COVERAGE = 0.80
SAGE_LAYERS = 2
BUNDLE_KEYS = ("mmd", "cmean_x", "cmean_y", "po")
METRIC_FLOORS = {"mmd": 0.02, "cmean_x": 0.15, "cmean_y": 0.05, "po": 1e-3}


def _as_2d(X) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        return X.reshape(-1, 1)
    return X


def merchant_bag_x(window: Mapping) -> tuple[np.ndarray, np.ndarray]:
    """One row per merchant: native merchant X plus mean order X. No Y, no region."""
    mids = np.asarray(window["merchant_id"])
    Xo = _as_2d(window["X_order"])
    Xm = _as_2d(window["X_merchant"])
    uniq = np.unique(mids)
    rows = []
    for m in uniq:
        mask = mids == m
        rows.append(np.concatenate([Xm[mask][0], Xo[mask].mean(axis=0)]))
    return uniq.astype(int), np.vstack(rows)


def build_order_graph(window: Mapping, *, k_nn: int = 4) -> nx.Graph:
    """networkx Graph of this window. Node attrs: ntype, x. No Y."""
    G = nx.Graph()
    G.graph["package"] = GRAPH_PACKAGE
    G.graph["package_version"] = GRAPH_PACKAGE_VERSION
    mids = np.asarray(window["merchant_id"]).astype(int)
    uids = np.asarray(window["user_id"]).astype(int)
    mer_ids, mer_x = merchant_bag_x(window)
    for m, x in zip(mer_ids, mer_x):
        G.add_node(("m", int(m)), ntype="merchant", x=np.asarray(x, dtype=float))
    u_x = _as_2d(window["X_user"])
    for u in np.unique(uids):
        mask = uids == u
        G.add_node(("u", int(u)), ntype="user", x=u_x[mask][0].astype(float))
    pair_w: dict[tuple[int, int], int] = defaultdict(int)
    for m, u in zip(mids, uids):
        pair_w[(int(m), int(u))] += 1
    for (m, u), w in pair_w.items():
        G.add_edge(("m", m), ("u", u), weight=float(w), etype="co_order")
    _add_merchant_knn(G, mer_ids, mer_x, k_nn=k_nn)
    return G


def _add_merchant_knn(G: nx.Graph, mer_ids, mer_x, k_nn: int) -> None:
    mer_x = _as_2d(mer_x)
    n = len(mer_ids)
    if n <= 1:
        return
    k = max(1, min(int(k_nn), n - 1))
    # z-score so amount and gmv share a scale; still no Y.
    sd = mer_x.std(axis=0)
    sd = np.where(sd < 1e-8, 1.0, sd)
    Z = (mer_x - mer_x.mean(axis=0)) / sd
    d2 = ((Z[:, None, :] - Z[None, :, :]) ** 2).sum(axis=2)
    np.fill_diagonal(d2, np.inf)
    for i, m in enumerate(mer_ids):
        nn = np.argpartition(d2[i], kth=k - 1)[:k]
        for j in nn:
            a, b = ("m", int(m)), ("m", int(mer_ids[j]))
            w = float(1.0 / (1.0 + np.sqrt(d2[i, j])))
            if G.has_edge(a, b):
                G[a][b]["weight"] = max(float(G[a][b].get("weight", 0.0)), w)
            else:
                G.add_edge(a, b, weight=w, etype="knn_x")


def merchant_subgraph(G: nx.Graph) -> nx.Graph:
    nodes = [n for n, d in G.nodes(data=True) if d.get("ntype") == "merchant"]
    return G.subgraph(nodes).copy()


def louvain_merchant_cut(G: nx.Graph, *, seed: int = 2026, resolution: float = 1.0) -> dict[int, str]:
    """Merchant id → community label C0, C1, … via networkx Louvain."""
    H = merchant_subgraph(G)
    if H.number_of_nodes() == 0:
        return {}
    parts = nx.community.louvain_communities(H, weight="weight", seed=int(seed), resolution=float(resolution))
    parts = sorted(parts, key=lambda s: (-len(s), min(n[1] for n in s)))
    out = {}
    for i, comm in enumerate(parts):
        lab = f"C{i}"
        for n in comm:
            if n[0] == "m":
                out[int(n[1])] = lab
    return out


def graphsage_mean_frozen(G: nx.Graph, *, layers: int = SAGE_LAYERS) -> dict:
    """Frozen GraphSAGE-mean on whatever node feature dim is present.

    Nodes with a missing neighbor mean get a self-loop. No torch, no W.
    """
    nodes = list(G.nodes())
    if not nodes:
        return {}
    x0 = {n: np.asarray(G.nodes[n]["x"], dtype=float).ravel() for n in nodes}
    dim0 = max(v.size for v in x0.values())
    h = {}
    for n, v in x0.items():
        if v.size < dim0:
            v = np.pad(v, (0, dim0 - v.size))
        h[n] = v
    for _ in range(int(layers)):
        nxt = {}
        for n in nodes:
            nbrs = [h[m] for m in G.neighbors(n)]
            mean = h[n] if not nbrs else np.mean(np.stack(nbrs, axis=0), axis=0)
            nxt[n] = np.concatenate([h[n], mean])
        h = nxt
    return h


def sage_merchant_matrix(G: nx.Graph) -> tuple[np.ndarray, np.ndarray]:
    emb = graphsage_mean_frozen(merchant_subgraph(G))
    mids = sorted(n[1] for n in emb)
    X = np.vstack([emb[("m", m)] for m in mids]) if mids else np.zeros((0, 1))
    return np.asarray(mids, dtype=int), X


def order_community_labels(merchant_id, cut: Mapping[int, str], default: str = "C?") -> np.ndarray:
    merchant_id = np.asarray(merchant_id).astype(int)
    return np.asarray([cut.get(int(m), default) for m in merchant_id], dtype=object)


def order_level_set_labels(
    entity_ids,
    loud_ids: Sequence[int],
    *,
    loud_name: str = "loud",
    other_name: str = "other",
) -> np.ndarray:
    """Binary order labels from a node level-set. Not a community partition."""
    loud = {int(i) for i in loud_ids}
    ids = np.asarray(entity_ids).astype(int)
    return np.asarray([loud_name if int(e) in loud else other_name for e in ids], dtype=object)


def community_south_frac(cut: Mapping[int, str], region_by_merchant: Mapping[int, str]) -> dict[str, float]:
    buckets: dict[str, list[str]] = defaultdict(list)
    for m, lab in cut.items():
        buckets[lab].append(str(region_by_merchant.get(int(m), "north")))
    return {lab: float(np.mean([r == "south" for r in regs])) for lab, regs in buckets.items()}


def graph_localize(window: Mapping, *, k_nn: int = 4, seed: int = 2026) -> dict:
    """Build the networkx graph, Louvain-cut merchants, freeze SAGE-mean."""
    G = build_order_graph(window, k_nn=k_nn)
    cut = louvain_merchant_cut(G, seed=seed)
    labels = order_community_labels(window["merchant_id"], cut)
    mids, sage = sage_merchant_matrix(G)
    region_by_m = {}
    m_arr = np.asarray(window["merchant_id"]).astype(int)
    r_arr = np.asarray(window["region"])
    for m, r in zip(m_arr, r_arr):
        region_by_m.setdefault(int(m), str(r))
    return {
        "package": GRAPH_PACKAGE,
        "package_version": GRAPH_PACKAGE_VERSION,
        "encoder": ENCODER,
        "cut": CUT_STRUCTURAL,
        "graph": G,
        "n_nodes": int(G.number_of_nodes()),
        "n_edges": int(G.number_of_edges()),
        "merchant_cut": cut,
        "order_community": labels,
        "n_communities": int(len(set(cut.values()))),
        "south_frac": community_south_frac(cut, region_by_m),
        "sage_merchant_ids": mids,
        "sage": sage,
        "sage_dim": int(sage.shape[1]) if sage.size else 0,
    }


def bundle_vec(row: Mapping, *, signed: bool = True) -> np.ndarray:
    """4-vector. CMean_Y keeps sign so same-direction hops attract. PO missing = 0."""
    cy = float(row.get("cmean_y") or 0.0)
    return np.array(
        [
            max(float(row.get("mmd") or 0.0), 0.0),
            max(float(row.get("cmean_x") or 0.0), 0.0),
            cy if signed else abs(cy),
            max(float(row.get("po") or 0.0), 0.0),
        ],
        dtype=float,
    )


def bundle_phi(row: Mapping, scale: np.ndarray | None = None) -> float:
    """Loudness = L1 of scaled |v|. Quiet metrics stay 0."""
    v = bundle_vec(row, signed=True)
    if scale is None:
        scale = np.ones(4)
    z = v / np.maximum(scale, 1e-8)
    return float(np.abs(z).sum())


def _column_scale(vectors: Sequence[np.ndarray]) -> np.ndarray:
    floors = np.array([0.02, 0.15, 0.05, 1e-3], dtype=float)
    if not vectors:
        return floors
    V = np.vstack(vectors)
    sd = V.std(axis=0)
    return np.maximum(np.where(sd < 1e-8, 0.0, sd), floors)


def bundled_shift_graph(
    scores: Mapping[int, Mapping],
    *,
    min_sim: float = 0.15,
    phi_floor: float = 1.0,
) -> nx.Graph:
    """Merchant graph whose edges are same-direction bundled-shift affinity.

    Isolated quiet nodes stay isolated. Louvain on this graph groups
    merchants that are loud *together* — a proxy for the changing subset,
    not modularity of co-purchase.
    """
    G = nx.Graph()
    G.graph["package"] = GRAPH_PACKAGE
    G.graph["objective"] = "bundled-shift"
    mids = [int(m) for m in scores]
    vecs = [bundle_vec(scores[m]) for m in mids]
    scale = _column_scale(vecs)
    phis = {}
    for m, v in zip(mids, vecs):
        phi = bundle_phi(scores[m], scale=scale)
        phis[m] = phi
        G.add_node(int(m), ntype="merchant", phi=phi, v=v.tolist())
    tau = float(phi_floor)
    if phis:
        tau = max(tau, 0.30 * max(phis.values()))
    Z = np.vstack([bundle_vec(scores[m]) / scale for m in mids]) if mids else np.zeros((0, 4))
    for i, mi in enumerate(mids):
        if phis[mi] < tau:
            continue
        zi = Z[i]
        ni = float(np.linalg.norm(zi))
        if ni <= 1e-12:
            continue
        for j in range(i + 1, len(mids)):
            mj = mids[j]
            if phis[mj] < tau:
                continue
            zj = Z[j]
            nj = float(np.linalg.norm(zj))
            if nj <= 1e-12:
                continue
            sim = float(np.dot(zi, zj) / (ni * nj))
            if sim >= float(min_sim):
                G.add_edge(int(mi), int(mj), weight=float(sim), etype="bundle")
    return G


def louvain_int_cut(G: nx.Graph, *, seed: int = 2026, resolution: float = 1.0) -> dict[int, str]:
    """int node id → C0, C1, … . Empty / edgeless graph → one community per node."""
    if G.number_of_nodes() == 0:
        return {}
    parts = nx.community.louvain_communities(G, weight="weight", seed=int(seed), resolution=float(resolution))
    parts = sorted(parts, key=lambda s: (-len(s), min(int(n) for n in s)))
    out = {}
    for i, comm in enumerate(parts):
        lab = f"C{i}"
        for n in comm:
            out[int(n)] = lab
    return out


def community_phi(cut: Mapping[int, str], scores: Mapping[int, Mapping], scale=None) -> dict[str, float]:
    """Sum of node potentials in each community — which community looks shifted."""
    acc: dict[str, float] = defaultdict(float)
    vecs = [bundle_vec(scores[m]) for m in scores] if scale is None else None
    sc = _column_scale(vecs) if scale is None else scale
    for m, lab in cut.items():
        row = scores.get(int(m))
        if row is None:
            continue
        acc[str(lab)] += bundle_phi(row, scale=sc)
    return dict(acc)


def loud_community(cut: Mapping[int, str], scores: Mapping[int, Mapping]) -> str | None:
    phis = community_phi(cut, scores)
    if not phis:
        return None
    if max(phis.values()) <= 1e-12:
        return None
    return max(phis, key=phis.get)


def node_phis(scores: Mapping[int, Mapping]) -> dict[int, float]:
    vecs = [bundle_vec(scores[m]) for m in scores]
    scale = _column_scale(vecs)
    return {int(m): bundle_phi(scores[m], scale=scale) for m in scores}


def phi_tau(phis: Mapping[int, float], floor: float = 1.0, frac: float = 0.30) -> float:
    if not phis:
        return float(floor)
    return max(float(floor), float(frac) * max(phis.values()))


def node_mass(scores: Mapping[int, Mapping]) -> dict[int, float]:
    """Bag size per node. Missing n → 1 so mass weight does not drop the node."""
    out = {}
    for i, row in scores.items():
        n = (row or {}).get("n")
        out[int(i)] = float(n) if n is not None and float(n) > 0 else 1.0
    return out


def _scan_priority(
    scores: Mapping[int, Mapping],
    *,
    weight: str = "phi",
) -> tuple[dict[int, float], dict[int, float]]:
    """ψ used to rank nodes. weight='phi' is intensity; 'mass' is φ·n."""
    phis = node_phis(scores)
    if weight == "phi":
        return phis, dict(phis)
    if weight != "mass":
        raise ValueError(f"weight must be phi or mass, got {weight!r}")
    ns = node_mass(scores)
    psi = {int(i): float(phis[i]) * float(ns.get(i, 1.0)) for i in phis}
    return phis, psi


def subset_scan(
    scores: Mapping[int, Mapping],
    *,
    rule: str = "level_set",
    floor: float = 1.0,
    frac: float = 0.30,
    coverage: float = SCAN_COVERAGE,
    weight: str = "phi",
) -> dict:
    """Fast subset scan: rank nodes by ψ, take a prefix. Not a partition.

    rule='level_set'  → {i : ψ_i ≥ max(floor, frac · max ψ)}
    rule='coverage'   → shortest prefix of the ranking with Σψ ≥ coverage · Σψ,
                        after dropping nodes with φ < floor (quiet).

    Graph is not used. Weight 'mass' scores φ·n so a tiny noisy bag cannot
    outrank a large shifted bag. This is LTSS (Neill): an additive score's
    maximizer is a prefix of the priority ranking — not modularity.
    """
    if rule not in ("level_set", "coverage"):
        raise ValueError(f"rule must be level_set or coverage, got {rule!r}")
    phis, psi = _scan_priority(scores, weight=weight)
    ranked = sorted(psi, key=lambda i: (-float(psi[i]), int(i)))
    if rule == "level_set":
        tau = phi_tau(psi, floor=floor, frac=frac)
        loud = [int(i) for i in ranked if float(psi[i]) >= tau]
        cut_name = CUT_LEVEL_SET if weight == "phi" else CUT_SCAN_MASS
    else:
        # Quiet nodes (φ below floor) never enter the prefix.
        live = [int(i) for i in ranked if float(phis[i]) >= float(floor)]
        total = float(sum(psi[i] for i in live))
        loud = []
        acc = 0.0
        target = float(coverage) * total
        if total > 0.0:
            for i in live:
                loud.append(int(i))
                acc += float(psi[i])
                if acc >= target:
                    break
        tau = float(psi[loud[-1]]) if loud else float(floor)
        cut_name = CUT_SCAN_MASS if weight == "mass" else CUT_SCAN_COVERAGE
    return {
        "loud_ids": loud,
        "ranked": [int(i) for i in ranked],
        "tau": float(tau),
        "phis": {str(k): float(v) for k, v in phis.items()},
        "psi": {str(k): float(v) for k, v in psi.items()},
        "n_loud": int(len(loud)),
        "rule": rule,
        "weight": weight,
        "coverage": float(coverage) if rule == "coverage" else None,
        "cut_name": cut_name,
    }


def level_set_ids(
    scores: Mapping[int, Mapping],
    *,
    floor: float = 1.0,
    frac: float = 0.30,
) -> dict:
    """Changing subset = {i : φ_i ≥ τ}. Prefix of a subset scan, not a community.

    Graph is not used. τ = max(floor, frac * max φ).
    """
    out = subset_scan(scores, rule="level_set", floor=floor, frac=frac, weight="phi")
    return {
        "loud_ids": sorted(out["loud_ids"]),
        "tau": out["tau"],
        "phis": out["phis"],
        "n_loud": out["n_loud"],
    }


def metric_level_set_ids(
    scores: Mapping[int, Mapping],
    metric: str,
    *,
    floor: float | None = None,
    frac: float = 0.30,
) -> dict:
    """Per-metric level set {i : |metric_i| ≥ τ}. Diagnostic, not a second cut.

    Covariate should light MMD / CMean_X; concept should light CMean_Y / PO.
    Union of metric slices is allowed to be larger than the bundled Ŝ.
    """
    if metric not in BUNDLE_KEYS:
        raise ValueError(f"metric must be one of {BUNDLE_KEYS}, got {metric!r}")
    if floor is None:
        floor = float(METRIC_FLOORS[metric])
    vals: dict[int, float] = {}
    for m, row in scores.items():
        raw = row.get(metric)
        if raw is None:
            continue
        v = abs(float(raw)) if metric == "cmean_y" else max(float(raw), 0.0)
        vals[int(m)] = v
    if not vals:
        return {"metric": metric, "loud_ids": [], "tau": float(floor), "n_loud": 0}
    tau = max(float(floor), float(frac) * max(vals.values()))
    loud = sorted(int(i) for i, v in vals.items() if v >= tau)
    return {
        "metric": metric,
        "loud_ids": loud,
        "tau": float(tau),
        "n_loud": int(len(loud)),
        "vals": {str(k): float(v) for k, v in vals.items()},
    }


def y_in_graph_attrs(G: nx.Graph) -> bool:
    """Leakage check: Y must not be a node or edge attribute."""
    banned = {"y", "Y", "label", "target", "outcome"}
    for _, data in G.nodes(data=True):
        if banned.intersection(data.keys()):
            return True
        if any(str(v) in banned for v in data.values() if isinstance(v, str)):
            return True
    for _, _, data in G.edges(data=True):
        if banned.intersection(data.keys()):
            return True
    return False
