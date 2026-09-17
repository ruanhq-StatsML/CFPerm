"""Graph localization → FSDS grain-unify → two-layer subset attribution.

No unique causal decomp. Shares are ReLU localization proxies, not Shapley
and not CATE. T is the batch label. Y is the outcome, never a feature,
never the subset key.

Pipeline (leakage-safe):
  1. Score every native grain (order / merchant-node / user-node) with
     MMD² vs D_ref, PO-risk, and Conditional Mean — before any join.
  2. Changing subset on each layer: level set {i : φ_i ≥ τ} vs own-ref.
     Graph is incidence only. Louvain is an optional contrast.
  3. FSDS: rank features inside each grain; select the loud ones.
  4. Unify selected features onto one grain (order, or merchant).
  5. Two-layer on that grain: loud vs other, then LOGO on feature blocks.
  6. Characterize the loud subset with the same three readouts.

D_ref freezes σ, quantile bins, entity profiles used as serving features,
and the subset vocabulary. New-batch Y is only an outcome for PO / CMean_Y.
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from logo_modality import as_groups, drop_group, logo_batch, two_layer_ratios
from order_graph_nx import (
    GRAPH_PACKAGE,
    GRAPH_PACKAGE_VERSION,
    BUNDLE_KEYS,
    CUT_BUNDLED,
    CUT_LEVEL_SET,
    CUT_STRUCTURAL,
    bundled_shift_graph,
    community_south_frac,
    graph_localize,
    level_set_ids,
    loud_community,
    louvain_int_cut,
    metric_level_set_ids,
    order_community_labels,
    order_level_set_labels,
    y_in_graph_attrs,
)
from streaming_po_risk import (
    large_deviation,
    mmd_vs_reference,
    rbf_bandwidth,
    streaming_po_and_mse,
)

OUTCOME_NAMES = frozenset({"y", "Y", "label", "target", "outcome"})
PO_MIN_N = 40
NODE_MIN_N = 8
DEVIATION_RATIO = 2.0


def _as_2d(X) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        return X.reshape(-1, 1)
    return X


def _names(names, p: int, prefix: str) -> tuple[str, ...]:
    if names is None:
        return tuple(f"{prefix}{j}" for j in range(p))
    out = tuple(str(n) for n in names)
    if len(out) != p:
        raise ValueError(f"{prefix} names ({len(out)}) != columns ({p})")
    return out


def assert_not_outcome(names: Sequence[str]) -> None:
    """HH / subset key / selected features must not be Y."""
    bad = [n for n in names if str(n) in OUTCOME_NAMES]
    if bad:
        raise ValueError(f"outcome leaked into features or keys: {bad}")


def pos_share(deltas: Mapping[str, float]) -> dict[str, float]:
    names = list(deltas)
    raw = np.array([max(float(deltas[n]), 0.0) for n in names], dtype=float)
    s = float(raw.sum())
    if s <= 1e-15:
        return {n: 0.0 for n in names}
    return {n: float(raw[i] / s) for i, n in enumerate(names)}


def gated_share(
    deltas: Mapping[str, float],
    *,
    abs_floor: float = 0.0,
    rel_gap: float = 1.4,
) -> dict[str, float]:
    """Shares only if some subset is actually louder. Tiny MMD jitter must not vote."""
    names = list(deltas)
    raw = np.array([max(float(deltas[n]), 0.0) for n in names], dtype=float)
    if raw.size == 0:
        return {}
    mx = float(raw.max())
    mn = float(raw.min())
    if mx <= max(float(abs_floor), 1e-15):
        return {n: 0.0 for n in names}
    if raw.size > 1 and mx <= rel_gap * max(mn, 1e-15) and (mx - mn) < 0.25 * max(mx, 1e-12):
        return {n: 0.0 for n in names}
    return pos_share(deltas)


def slice_stream(tables: Mapping, t: int) -> dict:
    """Reference window plus one incoming batch. No future rows."""
    meta = tables["meta"]
    n_ref = int(meta["n_ref"])
    n_new = int(meta["n_new"])
    lo = n_ref + int(t) * n_new
    hi = lo + n_new
    if hi > len(tables["Y"]):
        raise ValueError(f"batch {t} exceeds stream")
    ref = _cut(tables, 0, n_ref)
    new = _cut(tables, lo, hi)
    return {"ref": ref, "new": new, "meta": dict(meta), "t": int(t)}


def _cut(tables: Mapping, lo: int, hi: int) -> dict:
    sl = slice(int(lo), int(hi))
    return {
        "order_id": np.asarray(tables["order_id"])[sl],
        "merchant_id": np.asarray(tables["merchant_id"])[sl],
        "user_id": np.asarray(tables["user_id"])[sl],
        "region": np.asarray(tables["region"])[sl],
        "X_order": np.asarray(tables["X_order"], dtype=float)[sl],
        "X_merchant": np.asarray(tables["X_merchant"], dtype=float)[sl],
        "X_user": np.asarray(tables["X_user"], dtype=float)[sl],
        "Y": np.asarray(tables["Y"], dtype=float).ravel()[sl],
        "names_order": tuple(tables["names_order"]),
        "names_merchant": tuple(tables["names_merchant"]),
        "names_user": tuple(tables["names_user"]),
    }


def freeze_ref_stats(ref: Mapping, seed: int = 2026) -> dict:
    """Everything that must not be refit on the new batch."""
    Xo = _as_2d(ref["X_order"])
    assert_not_outcome(ref["names_order"])
    assert_not_outcome(ref["names_merchant"])
    assert_not_outcome(ref["names_user"])
    bins = {}
    for j, name in enumerate(ref["names_order"]):
        col = Xo[:, j]
        bins[name] = {
            "median": float(np.median(col)),
            "mean": float(np.mean(col)),
            "std": float(max(np.std(col), 1e-8)),
        }
    return {
        "sigma_order": rbf_bandwidth(Xo, seed=seed),
        "order_bins": bins,
        "merchant_profile": entity_profile(
            ref["merchant_id"], ref["X_merchant"], ref["names_merchant"]
        ),
        "user_profile": entity_profile(ref["user_id"], ref["X_user"], ref["names_user"]),
        "y_mean": float(np.mean(ref["Y"])),
        "seed": int(seed),
    }


def entity_profile(ids, X, names) -> dict:
    ids = np.asarray(ids)
    X = _as_2d(X)
    names = _names(names, X.shape[1], "e")
    assert_not_outcome(names)
    out = {"names": names, "mean": {}, "n": {}}
    for e in np.unique(ids):
        mask = ids == e
        out["mean"][int(e) if np.issubdtype(type(e), np.integer) else e] = X[mask].mean(axis=0)
        out["n"][int(e) if np.issubdtype(type(e), np.integer) else e] = int(mask.sum())
    out["global_mean"] = X.mean(axis=0)
    return out


def join_profile(ids, profile: Mapping) -> np.ndarray:
    """Left-join a frozen D_ref entity profile. Unseen ids get the ref mean."""
    ids = np.asarray(ids)
    g = np.asarray(profile["global_mean"], dtype=float)
    X = np.repeat(g.reshape(1, -1), len(ids), axis=0)
    means = profile["mean"]
    for i, e in enumerate(ids):
        key = int(e) if np.issubdtype(type(e), np.integer) else e
        if key in means:
            X[i] = means[key]
        elif e in means:
            X[i] = means[e]
    return X


def current_entity_matrix(ids, X) -> tuple[np.ndarray, np.ndarray]:
    """One row per entity in this window: current X, not Y."""
    ids = np.asarray(ids)
    X = _as_2d(X)
    uniq = []
    rows = []
    for e in np.unique(ids):
        mask = ids == e
        uniq.append(e)
        rows.append(X[mask].mean(axis=0))
    return np.asarray(uniq), np.vstack(rows)


def cmean_x(X_a, X_b) -> float:
    a = _as_2d(X_a).mean(axis=0)
    b = _as_2d(X_b).mean(axis=0)
    return float(np.linalg.norm(a - b))


def cmean_y(Y_a, Y_b) -> float:
    return float(np.mean(Y_a) - np.mean(Y_b))


def cmean_y_by_ref_bin(x_ref, y_ref, x_new, y_new, median: float) -> float:
    """Change in E[Y|high]−E[Y|low] with bins frozen on D_ref.

    A global E[Y] walk must not make every independent column look loud.
    """
    def _contrast(x, y):
        x = np.asarray(x, dtype=float).ravel()
        y = np.asarray(y, dtype=float).ravel()
        hi = x >= float(median)
        if hi.sum() < 3 or (~hi).sum() < 3:
            return 0.0
        return float(y[hi].mean() - y[~hi].mean())

    return abs(_contrast(x_new, y_new) - _contrast(x_ref, y_ref))


def three_metrics(
    X_ref,
    Y_ref,
    X_new,
    Y_new,
    *,
    sigma=None,
    seed: int = 2026,
    with_po: bool = True,
) -> dict:
    """MMD vs D_ref, PO-risk, Conditional Mean. Same clock as the freeze board."""
    X_ref = _as_2d(X_ref)
    X_new = _as_2d(X_new)
    Y_ref = np.asarray(Y_ref, dtype=float).ravel()
    Y_new = np.asarray(Y_new, dtype=float).ravel()
    if sigma is None:
        sigma = rbf_bandwidth(X_ref, seed=seed)
    mmd = mmd_vs_reference(X_ref, X_new, sigma=sigma, seed=seed)
    po = None
    mse = None
    if with_po and len(X_new) >= PO_MIN_N and len(X_ref) >= PO_MIN_N:
        po, mse = streaming_po_and_mse(X_ref, Y_ref, X_new, Y_new, seed=seed)
    return {
        "mmd": float(mmd),
        "po": None if po is None else float(po),
        "mse": None if mse is None else float(mse),
        "cmean_x": cmean_x(X_new, X_ref),
        "cmean_y": cmean_y(Y_new, Y_ref),
        "n_ref": int(len(X_ref)),
        "n_new": int(len(X_new)),
    }


def _metric_or_zero(row: Mapping, key: str) -> float:
    v = row.get(key)
    if v is None:
        return 0.0
    return float(v)


def grain_localization(
    ref: Mapping,
    new: Mapping,
    stats: Mapping,
    seed: int = 2026,
    with_po: bool = True,
) -> dict:
    """Step 1: every native grain, all three readouts, no join."""
    out = {
        "order": three_metrics(
            ref["X_order"],
            ref["Y"],
            new["X_order"],
            new["Y"],
            sigma=stats["sigma_order"],
            seed=seed,
            with_po=with_po,
        )
    }
    out["merchant_nodes"] = node_localization(
        ref["X_order"],
        ref["Y"],
        new["X_order"],
        new["Y"],
        ref["merchant_id"],
        new["merchant_id"],
        sigma=stats["sigma_order"],
        seed=seed,
        with_po=False,
    )
    out["user_nodes"] = node_localization(
        ref["X_order"],
        ref["Y"],
        new["X_order"],
        new["Y"],
        ref["user_id"],
        new["user_id"],
        sigma=stats["sigma_order"],
        seed=seed,
        with_po=False,
    )
    xm_r_id, xm_r = current_entity_matrix(ref["merchant_id"], ref["X_merchant"])
    xm_n_id, xm_n = current_entity_matrix(new["merchant_id"], new["X_merchant"])
    aligned_r, aligned_n = _align_entities(xm_r_id, xm_r, xm_n_id, xm_n)
    out["merchant_table"] = {
        "cmean_x": cmean_x(aligned_n, aligned_r) if len(aligned_r) else 0.0,
        "n": int(len(aligned_n)),
        "names": list(ref["names_merchant"]),
    }
    return out


def _align_entities(id_a, X_a, id_b, X_b):
    id_a = np.asarray(id_a)
    id_b = np.asarray(id_b)
    common = [e for e in id_a if e in set(id_b.tolist())]
    if not common:
        return np.zeros((0, X_a.shape[1])), np.zeros((0, X_b.shape[1]))
    ia = {e: i for i, e in enumerate(id_a)}
    ib = {e: i for i, e in enumerate(id_b)}
    return X_a[[ia[e] for e in common]], X_b[[ib[e] for e in common]]


def node_localization(
    X_ref,
    Y_ref,
    X_new,
    Y_new,
    ids_ref,
    ids_new,
    *,
    sigma,
    seed: int = 2026,
    min_n: int = NODE_MIN_N,
    with_po: bool = False,
) -> list[dict]:
    """Graph-node scores: orders attached to each entity vs D_ref (and vs own-ref)."""
    ids_ref = np.asarray(ids_ref)
    ids_new = np.asarray(ids_new)
    rows = []
    for e in np.unique(ids_new):
        m_n = ids_new == e
        n = int(m_n.sum())
        if n < int(min_n):
            continue
        m_r = ids_ref == e
        vs_full = three_metrics(
            X_ref,
            Y_ref,
            X_new[m_n],
            Y_new[m_n],
            sigma=sigma,
            seed=seed,
            with_po=with_po and n >= PO_MIN_N,
        )
        own = None
        if int(m_r.sum()) >= int(min_n):
            own = three_metrics(
                X_ref[m_r],
                Y_ref[m_r],
                X_new[m_n],
                Y_new[m_n],
                sigma=rbf_bandwidth(_as_2d(X_ref)[m_r], seed=seed),
                seed=seed,
                with_po=False,
            )
        rows.append(
            {
                "node": int(e) if np.issubdtype(type(e), np.integer) else e,
                "n": n,
                "n_ref_node": int(m_r.sum()),
                "vs_full_ref": vs_full,
                "vs_own_ref": own,
            }
        )
    rows.sort(key=lambda r: r["vs_full_ref"]["mmd"], reverse=True)
    return rows


def node_bundle_scores(node_rows: Sequence[Mapping], *, prefer: str = "vs_own_ref") -> dict[int, dict]:
    """Flatten node rows into merchant → bundled metrics.

    Default vs_own_ref: did *this* merchant move. vs_full_ref would just
    say every merchant differs from the global mean (heterogeneity, not shift).
    """
    out = {}
    for r in node_rows:
        if r.get(prefer) is None and prefer == "vs_own_ref":
            continue
        src = dict(r.get(prefer) or r.get("vs_full_ref") or {})
        full = r.get("vs_full_ref") or {}
        if src.get("po") is None and full.get("po") is not None:
            src["po"] = full["po"]
        nid = int(r["node"])
        out[nid] = {
            "mmd": src.get("mmd"),
            "po": src.get("po"),
            "cmean_x": src.get("cmean_x"),
            "cmean_y": src.get("cmean_y"),
            "n": r.get("n"),
        }
    return out


def graph_shift_cuts(
    ref: Mapping,
    new: Mapping,
    stats: Mapping,
    *,
    seed: int = 2026,
    k_nn: int = 4,
    min_n: int = NODE_MIN_N,
) -> dict:
    """Structural / bundled Louvain (contrast) plus the level-set cut (primary)."""
    structural = graph_localize(new, k_nn=k_nn, seed=seed)
    nodes = node_localization(
        ref["X_order"],
        ref["Y"],
        new["X_order"],
        new["Y"],
        ref["merchant_id"],
        new["merchant_id"],
        sigma=stats["sigma_order"],
        seed=seed,
        min_n=min_n,
        with_po=False,
    )
    scores = node_bundle_scores(nodes)
    G_b = bundled_shift_graph(scores)
    bundled_cut = louvain_int_cut(G_b, seed=seed)
    region_by_m = {}
    for m, r in zip(np.asarray(new["merchant_id"]).astype(int), np.asarray(new["region"])):
        region_by_m.setdefault(int(m), str(r))
    hot = loud_community(bundled_cut, scores)
    return {
        "package": GRAPH_PACKAGE,
        "package_version": GRAPH_PACKAGE_VERSION,
        "structural": structural,
        "bundled": {
            "cut_name": CUT_BUNDLED,
            "merchant_cut": bundled_cut,
            "n_communities": int(len(set(bundled_cut.values()))) if bundled_cut else 0,
            "n_nodes": int(G_b.number_of_nodes()),
            "n_edges": int(G_b.number_of_edges()),
            "south_frac": community_south_frac(bundled_cut, region_by_m),
            "loud_community": hot,
            "loud_south_frac": float(community_south_frac(bundled_cut, region_by_m).get(hot, 0.0))
            if hot
            else 0.0,
            "y_in_graph": y_in_graph_attrs(G_b),
        },
        "scores": {str(k): v for k, v in scores.items()},
        "region_by_merchant": {str(k): v for k, v in region_by_m.items()},
        "y_in_structural": y_in_graph_attrs(structural["graph"]),
        "own_vs_full": own_vs_full_summary(nodes, region_by_m),
        "layers": layer_eval_on_orders(ref, new, stats, merchant_nodes=nodes, seed=seed, min_n=min_n),
        "level_set": multilayer_level_set(
            ref, new, stats, merchant_nodes=nodes, seed=seed, min_n=min_n
        ),
    }


def jaccard_masks(a, b) -> float:
    a = np.asarray(a, dtype=bool).ravel()
    b = np.asarray(b, dtype=bool).ravel()
    if a.size != b.size:
        raise ValueError("masks must align on the order grain")
    inter = int(np.logical_and(a, b).sum())
    union = int(np.logical_or(a, b).sum())
    return float(inter / union) if union else 0.0


def lift_cut_to_orders(entity_ids, cut: Mapping[int, str], default: str = "C_quiet") -> np.ndarray:
    """Entity-layer labels → one label per order. Evaluation grain is always the order."""
    entity_ids = np.asarray(entity_ids).astype(int)
    return np.asarray([cut.get(int(e), default) for e in entity_ids], dtype=object)


def loud_order_mask(entity_ids, cut: Mapping[int, str], scores: Mapping[int, Mapping]):
    hot = loud_community(cut, scores)
    labs = lift_cut_to_orders(entity_ids, cut)
    if hot is None:
        return np.zeros(len(labs), dtype=bool), None, labs
    return labs == hot, hot, labs


def own_vs_full_summary(node_rows: Sequence[Mapping], region_by_id: Mapping[int, str]) -> dict:
    """Why own-ref is the changing-subset clock, and full-ref is heterogeneity.

    own:  this node's new bag vs this node's D_ref bag
    full: this node's new bag vs the global D_ref
    A niche merchant who did not move is loud on full, quiet on own.
    A planted shifter is loud on own.
    """
    south_own, north_own, south_full, north_full = [], [], [], []
    n_own, n_skip = 0, 0
    for r in node_rows:
        own, full = r.get("vs_own_ref"), r.get("vs_full_ref") or {}
        lab = str(region_by_id.get(int(r["node"]), ""))
        if own is None:
            n_skip += 1
            continue
        n_own += 1
        om, fm = float(own.get("mmd") or 0.0), float(full.get("mmd") or 0.0)
        oc, fc = float(own.get("cmean_x") or 0.0), float(full.get("cmean_x") or 0.0)
        if lab == "south":
            south_own.append((om, oc, abs(float(own.get("cmean_y") or 0.0))))
            south_full.append((fm, fc))
        else:
            north_own.append((om, oc, abs(float(own.get("cmean_y") or 0.0))))
            north_full.append((fm, fc))

    def _mean(rows, j):
        if not rows:
            return None
        return float(np.mean([r[j] for r in rows]))

    return {
        "n_with_own_ref": n_own,
        "n_cold_start_skipped": n_skip,
        "south_own_mmd": _mean(south_own, 0),
        "north_own_mmd": _mean(north_own, 0),
        "south_own_cmean_x": _mean(south_own, 1),
        "north_own_cmean_x": _mean(north_own, 1),
        "south_own_cmean_y": _mean(south_own, 2),
        "north_own_cmean_y": _mean(north_own, 2),
        "south_full_mmd": _mean(south_full, 0),
        "north_full_mmd": _mean(north_full, 0),
        "south_full_cmean_x": _mean(south_full, 1),
        "north_full_cmean_x": _mean(north_full, 1),
        "read": (
            "own-ref = did this merchant move; full-ref = does this merchant "
            "differ from the global mix (heterogeneity, not shift)"
        ),
    }


def _layer_cut(ref, new, stats, id_key: str, min_n: int, seed: int) -> dict:
    nodes = node_localization(
        ref["X_order"],
        ref["Y"],
        new["X_order"],
        new["Y"],
        ref[id_key],
        new[id_key],
        sigma=stats["sigma_order"],
        seed=seed,
        min_n=min_n,
        with_po=False,
    )
    scores = node_bundle_scores(nodes, prefer="vs_own_ref")
    G = bundled_shift_graph(scores)
    cut = louvain_int_cut(G, seed=seed)
    mask, hot, labs = loud_order_mask(new[id_key], cut, scores)
    planted = np.asarray(new["region"]) == "south"
    return {
        "layer": id_key.replace("_id", ""),
        "n_nodes": int(G.number_of_nodes()),
        "n_edges": int(G.number_of_edges()),
        "n_communities": int(len(set(cut.values()))) if cut else 0,
        "loud_community": hot,
        "n_loud_orders": int(mask.sum()),
        "jaccard_vs_planted_south": jaccard_masks(mask, planted),
        "south_frac_loud_orders": float(planted[mask].mean()) if int(mask.sum()) else 0.0,
        "cut": {str(k): v for k, v in cut.items()},
        "y_in_graph": y_in_graph_attrs(G),
        "_mask": mask,
    }


def layer_eval_on_orders(
    ref: Mapping,
    new: Mapping,
    stats: Mapping,
    *,
    merchant_nodes=None,
    seed: int = 2026,
    min_n: int = NODE_MIN_N,
) -> dict:
    """Cut each layer on its own nodes, then evaluate only after lifting to orders.

    Same-dimension localization stays on one grain. Multi-layer graphs are
    compared by Jaccard of loud order-sets — never by mixing merchant-MMD
    with user-MMD in one simplex.
    """
    mer = _layer_cut(ref, new, stats, "merchant_id", min_n=min_n, seed=seed)
    user_min = max(int(min_n), 4)
    usr = _layer_cut(ref, new, stats, "user_id", min_n=user_min, seed=seed)
    planted = np.asarray(new["region"]) == "south"
    mer_mask = mer.pop("_mask")
    usr_mask = usr.pop("_mask")
    return {
        "merchant": mer,
        "user": usr,
        "order_oracle_south_n": int(planted.sum()),
        "jaccard_merchant_vs_user": jaccard_masks(mer_mask, usr_mask),
        "jaccard_merchant_vs_south": float(mer["jaccard_vs_planted_south"]),
        "jaccard_user_vs_south": float(usr["jaccard_vs_planted_south"]),
        "read": (
            "evaluate layers only after lift-to-order; merchant layer should "
            "recover planted south when the DGP is merchant-local; user layer "
            "need not (users are mixed across merchants)"
        ),
    }


def lift_level_set_to_orders(entity_ids, loud_ids: Sequence[int]) -> np.ndarray:
    """Entity-layer level-set → one bool per order. Evaluation grain is always the order."""
    loud = {int(i) for i in loud_ids}
    ids = np.asarray(entity_ids).astype(int)
    return np.asarray([int(e) in loud for e in ids], dtype=bool)


def _south_frac_ids(ids: Sequence[int], region_by_id: Mapping[int, str]) -> float:
    if not ids:
        return 0.0
    return float(np.mean([str(region_by_id.get(int(i), "")) == "south" for i in ids]))


def _level_set_layer(
    ref,
    new,
    stats,
    id_key: str,
    min_n: int,
    seed: int,
    *,
    scores: Mapping[int, Mapping] | None = None,
    floor: float = 1.0,
    frac: float = 0.30,
) -> dict:
    """Score this layer vs own-ref, take {φ ≥ τ}, lift to orders. No Louvain."""
    if scores is None:
        nodes = node_localization(
            ref["X_order"],
            ref["Y"],
            new["X_order"],
            new["Y"],
            ref[id_key],
            new[id_key],
            sigma=stats["sigma_order"],
            seed=seed,
            min_n=min_n,
            with_po=False,
        )
        scores = node_bundle_scores(nodes, prefer="vs_own_ref")
    ls = level_set_ids(scores, floor=floor, frac=frac)
    mask = lift_level_set_to_orders(new[id_key], ls["loud_ids"])
    planted = np.asarray(new["region"]) == "south"
    region_by: dict[int, str] = {}
    for e, r in zip(np.asarray(new[id_key]).astype(int), np.asarray(new["region"])):
        region_by.setdefault(int(e), str(r))
    slices = {}
    for metric in BUNDLE_KEYS:
        sl = metric_level_set_ids(scores, metric, frac=frac)
        slices[metric] = {
            "n_loud": sl["n_loud"],
            "tau": sl["tau"],
            "south_frac": _south_frac_ids(sl["loud_ids"], region_by),
            "jaccard_vs_planted_south": jaccard_masks(
                lift_level_set_to_orders(new[id_key], sl["loud_ids"]), planted
            ),
        }
    return {
        "layer": id_key.replace("_id", ""),
        "cut_name": CUT_LEVEL_SET,
        "loud_ids": [int(i) for i in ls["loud_ids"]],
        "n_loud_nodes": int(ls["n_loud"]),
        "n_scored": int(len(scores)),
        "tau": float(ls["tau"]),
        "n_loud_orders": int(mask.sum()),
        "jaccard_vs_planted_south": jaccard_masks(mask, planted),
        "south_frac_loud_nodes": _south_frac_ids(ls["loud_ids"], region_by),
        "south_frac_loud_orders": float(planted[mask].mean()) if int(mask.sum()) else 0.0,
        "slices": slices,
        "phis": ls["phis"],
        "_mask": mask,
    }


def multilayer_level_set(
    ref: Mapping,
    new: Mapping,
    stats: Mapping,
    *,
    merchant_nodes=None,
    seed: int = 2026,
    min_n: int = NODE_MIN_N,
) -> dict:
    """Per-layer level sets, evaluated only after lifting to orders.

    Graph = incidence (order→merchant, order→user). The cut on each layer is
    Ŝ = {i : φ_i ≥ τ} against that node's own D_ref. No community detection.
    """
    mer_scores = None
    if merchant_nodes is not None:
        mer_scores = node_bundle_scores(merchant_nodes, prefer="vs_own_ref")
    mer = _level_set_layer(
        ref, new, stats, "merchant_id", min_n=min_n, seed=seed, scores=mer_scores
    )
    user_min = max(int(min_n), 4)
    usr = _level_set_layer(ref, new, stats, "user_id", min_n=user_min, seed=seed)
    planted = np.asarray(new["region"]) == "south"
    mer_mask = mer.pop("_mask")
    usr_mask = usr.pop("_mask")
    return {
        "cut_name": CUT_LEVEL_SET,
        "merchant": mer,
        "user": usr,
        "order_oracle_south_n": int(planted.sum()),
        "jaccard_merchant_vs_user": jaccard_masks(mer_mask, usr_mask),
        "jaccard_merchant_vs_south": float(mer["jaccard_vs_planted_south"]),
        "jaccard_user_vs_south": float(usr["jaccard_vs_planted_south"]),
        "read": (
            "changing subset = level set of own-ref φ, not Louvain; "
            "evaluate layers only after lift-to-order; merchant layer should "
            "recover planted south when the DGP is merchant-local; user layer "
            "need not (users are mixed across merchants)"
        ),
    }


def fsds_rank_columns(
    X_ref,
    X_new,
    Y_ref,
    Y_new,
    names,
    *,
    bins: Mapping | None = None,
    seed: int = 2026,
) -> list[dict]:
    """FSDS inside one grain: univariate MMD + CMean_X + binned CMean_Y.

    Bins / σ come from D_ref. Y is the outcome of CMean_Y, not a column.
    """
    X_ref = _as_2d(X_ref)
    X_new = _as_2d(X_new)
    Y_ref = np.asarray(Y_ref, dtype=float).ravel()
    Y_new = np.asarray(Y_new, dtype=float).ravel()
    names = _names(names, X_ref.shape[1], "x")
    assert_not_outcome(names)
    use_cy = len(X_ref) >= 30
    rows = []
    for j, name in enumerate(names):
        xr, xn = X_ref[:, [j]], X_new[:, [j]]
        sig = rbf_bandwidth(xr, seed=seed)
        mmd = mmd_vs_reference(xr, xn, sigma=sig, seed=seed)
        std = 1.0
        med = float(np.median(X_ref[:, j]))
        if bins and name in bins:
            std = float(bins[name]["std"])
            med = float(bins[name]["median"])
        dx = abs(float(X_new[:, j].mean()) - float(X_ref[:, j].mean())) / max(std, 1e-8)
        cy = (
            cmean_y_by_ref_bin(X_ref[:, j], Y_ref, X_new[:, j], Y_new, med)
            if use_cy
            else 0.0
        )
        quiet = _column_quiet(xr, Y_ref, seed=seed)
        score = max(
            mmd / max(quiet["mmd"], 1e-8),
            dx / max(quiet["cmean_x"], 1e-8),
            cy / max(quiet["cmean_y"], 1e-8),
        )
        rows.append(
            {
                "feature": name,
                "j": int(j),
                "mmd": float(mmd),
                "cmean_x": float(dx),
                "cmean_y": float(cy),
                "score": float(score),
                "loud": bool(
                    large_deviation(mmd, quiet["mmd"], ratio=DEVIATION_RATIO)
                    or large_deviation(dx, quiet["cmean_x"], ratio=DEVIATION_RATIO)
                    or (cy >= DEVIATION_RATIO * max(quiet["cmean_y"], 1e-8) and cy > 0.12)
                ),
            }
        )
    rows.sort(key=lambda r: r["score"], reverse=True)
    return rows


def _column_quiet(x_ref, y_ref, seed: int = 2026) -> dict:
    x_ref = _as_2d(x_ref)
    y_ref = np.asarray(y_ref, dtype=float).ravel()
    n = len(x_ref)
    half = max(n // 2, 1)
    sig = rbf_bandwidth(x_ref[:half], seed=seed)
    mmd = mmd_vs_reference(x_ref[:half], x_ref[half:], sigma=sig, seed=seed)
    dx = abs(float(x_ref[:half].mean()) - float(x_ref[half:].mean())) / max(
        float(np.std(x_ref)), 1e-8
    )
    med = float(np.median(x_ref))
    cy = cmean_y_by_ref_bin(x_ref[:half, 0], y_ref[:half], x_ref[half:, 0], y_ref[half:], med)
    return {
        "mmd": max(float(mmd), 0.015),
        "cmean_x": max(float(dx), 0.25),
        "cmean_y": max(float(cy), 0.10),
    }


def fsds_select(
    ref: Mapping,
    new: Mapping,
    stats: Mapping,
    *,
    top_k: int = 4,
    seed: int = 2026,
) -> dict:
    """Step 2: per-grain ranking, then the union of loud / top-k columns."""
    order_rank = fsds_rank_columns(
        ref["X_order"],
        new["X_order"],
        ref["Y"],
        new["Y"],
        ref["names_order"],
        bins=stats["order_bins"],
        seed=seed,
    )
    mer_r_id, mer_r = current_entity_matrix(ref["merchant_id"], ref["X_merchant"])
    mer_n_id, mer_n = current_entity_matrix(new["merchant_id"], new["X_merchant"])
    ar, an = _align_entities(mer_r_id, mer_r, mer_n_id, mer_n)
    y_mer_r = _entity_mean_y(ref["merchant_id"], ref["Y"], mer_r_id)
    y_mer_n = _entity_mean_y(new["merchant_id"], new["Y"], mer_n_id)
    y_r_al, y_n_al = _align_y(mer_r_id, y_mer_r, mer_n_id, y_mer_n)
    if len(ar):
        merchant_rank = fsds_rank_columns(
            ar, an, y_r_al, y_n_al, ref["names_merchant"], seed=seed
        )
    else:
        merchant_rank = []
    user_rank = fsds_rank_columns(
        ref["X_user"],
        new["X_user"],
        ref["Y"],
        new["Y"],
        ref["names_user"],
        seed=seed,
    )
    selected = {
        "order": _pick(order_rank, top_k, fill_if_quiet=True),
        "merchant": _pick(merchant_rank, max(1, top_k // 2), fill_if_quiet=False),
        "user": _pick(user_rank, max(1, top_k // 2), fill_if_quiet=False),
    }
    names = [r["feature"] for g in selected.values() for r in g]
    assert_not_outcome(names)
    return {
        "rank": {"order": order_rank, "merchant": merchant_rank, "user": user_rank},
        "selected": selected,
        "selected_names": names,
    }


def _entity_mean_y(ids, Y, uniq) -> np.ndarray:
    ids = np.asarray(ids)
    Y = np.asarray(Y, dtype=float).ravel()
    return np.array([float(Y[ids == e].mean()) for e in uniq], dtype=float)


def _align_y(id_a, y_a, id_b, y_b):
    id_a = np.asarray(id_a)
    id_b = np.asarray(id_b)
    common = [e for e in id_a if e in set(id_b.tolist())]
    ia = {e: i for i, e in enumerate(id_a)}
    ib = {e: i for i, e in enumerate(id_b)}
    return (
        np.array([y_a[ia[e]] for e in common], dtype=float),
        np.array([y_b[ib[e]] for e in common], dtype=float),
    )


def _pick(rank: Sequence[Mapping], k: int, fill_if_quiet: bool = False) -> list[dict]:
    """Prefer loud columns. Quiet merchant/user blocks stay empty."""
    if not rank:
        return []
    loud = [r for r in rank if r.get("loud")]
    if loud:
        return loud[: int(k)]
    if fill_if_quiet:
        return list(rank[: max(1, min(int(k), 2))])
    return []


def unify_to_order(
    window: Mapping,
    selected: Mapping,
    stats: Mapping,
    *,
    mode: str = "localize",
) -> dict:
    """Step 3: selected features onto the order grain.

    localize — current X of selected columns (object of study; still no Y).
    serve    — order X now, entity columns from frozen D_ref profiles.
    """
    names_o = list(window["names_order"])
    names_m = list(window["names_merchant"])
    names_u = list(window["names_user"])
    take_o = [r["feature"] for r in selected.get("order", [])]
    take_m = [r["feature"] for r in selected.get("merchant", [])]
    take_u = [r["feature"] for r in selected.get("user", [])]
    assert_not_outcome(take_o + take_m + take_u)
    cols = []
    groups = {"order": [], "merchant": [], "user": []}
    names = []

    def _append(block, X, take, source_names):
        idx = [source_names.index(n) for n in take if n in source_names]
        if not idx:
            return
        start = len(names)
        cols.append(_as_2d(X)[:, idx])
        names.extend([source_names[j] for j in idx])
        groups[block] = list(range(start, len(names)))

    _append("order", window["X_order"], take_o, names_o)
    if mode == "serve":
        Xm = join_profile(window["merchant_id"], stats["merchant_profile"])
        Xu = join_profile(window["user_id"], stats["user_profile"])
    elif mode == "localize":
        Xm = _as_2d(window["X_merchant"])
        Xu = _as_2d(window["X_user"])
    else:
        raise ValueError(f"unknown unify mode {mode!r}")
    _append("merchant", Xm, take_m, names_m)
    _append("user", Xu, take_u, names_u)
    if not cols:
        raise ValueError("FSDS selected an empty feature set")
    Z = np.hstack(cols)
    assert_not_outcome(names)
    return {
        "Z": Z,
        "Y": np.asarray(window["Y"], dtype=float).ravel(),
        "names": tuple(names),
        "groups": {k: v for k, v in groups.items() if v},
        "mode": mode,
        "region": np.asarray(window["region"]),
        "merchant_id": np.asarray(window["merchant_id"]),
        "user_id": np.asarray(window["user_id"]),
    }


def unify_to_merchant(order_unified: Mapping) -> dict:
    """Same selected columns, one row per merchant (order-mean). Still no Y in Z."""
    ids = np.asarray(order_unified["merchant_id"])
    Z = _as_2d(order_unified["Z"])
    Y = np.asarray(order_unified["Y"], dtype=float).ravel()
    region = np.asarray(order_unified["region"])
    uniq, Z_m = current_entity_matrix(ids, Z)
    y_m = _entity_mean_y(ids, Y, uniq)
    reg = []
    for e in uniq:
        r = region[ids == e]
        reg.append(str(r[0]) if len(r) else "unknown")
    return {
        "Z": Z_m,
        "Y": y_m,
        "names": order_unified["names"],
        "groups": order_unified["groups"],
        "mode": order_unified["mode"],
        "region": np.asarray(reg),
        "merchant_id": uniq,
        "user_id": None,
        "grain": "merchant",
    }


def subset_three_metrics(
    Z_ref,
    Y_ref,
    Z_new,
    Y_new,
    labels_ref,
    labels_new,
    *,
    sigma,
    seed: int = 2026,
    min_n: int = 20,
    with_po: bool = True,
) -> list[dict]:
    """Subset localization on the unified grain: MMD, PO, Conditional Mean."""
    labels_new = np.asarray(labels_new)
    labels_ref = np.asarray(labels_ref)
    Z_ref, Z_new = _as_2d(Z_ref), _as_2d(Z_new)
    Y_ref = np.asarray(Y_ref, dtype=float).ravel()
    Y_new = np.asarray(Y_new, dtype=float).ravel()
    rows = []
    for lab in sorted(set(labels_new.tolist()), key=str):
        m_n = labels_new == lab
        n = int(m_n.sum())
        if n < int(min_n):
            continue
        m_r = labels_ref == lab
        vs_full = three_metrics(
            Z_ref,
            Y_ref,
            Z_new[m_n],
            Y_new[m_n],
            sigma=sigma,
            seed=seed,
            with_po=bool(with_po) and n >= PO_MIN_N,
        )
        own = None
        if int(m_r.sum()) >= int(min_n):
            own = three_metrics(
                Z_ref[m_r],
                Y_ref[m_r],
                Z_new[m_n],
                Y_new[m_n],
                sigma=rbf_bandwidth(Z_ref[m_r], seed=seed),
                seed=seed,
                with_po=bool(with_po)
                and n >= PO_MIN_N
                and int(m_r.sum()) >= PO_MIN_N,
            )
        other = ~m_n
        vs_other = None
        if int(other.sum()) >= int(min_n):
            vs_other = three_metrics(
                Z_ref,
                Y_ref,
                Z_new[other],
                Y_new[other],
                sigma=sigma,
                seed=seed,
                with_po=bool(with_po) and int(other.sum()) >= PO_MIN_N,
            )
        rows.append(
            {
                "subset": lab,
                "n": n,
                "share_of_batch": float(n / max(len(Y_new), 1)),
                "vs_full_ref": vs_full,
                "vs_own_ref": own,
                "other": vs_other,
            }
        )
    return rows


def subset_shares(rows: Sequence[Mapping], key: str = "vs_own_ref") -> dict[str, dict]:
    """Layer-1 shares across subsets on each of the three readouts."""

    def _val(row, metric):
        src = row.get(key) or {}
        v = src.get(metric)
        if v is None and metric == "po":
            v = (row.get("vs_full_ref") or {}).get("po")
        if v is None and key != "vs_full_ref":
            v = (row.get("vs_full_ref") or {}).get(metric)
        if v is None:
            return 0.0
        if metric == "cmean_y":
            return abs(float(v))
        return max(float(v), 0.0)

    mmd = {r["subset"]: _val(r, "mmd") for r in rows}
    po = {r["subset"]: _val(r, "po") for r in rows}
    cx = {r["subset"]: _val(r, "cmean_x") for r in rows}
    cy = {r["subset"]: _val(r, "cmean_y") for r in rows}
    pi_mmd = gated_share(mmd, abs_floor=0.02)
    pi_po = gated_share(po, abs_floor=1e-3)
    pi_cx = gated_share(cx, abs_floor=0.15)
    pi_cy = gated_share(cy, abs_floor=0.05)
    mixed = {}
    for s in mmd:
        mixed[s] = {
            "pi_mmd": pi_mmd[s],
            "pi_po": pi_po[s],
            "pi_cmean_x": pi_cx[s],
            "pi_cmean_y": pi_cy[s],
            "pi_cmean": pos_share({k: 0.5 * pi_cy[k] + 0.5 * pi_cx[k] for k in cy})[s]
            if any(pi_cy[k] + pi_cx[k] > 0 for k in cy)
            else 0.0,
        }
        den = mixed[s]["pi_mmd"] + mixed[s]["pi_po"] + mixed[s]["pi_cmean"]
        mixed[s]["mix_mmd"] = 0.0 if den <= 1e-12 else mixed[s]["pi_mmd"] / den
        mixed[s]["mix_po"] = 0.0 if den <= 1e-12 else mixed[s]["pi_po"] / den
        mixed[s]["mix_cmean"] = 0.0 if den <= 1e-12 else mixed[s]["pi_cmean"] / den
        mixed[s]["contribution"] = mixed[s]["pi_mmd"] + mixed[s]["pi_po"] + mixed[s]["pi_cmean"]
    return mixed


def loud_subset(shares: Mapping[str, Mapping]) -> str | None:
    if not shares:
        return None
    if max(float(v["contribution"]) for v in shares.values()) <= 1e-12:
        return None
    return max(shares, key=lambda s: float(shares[s]["contribution"]))


def complement_gap(row: Mapping) -> dict:
    """Subset vs other, both scored against the same D_ref."""
    s = row.get("vs_full_ref") or {}
    o = row.get("other") or {}

    def gap(metric):
        a, b = s.get(metric), o.get(metric)
        if a is None or b is None:
            return None
        if metric == "cmean_y":
            return abs(float(a)) - abs(float(b))
        return float(a) - float(b)

    return {"mmd": gap("mmd"), "po": gap("po"), "cmean_x": gap("cmean_x"), "cmean_y": gap("cmean_y")}


def fingerprint(mix: Mapping, own: Mapping | None) -> str:
    """Mechanism tag from mix(PO, MMD, CMean). Not a unique decomp."""
    mix_mmd = float(mix.get("mix_mmd", 0.0))
    mix_po = float(mix.get("mix_po", 0.0))
    own = own or {}
    mmd = _metric_or_zero(own, "mmd")
    po = _metric_or_zero(own, "po")
    cx = _metric_or_zero(own, "cmean_x")
    cy = abs(_metric_or_zero(own, "cmean_y"))
    if max(mmd, po, cx, cy) <= 1e-12:
        return "quiet"
    if mix_mmd >= 0.5 and mix_po < 0.3:
        return "x_shift"
    if mix_mmd < 0.2 and mix_po >= 0.35:
        return "concept"
    if mix_mmd >= 0.3 and mix_po >= 0.25:
        return "both"
    if mix_mmd >= 0.5:
        return "x_shift"
    if cx > 1e-8 and cy > 1e-8 and mix_po < 0.3:
        return "x_shift"
    return "both"


def characterize_subset(row: Mapping, mix: Mapping) -> dict:
    """Concrete portrait of one subset on the unified grain."""
    own = dict(row.get("vs_own_ref") or {})
    full = row.get("vs_full_ref") or {}
    if own.get("po") is None:
        own["po"] = full.get("po")
    if not own:
        own = dict(full)
    tag = fingerprint(mix, own)
    cy = float((own or {}).get("cmean_y") or 0.0)
    cx = float((own or {}).get("cmean_x") or 0.0)
    reads = {
        "x_shift": (
            "this subset's P(X) moved vs its own D_ref; "
            "E[Y|S] may follow because f is not constant"
        ),
        "concept": (
            "this subset's P(Y|X) hopped; E[Y|S] moved while selected X stayed closer"
        ),
        "both": "P(X) and P(Y|X) both look loud on this subset — localization, not a split",
        "quiet": "three readouts are quiet on this subset",
    }
    return {
        "subset": row["subset"],
        "n": row["n"],
        "share_of_batch": row["share_of_batch"],
        "mmd": own.get("mmd"),
        "po": own.get("po"),
        "cmean_x": cx,
        "cmean_y": cy,
        "pi_mmd": mix["pi_mmd"],
        "pi_po": mix["pi_po"],
        "pi_cmean": mix["pi_cmean"],
        "mix_mmd": mix["mix_mmd"],
        "mix_po": mix["mix_po"],
        "mix_cmean": mix["mix_cmean"],
        "contribution": mix["contribution"],
        "gap_vs_other": complement_gap(row),
        "fingerprint": tag,
        "read": reads[tag],
    }


def layer2_logo(Z_ref, Y_ref, Z_new, Y_new, groups, seed: int = 2026) -> dict | None:
    """LOGO on selected blocks (order / merchant / user) of the unified table."""
    groups = {k: v for k, v in as_groups(groups).items() if len(v)}
    if len(groups) < 2:
        return None
    return logo_batch(Z_ref, Y_ref, Z_new, Y_new, groups, seed=seed)


def run_pipeline(
    tables: Mapping,
    t: int,
    *,
    grain: str = "order",
    mode: str = "localize",
    subset_by: str = "level_set",
    top_k: int = 4,
    seed: int = 2026,
    min_n: int = 20,
    with_logo: bool = True,
    with_po: bool = True,
    k_nn: int = 4,
) -> dict:
    """Level-set localization first (changing subset = {φ ≥ τ}), then FSDS, then two-layer.

    subset_by='level_set' labels orders loud vs other from the merchant-layer
    own-ref level set. community / structural Louvain and region are contrasts.
    """
    cut = slice_stream(tables, t)
    ref, new, meta = cut["ref"], cut["new"], cut["meta"]
    stats = freeze_ref_stats(ref, seed=seed)
    graph_pack = graph_shift_cuts(ref, new, stats, seed=seed, k_nn=k_nn)
    grains = grain_localization(ref, new, stats, seed=seed, with_po=with_po)
    fsds = fsds_select(ref, new, stats, top_k=top_k, seed=seed)
    uni_ref = unify_to_order(ref, fsds["selected"], stats, mode=mode)
    uni_new = unify_to_order(new, fsds["selected"], stats, mode=mode)
    if grain == "merchant":
        uni_ref = unify_to_merchant(uni_ref)
        uni_new = unify_to_merchant(uni_new)
        min_n = min(int(min_n), 2)
    elif grain != "order":
        raise ValueError(f"grain must be order or merchant, got {grain!r}")
    if subset_by == "level_set":
        loud_ids = graph_pack["level_set"]["merchant"]["loud_ids"]
        labels_new = order_level_set_labels(uni_new["merchant_id"], loud_ids)
        labels_ref = order_level_set_labels(uni_ref["merchant_id"], loud_ids)
        planted_new = np.asarray(uni_new["region"]) == "south"
        loud_mask = np.asarray(labels_new) == "loud"
        other_mask = np.asarray(labels_new) == "other"
        south_frac = {
            "loud": float(planted_new[loud_mask].mean()) if int(loud_mask.sum()) else 0.0,
            "other": float(planted_new[other_mask].mean()) if int(other_mask.sum()) else 0.0,
        }
        cut_name = CUT_LEVEL_SET
    elif subset_by == "community":
        merchant_cut = graph_pack["bundled"]["merchant_cut"]
        labels_new = order_community_labels(uni_new["merchant_id"], merchant_cut)
        labels_ref = order_community_labels(uni_ref["merchant_id"], merchant_cut)
        south_frac = graph_pack["bundled"]["south_frac"]
        cut_name = CUT_BUNDLED
    elif subset_by == "structural":
        merchant_cut = graph_pack["structural"]["merchant_cut"]
        labels_new = order_community_labels(uni_new["merchant_id"], merchant_cut)
        labels_ref = order_community_labels(uni_ref["merchant_id"], merchant_cut)
        south_frac = graph_pack["structural"]["south_frac"]
        cut_name = CUT_STRUCTURAL
    elif subset_by == "region":
        labels_ref, labels_new = uni_ref["region"], uni_new["region"]
        south_frac = {}
        cut_name = "region"
    else:
        raise ValueError(
            f"subset_by must be level_set, community, structural, or region, got {subset_by!r}"
        )
    sigma_z = rbf_bandwidth(uni_ref["Z"], seed=seed)
    subsets = subset_three_metrics(
        uni_ref["Z"],
        uni_ref["Y"],
        uni_new["Z"],
        uni_new["Y"],
        labels_ref,
        labels_new,
        sigma=sigma_z,
        seed=seed,
        min_n=min_n,
        with_po=with_po,
    )
    shares = subset_shares(subsets, key="vs_own_ref")
    hot = loud_subset(shares)
    portraits = []
    for row in subsets:
        portraits.append(characterize_subset(row, shares[row["subset"]]))
    portraits.sort(key=lambda r: r["contribution"], reverse=True)
    if subset_by in ("community", "structural", "level_set"):
        for p in portraits:
            p["south_frac"] = float(south_frac.get(str(p["subset"]), 0.0))
    logo_full = None
    logo_sub = None
    if with_logo:
        logo_full = layer2_logo(
            uni_ref["Z"], uni_ref["Y"], uni_new["Z"], uni_new["Y"], uni_new["groups"], seed=seed
        )
        if hot is not None:
            mask = np.asarray(labels_new) == hot
            if int(mask.sum()) >= PO_MIN_N:
                logo_sub = layer2_logo(
                    uni_ref["Z"],
                    uni_ref["Y"],
                    uni_new["Z"][mask],
                    uni_new["Y"][mask],
                    uni_new["groups"],
                    seed=seed,
                )
    leakage = {
        "y_in_Z": any(n in OUTCOME_NAMES for n in uni_new["names"]),
        "subset_is_y": False,
        "mode": mode,
        "sigma_from_ref": True,
        "bins_from_ref": True,
        "serve_profiles_from_ref": mode == "serve",
        "graph_package": GRAPH_PACKAGE,
        "graph_package_version": GRAPH_PACKAGE_VERSION,
        "graph_cut": cut_name,
        "y_in_bundled_graph": graph_pack["bundled"]["y_in_graph"],
        "y_in_structural_graph": graph_pack["y_in_structural"],
    }
    return {
        "t": int(t),
        "grain": grain,
        "subset_by": subset_by,
        "meta": meta,
        "graph": {
            "package": GRAPH_PACKAGE,
            "package_version": GRAPH_PACKAGE_VERSION,
            "encoder": graph_pack["structural"]["encoder"],
            "cut": cut_name,
            "bundled": {
                "n_communities": graph_pack["bundled"]["n_communities"],
                "n_nodes": graph_pack["bundled"]["n_nodes"],
                "n_edges": graph_pack["bundled"]["n_edges"],
                "south_frac": graph_pack["bundled"]["south_frac"],
                "loud_community": graph_pack["bundled"]["loud_community"],
                "loud_south_frac": graph_pack["bundled"]["loud_south_frac"],
                "merchant_cut": {str(k): v for k, v in graph_pack["bundled"]["merchant_cut"].items()},
            },
            "structural": {
                "n_communities": graph_pack["structural"]["n_communities"],
                "n_nodes": graph_pack["structural"]["n_nodes"],
                "n_edges": graph_pack["structural"]["n_edges"],
                "south_frac": graph_pack["structural"]["south_frac"],
                "merchant_cut": {str(k): v for k, v in graph_pack["structural"]["merchant_cut"].items()},
            },
            "own_vs_full": graph_pack["own_vs_full"],
            "level_set": {
                "cut": graph_pack["level_set"]["cut_name"],
                "merchant": {
                    k: v
                    for k, v in graph_pack["level_set"]["merchant"].items()
                    if k not in ("phis",)
                },
                "user": {
                    k: v
                    for k, v in graph_pack["level_set"]["user"].items()
                    if k not in ("phis",)
                },
                "jaccard_merchant_vs_user": graph_pack["level_set"]["jaccard_merchant_vs_user"],
                "jaccard_merchant_vs_south": graph_pack["level_set"]["jaccard_merchant_vs_south"],
                "jaccard_user_vs_south": graph_pack["level_set"]["jaccard_user_vs_south"],
                "order_oracle_south_n": graph_pack["level_set"]["order_oracle_south_n"],
                "read": graph_pack["level_set"]["read"],
            },
            "layers": {
                "merchant": {k: v for k, v in graph_pack["layers"]["merchant"].items() if k != "cut"},
                "user": {k: v for k, v in graph_pack["layers"]["user"].items() if k != "cut"},
                "jaccard_merchant_vs_user": graph_pack["layers"]["jaccard_merchant_vs_user"],
                "jaccard_merchant_vs_south": graph_pack["layers"]["jaccard_merchant_vs_south"],
                "jaccard_user_vs_south": graph_pack["layers"]["jaccard_user_vs_south"],
                "order_oracle_south_n": graph_pack["layers"]["order_oracle_south_n"],
                "read": graph_pack["layers"]["read"],
            },
        },
        "grains": grains,
        "fsds": {
            "rank": fsds["rank"],
            "selected_names": fsds["selected_names"],
            "selected": {g: [r["feature"] for r in rows] for g, rows in fsds["selected"].items()},
        },
        "unified_names": list(uni_new["names"]),
        "unified_groups": uni_new["groups"],
        "subsets": subsets,
        "shares": shares,
        "portraits": portraits,
        "loud_subset": hot,
        "logo_full": None
        if logo_full is None
        else {
            "pi_loss": logo_full["pi_loss"],
            "pi_po": logo_full["pi_po"],
            "pi_mmd": logo_full["pi_mmd"],
            "ratios": logo_full["ratios"],
            "global_action": logo_full["global_action"],
        },
        "logo_subset": None
        if logo_sub is None
        else {
            "pi_loss": logo_sub["pi_loss"],
            "pi_po": logo_sub["pi_po"],
            "pi_mmd": logo_sub["pi_mmd"],
            "ratios": logo_sub["ratios"],
        },
        "leakage": leakage,
        "n_selected": len(uni_new["names"]),
    }


# re-export for tests that want the LOGO mix helper
__all__ = [
    "assert_not_outcome",
    "characterize_subset",
    "cmean_x",
    "cmean_y",
    "complement_gap",
    "entity_profile",
    "fingerprint",
    "freeze_ref_stats",
    "fsds_rank_columns",
    "fsds_select",
    "gated_share",
    "CUT_LEVEL_SET",
    "CUT_BUNDLED",
    "CUT_STRUCTURAL",
    "GRAPH_PACKAGE",
    "graph_shift_cuts",
    "jaccard_masks",
    "layer_eval_on_orders",
    "lift_level_set_to_orders",
    "multilayer_level_set",
    "own_vs_full_summary",
    "node_bundle_scores",
    "grain_localization",
    "join_profile",
    "layer2_logo",
    "loud_subset",
    "node_localization",
    "pos_share",
    "run_pipeline",
    "slice_stream",
    "subset_shares",
    "subset_three_metrics",
    "three_metrics",
    "two_layer_ratios",
    "unify_to_merchant",
    "unify_to_order",
    "drop_group",
]
