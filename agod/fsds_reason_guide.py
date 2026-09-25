"""FSDS Reasoning Empowerment — adjust *next step* with ROI ledger.

Skill contract (fsds_reasoning_empowerment v1):
  X_node → FSDS attribution → importance(i) ∝ P(Y_i=1)
  → Guide {prune, expand, reorder, backtrack}
  → NetValue = Σ Y_i·V_task − Σ (1−Y_i)·C_node

Y_node(i) = 1 iff node i lies on the oracle optimal path.
This module is a *closed-loop scorecard*: synthetic RAP/ToT trees with
known oracle paths, FSDS-guided search vs score-greedy baseline.

Claims stay ledger-level (simulation Δ), not causal $ ROI.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold


# ---------------------------------------------------------------------------
# Feature schema (skill §4) — compact numeric blocks
# ---------------------------------------------------------------------------

STRUCT = ["depth", "num_children", "num_siblings", "path_length", "is_leaf", "is_root"]
SEMANTIC = ["state_length", "has_numbers", "has_operators", "has_keywords", "state_entropy"]
SCORE = ["score", "score_gap", "score_rank", "score_trend"]
SEARCH = ["num_expanded", "num_pruned", "num_backtracks", "prune_rate", "backtrack_rate"]
COST = ["llm_calls", "tokens", "latency", "cost"]
TIME = ["timestamp", "hour", "day_of_week", "batch_id"]
ENV = ["version", "model", "strategy", "max_depth", "beam_width"]

FEATURE_NAMES: List[str] = STRUCT + SEMANTIC + SCORE + SEARCH + COST + TIME + ENV
BLOCKS = {
    "struct": STRUCT,
    "semantic": SEMANTIC,
    "score": SCORE,
    "search": SEARCH,
    "cost": COST,
    "time": TIME,
    "env": ENV,
}


@dataclass
class ReasonNode:
    """One search node with skill-aligned features + oracle label."""

    node_id: int
    parent_id: Optional[int]
    depth: int
    score: float
    x: np.ndarray
    y_opt: int  # 1 if on oracle optimal path
    children: List[int] = field(default_factory=list)
    expanded: bool = False
    pruned: bool = False


@dataclass
class Trajectory:
    """One RAP/ToT episode."""

    nodes: Dict[int, ReasonNode]
    root_id: int
    opt_path: List[int]
    task_id: int


# ---------------------------------------------------------------------------
# Synthetic trajectories (known Y)
# ---------------------------------------------------------------------------

def _node_features(
    *,
    depth: int,
    n_sib: int,
    n_child: int,
    is_leaf: int,
    is_root: int,
    score: float,
    score_gap: float,
    score_rank: float,
    score_trend: float,
    path_length: int,
    state_length: float,
    has_numbers: float,
    has_operators: float,
    has_keywords: float,
    state_entropy: float,
    n_exp: float,
    n_prn: float,
    n_bt: float,
    llm_calls: float,
    tokens: float,
    latency: float,
    cost: float,
    timestamp: float,
    hour: float,
    dow: float,
    batch_id: float,
    version: float,
    model: float,
    strategy: float,
    max_depth: float,
    beam: float,
) -> np.ndarray:
    prune_rate = n_prn / max(n_exp + n_prn, 1.0)
    bt_rate = n_bt / max(n_exp + n_bt, 1.0)
    vals = [
        float(depth),
        float(n_child),
        float(n_sib),
        float(path_length),
        float(is_leaf),
        float(is_root),
        state_length,
        has_numbers,
        has_operators,
        has_keywords,
        state_entropy,
        score,
        score_gap,
        score_rank,
        score_trend,
        n_exp,
        n_prn,
        n_bt,
        prune_rate,
        bt_rate,
        llm_calls,
        tokens,
        latency,
        cost,
        timestamp,
        hour,
        dow,
        batch_id,
        version,
        model,
        strategy,
        max_depth,
        beam,
    ]
    assert len(vals) == len(FEATURE_NAMES)
    return np.asarray(vals, dtype=float)


def simulate_trajectory(
    task_id: int,
    *,
    seed: int,
    max_depth: int = 4,
    branching: int = 3,
    beam_width: int = 3,
) -> Trajectory:
    """Build a tree; one random root→leaf path is oracle-optimal.

    Optimal nodes get a *latent* quality bump that leaks into score +
    semantic features (so FSDS can learn X→Y); baseline score-greedy
    still sees noise and decoys.
    """
    rng = np.random.default_rng(seed + 17 * task_id)
    nodes: Dict[int, ReasonNode] = {}
    nid = 0

    def make(
        parent: Optional[int],
        depth: int,
        on_opt: bool,
        parent_score: float,
        n_sib: int,
    ) -> int:
        nonlocal nid
        mid = nid
        nid += 1
        # latent quality: opt path higher; decoys can look good by chance
        latent = (1.2 if on_opt else 0.0) + rng.normal(0, 0.55)
        score = 0.55 * latent + 0.45 * rng.normal(0, 1.0)
        score_gap = score - parent_score
        # semantic leakage toward opt
        has_kw = float(rng.random() < (0.75 if on_opt else 0.25))
        has_op = float(rng.random() < (0.65 if on_opt else 0.35))
        has_num = float(rng.random() < (0.55 if on_opt else 0.40))
        state_len = float(8 + depth * 3 + (4 if on_opt else 0) + rng.integers(0, 5))
        entropy = float(np.clip(1.4 - 0.25 * depth + (0.3 if on_opt else -0.1) + rng.normal(0, 0.2), 0.1, 2.5))
        trend = float(np.tanh(score_gap + (0.4 if on_opt else -0.15)))
        tokens = float(40 + 18 * depth + rng.integers(0, 25))
        calls = float(1 + (depth > 0))
        latency = float(0.2 + 0.08 * depth + rng.random() * 0.1)
        cost = tokens * 0.002 + calls * 0.01
        x = _node_features(
            depth=depth,
            n_sib=n_sib,
            n_child=0,
            is_leaf=0,
            is_root=int(parent is None),
            score=float(score),
            score_gap=float(score_gap),
            score_rank=0.0,  # filled after siblings known
            score_trend=trend,
            path_length=depth,
            state_length=state_len,
            has_numbers=has_num,
            has_operators=has_op,
            has_keywords=has_kw,
            state_entropy=entropy,
            n_exp=0.0,
            n_prn=0.0,
            n_bt=0.0,
            llm_calls=calls,
            tokens=tokens,
            latency=latency,
            cost=cost,
            timestamp=float(task_id * 100 + mid),
            hour=float(rng.integers(0, 24)),
            dow=float(rng.integers(0, 7)),
            batch_id=float(task_id % 8),
            version=1.0,
            model=float(rng.integers(0, 3)),
            strategy=float(rng.integers(0, 2)),
            max_depth=float(max_depth),
            beam=float(beam_width),
        )
        nodes[mid] = ReasonNode(
            node_id=mid,
            parent_id=parent,
            depth=depth,
            score=float(score),
            x=x,
            y_opt=int(on_opt),
        )
        if parent is not None:
            nodes[parent].children.append(mid)
        return mid

    root = make(None, 0, True, 0.0, 0)
    # choose which child index is on the opt path at each depth
    opt_choices = [int(rng.integers(0, branching)) for _ in range(max_depth)]
    frontier = [(root, True)]
    for d in range(1, max_depth + 1):
        nxt: List[Tuple[int, bool]] = []
        for pid, parent_on_opt in frontier:
            kids: List[int] = []
            for b in range(branching):
                on_opt = bool(parent_on_opt and b == opt_choices[d - 1])
                kid = make(pid, d, on_opt, nodes[pid].score, branching - 1)
                kids.append(kid)
            # score ranks among siblings
            order = sorted(kids, key=lambda i: nodes[i].score, reverse=True)
            for rank, i in enumerate(order):
                nodes[i].x[FEATURE_NAMES.index("score_rank")] = float(rank)
            nodes[pid].x[FEATURE_NAMES.index("num_children")] = float(len(kids))
            if d == max_depth:
                for i in kids:
                    nodes[i].x[FEATURE_NAMES.index("is_leaf")] = 1.0
            else:
                nxt.extend((i, bool(nodes[i].y_opt)) for i in kids)
        frontier = nxt

    opt_path = [i for i, n in nodes.items() if n.y_opt]
    opt_path.sort(key=lambda i: nodes[i].depth)
    return Trajectory(nodes=nodes, root_id=root, opt_path=opt_path, task_id=task_id)


def stack_xy(trajs: Sequence[Trajectory]) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    xs, ys = [], []
    for tr in trajs:
        for n in tr.nodes.values():
            xs.append(n.x)
            ys.append(n.y_opt)
    X = np.vstack(xs)
    y = np.asarray(ys, dtype=int)
    return X, y, list(FEATURE_NAMES)


# ---------------------------------------------------------------------------
# FSDS four components (+ overlap router)
# ---------------------------------------------------------------------------

def _rbf_mmd2(X0: np.ndarray, X1: np.ndarray, *, max_n: int, rng: np.random.Generator) -> float:
    def _sub(X: np.ndarray) -> np.ndarray:
        if len(X) <= max_n:
            return X
        idx = rng.choice(len(X), size=max_n, replace=False)
        return X[idx]

    a, b = _sub(X0), _sub(X1)
    if len(a) < 2 or len(b) < 2:
        return 0.0
    Z = np.vstack([a, b])
    # median heuristic on pairwise dist sample
    n = min(80, len(Z))
    s = Z[rng.choice(len(Z), size=n, replace=False)]
    d2 = ((s[:, None, :] - s[None, :, :]) ** 2).sum(-1)
    med = float(np.median(d2[d2 > 0])) if np.any(d2 > 0) else 1.0
    gamma = 1.0 / max(med, 1e-8)

    def k(U, V):
        return np.exp(-gamma * ((U[:, None, :] - V[None, :, :]) ** 2).sum(-1))

    kaa = k(a, a)
    kbb = k(b, b)
    kab = k(a, b)
    na, nb = len(a), len(b)
    # unbiased-ish diagonal-zeroed
    np.fill_diagonal(kaa, 0.0)
    np.fill_diagonal(kbb, 0.0)
    return float(kaa.sum() / (na * (na - 1)) + kbb.sum() / (nb * (nb - 1)) - 2.0 * kab.mean())


def overlap_score(X: np.ndarray, y: np.ndarray, *, seed: int = 0) -> float:
    """Propensity-style overlap proxy: 2·min(p,1−p) mean on RF P(Y=1|X).

    Low overlap → sparse support for Y=1 under X → treat as covariate-like
    separation (MMD/RF-VIMP). High overlap → concept-style PO/PermuCATE.
    """
    if len(np.unique(y)) < 2:
        return 0.0
    clf = RandomForestClassifier(
        n_estimators=40, max_depth=6, min_samples_leaf=3, random_state=seed, n_jobs=1
    )
    clf.fit(X, y)
    p = clf.predict_proba(X)[:, 1]
    return float(np.mean(2.0 * np.minimum(p, 1.0 - p)))


def mmd_loco_importance(
    X: np.ndarray, y: np.ndarray, names: Sequence[str], *, seed: int = 0, max_n: int = 120
) -> Dict[str, float]:
    """MMD-LOCO: Importance(j) = MMD_full − MMD_{−j} between Y=0 and Y=1 supports."""
    rng = np.random.default_rng(seed)
    X0, X1 = X[y == 0], X[y == 1]
    if len(X0) < 2 or len(X1) < 2:
        return {n: 0.0 for n in names}
    full = _rbf_mmd2(X0, X1, max_n=max_n, rng=rng)
    out: Dict[str, float] = {}
    for j, n in enumerate(names):
        keep = [i for i in range(X.shape[1]) if i != j]
        d = _rbf_mmd2(X0[:, keep], X1[:, keep], max_n=max_n, rng=rng) if keep else 0.0
        out[n] = float(full - d)
    return out


def rf_binary_vimp(
    X: np.ndarray, y: np.ndarray, names: Sequence[str], *, seed: int = 0
) -> Dict[str, Any]:
    """RF-Binary-VIMP: classify optimal-path vs not; Gini importance + AUC."""
    if len(np.unique(y)) < 2:
        return {"vimp": {n: 0.0 for n in names}, "auc": 0.5, "proba_model": None}
    clf = RandomForestClassifier(
        n_estimators=80, max_depth=8, min_samples_leaf=2, random_state=seed, n_jobs=1
    )
    # OOF AUC
    n_splits = 3 if len(y) >= 60 else 2
    oof = np.zeros(len(y), dtype=float)
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for fold, (tr, te) in enumerate(cv.split(X, y)):
        c = RandomForestClassifier(
            n_estimators=60, max_depth=8, min_samples_leaf=2, random_state=seed + fold, n_jobs=1
        )
        c.fit(X[tr], y[tr])
        oof[te] = c.predict_proba(X[te])[:, 1]
    auc = float(roc_auc_score(y, oof)) if len(np.unique(y)) > 1 else 0.5
    clf.fit(X, y)
    vimp = {n: float(v) for n, v in zip(names, clf.feature_importances_)}
    return {"vimp": vimp, "auc": auc, "proba_model": clf, "oof": oof}


def po_risk_loco(
    X: np.ndarray,
    y: np.ndarray,
    names: Sequence[str],
    *,
    seed: int = 0,
    n_trees: int = 30,
) -> Dict[str, float]:
    """PO-risk-LOCO: Importance(j) = PO-risk_full − PO-risk_{−j}.

    Treatment proxy T = 1{score_gap > median} (expand-worthy local signal);
    outcome Y = on optimal path. Matches skill §5.3 spirit on search nodes.
    """
    gap_i = names.index("score_gap") if "score_gap" in names else 0
    T = (X[:, gap_i] > np.median(X[:, gap_i])).astype(int)
    if len(np.unique(T)) < 2 or len(np.unique(y)) < 2:
        return {n: 0.0 for n in names}

    def _risk(Xc: np.ndarray) -> float:
        Yf = y.astype(float)
        m_hat = np.zeros(len(Yf))
        e_hat = np.zeros(len(Yf))
        n_splits = 3 if len(Yf) >= 60 else 2
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        for fold, (tr, te) in enumerate(cv.split(Xc, T)):
            m = RandomForestRegressor(
                n_estimators=n_trees, max_depth=5, min_samples_leaf=4, random_state=seed + fold, n_jobs=1
            )
            e = RandomForestClassifier(
                n_estimators=n_trees, max_depth=5, min_samples_leaf=4, random_state=seed + 40 + fold, n_jobs=1
            )
            m.fit(Xc[tr], Yf[tr])
            e.fit(Xc[tr], T[tr])
            m_hat[te] = m.predict(Xc[te])
            e_hat[te] = e.predict_proba(Xc[te])[:, 1]
        e_hat = np.clip(e_hat, 0.05, 0.95)
        po = (Yf - m_hat) * (T - e_hat)
        tau = RandomForestRegressor(
            n_estimators=max(40, n_trees), max_depth=6, min_samples_leaf=4, random_state=seed + 7, n_jobs=1
        )
        tau.fit(Xc, po)
        return float(np.mean(tau.predict(Xc) ** 2))

    full = _risk(X)
    out: Dict[str, float] = {}
    # cheap LOCO: only top-score + semantic + struct dims to keep scorecard fast
    focus = [
        i
        for i, n in enumerate(names)
        if n in STRUCT + SEMANTIC + SCORE or i < 12
    ]
    focus = sorted(set(focus))
    for j in range(len(names)):
        if j not in focus:
            out[names[j]] = 0.0
            continue
        keep = [i for i in range(X.shape[1]) if i != j]
        out[names[j]] = float(full - _risk(X[:, keep]))
    return out


def permu_cate(
    X: np.ndarray,
    y: np.ndarray,
    names: Sequence[str],
    *,
    seed: int = 0,
    n_perm: int = 3,
) -> Dict[str, float]:
    """PermuCATE-lite: E[(τ(X)−τ(X_perm_j))²] with T=1{score_gap>med}."""
    gap_i = names.index("score_gap") if "score_gap" in names else 0
    T = (X[:, gap_i] > np.median(X[:, gap_i])).astype(int)
    if len(np.unique(T)) < 2:
        return {n: 0.0 for n in names}
    Yf = y.astype(float)
    # simple DR residual then regress
    m = RandomForestRegressor(n_estimators=40, max_depth=5, min_samples_leaf=4, random_state=seed, n_jobs=1)
    e = RandomForestClassifier(n_estimators=40, max_depth=5, min_samples_leaf=4, random_state=seed + 1, n_jobs=1)
    m.fit(X, Yf)
    e.fit(X, T)
    e_hat = np.clip(e.predict_proba(X)[:, 1], 0.05, 0.95)
    po = (Yf - m.predict(X)) * (T - e_hat)
    tau = RandomForestRegressor(n_estimators=50, max_depth=6, min_samples_leaf=3, random_state=seed + 2, n_jobs=1)
    tau.fit(X, po)
    base = tau.predict(X)
    rng = np.random.default_rng(seed)
    out: Dict[str, float] = {}
    focus = [i for i, n in enumerate(names) if n in STRUCT + SEMANTIC + SCORE]
    for j, n in enumerate(names):
        if j not in focus:
            out[n] = 0.0
            continue
        diffs = []
        for _ in range(n_perm):
            Xp = X.copy()
            Xp[:, j] = rng.permutation(Xp[:, j])
            diffs.append(float(np.mean((base - tau.predict(Xp)) ** 2)))
        out[n] = float(np.mean(diffs))
    return out


def fsds_attribute(
    X: np.ndarray, y: np.ndarray, names: Sequence[str], *, seed: int = 0
) -> Dict[str, Any]:
    """Run overlap router + selected FSDS pair; fuse to feature importance."""
    ov = overlap_score(X, y, seed=seed)
    rf = rf_binary_vimp(X, y, names, seed=seed)
    regime = "covariate_shift" if ov < 0.3 else "concept_drift"
    if regime == "covariate_shift":
        mmd = mmd_loco_importance(X, y, names, seed=seed)
        po = {n: 0.0 for n in names}
        pc = {n: 0.0 for n in names}
        pair = ("MMD-LOCO", "RF-Binary-VIMP")
        # fuse: z-score then mean
        fuse_src = [mmd, rf["vimp"]]
    else:
        mmd = {n: 0.0 for n in names}
        po = po_risk_loco(X, y, names, seed=seed)
        pc = permu_cate(X, y, names, seed=seed)
        pair = ("PO-risk-LOCO", "PermuCATE")
        fuse_src = [po, pc, rf["vimp"]]

    def _z(d: Dict[str, float]) -> Dict[str, float]:
        v = np.asarray([d[n] for n in names], dtype=float)
        s = float(v.std()) + 1e-8
        return {n: float((d[n] - v.mean()) / s) for n in names}

    zs = [_z(d) for d in fuse_src]
    fused = {n: float(np.mean([z[n] for z in zs])) for n in names}
    # shift to non-negative mass
    m = min(fused.values())
    fused = {n: fused[n] - m for n in names}
    s = sum(fused.values()) + 1e-12
    fused = {n: fused[n] / s for n in names}
    top = sorted(fused.items(), key=lambda kv: -kv[1])[:8]
    return {
        "overlap": ov,
        "regime": regime,
        "components": list(pair) + (["RF-Binary-VIMP"] if "RF-Binary-VIMP" not in pair else []),
        "mmd_loco": mmd,
        "rf_vimp": rf["vimp"],
        "po_loco": po,
        "permu_cate": pc,
        "fused_importance": fused,
        "top_features": top,
        "auc_y": rf["auc"],
        "proba_model": rf["proba_model"],
    }


def node_importance(model: Any, x: np.ndarray) -> float:
    """importance(i) ∝ P(Y_i=1 | X_i)."""
    if model is None:
        return 0.5
    return float(model.predict_proba(x.reshape(1, -1))[0, 1])


# ---------------------------------------------------------------------------
# Guide: prune / expand / reorder / backtrack  (skill §6)
# ---------------------------------------------------------------------------

@dataclass
class GuideConfig:
    """Thresholds relative to train importance distribution when possible."""

    imp_expand: float = 0.35
    imp_prune: float = 0.12
    gap_expand: float = -0.25
    gap_prune: float = -0.85
    max_depth: int = 4
    anomaly_z: float = -2.0
    consec_fail: int = 3
    beam: int = 2
    # never prune the best score-rank sibling if imp above floor
    keep_top_rank_floor: float = 0.08


def guide_decision(
    node: ReasonNode,
    *,
    importance: float,
    cfg: GuideConfig,
    consec_failures: int = 0,
    imp_mean: float = 0.5,
    imp_std: float = 0.2,
) -> Dict[str, Any]:
    """Decide next-step action for one frontier node."""
    gap = float(node.x[FEATURE_NAMES.index("score_gap")])
    trend = float(node.x[FEATURE_NAMES.index("score_trend")])
    rank = float(node.x[FEATURE_NAMES.index("score_rank")])
    depth = node.depth
    z = (importance - imp_mean) / (imp_std + 1e-8)
    anomaly = z < cfg.anomaly_z and trend < -0.7

    # hard depth cap
    if depth > cfg.max_depth:
        return {
            "action": "prune",
            "priority": 0.9,
            "reason": f"depth>{cfg.max_depth}",
            "expected_gain": "cost↓",
        }

    # prune only clear losers (low imp AND bad gap/trend); keep top-rank if not terrible
    clear_loser = importance < cfg.imp_prune and (gap < cfg.gap_prune or trend < -0.5)
    if clear_loser and not (rank <= 0.0 and importance >= cfg.keep_top_rank_floor):
        return {
            "action": "prune",
            "priority": 0.9,
            "reason": f"imp={importance:.2f} gap={gap:.2f} trend={trend:.2f} depth={depth}",
            "expected_gain": "cost↓",
        }
    if anomaly or (consec_failures >= cfg.consec_fail and importance < imp_mean):
        return {
            "action": "backtrack",
            "priority": 0.85,
            "reason": f"anomaly_z={z:.2f} fails={consec_failures}",
            "expected_gain": "fail↓",
        }
    if importance >= cfg.imp_expand and gap >= cfg.gap_expand and depth < cfg.max_depth and trend > -0.25:
        return {
            "action": "expand",
            "priority": 0.8,
            "reason": f"imp={importance:.2f} gap={gap:.2f} trend={trend:.2f}",
            "expected_gain": "revenue↑",
        }
    # reorder default — rank by importance-weighted gap
    pri = importance * (1.0 + max(gap, -0.2)) / max(depth, 1)
    return {
        "action": "reorder",
        "priority": 0.7,
        "reason": f"priority={pri:.3f}",
        "expected_gain": "efficiency↑",
        "sort_key": float(pri + 0.15 * importance),
    }


# ---------------------------------------------------------------------------
# Search policies: baseline score-greedy vs FSDS-guided
# ---------------------------------------------------------------------------

@dataclass
class SearchResult:
    success: bool
    path: List[int]
    nodes_visited: int
    tokens: float
    cost: float
    llm_calls: float
    n_expand: int
    n_prune: int
    n_backtrack: int
    n_reorder: int
    actions: List[Dict[str, Any]] = field(default_factory=list)


def _leaf_ids(tr: Trajectory) -> List[int]:
    return [i for i, n in tr.nodes.items() if not n.children]


def run_score_greedy(tr: Trajectory, *, beam: int = 2) -> SearchResult:
    """Baseline: always expand top-`beam` children by raw score."""
    cur = [tr.root_id]
    path = [tr.root_id]
    tokens = float(tr.nodes[tr.root_id].x[FEATURE_NAMES.index("tokens")])
    cost = float(tr.nodes[tr.root_id].x[FEATURE_NAMES.index("cost")])
    calls = float(tr.nodes[tr.root_id].x[FEATURE_NAMES.index("llm_calls")])
    visited = 1
    n_exp = 0
    while True:
        kids = []
        for p in cur:
            kids.extend(tr.nodes[p].children)
        if not kids:
            break
        kids = sorted(kids, key=lambda i: tr.nodes[i].score, reverse=True)[:beam]
        n_exp += len(kids)
        visited += len(kids)
        for i in kids:
            tokens += float(tr.nodes[i].x[FEATURE_NAMES.index("tokens")])
            cost += float(tr.nodes[i].x[FEATURE_NAMES.index("cost")])
            calls += float(tr.nodes[i].x[FEATURE_NAMES.index("llm_calls")])
        # pick best among beam as path tip
        best = max(kids, key=lambda i: tr.nodes[i].score)
        path.append(best)
        cur = kids
        if not tr.nodes[best].children:
            break
    success = all(tr.nodes[i].y_opt for i in path) and path[-1] in _leaf_ids(tr) and tr.nodes[path[-1]].y_opt == 1
    # stricter: final leaf on opt path
    success = bool(tr.nodes[path[-1]].y_opt == 1 and not tr.nodes[path[-1]].children)
    return SearchResult(
        success=success,
        path=path,
        nodes_visited=visited,
        tokens=tokens,
        cost=cost,
        llm_calls=calls,
        n_expand=n_exp,
        n_prune=0,
        n_backtrack=0,
        n_reorder=0,
    )


def run_fsds_guided(
    tr: Trajectory,
    model: Any,
    *,
    cfg: GuideConfig,
    imp_mean: float,
    imp_std: float,
) -> SearchResult:
    """FSDS-guided next-step: importance → prune/expand/reorder/backtrack.

    Beam is ranked by P(Y=1|X), not raw score. Prune drops clear losers;
    tip advances along max-importance survivor; if the tip's children are
    all pruned, backtrack to the next beam survivor at the same depth.
    """
    path = [tr.root_id]
    tokens = float(tr.nodes[tr.root_id].x[FEATURE_NAMES.index("tokens")])
    cost = float(tr.nodes[tr.root_id].x[FEATURE_NAMES.index("cost")])
    calls = float(tr.nodes[tr.root_id].x[FEATURE_NAMES.index("llm_calls")])
    visited = 1
    n_exp = n_prn = n_bt = n_ord = 0
    actions: List[Dict[str, Any]] = []
    consec_fail = 0
    # beam state: list of node ids at current depth
    frontier = list(tr.nodes[tr.root_id].children)

    while frontier:
        scored = []
        for i in frontier:
            imp = node_importance(model, tr.nodes[i].x)
            dec = guide_decision(
                tr.nodes[i],
                importance=imp,
                cfg=cfg,
                consec_failures=consec_fail,
                imp_mean=imp_mean,
                imp_std=imp_std,
            )
            scored.append((i, imp, dec))
            actions.append({"node": i, "imp": imp, **dec})

        kept = []
        for i, imp, dec in scored:
            if dec["action"] == "prune":
                n_prn += 1
                continue
            if dec["action"] == "backtrack":
                n_bt += 1
                continue
            kept.append((i, imp, dec))

        # safety: if everything pruned, keep top-1 by importance (avoid empty search)
        if not kept:
            n_bt += 1
            scored.sort(key=lambda t: t[1], reverse=True)
            kept = [scored[0]]
            kept[0] = (kept[0][0], kept[0][1], {**kept[0][2], "action": "backtrack", "reason": "rescue_top_imp"})

        def _key(t):
            i, imp, dec = t
            return float(dec.get("sort_key", imp) if dec["action"] == "reorder" else imp)

        kept.sort(key=_key, reverse=True)
        n_ord += 1
        beam = kept[: cfg.beam]
        n_exp += len(beam)
        visited += len(beam)
        for i, _, _ in beam:
            tokens += float(tr.nodes[i].x[FEATURE_NAMES.index("tokens")])
            cost += float(tr.nodes[i].x[FEATURE_NAMES.index("cost")])
            calls += float(tr.nodes[i].x[FEATURE_NAMES.index("llm_calls")])

        tip = beam[0][0]
        path.append(tip)
        if tr.nodes[tip].y_opt:
            consec_fail = 0
        else:
            consec_fail += 1

        # next frontier = children of all beam nodes (importance-ranked widen)
        nxt: List[int] = []
        for i, _, _ in beam:
            nxt.extend(tr.nodes[i].children)
        if not nxt:
            break
        # de-dupe preserve order
        seen = set()
        frontier = []
        for i in nxt:
            if i not in seen:
                seen.add(i)
                frontier.append(i)

    success = bool(tr.nodes[path[-1]].y_opt == 1 and not tr.nodes[path[-1]].children)
    return SearchResult(
        success=success,
        path=path,
        nodes_visited=visited,
        tokens=tokens,
        cost=cost,
        llm_calls=calls,
        n_expand=n_exp,
        n_prune=n_prn,
        n_backtrack=n_bt,
        n_reorder=n_ord,
        actions=actions,
    )


# ---------------------------------------------------------------------------
# Economic ledger (skill §7)
# ---------------------------------------------------------------------------

@dataclass
class EconParams:
    v_task: float = 5.0  # value of a successful task
    c_token: float = 0.002
    c_call: float = 0.01
    v_auc: float = 2.0  # valuation of +1.0 AUC in selection model
    cost_fsds: float = 0.15  # amortized attribution cost per task batch share


def net_value_path(tr: Trajectory, res: SearchResult, *, p: EconParams) -> Dict[str, float]:
    """NetValue = Σ Y_i·V_task/|path| style + task success − node costs.

    Ledger decomposition (simulation, not causal $):
      revenue = 1{success} · V_task
      cost_nodes = Σ_{visited, y=0} C_node  +  token/call tariffs
    """
    rev = float(res.success) * p.v_task
    waste = 0.0
    for i in res.path:
        if tr.nodes[i].y_opt == 0:
            waste += float(tr.nodes[i].x[FEATURE_NAMES.index("cost")])
    # also charge full token/call stream
    stream = res.tokens * p.c_token + res.llm_calls * p.c_call
    # prefer not double-count: stream already includes path costs roughly;
    # keep waste as *extra* penalty for off-opt steps
    cost = stream + 0.5 * waste
    return {
        "revenue": rev,
        "cost": cost,
        "waste_off_opt": waste,
        "net": rev - cost,
        "success": float(res.success),
        "tokens": res.tokens,
        "llm_calls": res.llm_calls,
    }


def incremental_roi(
    base: Dict[str, float],
    guided: Dict[str, float],
    *,
    auc_base: float,
    auc_fsds: float,
    p: EconParams,
    n_tasks: int,
) -> Dict[str, float]:
    """Δ ledger + ROI vs score-greedy baseline."""
    d_rev = guided["revenue"] - base["revenue"]
    d_cost = base["cost"] - guided["cost"]  # positive = savings
    d_auc = (auc_fsds - auc_base) * p.v_auc
    fsds_cost = p.cost_fsds * max(n_tasks, 1) / max(n_tasks, 1)  # per-task share already
    # when aggregating outside, caller multiplies; here per-task mean inputs
    delta = d_rev + d_cost + d_auc - p.cost_fsds
    invest = p.cost_fsds + max(guided["cost"] - base["cost"], 0.0)
    roi = delta / invest if invest > 1e-9 else float("inf")
    return {
        "delta_revenue": float(d_rev),
        "delta_cost_savings": float(d_cost),
        "delta_auc_value": float(d_auc),
        "cost_fsds": float(p.cost_fsds),
        "net_incremental": float(delta),
        "roi": float(roi),
        "auc_base": float(auc_base),
        "auc_fsds": float(auc_fsds),
    }


def run_suite(
    *,
    n_train: int = 40,
    n_test: int = 30,
    seed: int = 0,
    max_depth: int = 4,
    branching: int = 3,
    econ: Optional[EconParams] = None,
) -> Dict[str, Any]:
    """Train FSDS on train trajectories; compare guided vs greedy on test."""
    econ = econ or EconParams()
    train = [simulate_trajectory(i, seed=seed, max_depth=max_depth, branching=branching) for i in range(n_train)]
    test = [
        simulate_trajectory(10_000 + i, seed=seed + 1, max_depth=max_depth, branching=branching)
        for i in range(n_test)
    ]
    X, y, names = stack_xy(train)
    attr = fsds_attribute(X, y, names, seed=seed)
    model = attr["proba_model"]
    imps = [node_importance(model, n.x) for tr in train for n in tr.nodes.values()]
    imp_mean = float(np.mean(imps)) if imps else 0.5
    imp_std = float(np.std(imps)) if imps else 0.2
    # calibrate prune/expand to train importance percentiles (skill: relative thresholds)
    if imps:
        imp_prune = float(np.percentile(imps, 20))
        imp_expand = float(np.percentile(imps, 55))
    else:
        imp_prune, imp_expand = 0.12, 0.35
    cfg = GuideConfig(
        max_depth=max_depth,
        beam=2,
        imp_prune=imp_prune,
        imp_expand=imp_expand,
        keep_top_rank_floor=float(np.percentile(imps, 10)) if imps else 0.08,
    )

    # AUC of score-only baseline for Y (for ΔAUC ledger)
    scores = np.asarray([n.score for tr in train for n in tr.nodes.values()])
    try:
        auc_base = float(roc_auc_score(y, scores))
    except ValueError:
        auc_base = 0.5

    base_nets, guid_nets = [], []
    base_succ, guid_succ = [], []
    action_counts = {"prune": 0, "expand": 0, "reorder": 0, "backtrack": 0}
    for tr in test:
        b = run_score_greedy(tr, beam=cfg.beam)
        g = run_fsds_guided(tr, model, cfg=cfg, imp_mean=imp_mean, imp_std=imp_std)
        base_nets.append(net_value_path(tr, b, p=econ))
        guid_nets.append(net_value_path(tr, g, p=econ))
        base_succ.append(b.success)
        guid_succ.append(g.success)
        for a in g.actions:
            action_counts[a["action"]] = action_counts.get(a["action"], 0) + 1

    def _mean(rows: List[Dict[str, float]], k: str) -> float:
        return float(np.mean([r[k] for r in rows])) if rows else 0.0

    base_agg = {k: _mean(base_nets, k) for k in ("revenue", "cost", "net", "tokens", "llm_calls", "waste_off_opt")}
    guid_agg = {k: _mean(guid_nets, k) for k in ("revenue", "cost", "net", "tokens", "llm_calls", "waste_off_opt")}
    base_agg["success_rate"] = float(np.mean(base_succ))
    guid_agg["success_rate"] = float(np.mean(guid_succ))

    inc = incremental_roi(
        base_agg,
        guid_agg,
        auc_base=auc_base,
        auc_fsds=float(attr["auc_y"]),
        p=econ,
        n_tasks=n_test,
    )
    # scale revenue delta using success rates explicitly for clarity
    inc["delta_success"] = guid_agg["success_rate"] - base_agg["success_rate"]
    inc["delta_tokens"] = base_agg["tokens"] - guid_agg["tokens"]
    inc["delta_net_per_task"] = guid_agg["net"] - base_agg["net"]

    return {
        "attr": {k: v for k, v in attr.items() if k != "proba_model"},
        "baseline": base_agg,
        "guided": guid_agg,
        "incremental": inc,
        "action_counts": action_counts,
        "n_train": n_train,
        "n_test": n_test,
        "econ": {
            "v_task": econ.v_task,
            "c_token": econ.c_token,
            "c_call": econ.c_call,
            "v_auc": econ.v_auc,
            "cost_fsds": econ.cost_fsds,
        },
        "feature_names": names,
        "top_features": attr["top_features"],
        "regime": attr["regime"],
        "overlap": attr["overlap"],
        "auc_y": attr["auc_y"],
    }


def justify_bullets(suite: Dict[str, Any]) -> List[str]:
    """Short justification lines for the scorecard (no overclaim)."""
    inc = suite["incremental"]
    b, g = suite["baseline"], suite["guided"]
    lines = [
        f"Overlap={suite['overlap']:.3f} → regime={suite['regime']} (skill §5.5 router).",
        f"FSDS P(Y=1|X) AUC={suite['auc_y']:.3f} vs score-only AUC={inc['auc_base']:.3f} on train nodes.",
        f"Test success: greedy {b['success_rate']:.1%} → guided {g['success_rate']:.1%} (Δ={inc['delta_success']:+.1%}).",
        f"Tokens/task: {b['tokens']:.0f} → {g['tokens']:.0f} (Δ saved={inc['delta_tokens']:+.0f}).",
        f"Net/task ledger: {b['net']:.3f} → {g['net']:.3f} (Δ={inc['delta_net_per_task']:+.3f}).",
        f"Incremental NetValue≈{inc['net_incremental']:.3f}/task; ROI≈{inc['roi']:.2f} (sim ledger, not causal $).",
        "Y=1{node on oracle path} is the direct driver of task success → revenue; Guide maps importance→prune/expand.",
    ]
    return lines
