"""PO-risk Var(u): CF leaf-debiased / BCF / bootstrap. ~150 lines."""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Optional

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LogisticRegression

from as_dta_aigc import SEED
from as_dta_permucate_magvar import permucate_unit_po

MIN_LEAF_MERGE_SMALL, MIN_LEAF_MERGE_LARGE, MIN_LEAF_MERGE_N_SWITCH = 30, 50, 200


def resolve_min_leaf_merge_n(n: int, override: Optional[int] = None) -> int:
    if override is not None:
        return int(max(2, override))
    return MIN_LEAF_MERGE_LARGE if int(n) >= MIN_LEAF_MERGE_N_SWITCH else MIN_LEAF_MERGE_SMALL


def tree_parents(tree) -> np.ndarray:
    """parent[v]; parent[root] = -1."""
    tr = tree.tree_ if hasattr(tree, "tree_") else tree
    left, right = np.asarray(tr.children_left), np.asarray(tr.children_right)
    parent = np.full(len(left), -1, dtype=int)
    for p in range(len(left)):
        if left[p] >= 0:
            parent[left[p]] = p
        if right[p] >= 0:
            parent[right[p]] = p
    return parent


def merge_leaf_u(tree, u: np.ndarray, leaf_ids: np.ndarray, min_n: int = 30):
    """Climb parent until n_node >= min_n; mean(u) in cluster. → (u_merged, merge_ids)."""
    u, leaf_ids = np.asarray(u, float).ravel(), np.asarray(leaf_ids, int).ravel()
    tr = tree.tree_ if hasattr(tree, "tree_") else tree
    parent, n_node, min_n = tree_parents(tree), np.asarray(tr.n_node_samples, float), int(max(1, min_n))
    merge_ids = np.empty(len(leaf_ids), dtype=int)
    for i, lid in enumerate(leaf_ids):
        v = int(lid)
        while parent[v] >= 0 and n_node[v] < min_n:
            v = int(parent[v])
        merge_ids[i] = v
    u_merged = np.empty_like(u)
    for gid in np.unique(merge_ids):
        m = merge_ids == gid
        u_merged[m] = u[m].mean()
    return u_merged, merge_ids


def causal_forest_var_u(
    X, Y, W, phi, *, n_estimators=80, max_depth=5, min_samples_leaf=None,
    seed=SEED, leaf_merge=True, min_leaf_merge_n=None,
):
    """V_debiased = V_naive - B_leaf; Var(u) across trees after leaf merge."""
    from econml.dml import CausalForestDML

    X, Y, W, phi = np.asarray(X, float), np.asarray(Y, float).ravel(), np.asarray(W, int).ravel(), np.asarray(phi, float).ravel()
    n = len(phi)
    merge_n = resolve_min_leaf_merge_n(n, min_leaf_merge_n)
    if min_samples_leaf is None:
        min_samples_leaf = merge_n
    n_estimators = max(4, int(n_estimators) - int(n_estimators) % 4)
    cf = CausalForestDML(
        n_estimators=n_estimators, max_depth=int(max_depth), min_samples_leaf=int(min_samples_leaf),
        random_state=int(seed), discrete_treatment=True, honest=True, n_jobs=1,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cf.fit(Y, W, X=X)
    forest, trees = cf.model_final_.estimators_[0], cf.model_final_.estimators_[0].estimators_
    B, leaves = len(trees), np.asarray(forest.apply(X), int)
    tau = np.column_stack([np.asarray(t.predict(X, check_input=False)).ravel() for t in trees])
    tau_bar = tau.mean(1)
    V_naive = ((tau - tau_bar[:, None]) ** 2).sum(1) / (B * (B - 1))
    B_leaf = np.zeros(n)
    for b, tree in enumerate(trees):
        tr, lid = tree.tree_, leaves[:, b]
        B_leaf += np.maximum(tr.impurity[lid], 0.0) / np.maximum(tr.n_node_samples[lid], 1.0)
    B_leaf *= (B - 1) / (B ** 2)
    V_debiased = np.maximum(V_naive - B_leaf, 0.0)
    u_trees = (phi[:, None] - tau) ** 2
    u_for_var, merged = u_trees.copy(), leaves.copy()
    if leaf_merge:
        for b, tree in enumerate(trees):
            u_for_var[:, b], merged[:, b] = merge_leaf_u(tree, u_trees[:, b], leaves[:, b], min_n=merge_n)
    r = phi - tau_bar
    u_hat = r ** 2
    u_leaf = merge_leaf_u(trees[0], u_hat, leaves[:, 0], min_n=merge_n)[0] if leaf_merge else u_hat
    Var_u_tree = u_for_var.var(1, ddof=1)
    return {
        "estimator": "causal_forest_leaf_debiased", "n_estimators": B,
        "min_samples_leaf": int(min_samples_leaf), "min_leaf_merge_n": int(merge_n),
        "leaf_assignment": leaves, "leaf_assignment_merged": merged,
        "tau_bar": tau_bar, "u_hat": u_hat, "u_leaf_merged": u_leaf,
        "V_naive": V_naive, "B_leaf": B_leaf, "V_debiased": V_debiased,
        "Var_u_tree": Var_u_tree, "Var_u_delta": 4.0 * V_debiased * (r ** 2),
        "u_draws": u_for_var.T.copy(), "model": cf, "forest": forest,
        "mean_V_naive": float(V_naive.mean()), "mean_V_debiased": float(V_debiased.mean()),
        "mean_Var_u_tree": float(Var_u_tree.mean()), "mean_u": float(u_hat.mean()),
        "bias_ratio": float(B_leaf.mean() / max(V_naive.mean(), 1e-12)),
    }


def bcf_var_u(X, Y, W, phi, *, n_burn=40, n_sim=60, seed=SEED, pihat=None):
    """rpy2 bcf::bcf → tau (n_sim, n) → Var_s (φ - τ_s)²."""
    import rpy2.robjects as ro
    from rpy2.rinterface_lib import callbacks as rpy_cb
    from rpy2.robjects.packages import importr

    X, Y, W, phi = np.asarray(X, float), np.asarray(Y, float).ravel(), np.asarray(W, float).ravel(), np.asarray(phi, float).ravel()
    n, p = X.shape
    if pihat is None:
        pihat = np.clip(LogisticRegression(max_iter=200, random_state=int(seed)).fit(X, W.astype(int)).predict_proba(X)[:, 1], 0.05, 0.95)
    user = str(Path.home() / "R" / "library")
    if Path(user).is_dir():
        ro.r(".libPaths")(ro.StrVector([user] + [str(x) for x in ro.r(".libPaths()")]))
    prev, rpy_cb.consolewrite_print = rpy_cb.consolewrite_print, (lambda *_: None)
    try:
        fit = importr("bcf").bcf(
            y=ro.FloatVector(Y.tolist()), z=ro.FloatVector(W.tolist()),
            x_control=ro.r.matrix(ro.FloatVector(X.ravel().tolist()), nrow=n, ncol=p, byrow=True),
            pihat=ro.FloatVector(np.asarray(pihat, float).ravel().tolist()),
            nburn=int(n_burn), nsim=int(n_sim), n_chains=1, n_threads=1,
            include_pi="control", verbose=False, no_output=True, random_seed=int(seed),
        )
        tau = np.asarray(fit.rx2("tau"), float)
    finally:
        rpy_cb.consolewrite_print = prev
    if tau.shape[1] != n:
        tau = tau.T
    u, tb = (phi - tau) ** 2, tau.mean(0)
    V, uh = u.var(0, ddof=1), (phi - tb) ** 2
    return {"estimator": "bcf_posterior", "backend": "rpy2_bcf", "n_draws": int(tau.shape[0]),
            "tau_bar": tb, "u_hat": uh, "Var_u_posterior": V, "u_draws": u,
            "mean_Var_u": float(V.mean()), "mean_u": float(uh.mean())}


def bootstrap_var_u(X, Y, W, *, B=40, n_estimators=40, seed=SEED):
    """permucate φ → boot RF(φ~X) → Var_b (φ-τ̂_b)²."""
    X, Y, W = np.asarray(X, float), np.asarray(Y, float).ravel(), np.asarray(W, int).ravel()
    n = len(Y)
    phi = permucate_unit_po(X, Y, W, seed=seed, n_estimators=n_estimators)["phi"]
    rng, u = np.random.default_rng(int(seed)), np.empty((B, n))
    rf = lambda b: RandomForestRegressor(n_estimators=max(20, n_estimators // 2), max_depth=4, random_state=seed + b, n_jobs=1)
    for b in range(B):
        idx = rng.integers(0, n, n)
        u[b] = (phi - rf(b).fit(X[idx], phi[idx]).predict(X)) ** 2
    V = u.var(0, ddof=1)
    return {"estimator": "bootstrap", "B": B, "u_hat": u.mean(0), "Var_u_boot": V, "u_draws": u,
            "mean_Var_u": float(V.mean()), "mean_u": float(u.mean())}
