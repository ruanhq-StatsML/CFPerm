"""PO-risk variance estimators (legacy import path).

Prefer the packaged API:
  pip install -e Python/po_risk_var
  from po_risk_var import monitor_batch, triage, run_justify_suite
  po-risk-var justify --quick
  po-risk-var products

This module keeps the CF / BCF / bootstrap implementations used by the package.

Defaults
--------
  CF leaf merge: min n=30 (n<200) or 50 (n>=200)
  BCF: rpy2 -> R package bcf with pihat (Rscript bridge fallback)
  Bootstrap: intuitive sanity check
"""
from __future__ import annotations

import json
import subprocess
import tempfile
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LogisticRegression

from as_dta_aigc import ROOT, SEED
from as_dta_permucate_magvar import (
    build_shift_pair,
    mag_var_from_unit,
    merge_J,
    permucate_unit_po,
    pred_shift,
)

OUT_DEFAULT = ROOT / "results" / "as_dta_po_risk_var"
CLIP = 1e-2
MIN_LEAF_MERGE_SMALL = 30
MIN_LEAF_MERGE_LARGE = 50
MIN_LEAF_MERGE_N_SWITCH = 200


def resolve_min_leaf_merge_n(n: int, override: Optional[int] = None) -> int:
    """Product leaf-merge floor: 30 if n<200 else 50 (override wins)."""
    if override is not None:
        return int(max(2, override))
    return MIN_LEAF_MERGE_LARGE if int(n) >= MIN_LEAF_MERGE_N_SWITCH else MIN_LEAF_MERGE_SMALL


# ---------------------------------------------------------------------------
# Leaf merge helpers. econml/sklearn Tree only stores children_*; parent_
# is not kept. We build parent once: parent[root] = -1 (null).
# ---------------------------------------------------------------------------


def tree_parents(tree) -> np.ndarray:
    """parent[v] for every node; parent[root] = -1."""
    tr = tree.tree_ if hasattr(tree, "tree_") else tree
    left = np.asarray(tr.children_left)
    right = np.asarray(tr.children_right)
    parent = np.full(len(left), -1, dtype=int)  # root stays -1
    for p in range(len(left)):
        if left[p] >= 0:
            parent[left[p]] = p
        if right[p] >= 0:
            parent[right[p]] = p
    return parent


def merge_leaf_u(tree, u: np.ndarray, leaf_ids: np.ndarray, min_n: int = 30):
    """Climb parent until n_node >= min_n (stop at root), then mean(u) in cluster.

    Returns (u_merged, merge_ids), both shape (n,).
    """
    u = np.asarray(u, dtype=float).ravel()
    leaf_ids = np.asarray(leaf_ids, dtype=int).ravel()
    tr = tree.tree_ if hasattr(tree, "tree_") else tree
    parent = tree_parents(tree)
    n_node = np.asarray(tr.n_node_samples, dtype=float)
    min_n = int(max(1, min_n))

    merge_ids = np.empty(len(leaf_ids), dtype=int)
    for i, lid in enumerate(leaf_ids):
        v = int(lid)
        while parent[v] >= 0 and n_node[v] < min_n:  # root parent==-1 → stop
            v = int(parent[v])
        merge_ids[i] = v

    u_merged = np.empty_like(u)
    for gid in np.unique(merge_ids):
        m = merge_ids == gid
        u_merged[m] = u[m].mean()
    return u_merged, merge_ids


# ---------------------------------------------------------------------------
# Causal Forest: leaf-size-debiased Var(τ̂) → Var(u)
# ---------------------------------------------------------------------------


def _fit_causal_forest(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    *,
    n_estimators: int = 80,
    max_depth: int = 5,
    min_samples_leaf: int = 30,
    seed: int = SEED,
):
    from econml.dml import CausalForestDML
    # n_estimators must be divisible by subforest_size (default 4)
    n_estimators = int(max(4, n_estimators))
    n_estimators -= n_estimators % 4
    if n_estimators < 4:
        n_estimators = 4
    cf = CausalForestDML(
        n_estimators=n_estimators,
        max_depth=int(max_depth),
        min_samples_leaf=int(min_samples_leaf),
        random_state=int(seed),
        discrete_treatment=True,
        honest=True,
        n_jobs=1,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Y, T=treatment, X=features; no extra confounders W
        cf.fit(Y, W, X=X)
    forest = cf.model_final_.estimators_[0]
    return cf, forest


def causal_forest_var_u(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    phi: np.ndarray,
    *,
    n_estimators: int = 80,
    max_depth: int = 5,
    min_samples_leaf: Optional[int] = None,
    seed: int = SEED,
    leaf_merge: bool = True,
    min_leaf_merge_n: Optional[int] = None,
) -> dict:
    """Leaf-size-debiased tree-to-tree Var(tau) + Var(u) via u^{(b)}.

    Naive:
      V_naive = 1/(B(B-1)) Sum_b (tau_b - tau_bar)^2
    Debiased:
      V_debiased = V_naive - B_leaf,
      B_leaf = ((B-1)/B^2) Sum_b sigma^2_{b,i} / n_{b,i}

    Leaf merge floor: min_n=30 (n<200) or 50 (n>=200). Climb ancestors until
    the merge cluster has at least that many training samples.

    Then for each tree b: u_i^{(b)} = (phi_i - tau_b(X_i))^2,
    optionally leaf-merge u^{(b)}, and
      Var(u)_tree = sample var across b of u^{(b)}
      Var(u)_delta ~ 4 V_debiased * (phi - tau_bar)^2
    """
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    W = np.asarray(W, dtype=int).ravel()
    phi = np.asarray(phi, dtype=float).ravel()
    n = len(phi)
    merge_n = resolve_min_leaf_merge_n(n, min_leaf_merge_n)
    if min_samples_leaf is None:
        min_samples_leaf = merge_n

    cf, forest = _fit_causal_forest(
        X,
        Y,
        W,
        n_estimators=n_estimators,
        max_depth=max_depth,
        min_samples_leaf=min_samples_leaf,
        seed=seed,
    )
    B = int(forest.n_estimators)
    trees = forest.estimators_
    leaves = np.asarray(forest.apply(X), dtype=int)  # (n, B)

    tau_trees = np.column_stack(
        [np.asarray(t.predict(X, check_input=False)).ravel() for t in trees]
    )  # (n, B)
    tau_bar = tau_trees.mean(axis=1)
    # V_naive = 1/(B(B-1)) Sum (tau_b - tau_bar)^2  == var(ddof=1)/B
    sq = (tau_trees - tau_bar[:, None]) ** 2
    V_naive = sq.sum(axis=1) / (B * (B - 1))

    # leaf-size bias term
    B_leaf = np.zeros(n, dtype=float)
    sigma_leaf = np.zeros(n, dtype=float)
    n_leaf = np.zeros(n, dtype=float)
    for b, tree in enumerate(trees):
        tr = tree.tree_
        n_node = np.asarray(tr.n_node_samples, dtype=float)
        # honest estimation-set impurity as sigma^2 proxy (bounded)
        impurity = np.asarray(tr.impurity, dtype=float)
        lid = leaves[:, b]
        n_b = np.maximum(n_node[lid], 1.0)
        sig2 = np.maximum(impurity[lid], 0.0)
        B_leaf += sig2 / n_b
        sigma_leaf += sig2
        n_leaf += n_b
    B_leaf *= (B - 1) / (B ** 2)
    sigma_leaf /= B
    n_leaf /= B

    V_debiased = np.maximum(V_naive - B_leaf, 0.0)

    # leaf assignment → merge_leaf_u per tree → per-sample Var across trees
    u_trees = (phi[:, None] - tau_trees) ** 2  # (n, B)
    leaf_assignment = leaves
    leaf_assignment_merged = leaves.copy()
    u_for_var = u_trees.copy()
    if leaf_merge:
        for b, tree in enumerate(trees):
            u_for_var[:, b], leaf_assignment_merged[:, b] = merge_leaf_u(
                tree, u_trees[:, b], leaves[:, b], min_n=merge_n
            )

    Var_u_tree = u_for_var.var(axis=1, ddof=1)
    u_draws = u_for_var.T.copy()
    r_hat = phi - tau_bar
    Var_u_delta = 4.0 * V_debiased * (r_hat ** 2)

    u_hat = r_hat ** 2
    if leaf_merge:
        u_leaf, _ = merge_leaf_u(trees[0], u_hat, leaves[:, 0], min_n=merge_n)
    else:
        u_leaf = u_hat

    # also expose econml built-in predict_and_var for comparison
    try:
        _tau_b, V_builtin = forest.predict_and_var(X)
        V_builtin = np.asarray(V_builtin).ravel()
    except Exception:  # noqa: BLE001
        V_builtin = np.full(n, np.nan)

    return {
        "estimator": "causal_forest_leaf_debiased",
        "n_estimators": B,
        "min_samples_leaf": int(min_samples_leaf),
        "min_leaf_merge_n": int(merge_n),
        "leaf_assignment": leaf_assignment,  # (n, B) raw leaf id per sample
        "leaf_assignment_merged": leaf_assignment_merged,  # (n, B) after min_n climb
        "tau_bar": tau_bar,
        "u_hat": u_hat,
        "u_leaf_merged": u_leaf,
        "V_naive": V_naive,
        "B_leaf": B_leaf,
        "V_debiased": V_debiased,
        "V_builtin_econml": V_builtin,
        "Var_u_tree": Var_u_tree,
        "Var_u_delta": Var_u_delta,
        "u_draws": u_draws,
        "mean_V_naive": float(np.mean(V_naive)),
        "mean_V_debiased": float(np.mean(V_debiased)),
        "mean_Var_u_tree": float(np.mean(Var_u_tree)),
        "mean_Var_u_delta": float(np.mean(Var_u_delta)),
        "mean_u": float(np.mean(u_hat)),
        "mean_u_leaf": float(np.mean(u_leaf)),
        "mean_sigma_leaf": float(np.mean(sigma_leaf)),
        "mean_n_leaf": float(np.mean(n_leaf)),
        "bias_ratio": float(np.mean(B_leaf) / max(np.mean(V_naive), 1e-12)),
        "model": cf,
        "forest": forest,
    }


# ---------------------------------------------------------------------------
# BCF posterior Var(u) — rpy2 / Rscript / Python bagged posterior
# ---------------------------------------------------------------------------


def _estimate_pihat(X: np.ndarray, W: np.ndarray, seed: int = SEED) -> np.ndarray:
    """Length-n propensity scores for BCF (clipped away from 0/1)."""
    X = np.asarray(X, dtype=float)
    W = np.asarray(W, dtype=int).ravel()
    try:
        lr = LogisticRegression(max_iter=200, random_state=int(seed))
        lr.fit(X, W)
        pi = lr.predict_proba(X)[:, 1]
    except Exception:  # noqa: BLE001
        pi = np.full(len(W), float(np.mean(W)))
    return np.clip(np.asarray(pi, dtype=float), 0.05, 0.95)


def _ensure_r_user_lib() -> None:
    """Prefer ~/R/library so rpy2 / Rscript see user-installed bcf."""
    import os

    user_lib = str(Path.home() / "R" / "library")
    if Path(user_lib).is_dir():
        cur = os.environ.get("R_LIBS_USER", "")
        if user_lib not in cur.split(os.pathsep):
            os.environ["R_LIBS_USER"] = (
                user_lib if not cur else user_lib + os.pathsep + cur
            )
        # also prepend for this process's R session if already started
        try:
            import rpy2.robjects as ro

            paths = [str(p) for p in ro.r(".libPaths()")]
            if user_lib not in paths:
                ro.r(".libPaths")(ro.StrVector([user_lib] + paths))
        except Exception:  # noqa: BLE001
            pass


def _bcf_via_rpy2(
    X, Y, W, phi, n_burn=50, n_sim=80, seed=SEED, pihat: Optional[np.ndarray] = None
) -> Optional[dict]:
    """BCF via rpy2 + R package bcf (requires pihat)."""
    try:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr

        _ensure_r_user_lib()
        # requireNamespace via ro.r can return None under some rpy2 builds;
        # importr is the reliable probe.
        try:
            bcf = importr("bcf")
        except Exception:
            return None
        X = np.asarray(X, dtype=float)
        Y = np.asarray(Y, dtype=float).ravel()
        W = np.asarray(W, dtype=float).ravel()
        if pihat is None:
            pihat = _estimate_pihat(X, W.astype(int), seed=seed)
        pihat = np.asarray(pihat, dtype=float).ravel()
        n, p = X.shape
        x_r = ro.r.matrix(
            ro.FloatVector(X.ravel().tolist()), nrow=n, ncol=p, byrow=True
        )
        # Mute BCF progress bar (binary bytes break rpy2 utf-8 console callback)
        try:
            from rpy2.rinterface_lib import callbacks as rpy_cb

            _prev_write = rpy_cb.consolewrite_print

            def _quiet_console(x):  # noqa: ANN001
                return None

            rpy_cb.consolewrite_print = _quiet_console
        except Exception:  # noqa: BLE001
            _prev_write = None
        try:
            # Keep MCMC light for product smoke; n_chains=1 avoids thread storms
            fit = bcf.bcf(
                y=ro.FloatVector(Y.tolist()),
                z=ro.FloatVector(W.tolist()),
                x_control=x_r,
                pihat=ro.FloatVector(pihat.tolist()),
                nburn=int(n_burn),
                nsim=int(n_sim),
                n_chains=1,
                n_threads=1,
                ntree_control=min(50, max(20, 2 * p)),
                ntree_moderate=min(20, max(10, p)),
                include_pi="control",
                verbose=False,
                no_output=True,
                random_seed=int(seed),
            )
        finally:
            if _prev_write is not None:
                try:
                    rpy_cb.consolewrite_print = _prev_write
                except Exception:  # noqa: BLE001
                    pass
        tau_draws = None
        names = list(fit.names) if getattr(fit, "names", None) is not None else []
        for key in ("tau", "tau_post", "tauhat"):
            if key in names:
                tau_draws = np.asarray(fit.rx2(key))
                break
        if tau_draws is None:
            # last resort: positional / attribute
            try:
                tau_draws = np.asarray(fit.rx2("tau"))
            except Exception:  # noqa: BLE001
                return {"_error": "rpy2_bcf: no tau draws in fit"}
        if tau_draws.ndim == 1:
            tau_draws = tau_draws.reshape(1, -1)
        if tau_draws.shape[0] == len(Y) and tau_draws.shape[1] != len(Y):
            tau_draws = tau_draws.T
        return _posterior_var_from_tau_draws(phi, tau_draws, backend="rpy2_bcf")
    except Exception as exc:  # noqa: BLE001
        return {"_error": "rpy2_bcf: %s" % exc}


def _bcf_via_rscript(
    X, Y, W, phi, n_burn=40, n_sim=60, seed=SEED, pihat: Optional[np.ndarray] = None
) -> Optional[dict]:
    """Call R bcf via Rscript CSV bridge when rpy2 path fails."""
    import os

    _ensure_r_user_lib()
    env = os.environ.copy()
    user_lib = str(Path.home() / "R" / "library")
    env["R_LIBS_USER"] = user_lib
    try:
        chk = subprocess.run(
            ["Rscript", "-e", 'cat(requireNamespace("bcf", quietly=TRUE))'],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=30,
            env=env,
        )
        if "TRUE" not in (chk.stdout or ""):
            return None
    except Exception:  # noqa: BLE001
        return None

    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    W = np.asarray(W, dtype=int).ravel()
    if pihat is None:
        pihat = _estimate_pihat(X, W, seed=seed)
    pihat = np.asarray(pihat, dtype=float).ravel()

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        data = td / "data.csv"
        out = td / "tau.csv"
        import pandas as pd

        cols = ["Y", "W", "pihat"] + ["x%d" % j for j in range(X.shape[1])]
        pd.DataFrame(
            np.column_stack([Y, W, pihat, X]),
            columns=cols,
        ).to_csv(data, index=False)
        ntree_c = min(50, max(20, 2 * X.shape[1]))
        ntree_m = min(20, max(10, X.shape[1]))
        rcode = f"""
        .libPaths(c("{user_lib}", .libPaths()))
        library(bcf)
        d <- read.csv("{data}")
        y <- d$Y; z <- as.numeric(d$W); pihat <- d$pihat
        x <- as.matrix(d[, -(1:3), drop=FALSE])
        fit <- bcf(y=y, z=z, x_control=x, pihat=pihat,
                   nburn={int(n_burn)}, nsim={int(n_sim)},
                   n_chains=1, n_threads=1,
                   ntree_control={ntree_c}, ntree_moderate={ntree_m},
                   include_pi="control", verbose=FALSE,
                   random_seed={int(seed)})
        tau <- fit$tau
        write.csv(tau, "{out}", row.names=FALSE)
        """
        try:
            proc = subprocess.run(
                ["Rscript", "-e", rcode],
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=300,
                cwd=str(td),
                env=env,
            )
            if proc.returncode != 0 and not out.exists():
                return {
                    "_error": "rscript_bcf: rc=%s %s"
                    % (proc.returncode, (proc.stderr or "")[-400:])
                }
        except Exception as exc:  # noqa: BLE001
            return {"_error": "rscript_bcf: %s" % exc}
        if not out.exists():
            return None
        tau_draws = np.loadtxt(out, delimiter=",", skiprows=1)
        if tau_draws.ndim == 1:
            tau_draws = tau_draws.reshape(1, -1)
        if tau_draws.shape[1] == len(Y):
            pass
        elif tau_draws.shape[0] == len(Y):
            tau_draws = tau_draws.T
        else:
            return {"_error": "rscript_bcf: bad tau shape %s" % (tau_draws.shape,)}
        return _posterior_var_from_tau_draws(phi, tau_draws, backend="rscript_bcf")


def _bcf_python_bagged(
    X: np.ndarray,
    phi: np.ndarray,
    *,
    n_draws: int = 60,
    n_estimators: int = 40,
    seed: int = SEED,
) -> dict:
    """Python BCF-style posterior: bagged honest RF on φ ~ X (regularized τ̂ draws).

    Not full BCF MCMC, but same *pathway*: draws of τ̂ → u^{(s)} → Var(u).
    Priors / honesty via min_samples_leaf + bootstrap subsample of rows.
    """
    X = np.asarray(X, dtype=float)
    phi = np.asarray(phi, dtype=float).ravel()
    n = len(phi)
    rng = np.random.default_rng(int(seed))
    tau_draws = np.zeros((n_draws, n), dtype=float)
    for s in range(int(n_draws)):
        idx = rng.choice(n, size=n, replace=True)
        oob = np.ones(n, dtype=bool)
        oob[idx] = False
        rf = RandomForestRegressor(
            n_estimators=max(20, int(n_estimators) // 2),
            max_depth=4,
            min_samples_leaf=max(5, n // 40),
            random_state=int(seed) + s,
            n_jobs=1,
        )
        rf.fit(X[idx], phi[idx])
        tau_draws[s] = rf.predict(X)
    return _posterior_var_from_tau_draws(phi, tau_draws, backend="python_bagged_bcf")


def _posterior_var_from_tau_draws(phi: np.ndarray, tau_draws: np.ndarray, backend: str) -> dict:
    phi = np.asarray(phi, dtype=float).ravel()
    tau_draws = np.asarray(tau_draws, dtype=float)
    if tau_draws.ndim != 2:
        raise ValueError("tau_draws must be 2d")
    if tau_draws.shape[1] != len(phi):
        if tau_draws.shape[0] == len(phi):
            tau_draws = tau_draws.T
        else:
            raise ValueError("tau_draws shape mismatch")
    u_draws = (phi[None, :] - tau_draws) ** 2
    Var_u = u_draws.var(axis=0, ddof=1)
    tau_bar = tau_draws.mean(axis=0)
    u_hat = (phi - tau_bar) ** 2
    return {
        "estimator": "bcf_posterior",
        "backend": backend,
        "n_draws": int(tau_draws.shape[0]),
        "tau_bar": tau_bar,
        "u_hat": u_hat,
        "Var_u_posterior": Var_u,
        "u_draws": u_draws,
        "mean_Var_u": float(np.mean(Var_u)),
        "mean_u": float(np.mean(u_hat)),
        "u_draws_std": float(np.std(u_draws)),
    }


def bcf_var_u(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    phi: np.ndarray,
    *,
    n_burn: int = 40,
    n_sim: int = 60,
    seed: int = SEED,
    prefer_python: bool = False,
    pihat: Optional[np.ndarray] = None,
) -> dict:
    """BCF posterior Var(u). Tries rpy2 (with pihat) -> Rscript -> Python bagged."""
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    W = np.asarray(W, dtype=int).ravel()
    phi = np.asarray(phi, dtype=float).ravel()
    if pihat is None:
        pihat = _estimate_pihat(X, W, seed=seed)
    errors = []
    if not prefer_python:
        for fn in (_bcf_via_rpy2, _bcf_via_rscript):
            out = fn(
                X, Y, W, phi, n_burn=n_burn, n_sim=n_sim, seed=seed, pihat=pihat
            )
            if out is None:
                continue
            if "_error" in out and "Var_u_posterior" not in out:
                errors.append(out["_error"])
                continue
            if "Var_u_posterior" in out:
                out["fallback_errors"] = errors
                return out
    out = _bcf_python_bagged(X, phi, n_draws=n_sim, seed=seed)
    out["fallback_errors"] = errors
    return out


# ---------------------------------------------------------------------------
# Bootstrap fallback
# ---------------------------------------------------------------------------


def bootstrap_var_u(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    *,
    B: int = 40,
    n_estimators: int = 40,
    seed: int = SEED,
) -> dict:
    """φ once; bootstrap RF(φ~X) → u=(φ-τ̂)² → Var across boots."""
    X = np.asarray(X, float)
    Y = np.asarray(Y, float).ravel()
    W = np.asarray(W, int).ravel()
    n = len(Y)
    phi = permucate_unit_po(X, Y, W, seed=seed, n_estimators=n_estimators)["phi"]
    rng = np.random.default_rng(int(seed))
    kw = dict(
        n_estimators=max(20, n_estimators // 2),
        max_depth=4,
        min_samples_leaf=max(3, n // 50),
        n_jobs=1,
    )
    u_draws = np.empty((int(B), n))
    for b in range(int(B)):
        idx = rng.integers(0, n, n)
        if W[idx].min() == W[idx].max():  # keep both arms
            idx[0], idx[1] = np.flatnonzero(W == 0)[0], np.flatnonzero(W == 1)[0]
        tau_b = RandomForestRegressor(random_state=int(seed) + b, **kw).fit(
            X[idx], phi[idx]
        ).predict(X)
        u_draws[b] = (phi - tau_b) ** 2
    Var_u = u_draws.var(axis=0, ddof=1)
    return {
        "estimator": "bootstrap",
        "B": int(B),
        "u_hat": u_draws.mean(axis=0),
        "Var_u_boot": Var_u,
        "u_draws": u_draws,
        "mean_Var_u": float(Var_u.mean()),
        "mean_u": float(u_draws.mean()),
    }


# ---------------------------------------------------------------------------
# Board eval: benign vs malignant — compare estimators
# ---------------------------------------------------------------------------


def _stack_ref_live(ref: dict, live: dict) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    X = np.vstack([ref["X"], live["X"]])
    Y = np.concatenate([ref["Y"], live["Y"]])
    W = np.array([0] * len(ref["Y"]) + [1] * len(live["Y"]), dtype=int)
    return X, Y, W


def eval_shift_var_estimators(
    kind: str,
    pair: dict,
    *,
    seed: int = SEED,
    n_estimators: int = 60,
    boot_B: int = 24,
    bcf_sim: int = 40,
    quick: bool = False,
) -> dict:
    """Compare CF-debiased / BCF / bootstrap / subenv Mag×Var on one shift."""
    live = pair[kind]
    X, Y, W = _stack_ref_live(pair["ref"], live)
    if quick:
        n_estimators, boot_B, bcf_sim = 40, 12, 24

    fit = permucate_unit_po(X, Y, W, seed=seed, n_estimators=n_estimators)
    phi = fit["phi"]
    unit = fit["unit_po"]
    live_mask = W == 1
    unit_live = unit[live_mask]
    mv = mag_var_from_unit(unit_live, K=24 if quick else 36, seed=seed)

    # empirical fourth-moment Var(u) on LIVE
    u_c = unit_live - unit_live.mean()
    emp_var_u = float(np.mean(u_c ** 2))

    cf = causal_forest_var_u(
        X,
        Y,
        W,
        phi,
        n_estimators=n_estimators,
        seed=seed,
        leaf_merge=True,
    )
    bcf = bcf_var_u(
        X, Y, W, phi, n_sim=bcf_sim, n_burn=max(20, bcf_sim // 2), seed=seed
    )
    boot = bootstrap_var_u(
        X, Y, W, B=boot_B, n_estimators=max(20, n_estimators // 2), seed=seed
    )

    # LIVE-only means of model-uncertainty Var(u)
    def _live_mean(arr):
        return float(np.mean(np.asarray(arr)[live_mask]))

    row = {
        "kind": kind,
        "n_live": int(live_mask.sum()),
        "subenv_mu": mv["mu"],
        "subenv_s": mv["s"],
        "emp_var_u": emp_var_u,
        "cf_mean_V_naive": float(np.mean(cf["V_naive"][live_mask])),
        "cf_mean_V_debiased": float(np.mean(cf["V_debiased"][live_mask])),
        "cf_bias_ratio": float(
            np.mean(cf["B_leaf"][live_mask])
            / max(np.mean(cf["V_naive"][live_mask]), 1e-12)
        ),
        "cf_mean_Var_u_tree": _live_mean(cf["Var_u_tree"]),
        "cf_mean_Var_u_delta": _live_mean(cf["Var_u_delta"]),
        "cf_mean_u_leaf": _live_mean(cf["u_leaf_merged"]),
        "bcf_backend": bcf.get("backend"),
        "bcf_mean_Var_u": _live_mean(bcf["Var_u_posterior"]),
        "bcf_mean_u": _live_mean(bcf["u_hat"]),
        "boot_mean_Var_u": _live_mean(boot["Var_u_boot"]),
        "boot_mean_u": _live_mean(boot["u_hat"]),
        "J_subenv": None,  # filled by caller vs REF
    }
    return row


def run_po_risk_var_suite(
    *,
    seed: int = SEED,
    out_dir: Optional[Path] = None,
    quick: bool = False,
    n: int = 280,
) -> dict:
    """Benign vs malignant tabular board — effect of variance estimators."""
    out_dir = Path(out_dir or OUT_DEFAULT)
    out_dir.mkdir(parents=True, exist_ok=True)
    if quick:
        n = 180

    pair = build_shift_pair(n=n, p=16, seed=seed)
    # REF baseline for Δ
    Xr, Yr, Wr = _stack_ref_live(pair["ref"], pair["ref"])
    # use half ref as pseudo-live for REF magvar baseline
    fit_ref = permucate_unit_po(
        np.vstack([pair["ref"]["X"], pair["ref"]["X"]]),
        np.concatenate([pair["ref"]["Y"], pair["ref"]["Y"]]),
        np.array([0] * n + [1] * n),
        seed=seed,
        n_estimators=40 if quick else 60,
    )
    mv_ref = mag_var_from_unit(
        fit_ref["unit_po"][n:], K=20 if quick else 32, seed=seed
    )

    rows = {}
    for kind in ("benign", "malignant"):
        print("[po-var] eval", kind, flush=True)
        rows[kind] = eval_shift_var_estimators(
            kind, pair, seed=seed, quick=quick
        )
        rows[kind]["d_mu"] = rows[kind]["subenv_mu"] - mv_ref["mu"]
        rows[kind]["d_s"] = rows[kind]["subenv_s"] - mv_ref["s"]
        rows[kind]["J_subenv"] = merge_J(
            rows[kind]["d_mu"], rows[kind]["d_s"], mv_ref["s"]
        )
        rows[kind]["pred_subenv"] = pred_shift(
            rows[kind]["d_mu"], rows[kind]["d_s"], rows[kind]["J_subenv"]
        )

    # Separation scores: malignant / benign ratio on each Var(u) estimator
    def _sep(key):
        b, m = rows["benign"][key], rows["malignant"][key]
        return {
            "benign": b,
            "malignant": m,
            "gap": m - b,
            "ratio": (m / b) if abs(b) > 1e-12 else float("inf"),
            "separates": bool(m > b * 1.05 + 1e-6),
        }

    separation = {
        "subenv_s": _sep("subenv_s"),
        "emp_var_u": _sep("emp_var_u"),
        "cf_V_naive": _sep("cf_mean_V_naive"),
        "cf_V_debiased": _sep("cf_mean_V_debiased"),
        "cf_Var_u_tree": _sep("cf_mean_Var_u_tree"),
        "cf_Var_u_delta": _sep("cf_mean_Var_u_delta"),
        "bcf_Var_u": _sep("bcf_mean_Var_u"),
        "boot_Var_u": _sep("boot_mean_Var_u"),
    }
    n_ok = sum(1 for v in separation.values() if v["separates"])

    payload = {
        "scope": (
            "PO-risk Var(u) estimators on PermuCATE φ: Causal Forest "
            "leaf-size-debiased (+ leaf merge), BCF posterior, bootstrap fallback."
        ),
        "unit_po": "(φ - τ̂)^2",
        "ref_magvar": mv_ref,
        "shifts": rows,
        "separation": separation,
        "n_estimators_separating": n_ok,
        "n_estimators_total": len(separation),
        "hit_rate": float(n_ok) / max(len(separation), 1),
        "recommend": (
            "Default: Causal Forest leaf-size-debiased + leaf-wise merge "
            "(min n=30 if n<200 else 50). BCF via rpy2+pihat when Bayesian "
            "path is required. Bootstrap only as fallback / sanity check."
        ),
        "bcf_backend_used": rows["benign"].get("bcf_backend"),
        "min_leaf_merge_n": resolve_min_leaf_merge_n(n),
    }
    (out_dir / "po_risk_var.json").write_text(json.dumps(payload, indent=2, default=str))
    _plot(payload, out_dir / "po_risk_var.png")
    _readme(out_dir / "README.md", payload)
    return payload


def _plot(payload: dict, path: Path):
    import matplotlib.pyplot as plt

    sep = payload["separation"]
    names = list(sep.keys())
    gaps = [sep[k]["gap"] for k in names]
    ratios = [min(sep[k]["ratio"], 20.0) for k in names]
    colors = ["#2F6B4F" if sep[k]["separates"] else "#A33B24" for k in names]

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.4))
    x = np.arange(len(names))

    ax = axes[0]
    ax.bar(x, gaps, color=colors, width=0.7)
    ax.axhline(0, color="#666", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("gap = malignant − benign")
    ax.set_title("Var(u) separation gap", loc="left", fontweight="bold")
    ax.grid(True, axis="y", alpha=0.3)

    ax = axes[1]
    ax.bar(x, ratios, color=colors, width=0.7)
    ax.axhline(1.0, ls="--", color="#999", lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("ratio malignant / benign (capped 20)")
    ax.set_title(
        "Hit-rate %.0f%% (%d/%d) — CF debiased is default"
        % (100 * payload["hit_rate"], payload["n_estimators_separating"], payload["n_estimators_total"]),
        loc="left",
        fontweight="bold",
    )
    ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle(
        "PO-risk Var(u) estimators — benign vs malignant",
        fontsize=12,
        fontweight="bold",
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def _readme(path: Path, payload: dict):
    lines = [
        "# PO-risk Var(u) estimators",
        "",
        payload["scope"],
        "",
        "**Recommend:** " + payload["recommend"],
        "",
        "BCF backend: `%s`" % payload.get("bcf_backend_used"),
        "",
        "## Separation (malignant vs benign)",
        "",
        "| estimator | benign | malignant | gap | ratio | OK? |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for k, v in payload["separation"].items():
        lines.append(
            "| %s | %.4g | %.4g | %.4g | %.3g | %s |"
            % (
                k,
                v["benign"],
                v["malignant"],
                v["gap"],
                v["ratio"] if np.isfinite(v["ratio"]) else float("inf"),
                "yes" if v["separates"] else "no",
            )
        )
    lines.extend(
        [
            "",
            "hit-rate: **%.2f**" % payload["hit_rate"],
            "",
            "```bash",
            "python3 scripts/run_as_dta_po_risk_var.py --quick",
            "```",
            "",
        ]
    )
    path.write_text("\n".join(lines))
