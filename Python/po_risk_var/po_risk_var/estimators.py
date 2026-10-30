"""Var(u) estimators — CF leaf-debiased (default) / BCF / bootstrap.

Public API
----------
  estimate_var_u(X, Y, W, method="cf")  → dict with V_debiased / Var_u / …
  monitor_batch(ref, live, …)           → Mag×Var + Var(u) + TriageDecision
"""
from __future__ import annotations

from typing import Dict, Optional, Sequence, Tuple

import numpy as np

from po_risk_var._path import ensure_cfperm_src
from po_risk_var.product import TriageDecision, estimator_priority, triage

ensure_cfperm_src()

from as_dta_permucate_magvar import (  # noqa: E402
    mag_var_from_unit,
    merge_J,
    permucate_unit_po,
)
from as_dta_po_risk_var import (  # noqa: E402
    bcf_var_u,
    bootstrap_var_u,
    causal_forest_var_u,
    leaf_wise_merge_u,
)


def estimate_var_u(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    *,
    method: str = "cf",
    phi: Optional[np.ndarray] = None,
    n_estimators: int = 60,
    seed: int = 2026,
    leaf_merge: bool = True,
    boot_B: int = 24,
    bcf_sim: int = 40,
) -> dict:
    """Estimate model-uncertainty Var(u) on stacked REF/LIVE.

    method:
      "cf"   — Causal Forest leaf-size-debiased + leaf merge  (product default)
      "bcf"  — BCF / bagged posterior pathway
      "boot" — bootstrap fallback
    """
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    W = np.asarray(W, dtype=int).ravel()
    if phi is None:
        fit = permucate_unit_po(X, Y, W, seed=seed, n_estimators=n_estimators)
        phi = fit["phi"]
        unit = fit["unit_po"]
    else:
        phi = np.asarray(phi, dtype=float).ravel()
        fit = None
        unit = None

    method = str(method).lower()
    if method in ("cf", "causal_forest", "default"):
        out = causal_forest_var_u(
            X,
            Y,
            W,
            phi,
            n_estimators=n_estimators,
            seed=seed,
            leaf_merge=leaf_merge,
        )
        out["method"] = "cf"
        out["Var_u"] = out["Var_u_tree"]
        out["V_for_gate"] = out["V_debiased"]
    elif method in ("bcf", "bayes"):
        out = bcf_var_u(X, Y, W, phi, n_sim=bcf_sim, seed=seed)
        out["method"] = "bcf"
        out["Var_u"] = out["Var_u_posterior"]
        out["V_for_gate"] = out["Var_u_posterior"]
        out["V_debiased"] = out["Var_u_posterior"]
        out["V_naive"] = out["Var_u_posterior"]
    elif method in ("boot", "bootstrap"):
        out = bootstrap_var_u(X, Y, W, B=boot_B, n_estimators=n_estimators, seed=seed)
        out["method"] = "boot"
        out["Var_u"] = out["Var_u_boot"]
        out["V_for_gate"] = out["Var_u_boot"]
        out["V_debiased"] = out["Var_u_boot"]
        out["V_naive"] = out["Var_u_boot"]
    else:
        raise ValueError("method must be cf|bcf|boot, got %r" % method)

    out["phi"] = phi
    out["unit_po"] = unit
    out["priority"] = estimator_priority()
    return out


def monitor_batch(
    X_ref: np.ndarray,
    Y_ref: np.ndarray,
    X_live: np.ndarray,
    Y_live: np.ndarray,
    *,
    method: str = "cf",
    y_flip: float = 0.0,
    n_estimators: int = 60,
    K: int = 28,
    seed: int = 2026,
    quick: bool = False,
) -> Dict:
    """Full product step: Mag×Var + Var(u) + triage decision.

    Returns dict with keys: magvar, var_u, decision, live_mask stats.
    """
    if quick:
        n_estimators, K = 36, 16
    X = np.vstack([np.asarray(X_ref, float), np.asarray(X_live, float)])
    Y = np.concatenate([np.asarray(Y_ref, float).ravel(), np.asarray(Y_live, float).ravel()])
    n0, n1 = len(Y_ref), len(Y_live)
    W = np.array([0] * n0 + [1] * n1, dtype=int)

    fit = permucate_unit_po(X, Y, W, seed=seed, n_estimators=n_estimators)
    phi, unit = fit["phi"], fit["unit_po"]
    mv_ref = mag_var_from_unit(unit[W == 0], K=K, seed=seed)
    mv_live = mag_var_from_unit(unit[W == 1], K=K, seed=seed + 3)
    d_mu = float(mv_live["mu"] - mv_ref["mu"])
    d_s = float(mv_live["s"] - mv_ref["s"])
    J = float(merge_J(d_mu, d_s, mv_ref["s"]))

    var = estimate_var_u(
        X,
        Y,
        W,
        method=method,
        phi=phi,
        n_estimators=n_estimators,
        seed=seed,
    )
    live = W == 1
    V_deb = float(np.mean(np.asarray(var["V_debiased"])[live]))
    V_nai = float(np.mean(np.asarray(var["V_naive"])[live]))
    V_ref = float(np.mean(np.asarray(var["V_debiased"])[~live]))

    decision = triage(
        mu=d_mu,
        s=d_s,
        V_debiased=V_deb,
        V_naive=V_nai,
        J=J,
        y_flip=float(y_flip),
        V_ref=V_ref,
    )
    return {
        "magvar": {
            "mu_ref": mv_ref["mu"],
            "mu_live": mv_live["mu"],
            "s_ref": mv_ref["s"],
            "s_live": mv_live["s"],
            "d_mu": d_mu,
            "d_s": d_s,
            "J": J,
        },
        "var_u": {
            "method": var["method"],
            "mean_V_debiased_live": V_deb,
            "mean_V_naive_live": V_nai,
            "mean_V_debiased_ref": V_ref,
            "mean_Var_u_live": float(np.mean(np.asarray(var["Var_u"])[live])),
            "bias_ratio_live": decision.bias_ratio,
            "backend": var.get("backend"),
        },
        "decision": decision.to_dict(),
        "n_ref": int(n0),
        "n_live": int(n1),
        "y_flip": float(y_flip),
    }


__all__ = [
    "estimate_var_u",
    "monitor_batch",
    "causal_forest_var_u",
    "bcf_var_u",
    "bootstrap_var_u",
    "leaf_wise_merge_u",
    "triage",
    "TriageDecision",
    "estimator_priority",
]
