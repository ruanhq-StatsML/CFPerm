#!/usr/bin/env python3
"""PO-risk Var(u) board: CF / BCF / bootstrap on benign vs malignant.

  python3 scripts/run_as_dta_po_risk_var.py --quick
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from as_dta_aigc import SEED  # noqa: E402
from as_dta_permucate_magvar import (  # noqa: E402
    build_shift_pair,
    mag_var_from_unit,
    merge_J,
    permucate_unit_po,
    pred_shift,
)
from as_dta_po_risk_var import (  # noqa: E402
    bootstrap_var_u,
    bcf_var_u,
    causal_forest_var_u,
    resolve_min_leaf_merge_n,
)

OUT = ROOT / "results" / "as_dta_po_risk_var"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")


def _stack(ref, live):
    X = np.vstack([ref["X"], live["X"]])
    Y = np.concatenate([ref["Y"], live["Y"]])
    W = np.array([0] * len(ref["Y"]) + [1] * len(live["Y"]), dtype=int)
    return X, Y, W


def _eval(kind, pair, *, seed, quick):
    n_est, boot_B, bcf_sim = (40, 12, 24) if quick else (60, 24, 40)
    X, Y, W = _stack(pair["ref"], pair[kind])
    fit = permucate_unit_po(X, Y, W, seed=seed, n_estimators=n_est)
    phi, live = fit["phi"], W == 1
    unit_live = fit["unit_po"][live]
    mv = mag_var_from_unit(unit_live, K=24 if quick else 36, seed=seed)
    u_c = unit_live - unit_live.mean()
    cf = causal_forest_var_u(X, Y, W, phi, n_estimators=n_est, seed=seed, leaf_merge=True)
    bcf = bcf_var_u(X, Y, W, phi, n_sim=bcf_sim, n_burn=max(20, bcf_sim // 2), seed=seed)
    boot = bootstrap_var_u(X, Y, W, B=boot_B, n_estimators=max(20, n_est // 2), seed=seed)
    lm = lambda a: float(np.mean(np.asarray(a)[live]))
    return {
        "kind": kind, "n_live": int(live.sum()),
        "subenv_mu": mv["mu"], "subenv_s": mv["s"], "emp_var_u": float(np.mean(u_c ** 2)),
        "cf_mean_V_naive": float(np.mean(cf["V_naive"][live])),
        "cf_mean_V_debiased": float(np.mean(cf["V_debiased"][live])),
        "cf_bias_ratio": float(np.mean(cf["B_leaf"][live]) / max(np.mean(cf["V_naive"][live]), 1e-12)),
        "cf_mean_Var_u_tree": lm(cf["Var_u_tree"]), "cf_mean_Var_u_delta": lm(cf["Var_u_delta"]),
        "cf_mean_u_leaf": lm(cf["u_leaf_merged"]),
        "bcf_backend": bcf.get("backend"), "bcf_mean_Var_u": lm(bcf["Var_u_posterior"]), "bcf_mean_u": lm(bcf["u_hat"]),
        "boot_mean_Var_u": lm(boot["Var_u_boot"]), "boot_mean_u": lm(boot["u_hat"]),
    }


def run_po_risk_var_suite(*, seed=SEED, out_dir=None, quick=False, n=280):
    out_dir = Path(out_dir or OUT)
    out_dir.mkdir(parents=True, exist_ok=True)
    if quick:
        n = 180
    pair = build_shift_pair(n=n, p=16, seed=seed)
    fit_ref = permucate_unit_po(
        np.vstack([pair["ref"]["X"], pair["ref"]["X"]]),
        np.concatenate([pair["ref"]["Y"], pair["ref"]["Y"]]),
        np.array([0] * n + [1] * n), seed=seed, n_estimators=40 if quick else 60,
    )
    mv_ref = mag_var_from_unit(fit_ref["unit_po"][n:], K=20 if quick else 32, seed=seed)
    rows = {}
    for kind in ("benign", "malignant"):
        print("[po-var] eval", kind, flush=True)
        rows[kind] = _eval(kind, pair, seed=seed, quick=quick)
        rows[kind]["d_mu"] = rows[kind]["subenv_mu"] - mv_ref["mu"]
        rows[kind]["d_s"] = rows[kind]["subenv_s"] - mv_ref["s"]
        rows[kind]["J_subenv"] = merge_J(rows[kind]["d_mu"], rows[kind]["d_s"], mv_ref["s"])
        rows[kind]["pred_subenv"] = pred_shift(rows[kind]["d_mu"], rows[kind]["d_s"], rows[kind]["J_subenv"])

    def _sep(key):
        b, m = rows["benign"][key], rows["malignant"][key]
        return {"benign": b, "malignant": m, "gap": m - b,
                "ratio": (m / b) if abs(b) > 1e-12 else float("inf"),
                "separates": bool(m > b * 1.05 + 1e-6)}

    separation = {k: _sep(k) for k in (
        "subenv_s", "emp_var_u", "cf_mean_V_naive", "cf_mean_V_debiased",
        "cf_mean_Var_u_tree", "cf_mean_Var_u_delta", "bcf_mean_Var_u", "boot_mean_Var_u",
    )}
    # alias keys expected by docs/print
    separation = {
        "subenv_s": separation["subenv_s"], "emp_var_u": separation["emp_var_u"],
        "cf_V_naive": separation["cf_mean_V_naive"], "cf_V_debiased": separation["cf_mean_V_debiased"],
        "cf_Var_u_tree": separation["cf_mean_Var_u_tree"], "cf_Var_u_delta": separation["cf_mean_Var_u_delta"],
        "bcf_Var_u": separation["bcf_mean_Var_u"], "boot_Var_u": separation["boot_mean_Var_u"],
    }
    n_ok = sum(1 for v in separation.values() if v["separates"])
    payload = {
        "scope": "PO-risk Var(u): CF leaf-debiased, BCF, bootstrap.",
        "unit_po": "(φ - τ̂)^2", "ref_magvar": mv_ref, "shifts": rows, "separation": separation,
        "n_estimators_separating": n_ok, "n_estimators_total": len(separation),
        "hit_rate": float(n_ok) / max(len(separation), 1),
        "recommend": "Default CF leaf-debiased + merge (min n=30/50). BCF via rpy2. Bootstrap fallback.",
        "bcf_backend_used": rows["benign"].get("bcf_backend"),
        "min_leaf_merge_n": resolve_min_leaf_merge_n(n),
    }
    (out_dir / "po_risk_var.json").write_text(json.dumps(payload, indent=2, default=str))
    # tiny bar plot
    import matplotlib.pyplot as plt
    names = list(separation)
    gaps = [separation[k]["gap"] for k in names]
    colors = ["#2F6B4F" if separation[k]["separates"] else "#A33B24" for k in names]
    fig, ax = plt.subplots(figsize=(10, 3.6))
    ax.bar(np.arange(len(names)), gaps, color=colors)
    ax.axhline(0, color="#666", lw=0.8)
    ax.set_xticks(np.arange(len(names)))
    ax.set_xticklabels(names, rotation=30, ha="right", fontsize=8)
    ax.set_title("Var(u) gap malignant−benign  hit=%.0f%%" % (100 * payload["hit_rate"]), loc="left")
    fig.tight_layout()
    fig.savefig(out_dir / "po_risk_var.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    lines = ["# PO-risk Var(u) estimators", "", payload["recommend"], "",
             "| est | ben | mal | gap | OK |", "| --- | --- | --- | --- | --- |"]
    for k, v in separation.items():
        lines.append("| %s | %.4g | %.4g | %.4g | %s |" % (k, v["benign"], v["malignant"], v["gap"], "yes" if v["separates"] else "no"))
    (out_dir / "README.md").write_text("\n".join(lines) + "\n")
    return payload


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()
    payload = run_po_risk_var_suite(seed=int(args.seed), out_dir=OUT, quick=bool(args.quick))
    print("[po-var] hit_rate=%.2f  bcf=%s" % (payload["hit_rate"], payload.get("bcf_backend_used")), flush=True)
    for k, v in payload["separation"].items():
        print("[po-var] %-16s  gap=%.4g  %s" % (k, v["gap"], "OK" if v["separates"] else "MISS"), flush=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "README.md", DOCS / "AS_DTA_po_risk_var_estimators.md")
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "po_risk_var.png", ARTIFACTS / "as_dta_po_risk_var.png")
    print("[po-var] wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
