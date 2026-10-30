"""Synthetic DGPs to justify Var(u) product logic.

Five boards, one disease model each:
  quiet            — no shift
  benign_covariate — uniform X drift, Y intact
  malignant_sparse — subset X attack + Y break
  malignant_dense  — large-fraction attack + Y break
  label_delay      — global Y flip, mild X (concept drift)

Success: product triage.call matches oracle family
  quiet|benign_covariate|label_delay → not VAI_freeze for quiet/benign;
  malignant_* → malignant; label_delay → malignant (outcome break).
"""
from __future__ import annotations

from typing import Dict, List

import numpy as np

from po_risk_var.estimators import monitor_batch


def _base(n: int, p: int, seed: int):
    rng = np.random.default_rng(int(seed))
    X = rng.normal(size=(n, p))
    beta = np.zeros(p)
    beta[:4] = rng.uniform(0.6, 1.2, size=4)
    Y = X @ beta + rng.normal(scale=0.5, size=n)
    return X, Y, beta, rng


def make_synthetic_board(
    name: str,
    *,
    n: int = 220,
    p: int = 14,
    seed: int = 2026,
) -> Dict:
    X0, Y0, beta, rng = _base(n, p, seed)
    ref = {"X": X0, "Y": Y0}
    y_flip = 0.0

    if name == "quiet":
        X1, Y1 = X0.copy(), Y0.copy()
        oracle = "quiet"
    elif name == "benign_covariate":
        X1 = X0 + rng.normal(scale=0.35, size=X0.shape)
        Y1 = X1 @ beta + rng.normal(scale=0.5, size=n)
        oracle = "benign"
    elif name == "malignant_sparse":
        X1 = X0.copy()
        S = rng.choice(n, size=max(40, n // 4), replace=False)
        cols = np.arange(0, min(6, p))
        X1[np.ix_(S, cols)] += rng.normal(loc=3.2, scale=0.5, size=(len(S), len(cols)))
        Y1 = X1 @ beta + rng.normal(scale=0.5, size=n)
        Y1[S] = -Y1[S] + rng.normal(scale=1.0, size=len(S))
        y_flip = float(len(S) / n)
        oracle = "malignant"
    elif name == "malignant_dense":
        X1 = X0.copy()
        S = rng.choice(n, size=max(80, n // 2), replace=False)
        cols = np.arange(0, min(8, p))
        X1[np.ix_(S, cols)] += rng.normal(loc=2.8, scale=0.5, size=(len(S), len(cols)))
        Y1 = X1 @ beta + rng.normal(scale=0.5, size=n)
        Y1[S] = -Y1[S] + rng.normal(scale=1.0, size=len(S))
        y_flip = float(len(S) / n)
        oracle = "malignant"
    elif name == "label_delay":
        X1 = X0 + rng.normal(scale=0.1, size=X0.shape)
        Y1 = -Y0 + rng.normal(scale=0.3, size=n)
        y_flip = 1.0
        oracle = "malignant"
    else:
        raise ValueError(name)

    return {
        "name": name,
        "oracle": oracle,
        "ref": ref,
        "live": {"X": X1, "Y": Y1},
        "y_flip": y_flip,
    }


SYNTH_BOARDS = (
    "quiet",
    "benign_covariate",
    "malignant_sparse",
    "malignant_dense",
    "label_delay",
)


def run_synthetic_apps(
    *,
    seed: int = 2026,
    quick: bool = True,
    method: str = "cf",
) -> dict:
    rows: List[dict] = []
    for name in SYNTH_BOARDS:
        board = make_synthetic_board(name, n=160 if quick else 240, seed=seed)
        print("[synth] ", name, flush=True)
        mon = monitor_batch(
            board["ref"]["X"],
            board["ref"]["Y"],
            board["live"]["X"],
            board["live"]["Y"],
            method=method,
            y_flip=board["y_flip"],
            quick=quick,
            seed=seed,
        )
        pred = mon["decision"]["call"]
        # quiet↔benign counted safe for non-malignant oracles
        if board["oracle"] == "malignant":
            ok = pred == "malignant"
        else:
            ok = pred in ("quiet", "benign")
        rows.append(
            {
                "board": name,
                "oracle": board["oracle"],
                "pred": pred,
                "actuator": mon["decision"]["actuator"],
                "correct": bool(ok),
                "y_flip": board["y_flip"],
                "d_mu": mon["magvar"]["d_mu"],
                "d_s": mon["magvar"]["d_s"],
                "J": mon["magvar"]["J"],
                "V_debiased": mon["var_u"]["mean_V_debiased_live"],
                "V_naive": mon["var_u"]["mean_V_naive_live"],
                "bias_ratio": mon["var_u"]["bias_ratio_live"],
                "reason": mon["decision"]["reason"],
            }
        )
    n_ok = sum(1 for r in rows if r["correct"])
    return {
        "scope": "Synthetic DGPs — product triage with CF-debiased Var(u)",
        "method": method,
        "boards": rows,
        "n_boards": len(rows),
        "n_correct": n_ok,
        "hit_rate": float(n_ok) / max(len(rows), 1),
    }
