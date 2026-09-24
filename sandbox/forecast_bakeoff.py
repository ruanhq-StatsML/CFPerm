"""Sandbox: classical next-step forecast bakeoff (NOT AGOD / PO / excess).

Scenario
--------
Time-ordered stream packs (metro, PM2.5, stocks, waymo_proxy).
Predict y[t+1] from lag features of X/y — plain supervised regression.

Methods (deliberately boring)
-----------------------------
- seasonal_naive: y[t+1] ≈ y[t]  (or y[t+1-period] if period set)
- ridge: Ridge on lagged y + last row of X
- hgb: HistGradientBoostingRegressor same features

Protocol
--------
70/30 chronological split. Report MAE / RMSE / R² on the holdout.
No null tests, no IPTW, no claim router — just numbers.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]


@dataclass
class Pack:
    name: str
    X: np.ndarray
    y: np.ndarray


def _load_metro(max_n: int = 12000) -> Pack:
    from agod.stream_packs import load_metro_interstate

    X, y, meta = load_metro_interstate(ROOT, max_n=max_n)
    return Pack(name=str(meta["name"]), X=X, y=y)


def _load_pm25(max_n: int = 12000) -> Pack:
    from agod.stream_packs import load_beijing_pm25

    X, y, meta = load_beijing_pm25(ROOT, max_n=max_n)
    return Pack(name=str(meta["name"]), X=X, y=y)


def _load_stock(ticker: str = "SPY", max_n: int = 8000) -> Optional[Pack]:
    from agod.stream_packs import load_stocks

    try:
        X, y, meta = load_stocks(ROOT, ticker=ticker, max_n=max_n)
    except Exception:
        return None
    return Pack(name=str(meta.get("name", ticker)), X=X, y=y)


def _load_waymo(max_n: int = 8000) -> Optional[Pack]:
    p = ROOT / "data/stream_packs/waymo_proxy/waymo_proxy_xy.npz"
    if not p.is_file():
        return None
    z = np.load(p, allow_pickle=True)
    X = np.asarray(z["X"], dtype=float)
    y = np.asarray(z["y"], dtype=float).ravel()
    n = min(max_n, len(X))
    return Pack(name="waymo_proxy", X=X[:n], y=y[:n])


def load_packs(max_n: int = 10000) -> List[Pack]:
    packs: List[Pack] = []
    for fn in (_load_metro, _load_pm25, _load_waymo):
        try:
            p = fn(max_n=max_n)
            if p is not None:
                packs.append(p)
        except Exception:
            continue
    for t in ("SPY", "AAPL"):
        p = _load_stock(t, max_n=max_n)
        if p is not None:
            packs.append(p)
    return packs


def make_supervised(
    X: np.ndarray,
    y: np.ndarray,
    *,
    n_lags: int = 5,
) -> Tuple[np.ndarray, np.ndarray]:
    """Rows i use y[i-n_lags:i] + X[i-1] → target y[i]."""
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float).ravel()
    n = min(len(X), len(y))
    X, y = X[:n], y[:n]
    rows, targs = [], []
    for i in range(n_lags, n):
        lag = y[i - n_lags : i]
        feat = np.concatenate([lag, X[i - 1]])
        rows.append(feat)
        targs.append(y[i])
    return np.asarray(rows, float), np.asarray(targs, float)


def _metrics(y_true: np.ndarray, y_hat: np.ndarray) -> Dict[str, float]:
    return {
        "mae": float(mean_absolute_error(y_true, y_hat)),
        "rmse": float(np.sqrt(mean_squared_error(y_true, y_hat))),
        "r2": float(r2_score(y_true, y_hat)),
        "n_test": float(len(y_true)),
    }


def run_bakeoff(
    pack: Pack,
    *,
    n_lags: int = 5,
    train_frac: float = 0.7,
    seed: int = 0,
) -> Dict[str, Any]:
    Z, t = make_supervised(pack.X, pack.y, n_lags=n_lags)
    if len(t) < 50:
        return {"dataset": pack.name, "ok": False, "reason": "too_short"}
    cut = int(train_frac * len(t))
    cut = max(20, min(cut, len(t) - 20))
    Z_tr, Z_te = Z[:cut], Z[cut:]
    y_tr, y_te = t[:cut], t[cut:]

    # seasonal naive: last lag (= y[t-1] which is last of lag block)
    naive = Z_te[:, n_lags - 1]
    results = {"naive_last": _metrics(y_te, naive)}

    scaler = StandardScaler()
    Z_trs = scaler.fit_transform(Z_tr)
    Z_tes = scaler.transform(Z_te)
    ridge = Ridge(alpha=1.0, random_state=seed)
    ridge.fit(Z_trs, y_tr)
    results["ridge"] = _metrics(y_te, ridge.predict(Z_tes))

    hgb = HistGradientBoostingRegressor(
        max_depth=4,
        learning_rate=0.08,
        max_iter=120,
        random_state=seed,
    )
    hgb.fit(Z_tr, y_tr)
    results["hgb"] = _metrics(y_te, hgb.predict(Z_te))

    # best by RMSE
    best = min(results, key=lambda k: results[k]["rmse"])
    # relative lift vs naive
    base = results["naive_last"]["rmse"]
    lifts = {
        k: float(1.0 - results[k]["rmse"] / base) if base > 0 else float("nan")
        for k in results
    }
    return {
        "dataset": pack.name,
        "ok": True,
        "n_train": cut,
        "n_test": len(y_te),
        "n_lags": n_lags,
        "d_feat": int(Z.shape[1]),
        "models": results,
        "rmse_lift_vs_naive": lifts,
        "best_rmse": best,
        "note": "chronological 70/30; classical regression only — not PO/excess",
    }


def bakeoff_all(
    *,
    max_n: int = 10000,
    n_lags: int = 5,
    seed: int = 0,
) -> Dict[str, Any]:
    packs = load_packs(max_n=max_n)
    cards = [run_bakeoff(p, n_lags=n_lags, seed=seed) for p in packs]
    ok = [c for c in cards if c.get("ok")]
    wins = {}
    for c in ok:
        wins[c["best_rmse"]] = wins.get(c["best_rmse"], 0) + 1
    # mean HGB lift vs naive
    hgb_lifts = [
        c["rmse_lift_vs_naive"]["hgb"]
        for c in ok
        if np.isfinite(c["rmse_lift_vs_naive"].get("hgb", np.nan))
    ]
    return {
        "n_packs": len(ok),
        "best_counts": wins,
        "mean_hgb_rmse_lift_vs_naive": float(np.mean(hgb_lifts)) if hgb_lifts else float("nan"),
        "cards": cards,
        "scenario": "next-step stream forecast bakeoff",
        "method": "naive_last / ridge / HGB — time split MAE·RMSE·R²",
        "distance_from_agod": (
            "no excess_auc, no IPTW, no soft_burn, no claim_router, no ToT"
        ),
    }
