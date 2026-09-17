"""OnlineRFPerm: shift-onset on a frozen RandomForestRegressor.

Fit once on D_ref. Each incoming batch is:

    pred = rf.predict(np.asarray(X_new))
    T_t  = MSE(Y_new, pred) − E_ref

Last-two hop_fires is WHEN the stream moved. Rank-p is the backup mark.
PO × MSE × MMD²(X_new, X_ref) says WHAT it was.
No online-bootstrap.
"""
from __future__ import annotations

import numpy as np
from sklearn.ensemble import RandomForestRegressor

GATE = 1.5
P_ALPHA = 0.05
NULL_B = 200


def fit_frozen_rf(X, y, *, seed: int = 0) -> RandomForestRegressor:
    """Shallow RF. The OOD probe is rf.predict(X_new)."""
    rf = RandomForestRegressor(
        n_estimators=20,
        max_depth=4,
        min_samples_leaf=3,
        random_state=int(seed),
        n_jobs=1,
    )
    rf.fit(np.asarray(X, dtype=float), np.asarray(y, dtype=float).ravel())
    return rf


def probe_mse(model, X, y) -> float:
    """mean((Y − rf.predict(X))²). Same one-liner as the board probe."""
    pred = np.asarray(model.predict(np.asarray(X, dtype=float)), dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    return float(np.mean((y - pred) ** 2))


def error_floor(n) -> float:
    return max(1.0 / max(int(n), 1), 1e-8)


def shift_ratio(e_now, e_prev, e_floor=0.0) -> float:
    denom = max(float(e_prev), float(e_floor or 0.0), 1e-8)
    return float(e_now) / denom


def hop_fires(e_now, e_prev, gate: float = GATE, e_floor: float = 0.0) -> bool:
    """Consecutive OOS gate. Vacuous e_prev → quiet. Floor only pads the denom."""
    if e_prev is None:
        return False
    ratio = shift_ratio(e_now, e_prev, e_floor=e_floor)
    return bool(np.isfinite(ratio) and ratio >= float(gate))


class FrozenRFPerm:
    """Fit RandomForestRegressor once on D_ref. Each batch is one T_t / p_t / hop.

    Rank-p for FDR is against a fixed ref-null pool (subsets of D_ref), not the
    short stream of previous T's — those ranks are too coarse for ADDIS/SAFFRON.
    """

    def __init__(self, X_ref, Y_ref, *, seed: int = 2026, gate: float = GATE, null_b: int = NULL_B):
        X_ref = np.asarray(X_ref, dtype=float)
        Y_ref = np.asarray(Y_ref, dtype=float).ravel()
        self.gate = float(gate)
        self.probe = fit_frozen_rf(X_ref, Y_ref, seed=seed)
        self.e_ref = probe_mse(self.probe, X_ref, Y_ref)
        self.null_pool = self._ref_null_pool(X_ref, Y_ref, seed=seed, b=int(null_b))
        self.pool: list[float] = []
        self.e_prev = None

    def _ref_null_pool(self, X_ref, Y_ref, *, seed: int, b: int) -> np.ndarray:
        rng = np.random.default_rng(int(seed) + 17)
        n = len(Y_ref)
        m = min(n, max(64, n // 5))
        out = np.empty(max(int(b), 1), dtype=float)
        for i in range(len(out)):
            sl = rng.choice(n, size=m, replace=False)
            out[i] = probe_mse(self.probe, X_ref[sl], Y_ref[sl]) - self.e_ref
        return out

    def _pval(self, T: float) -> float:
        pool = np.asarray(self.null_pool, dtype=float)
        if pool.size == 0:
            return 1.0
        return float(np.sum(float(T) <= pool) + 1.0) / float(pool.size + 1.0)

    def step(self, X_new, Y_new) -> dict:
        mse = probe_mse(self.probe, X_new, Y_new)
        T = float(mse - self.e_ref)
        p = self._pval(T)
        pool_arr = np.asarray(self.pool, dtype=float)
        p_seq = 1.0 if pool_arr.size == 0 else float(np.sum(T <= pool_arr) + 1) / float(len(self.pool) + 1)
        self.pool.append(T)
        n = len(np.asarray(Y_new).ravel())
        e_fl = error_floor(n)
        hop = hop_fires(mse, self.e_prev, gate=self.gate, e_floor=e_fl)
        ratio = 1.0 if self.e_prev is None else shift_ratio(mse, self.e_prev, e_floor=e_fl)
        rec = {
            "rfperm_mse": float(mse),
            "rfperm_T": float(T),
            "rfperm_p": float(p),
            "rfperm_p_seq": float(p_seq),
            "rfperm_hop": bool(hop),
            "rfperm_ratio": float(ratio),
            "rfperm_e_ref": float(self.e_ref),
        }
        self.e_prev = float(mse)
        return rec


def onset_from_rows(rows, *, p_alpha: float = P_ALPHA) -> dict:
    """Shift-onset point: first last-two hop. Rank-p is the backup mark."""
    hop_t = next((int(r["t"]) for r in rows if r.get("rfperm_hop")), None)
    rank_t = next(
        (
            int(r["t"])
            for r in rows
            if r.get("rfperm_p") is not None
            and float(r["rfperm_p"]) < float(p_alpha)
            and float(r.get("rfperm_T") or 0.0) > 0.0
        ),
        None,
    )
    return {
        "onset_hat": hop_t,
        "onset_rank": rank_t,
        "n_rfperm_hop": int(sum(1 for r in rows if r.get("rfperm_hop"))),
    }
