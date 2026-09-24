"""Classical-forecast × agent flywheel (sandbox; not AGOD).

Flywheel
--------
  observe pack window → classical forecast → residual/surprise
       ↑                                         │
       └──── log / optional light re-fit ← decide ┘

Opportunities (where value shows up)
------------------------------------
1. Pack routing: metro (HGB wins) ≠ pm25 (naive wins) → different agents
2. Surprise ticks: |residual|/RMSE_train large → escalate / retrain
3. Model switch: if rolling lift(HGB)<0 for K steps → fall back to naive
4. Theme agents (DiffDB): cluster id as context for generation/review agents
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

from sandbox.forecast_bakeoff import Pack, make_supervised


@dataclass
class FlywheelTick:
    t: int
    y_true: float
    y_hat: float
    residual: float
    surprise: float
    action: str
    model: str


@dataclass
class FlywheelState:
    model_name: str = "hgb"
    ticks: List[FlywheelTick] = field(default_factory=list)
    n_retrain: int = 0
    n_surprise: int = 0
    n_fallback: int = 0


def opportunity_map_from_bakeoff(bakeoff: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Translate bakeoff cards into concrete agent-flywheel opportunities."""
    opps: List[Dict[str, Any]] = []
    for c in bakeoff.get("cards") or []:
        if not c.get("ok"):
            continue
        name = c["dataset"]
        best = c["best_rmse"]
        lift = float(c["rmse_lift_vs_naive"].get("hgb", float("nan")))
        if best == "hgb" and lift > 0.15:
            opps.append(
                {
                    "id": f"FW-{name}-hgb",
                    "pack": name,
                    "opportunity": "specialize an HGB forecast agent",
                    "why": f"HGB RMSE lift vs naive={lift:.3f}; classical signal is real",
                    "flywheel": "surprise→retrain HGB on recent window; idle when calm",
                    "priority": "P0",
                }
            )
        elif best == "naive_last" or (np.isfinite(lift) and lift < 0.02):
            opps.append(
                {
                    "id": f"FW-{name}-naive",
                    "pack": name,
                    "opportunity": "do NOT spend an HGB agent here",
                    "why": "near random-walk — last-value is enough; HGB adds noise/cost",
                    "flywheel": "naive agent only; surprise still logged for regime change",
                    "priority": "P1",
                }
            )
        elif best == "ridge":
            opps.append(
                {
                    "id": f"FW-{name}-ridge",
                    "pack": name,
                    "opportunity": "cheap linear forecast agent",
                    "why": "ridge beats both naive and HGB on holdout RMSE",
                    "flywheel": "ridge online update / periodic refit; HGB as shadow only",
                    "priority": "P0",
                }
            )
    if len(opps) >= 2:
        opps.append(
            {
                "id": "FW-router",
                "pack": "*",
                "opportunity": "pack→model router agent",
                "why": "best model differs by pack (hgb/naive/ridge) — one global HGB is wrong",
                "flywheel": "router reads last bakeoff card → dispatches specialist",
                "priority": "P0",
            }
        )
    return opps


class _ModelBox:
    """Fit-once predictor; refit only when asked."""

    def __init__(self, name: str, n_lags: int, seed: int = 0):
        self.name = name
        self.n_lags = n_lags
        self.seed = seed
        self._model: Any = None
        self._scaler: Optional[StandardScaler] = None

    def fit(self, Z: np.ndarray, y: np.ndarray) -> None:
        if self.name == "naive_last":
            self._model = "naive"
            return
        if self.name == "ridge":
            self._scaler = StandardScaler()
            Zs = self._scaler.fit_transform(Z)
            m = Ridge(alpha=1.0, random_state=self.seed)
            m.fit(Zs, y)
            self._model = m
            return
        m = HistGradientBoostingRegressor(
            max_depth=3, learning_rate=0.1, max_iter=80, random_state=self.seed
        )
        m.fit(Z, y)
        self._model = m

    def predict_one(self, z_row: np.ndarray) -> float:
        if self.name == "naive_last" or self._model == "naive":
            return float(z_row[self.n_lags - 1])
        if self.name == "ridge":
            assert self._scaler is not None
            return float(
                self._model.predict(self._scaler.transform(z_row.reshape(1, -1)))[0]
            )
        return float(self._model.predict(z_row.reshape(1, -1))[0])


def run_flywheel(
    pack: Pack,
    *,
    model_name: str = "hgb",
    n_lags: int = 5,
    warm: int = 200,
    surprise_k: float = 2.5,
    retrain_every: int = 50,
    seed: int = 0,
    max_steps: int = 300,
) -> Dict[str, Any]:
    """Walk the stream: forecast → surprise → decide (idle/retrain/fallback)."""
    Z, y = make_supervised(pack.X, pack.y, n_lags=n_lags)
    if len(y) < warm + 30:
        return {"ok": False, "reason": "too_short", "dataset": pack.name}
    state = FlywheelState(model_name=model_name)
    tr_end = warm
    scale = float(np.std(y[:warm]) + 1e-6)
    fallback = False
    steps = min(max_steps, len(y) - warm)
    lift_window: List[float] = []

    box = _ModelBox(model_name, n_lags=n_lags, seed=seed)
    box.fit(Z[:tr_end], y[:tr_end])
    naive_box = _ModelBox("naive_last", n_lags=n_lags, seed=seed)
    naive_box.fit(Z[:tr_end], y[:tr_end])

    for step in range(steps):
        i = warm + step
        use = "naive_last" if fallback else state.model_name
        predictor = naive_box if fallback else box
        y_hat = predictor.predict_one(Z[i])
        resid = float(y[i] - y_hat)
        surprise = abs(resid) / scale
        action = "idle"
        if surprise >= surprise_k:
            action = "log_surprise"
            state.n_surprise += 1
        if (step + 1) % retrain_every == 0 and not fallback:
            action = "retrain"
            tr_end = i
            box.fit(Z[:tr_end], y[:tr_end])
            state.n_retrain += 1
            scale = float(np.std(y[:tr_end]) + 1e-6)
        naive_hat = naive_box.predict_one(Z[i])
        lift_window.append(abs(y[i] - naive_hat) - abs(y[i] - y_hat))
        if len(lift_window) >= 40:
            mean_lift = float(np.mean(lift_window[-40:]))
            if mean_lift < 0 and not fallback and state.model_name == "hgb":
                action = "fallback_naive"
                fallback = True
                state.n_fallback += 1
        state.ticks.append(
            FlywheelTick(
                t=i,
                y_true=float(y[i]),
                y_hat=y_hat,
                residual=resid,
                surprise=surprise,
                action=action,
                model=use,
            )
        )

    mae = float(np.mean([abs(t.residual) for t in state.ticks]))
    return {
        "ok": True,
        "dataset": pack.name,
        "model": model_name,
        "n_steps": len(state.ticks),
        "mae": mae,
        "n_surprise": state.n_surprise,
        "n_retrain": state.n_retrain,
        "n_fallback": state.n_fallback,
        "surprise_rate": state.n_surprise / max(len(state.ticks), 1),
        "fell_back_to_naive": fallback,
        "note": "classical forecast flywheel — observe→predict→surprise→retrain/fallback",
    }


def flywheel_suite(
    bakeoff: Dict[str, Any], packs: Sequence[Pack], seed: int = 0
) -> Dict[str, Any]:
    opps = opportunity_map_from_bakeoff(bakeoff)
    by_name = {c["dataset"]: c for c in bakeoff.get("cards") or [] if c.get("ok")}
    runs = []
    for p in packs:
        card = by_name.get(p.name)
        model = (card or {}).get("best_rmse") or "hgb"
        if model not in ("hgb", "ridge", "naive_last"):
            model = "hgb"
        runs.append(run_flywheel(p, model_name=model, seed=seed))
    return {
        "opportunities": opps,
        "runs": runs,
        "headline": _headline(opps, runs),
    }


def _headline(opps: List[Dict[str, Any]], runs: List[Dict[str, Any]]) -> str:
    p0 = sum(1 for o in opps if o.get("priority") == "P0")
    fb = sum(1 for r in runs if r.get("fell_back_to_naive"))
    sur = float(
        np.mean([r.get("surprise_rate", 0) for r in runs if r.get("ok")] or [0])
    )
    return (
        f"{p0} P0 flywheel opportunities; "
        f"mean surprise_rate={sur:.3f}; "
        f"{fb} packs fell back to naive online"
    )
