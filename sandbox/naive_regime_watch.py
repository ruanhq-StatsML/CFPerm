"""FW-pm25-naive: regime watcher for random-walk packs.

Opportunity
-----------
On packs where bakeoff says ``naive_last`` wins (near RW), do **not** run HGB.
Still watch surprises: a streak of large |y_t - y_{t-1}|/σ is a *regime alert*
(distribution may have changed — then re-bakeoff / reconsider routing).

Stays in sandbox — classical only.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence

import numpy as np

from sandbox.forecast_bakeoff import Pack, load_packs, make_supervised, run_bakeoff
from sandbox.pack_router import route_model_from_card


def naive_regime_watch(
    pack: Pack,
    *,
    warm: int = 150,
    max_steps: int = 250,
    surprise_k: float = 2.5,
    streak: int = 5,
    n_lags: int = 5,
) -> Dict[str, Any]:
    """Walk with naive_last only; alert if surprise streak ≥ ``streak``."""
    Z, y = make_supervised(pack.X, pack.y, n_lags=n_lags)
    if len(y) < warm + streak + 5:
        return {"ok": False, "reason": "too_short", "dataset": pack.name}
    scale = float(np.std(y[:warm]) + 1e-6)
    steps = min(max_steps, len(y) - warm)
    alerts: List[Dict[str, Any]] = []
    run = 0
    n_surprise = 0
    abs_resid = []
    for step in range(steps):
        i = warm + step
        y_hat = float(Z[i, n_lags - 1])  # last lag = y[t-1]
        resid = float(y[i] - y_hat)
        surprise = abs(resid) / scale
        abs_resid.append(abs(resid))
        action = "idle"
        if surprise >= surprise_k:
            n_surprise += 1
            run += 1
            action = "log_surprise"
            if run >= streak:
                action = "regime_alert"
                alerts.append(
                    {
                        "step": step,
                        "t": int(i),
                        "surprise": surprise,
                        "streak": run,
                        "decision_path": [
                            "observe",
                            "forecast:naive_last",
                            "flag_surprise",
                            "decide:regime_alert",
                            "suggest:rebakeoff_router",
                            "log",
                        ],
                    }
                )
                run = 0  # reset after alert
        else:
            run = 0
        # refresh scale slowly
        if (step + 1) % 50 == 0:
            scale = float(np.std(y[: i + 1]) + 1e-6)
    return {
        "ok": True,
        "dataset": pack.name,
        "model": "naive_last",
        "n_steps": steps,
        "mae": float(np.mean(abs_resid)) if abs_resid else float("nan"),
        "n_surprise": n_surprise,
        "surprise_rate": n_surprise / max(steps, 1),
        "n_regime_alerts": len(alerts),
        "alerts_head": alerts[:5],
        "note": "RW-pack watcher — no HGB; regime_alert ⇒ re-bakeoff router",
    }


def naive_watch_suite(
    packs: Sequence[Pack],
    *,
    bakeoff_cards: Optional[Sequence[Dict[str, Any]]] = None,
    seed: int = 0,
) -> Dict[str, Any]:
    """Run regime watch on packs routed to naive (or force all for contrast)."""
    by = {c["dataset"]: c for c in (bakeoff_cards or []) if c.get("ok")}
    rows = []
    for p in packs:
        card = by.get(p.name)
        routed = route_model_from_card(card) if card else "naive_last"
        # Always run naive watch; flag whether router agrees this is a naive pack
        watch = naive_regime_watch(p)
        if not watch.get("ok"):
            continue
        rows.append(
            {
                **{k: watch[k] for k in watch if k != "alerts_head"},
                "router_model": routed,
                "is_naive_route": routed == "naive_last",
                "alerts_head": watch.get("alerts_head"),
            }
        )
    naive_rows = [r for r in rows if r.get("is_naive_route")]
    return {
        "rows": rows,
        "n_naive_routed": len(naive_rows),
        "total_regime_alerts": int(sum(r.get("n_regime_alerts", 0) for r in rows)),
        "naive_regime_alerts": int(
            sum(r.get("n_regime_alerts", 0) for r in naive_rows)
        ),
        "headline": (
            f"FW-pm25-naive watch: {len(naive_rows)} naive-routed packs; "
            f"regime_alerts on naive="
            f"{sum(r.get('n_regime_alerts', 0) for r in naive_rows)}, "
            f"all_packs={sum(r.get('n_regime_alerts', 0) for r in rows)}"
        ),
        "opportunity": "FW-pm25-naive",
    }


def suite_from_loaders(*, max_n: int = 5000, seed: int = 0) -> Dict[str, Any]:
    packs = load_packs(max_n=max_n)
    cards = [run_bakeoff(p, seed=seed) for p in packs]
    return naive_watch_suite(packs, bakeoff_cards=cards, seed=seed)
