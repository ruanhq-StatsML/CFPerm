"""180-agent continuous-learning fleet: even force + specialty cohorts.

RSI pivot (away from soft-burn FLOPs polish)
------------------------------------------
Use-case: continuous learning that *empowers* a large agent fleet so
capacity is spent evenly across specialty lanes — e.g. ~10 antifraud,
~10 CUPED — instead of everyone piling onto one shiny metric.

Score (commercial, not math-Eval)
---------------------------------
  Score = 1{shipped} · (1{wired} + 1{used} + 1{roi_crude})

Fairness
--------
  Gini / HHI on WIP or completed scores across cohorts.
  Anti-pile tax if any cohort attracts > soft_cap concurrent agents.

Continuous learning
-------------------
  Tick t: observe feedback (fraud useful-rate, CUPED var-ratio),
  update cohort priors, reallocate idle agents toward under-served
  high-marginal-value lanes (not toward the currently loudest lane).
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

# Default specialty mix ≈ 180 agents (illustrative; captain can reweight).
DEFAULT_COHORTS: Dict[str, int] = {
    "antifraud": 10,  # wave→K*→clue→ticket / watchlist
    "cuped": 10,  # A/B variance reduction + experiment ops
    "review_accel": 20,  # human review cards / useful feedback
    "localize_fsds": 20,  # Drill / FSDS / tip+direction
    "board_stream": 15,  # adjacent-board / stream packs
    "eta_gate": 10,  # ETA behavior under R_t=1
    "gray_rollback": 15,  # L1/L2 flags + rollback
    "data_contract": 15,  # OpenAPI / schemas / audit
    "roi_ops": 15,  # cost dashboards / weekly ROI bots
    "flex_reserve": 50,  # float capacity; never >30% on one specialty
}

SOFT_CAP_FRAC = 0.12  # no single specialty should hold >12% concurrent WIP


@dataclass
class Agent:
    agent_id: int
    cohort: str
    wip: int = 0
    shipped: int = 0
    wired: int = 0
    used: int = 0
    roi: int = 0
    note: str = ""

    @property
    def score(self) -> int:
        if self.shipped <= 0:
            return 0
        return int(self.wired > 0) + int(self.used > 0) + int(self.roi > 0)


@dataclass
class FleetState:
    agents: List[Agent]
    tick: int = 0
    history: List[Dict[str, Any]] = field(default_factory=list)

    @property
    def n(self) -> int:
        return len(self.agents)


def build_fleet(
    cohorts: Optional[Mapping[str, int]] = None,
    *,
    seed: int = 0,
) -> FleetState:
    """Allocate agents to specialty cohorts (even-force prior)."""
    cohorts = dict(cohorts or DEFAULT_COHORTS)
    agents: List[Agent] = []
    i = 0
    for name, k in cohorts.items():
        for _ in range(int(k)):
            agents.append(Agent(agent_id=i, cohort=str(name)))
            i += 1
    return FleetState(agents=agents)


def cohort_counts(state: FleetState) -> Dict[str, int]:
    return dict(Counter(a.cohort for a in state.agents))


def gini(x: Sequence[float]) -> float:
    """Gini coefficient on non-negative loads (0=even, 1=all on one)."""
    a = np.asarray(x, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return float("nan")
    a = np.clip(a, 0.0, None)
    if float(a.sum()) <= 0:
        return 0.0
    a = np.sort(a)
    n = a.size
    idx = np.arange(1, n + 1, dtype=float)
    return float((2.0 * np.sum(idx * a) / (n * a.sum())) - (n + 1.0) / n)


def hhi(shares: Sequence[float]) -> float:
    """Herfindahl–Hirschman on share vector (1/n → even; 1 → monopoly)."""
    s = np.asarray(shares, dtype=float)
    s = s[np.isfinite(s)]
    if s.size == 0 or float(s.sum()) <= 0:
        return float("nan")
    p = s / s.sum()
    return float(np.sum(p * p))


def fairness_report(state: FleetState) -> Dict[str, Any]:
    """Even-force diagnostics across cohorts."""
    by = defaultdict(list)
    for a in state.agents:
        by[a.cohort].append(a)
    loads = []
    rows = []
    for c, members in sorted(by.items()):
        wip = sum(m.wip for m in members)
        scores = [m.score for m in members]
        loads.append(float(wip + 1e-9))
        rows.append(
            {
                "cohort": c,
                "n_agents": len(members),
                "wip": int(wip),
                "mean_score": float(np.mean(scores)) if scores else 0.0,
                "n_shipped": int(sum(m.shipped for m in members)),
                "share_agents": len(members) / max(state.n, 1),
            }
        )
    shares = [r["n_agents"] for r in rows]
    wip_shares = [float(r["wip"]) + 1e-9 for r in rows]
    return {
        "n_agents": state.n,
        "n_cohorts": len(rows),
        "gini_wip": gini(wip_shares),
        "hhi_headcount": hhi(shares),
        "hhi_wip": hhi(wip_shares),
        "even_hhi_ref": 1.0 / max(len(rows), 1),
        "soft_cap_frac": SOFT_CAP_FRAC,
        "cohorts": rows,
    }


def pile_on_baseline(
    n_agents: int = 180,
    *,
    magnet: str = "antifraud",
    seed: int = 0,
) -> FleetState:
    """Anti-pattern: almost everyone piles onto one loud specialty."""
    rng = np.random.default_rng(seed)
    agents = []
    for i in range(n_agents):
        # 70% magnet, rest sprinkle
        c = magnet if rng.random() < 0.70 else rng.choice(
            list(DEFAULT_COHORTS.keys())
        )
        agents.append(Agent(agent_id=i, cohort=str(c), wip=int(rng.integers(0, 3))))
    return FleetState(agents=agents)


def even_force_assign_wip(
    state: FleetState,
    *,
    tasks_per_cohort: Optional[Mapping[str, int]] = None,
    seed: int = 0,
) -> FleetState:
    """WIP=1 per agent inside soft-cap; excess tasks stay queued (not stolen)."""
    rng = np.random.default_rng(seed)
    tasks = dict(tasks_per_cohort or {c: max(1, n // 2) for c, n in cohort_counts(state).items()})
    by = defaultdict(list)
    for a in state.agents:
        a.wip = 0
        by[a.cohort].append(a)
    for c, members in by.items():
        soft = max(1, int(np.floor(SOFT_CAP_FRAC * state.n)))
        cap = min(len(members), soft, int(tasks.get(c, 0)))
        pick = list(members)
        rng.shuffle(pick)
        for a in pick[:cap]:
            a.wip = 1
    return state


# ---- Continuous learning signals -----------------------------------------


def cuped_variance_ratio(
    y: np.ndarray,
    x_pre: np.ndarray,
) -> Dict[str, float]:
    """Classic CUPED: Y_adj = Y - θ(X − Ē[X]); report Var(Y_adj)/Var(Y).

    θ = Cov(Y,X)/Var(X).  Ratio < 1 ⇒ variance reduced (experiment ops win).
    """
    y = np.asarray(y, dtype=float).ravel()
    x = np.asarray(x_pre, dtype=float).ravel()
    n = min(y.size, x.size)
    y, x = y[:n], x[:n]
    if n < 3:
        return {"var_y": float("nan"), "var_cuped": float("nan"), "ratio": float("nan"), "theta": float("nan")}
    vx = float(np.var(x))
    vy = float(np.var(y))
    if vx <= 1e-18:
        return {"var_y": vy, "var_cuped": vy, "ratio": 1.0, "theta": 0.0}
    theta = float(np.cov(y, x, ddof=0)[0, 1] / vx)
    y_adj = y - theta * (x - float(np.mean(x)))
    vadj = float(np.var(y_adj))
    return {
        "var_y": vy,
        "var_cuped": vadj,
        "ratio": float(vadj / vy) if vy > 0 else float("nan"),
        "theta": theta,
        "n": float(n),
    }


def antifraud_useful_rate(
    useful: Sequence[int],
    *,
    baseline: float = 0.35,
) -> Dict[str, float]:
    """Human 'useful' clicks on clue cards — commercial eval, not AUC."""
    u = np.asarray(list(useful), dtype=float)
    if u.size == 0:
        return {"rate": float("nan"), "lift_vs_baseline": float("nan"), "n": 0.0}
    rate = float(np.mean(u))
    return {
        "rate": rate,
        "lift_vs_baseline": float(rate - baseline),
        "n": float(u.size),
        "baseline": float(baseline),
    }


def continuous_learn_tick(
    state: FleetState,
    *,
    fraud_useful: Sequence[int],
    cuped_y: np.ndarray,
    cuped_x: np.ndarray,
    seed: int = 0,
) -> Dict[str, Any]:
    """One CL tick: observe specialty feedback, update scores, rebalance WIP.

    Reflection rule (RSI):
      - If fraud useful lift > 0 and cuped ratio < 1 → keep specialty mix.
      - If gini_wip high → peel WIP from crowded cohort into flex_reserve
        then re-home toward under-scored specialties (not the magnet).
    """
    rng = np.random.default_rng(seed + state.tick)
    fraud = antifraud_useful_rate(fraud_useful)
    cuped = cuped_variance_ratio(cuped_y, cuped_x)
    fair0 = fairness_report(state)

    # Credit shipped agents in specialty from feedback quality
    for a in state.agents:
        if a.wip <= 0:
            continue
        if a.cohort == "antifraud" and np.isfinite(fraud.get("lift_vs_baseline", np.nan)):
            a.shipped = 1
            a.wired = 1
            a.used = int(fraud["lift_vs_baseline"] > 0)
            a.roi = int(fraud["lift_vs_baseline"] > 0.05)
        elif a.cohort == "cuped" and np.isfinite(cuped.get("ratio", np.nan)):
            a.shipped = 1
            a.wired = 1
            a.used = int(cuped["ratio"] < 1.0)
            a.roi = int(cuped["ratio"] < 0.85)
        else:
            # other lanes: small random commercial progress (smoke)
            a.shipped = 1
            a.wired = int(rng.random() < 0.8)
            a.used = int(rng.random() < 0.5)
            a.roi = int(rng.random() < 0.3)
        a.wip = 0

    # Rebalance: if WIP concentration was high, nudge headcount via flex
    fair1 = fairness_report(state)
    reflection = _reflect(fair0, fair1, fraud, cuped)
    if reflection.get("action") == "rehome_flex_to_underserved":
        _rehome_flex(state, underserved=reflection.get("underserved") or [], rng=rng)

    state.tick += 1
    snap = {
        "tick": state.tick,
        "fraud": fraud,
        "cuped": cuped,
        "fairness_before": {k: fair0[k] for k in ("gini_wip", "hhi_wip", "hhi_headcount")},
        "fairness_after": {k: fair1[k] for k in ("gini_wip", "hhi_wip", "hhi_headcount")},
        "reflection": reflection,
        "mean_score": float(np.mean([a.score for a in state.agents])),
    }
    state.history.append(snap)
    return snap


def _reflect(
    fair0: Mapping[str, Any],
    fair1: Mapping[str, Any],
    fraud: Mapping[str, float],
    cuped: Mapping[str, float],
) -> Dict[str, Any]:
    cohorts = fair1.get("cohorts") or fair0.get("cohorts") or []
    scores = [(c["cohort"], c["mean_score"]) for c in cohorts if c["cohort"] != "flex_reserve"]
    scores_sorted = sorted(scores, key=lambda t: t[1])
    underserved = [c for c, s in scores_sorted[:3]]
    crowded = None
    for c in cohorts:
        if c["share_agents"] > SOFT_CAP_FRAC + 1e-9 and c["cohort"] != "flex_reserve":
            # headcount over soft cap — only flex should be large
            if c["cohort"] not in DEFAULT_COHORTS or DEFAULT_COHORTS.get(c["cohort"], 0) > int(
                SOFT_CAP_FRAC * fair0.get("n_agents", 180)
            ):
                pass
        if c["wip"] > 0 and c["n_agents"] > 0:
            if c["wip"] / max(c["n_agents"], 1) > 0.9 and c["cohort"] not in (
                "antifraud",
                "cuped",
            ):
                crowded = c["cohort"]
    gini0 = float(fair0.get("gini_wip", 0) or 0)
    keep_mix = bool(
        (fraud.get("lift_vs_baseline") or 0) >= 0
        and (cuped.get("ratio") or 1) <= 1.0
        and gini0 < 0.45
    )
    if keep_mix:
        return {
            "action": "keep_specialty_mix",
            "reason": "fraud lift≥0, CUPED ratio≤1, WIP gini ok — average force holds",
            "underserved": underserved,
        }
    return {
        "action": "rehome_flex_to_underserved",
        "reason": "feedback or concentration asks for rebalance toward low-score specialties",
        "underserved": underserved,
        "crowded": crowded,
    }


def _rehome_flex(
    state: FleetState,
    *,
    underserved: Sequence[str],
    rng: np.random.Generator,
    n_move: int = 5,
) -> None:
    flex = [a for a in state.agents if a.cohort == "flex_reserve"]
    if not flex or not underserved:
        return
    rng.shuffle(flex)
    targets = list(underserved)
    for i, a in enumerate(flex[:n_move]):
        a.cohort = str(targets[i % len(targets)])
        a.note = "rehomed_from_flex"


def compare_even_vs_pile(
    *,
    n_agents: int = 180,
    seed: int = 0,
) -> Dict[str, Any]:
    """Scorecard: even-force fleet vs pile-on magnet (antifraud)."""
    even = build_fleet(seed=seed)
    even_force_assign_wip(even, seed=seed)
    pile = pile_on_baseline(n_agents, seed=seed)
    fair_even_wip = fairness_report(even)
    fair_pile_wip = fairness_report(pile)
    # synthetic CL feedback favoring specialty competence
    rng = np.random.default_rng(seed)
    useful_even = (rng.random(80) < 0.55).astype(int)
    useful_pile = (rng.random(80) < 0.38).astype(int)
    x = rng.normal(size=500)
    y = 0.6 * x + rng.normal(scale=0.8, size=500)  # CUPED helps
    snap_even = continuous_learn_tick(
        even, fraud_useful=useful_even, cuped_y=y, cuped_x=x, seed=seed
    )
    snap_pile = continuous_learn_tick(
        pile, fraud_useful=useful_pile, cuped_y=y, cuped_x=x, seed=seed + 1
    )
    return {
        "even": {
            "fairness": fair_even_wip,
            "fairness_after_tick": fairness_report(even),
            "tick": snap_even,
            "mean_score": snap_even["mean_score"],
            "cohort_counts": cohort_counts(even),
        },
        "pile_on": {
            "fairness": fair_pile_wip,
            "fairness_after_tick": fairness_report(pile),
            "tick": snap_pile,
            "mean_score": snap_pile["mean_score"],
            "cohort_counts": cohort_counts(pile),
        },
        "verdict": _verdict(fair_even_wip, fair_pile_wip, snap_even, snap_pile),
    }


def _verdict(
    fe: Mapping[str, Any],
    fp: Mapping[str, Any],
    se: Mapping[str, Any],
    sp: Mapping[str, Any],
) -> str:
    parts = []
    if fe["hhi_headcount"] < fp["hhi_headcount"]:
        parts.append(
            f"even-force lower headcount HHI ({fe['hhi_headcount']:.3f}<{fp['hhi_headcount']:.3f})"
        )
    if np.isfinite(fe.get("gini_wip", np.nan)) and np.isfinite(fp.get("gini_wip", np.nan)):
        if fe["gini_wip"] + 1e-9 < fp["gini_wip"]:
            parts.append("even-force lower WIP Gini")
    if se["mean_score"] >= sp["mean_score"]:
        parts.append("even-force mean commercial score ≥ pile-on")
    fr = (se.get("fraud") or {}).get("lift_vs_baseline")
    if fr is not None and fr > 0:
        parts.append("antifraud useful lift>0 under specialty-10")
    cr = (se.get("cuped") or {}).get("ratio")
    if cr is not None and cr < 1:
        parts.append(f"CUPED var ratio={cr:.2f}<1")
    if not parts:
        return "inconclusive — check seeds / feedback"
    return "; ".join(parts)


def multi_tick_cl_demo(
    *,
    n_ticks: int = 5,
    seed: int = 0,
) -> Dict[str, Any]:
    """Reflective CL loop: starve cuped feedback early → rehome → recover ratio use."""
    st = build_fleet(seed=seed)
    rng = np.random.default_rng(seed)
    snaps = []
    for t in range(n_ticks):
        even_force_assign_wip(st, seed=seed + t)
        # Early ticks: weak fraud; later: recover. CUPED always informative.
        p_useful = 0.25 if t < 2 else 0.60
        useful = (rng.random(60) < p_useful).astype(int)
        x = rng.normal(size=400)
        # ticks 0-1: weak pre-period; later strong → ratio drops
        strength = 0.1 if t < 2 else 0.65
        y = strength * x + rng.normal(scale=0.9, size=400)
        snaps.append(
            continuous_learn_tick(
                st, fraud_useful=useful, cuped_y=y, cuped_x=x, seed=seed + t
            )
        )
    actions = [s["reflection"]["action"] for s in snaps]
    return {
        "n_ticks": n_ticks,
        "actions": actions,
        "final_cohort_counts": cohort_counts(st),
        "mean_scores": [s["mean_score"] for s in snaps],
        "cuped_ratios": [s["cuped"]["ratio"] for s in snaps],
        "fraud_lifts": [s["fraud"]["lift_vs_baseline"] for s in snaps],
        "reading": (
            "early weak feedback → rehome_flex_to_underserved; "
            "later fraud lift↑ / CUPED ratio↓ → keep_specialty_mix"
            if any(a == "rehome_flex_to_underserved" for a in actions)
            else "mix held; check seed if no rehome observed"
        ),
        "snaps": snaps,
    }


def fleet_observation_spec() -> Dict[str, str]:
    return {
        "X": "per-cohort feedback features (fraud useful, CUPED pre-period X, WIP)",
        "Y": "commercial Score = shipped·(wired+used+roi)",
        "intermediate": "fleet Thought = {cohort mix, gini, rehome action} — not Ŷ",
        "N": "≈180 agents; soft_cap≈12% per specialty (flex_reserve excepted)",
    }
