"""Tracking, OOF features, and π variation — protocol, not meta-loss.

Half-stream expert switch (Herbster & Warmuth 1998):
    static regret vs the best *fixed* expert is the wrong oracle when the
    leader changes at T/2. Switching regret vs a k-shift comparator is
    the right one. Discrete OSL and Hedge with share=0 lock on the
    first-half winner. Fixed-Share, sliding-window OSL, and sticky
    mixing are three ways to *spend* π-variation on the switch.

OOF features (Wolpert 1992; van der Laan Super Learner):
    the meta-learner may only see votes that were not trained on the
    scored points. Streaming analogue: probe vs holdout. Extra columns
    in the level-1 matrix each need their own OOF; stuffing in-sample
    stats is leaky stacking. Putting x or t into the gate is how this
    slides toward MoE — don't.

π variation:
    TV(π_t, π_{t+1}) = ½ ‖π_{t+1}−π_t‖₁ is the cost of tracking.
    Stationary stream + large TV_path = overfitting the last window
    (actuator chatter). Switch stream: want a *localized* TV spike
    after T/2, not uniform jitter. Two-timescale: fast π for the
    meta-loss, sticky π for the LR actuator (FWD stays on).
"""
from __future__ import annotations

from collections import deque
from typing import Mapping, Sequence

import numpy as np

from .online_stacking import EPS, OnlineStacker, _simplex, honest_step, project_simplex


def tv(p: Mapping[str, float] | np.ndarray, q: Mapping[str, float] | np.ndarray) -> float:
    """Total variation between two simplex weights: ½ ‖p−q‖₁."""
    if isinstance(p, Mapping) and isinstance(q, Mapping):
        keys = list(p.keys())
        a = np.array([float(p[k]) for k in keys], float)
        b = np.array([float(q[k]) for k in keys], float)
    else:
        a = np.asarray(p, float).ravel()
        b = np.asarray(q, float).ravel()
    return float(0.5 * np.abs(a - b).sum())


def path_tv(traj: Sequence[Mapping[str, float]]) -> float:
    if len(traj) < 2:
        return 0.0
    return float(sum(tv(traj[i], traj[i + 1]) for i in range(len(traj) - 1)))


def tv_after_frac(traj: Sequence[Mapping[str, float]], t_switch: int) -> float:
    """Share of path TV that occurs at or after the switch index."""
    total = path_tv(traj)
    if total <= EPS or len(traj) < 2:
        return float("nan")
    after = float(sum(tv(traj[i], traj[i + 1]) for i in range(t_switch, len(traj) - 1)))
    return after / total


def sticky_mix(
    pi_slow: Mapping[str, float],
    pi_fast: Mapping[str, float],
    experts: Sequence[str],
    *,
    lam: float,
) -> dict[str, float]:
    """Two-timescale π: π_slow ← (1−λ) π_slow + λ π_fast.

    λ=1 is the raw meta-update (chatters). λ→0 freezes the actuator.
    """
    experts = list(experts)
    lam = float(np.clip(lam, 0.0, 1.0))
    w = (1.0 - lam) * np.array([float(pi_slow[e]) for e in experts])
    w = w + lam * np.array([float(pi_fast[e]) for e in experts])
    w = _simplex(w)
    return {e: float(w[i]) for i, e in enumerate(experts)}


def switch_delay(
    traj: Sequence[Mapping[str, float]],
    *,
    new_leader: str,
    t_switch: int,
    thresh: float = 0.45,
) -> int | None:
    """First t ≥ t_switch with π[new_leader] ≥ thresh. None = lock-in."""
    for t in range(max(int(t_switch), 0), len(traj)):
        if float(traj[t].get(new_leader, 0.0)) >= float(thresh):
            return int(t)
    return None


def locked_on(
    traj: Sequence[Mapping[str, float]],
    *,
    old_leader: str,
    new_leader: str,
) -> bool:
    """End of stream still prefers the first-half winner."""
    if not traj:
        return False
    last = traj[-1]
    return float(last.get(old_leader, 0.0)) > float(last.get(new_leader, 0.0))


class RollingOSL:
    """Discrete OSL on a sliding window of expert losses (forgets first half)."""

    def __init__(self, experts: Sequence[str], *, window: int):
        self.experts = list(experts)
        self.window = int(max(window, 1))
        self.buf: deque[np.ndarray] = deque(maxlen=self.window)
        self.w = np.full(len(self.experts), 1.0 / max(len(self.experts), 1))
        self.history: list[dict[str, float]] = []

    def pi(self) -> dict[str, float]:
        return {e: float(self.w[i]) for i, e in enumerate(self.experts)}

    def update(self, expert_losses: Mapping[str, float]) -> dict[str, float]:
        loss = np.array([float(expert_losses[e]) for e in self.experts], float)
        self.buf.append(loss)
        tot = np.sum(np.stack(list(self.buf), axis=0), axis=0)
        self.w = np.zeros(len(self.experts))
        self.w[int(np.argmin(tot))] = 1.0
        self.history.append(self.pi())
        return self.pi()


def ridge_fit(x: np.ndarray, y: np.ndarray, *, lam: float = 0.3) -> np.ndarray:
    x = np.asarray(x, float)
    y = np.asarray(y, float).ravel()
    d = x.shape[1]
    a = x.T @ x + float(lam) * np.eye(d)
    return np.linalg.solve(a, x.T @ y)


def kfold_oof_predict(
    x: np.ndarray,
    y: np.ndarray,
    *,
    n_folds: int = 4,
    lam: float = 0.3,
) -> np.ndarray:
    """Level-1 OOF column: each row is predicted from a model that never saw it."""
    x = np.asarray(x, float)
    y = np.asarray(y, float).ravel()
    n = len(y)
    n_folds = int(max(min(n_folds, n), 2))
    pred = np.zeros(n, float)
    idx = np.arange(n)
    folds = np.array_split(idx, n_folds)
    for i, te in enumerate(folds):
        tr = np.concatenate([folds[j] for j in range(n_folds) if j != i])
        w = ridge_fit(x[tr], y[tr], lam=lam)
        pred[te] = x[te] @ w
    return pred


def leaky_in_sample_predict(x: np.ndarray, y: np.ndarray, *, lam: float = 0.3) -> np.ndarray:
    """In-sample votes — Wolpert's forbidden stacking features."""
    w = ridge_fit(x, y, lam=lam)
    return np.asarray(x, float) @ w


def oof_vs_leaky_stack(
    x_good: np.ndarray,
    x_noise: np.ndarray,
    y: np.ndarray,
    *,
    n_folds: int = 4,
    lam: float = 0.3,
) -> dict:
    """Fit π on OOF vs leaky level-1 features; score a held-out tail.

    Expert ``good`` sees the signal columns. Expert ``noise`` sees junk
    columns that can *memorize* in-sample (high-d noise) but fail OOF.
    """
    y = np.asarray(y, float).ravel()
    n = len(y)
    n_fit = int(0.75 * n)
    y_fit, y_te = y[:n_fit], y[n_fit:]
    xg, xn = np.asarray(x_good, float), np.asarray(x_noise, float)

    z_oof = np.column_stack(
        [
            kfold_oof_predict(xg[:n_fit], y_fit, n_folds=n_folds, lam=lam),
            kfold_oof_predict(xn[:n_fit], y_fit, n_folds=n_folds, lam=lam),
        ]
    )
    z_leak = np.column_stack(
        [
            leaky_in_sample_predict(xg[:n_fit], y_fit, lam=lam),
            leaky_in_sample_predict(xn[:n_fit], y_fit, lam=lam),
        ]
    )

    def _pi_ols(z, yy):
        # unconstrained LS then simplex — the stacking meta
        ones = np.ones((len(yy), 1))
        # drop intercept; π on the two votes
        g = z.T @ z + 1e-6 * np.eye(z.shape[1])
        raw = np.linalg.solve(g, z.T @ yy)
        w = project_simplex(raw)
        return w

    pi_oof = _pi_ols(z_oof, y_fit)
    pi_leak = _pi_ols(z_leak, y_fit)

    def _corr(a, b):
        a = np.asarray(a, float).ravel()
        b = np.asarray(b, float).ravel()
        if float(a.std()) < 1e-12 or float(b.std()) < 1e-12:
            return 0.0
        return float(np.corrcoef(a, b)[0, 1])

    w_g = ridge_fit(xg[:n_fit], y_fit, lam=lam)
    w_n = ridge_fit(xn[:n_fit], y_fit, lam=lam)
    z_te = np.column_stack([xg[n_fit:] @ w_g, xn[n_fit:] @ w_n])

    def _mse(z, pi, yy):
        return float(np.mean((z @ pi - yy) ** 2))

    in_oof = _mse(z_oof, pi_oof, y_fit)
    in_leak = _mse(z_leak, pi_leak, y_fit)
    te_oof = _mse(z_te, pi_oof, y_te)
    te_leak = _mse(z_te, pi_leak, y_te)
    return {
        "pi_oof": {"good": float(pi_oof[0]), "noise": float(pi_oof[1])},
        "pi_leak": {"good": float(pi_leak[0]), "noise": float(pi_leak[1])},
        "in_sample_oof": in_oof,
        "in_sample_leaky": in_leak,
        "holdout_oof": te_oof,
        "holdout_leaky": te_leak,
        "optimism_oof": float(te_oof - in_oof),
        "optimism_leaky": float(te_leak - in_leak),
        "corr_oof_good": _corr(z_oof[:, 0], y_fit),
        "corr_oof_noise": _corr(z_oof[:, 1], y_fit),
        "corr_leak_good": _corr(z_leak[:, 0], y_fit),
        "corr_leak_noise": _corr(z_leak[:, 1], y_fit),
    }


def disjoint_probe_holdout(
    n: int,
    *,
    probe_frac: float = 0.5,
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Honest gradient protocol: votes from P_t, target from disjoint H_t."""
    rng = rng or np.random.default_rng(0)
    idx = rng.permutation(int(n))
    k = int(max(1, min(n - 1, round(probe_frac * n))))
    return idx[:k], idx[k:]


def alignment_leak(
    g_probe: np.ndarray,
    g_hold: np.ndarray,
    g_same: np.ndarray,
) -> dict:
    """Using the same-batch gradient as both vote and target inflates cosine."""

    def _cos(a, b):
        a = np.asarray(a, float).ravel()
        b = np.asarray(b, float).ravel()
        na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
        if na < EPS or nb < EPS:
            return 0.0
        return float(np.dot(a, b) / (na * nb))

    return {
        "honest_cos": _cos(g_probe, g_hold),
        "leaky_cos": _cos(g_same, g_same),
        "gap": float(_cos(g_same, g_same) - _cos(g_probe, g_hold)),
    }


def run_switch_methods(
    votes: np.ndarray,
    y: np.ndarray,
    names: Sequence[str],
    *,
    t_switch: int,
    eta: float = 0.85,
) -> dict:
    """Honest paths for share=0 / Fixed-Share / window OSL / discrete OSL."""
    names = list(names)
    y = np.asarray(y, float).ravel()
    votes = np.asarray(votes, float)
    specs = {
        "hedge_share0": dict(method="hedge", eta=eta, share=0.0),
        "hedge_share": dict(method="hedge", eta=eta, share=0.08),
        "osl_disc": dict(method="osl_disc", eta=eta, share=0.0),
        "osl_sgd": dict(method="osl_sgd", eta=min(eta, 0.35), share=0.0),
    }
    out: dict = {}
    for key, kw in specs.items():
        st = OnlineStacker(names, **kw)
        traj = []
        preq = []
        for t in range(len(y)):
            v = {n: float(votes[t, i]) for i, n in enumerate(names)}
            scored = honest_step(st, v, float(y[t]))
            traj.append(scored["pi"])
            preq.append(scored["preq_loss"])
        out[key] = _pack_path(traj, preq, names, t_switch)

    roll = RollingOSL(names, window=max(8, len(y) // 5))
    traj, preq = [], []
    w0 = {n: 1.0 / len(names) for n in names}
    for t in range(len(y)):
        v = {n: float(votes[t, i]) for i, n in enumerate(names)}
        yhat = sum(w0[n] * v[n] for n in names)
        preq.append(float((yhat - float(y[t])) ** 2))
        traj.append(dict(w0))
        losses = {n: float((v[n] - float(y[t])) ** 2) for n in names}
        w0 = roll.update(losses)
    out["osl_window"] = _pack_path(traj, preq, names, t_switch)
    return out


def _n_eff(pi: Mapping[str, float], names: Sequence[str]) -> float:
    p = np.array([float(pi[n]) for n in names], float)
    p = np.clip(p, EPS, None)
    p = p / p.sum()
    return float(np.exp(-np.sum(p * np.log(p))))


def _pack_path(
    traj: Sequence[Mapping[str, float]],
    preq: Sequence[float],
    names: Sequence[str],
    t_switch: int,
) -> dict:
    names = list(names)
    new_leader = names[1] if len(names) > 1 else names[0]
    delay = switch_delay(traj, new_leader=new_leader, t_switch=t_switch, thresh=0.45)
    final = dict(traj[-1]) if traj else {}
    final_leader = max(names, key=lambda n: float(final.get(n, 0.0))) if final else ""
    return {
        "traj": [dict(p) for p in traj],
        "preq_mean": float(np.mean(preq)) if preq else float("nan"),
        "final_pi": final,
        "final_leader": final_leader,
        "path_tv": path_tv(traj),
        "tv_after_frac": tv_after_frac(traj, t_switch),
        "delay": delay,
        "delay_from_switch": None if delay is None else int(delay - t_switch),
        "locked": final_leader != new_leader,
        "n_eff_end": _n_eff(traj[-1], names) if traj else float("nan"),
    }


def sticky_path(
    traj_fast: Sequence[Mapping[str, float]],
    experts: Sequence[str],
    *,
    lam: float,
) -> list[dict[str, float]]:
    experts = list(experts)
    slow = {e: 1.0 / len(experts) for e in experts}
    out = []
    for fast in traj_fast:
        slow = sticky_mix(slow, fast, experts, lam=lam)
        out.append(dict(slow))
    return out
