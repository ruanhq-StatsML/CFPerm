"""Modality-agnostic online stacking of named experts.

The previous GLS snapshot ``π ∝ R^{-1} α`` is *static* stacking. Online
stacking is a meta-learner on the simplex that updates **after** scoring
the current window with frozen weights (prequential / one-step-ahead CV).

Experts are just names. They may be modalities, towers, seeds, or anything
that emits a vote and a scalar loss. Nothing here knows about video/text.

Honesty (the part that is not "easy code"):
    score with π_t *before* ``update``. Training the bases on the same
    window then voting is leaky stacking (in-sample meta-features).
    This is the streaming analogue of out-of-fold Super Learner.

Meta-learners (literature map in ``docs/agod/AGOD_online_stacking.md``):
    equal       — 1/|E| (forecast-combination baseline)
    hedge       — multiplicative weights (Freund & Schapire)
    osl_disc    — discrete Online Super Learner (best cumulative expert)
    osl_sgd     — convex OSL, projected SGD on stacked loss (Benkeser et al.)
    bg          — Bates–Granger inverse-MSE (ignore covariance)
    gls_ewma    — min-var combo from EWMA error covariance (Newbold–Granger)
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

EPS = 1e-12
METHODS = ("equal", "hedge", "osl_disc", "osl_sgd", "bg", "gls_ewma")


def _simplex(w: np.ndarray) -> np.ndarray:
    w = np.clip(np.asarray(w, float), 0.0, None)
    s = float(w.sum())
    if s <= EPS:
        return np.full(len(w), 1.0 / max(len(w), 1))
    return w / s


def _gls_weights_from_cov(cov: np.ndarray, *, shrink: float = 0.25) -> np.ndarray:
    """Min-var combo π ∝ Σ̃^{-1} 1 with eigenvalue shrinkage toward mean."""
    c = 0.5 * (np.asarray(cov, float) + np.asarray(cov, float).T)
    evals, evecs = np.linalg.eigh(c)
    evals = np.clip(evals, 1e-10, None)
    mu = float(np.mean(evals))
    evals = (1.0 - float(shrink)) * evals + float(shrink) * mu
    inv = (evecs * (1.0 / evals)) @ evecs.T
    raw = inv @ np.ones(c.shape[0])
    return _simplex(raw)


def project_simplex(v: np.ndarray) -> np.ndarray:
    """Euclidean projection onto {w ≥ 0, 1ᵀw = 1} (Duchi et al. 2008)."""
    v = np.asarray(v, float).ravel()
    n = v.size
    u = np.sort(v)[::-1]
    cssv = np.cumsum(u)
    rho = np.nonzero(u * np.arange(1, n + 1) > (cssv - 1.0))[0]
    if len(rho) == 0:
        return _simplex(np.ones(n))
    theta = (cssv[rho[-1]] - 1.0) / float(rho[-1] + 1)
    return np.clip(v - theta, 0.0, None)


class OnlineStacker:
    """Prequential meta-learner over named experts.

    Typical honest step::

        scored = stacker.score(votes, y)          # frozen π_t
        stacker.update(votes, expert_losses, y)   # π_{t+1}
    """

    def __init__(
        self,
        experts: Sequence[str],
        *,
        method: str = "hedge",
        eta: float = 0.35,
        ewma: float = 0.20,
        share: float = 0.06,
        prior: Mapping[str, float] | None = None,
    ):
        if method not in METHODS:
            raise ValueError(f"unknown method {method!r}; expected {METHODS}")
        self.experts = list(experts)
        self.method = str(method)
        self.eta = float(eta)
        self.ewma = float(ewma)
        self.share = float(share)
        n = len(self.experts)
        if prior:
            w0 = np.array([max(float(prior.get(e, 0.0)), 0.0) for e in self.experts])
            self.w = _simplex(w0)
        else:
            self.w = np.full(n, 1.0 / n)
        self.cum_loss = np.zeros(n, float)
        self.mse = np.ones(n, float)
        self.cov = np.eye(n, dtype=float)
        self.t = 0
        self.history: list[dict] = []

    def pi(self) -> dict[str, float]:
        return {e: float(self.w[i]) for i, e in enumerate(self.experts)}

    def n_eff(self) -> float:
        p = np.clip(self.w, EPS, None)
        p = p / p.sum()
        return float(np.exp(-np.sum(p * np.log(p))))

    def _vec(self, d: Mapping[str, float], default: float = 0.0) -> np.ndarray:
        return np.array([float(d.get(e, default)) for e in self.experts], float)

    def score(
        self,
        votes: Mapping[str, float],
        y: float | None = None,
    ) -> dict:
        """Evaluate frozen π_t. Does **not** write weights."""
        v = self._vec(votes)
        yhat = float(self.w @ v)
        out = {
            "t": self.t,
            "pi": self.pi(),
            "yhat": yhat,
            "n_eff": self.n_eff(),
            "method": self.method,
        }
        if y is not None:
            out["y"] = float(y)
            out["preq_loss"] = float((yhat - float(y)) ** 2)
            err = v - float(y)
            out["expert_se"] = {e: float(err[i] ** 2) for i, e in enumerate(self.experts)}
        return out

    def update(
        self,
        votes: Mapping[str, float],
        expert_losses: Mapping[str, float] | None = None,
        y: float | None = None,
    ) -> dict[str, float]:
        """Write π_{t+1} from the window that was *already* scored."""
        v = self._vec(votes)
        if expert_losses is None:
            if y is None:
                raise ValueError("update needs expert_losses or y")
            loss = (v - float(y)) ** 2
        else:
            loss = self._vec(expert_losses, default=1.0)
        loss = np.clip(loss, 0.0, None)
        self.cum_loss += loss
        lam = self.ewma
        self.mse = (1.0 - lam) * self.mse + lam * loss

        if y is not None:
            err = (v - float(y)).reshape(-1, 1)
            ee = err @ err.T
            self.cov = (1.0 - lam) * self.cov + lam * ee

        if self.method == "equal":
            self.w = np.full(len(self.experts), 1.0 / len(self.experts))
        elif self.method == "hedge":
            # multiplicative weights + Fixed-Share mix (Herbster & Warmuth 1998)
            m = float(loss.max()) if len(loss) else 1.0
            z = loss / max(m, EPS)
            w_hat = _simplex(self.w * np.exp(-self.eta * z))
            n = len(self.experts)
            self.w = (1.0 - self.share) * w_hat + self.share / n
        elif self.method == "osl_disc":
            k = int(np.argmin(self.cum_loss))
            self.w = np.zeros(len(self.experts))
            self.w[k] = 1.0
        elif self.method == "osl_sgd":
            if y is None:
                # fall back to loss-weighted hedge-like step
                g = loss - float(loss.mean())
            else:
                yhat = float(self.w @ v)
                g = 2.0 * (yhat - float(y)) * v
            self.w = project_simplex(self.w - self.eta * g)
        elif self.method == "bg":
            # Bates–Granger: inverse MSE, ignore off-diagonal (estimation)
            self.w = _simplex(1.0 / np.clip(self.mse, EPS, None))
        else:  # gls_ewma
            self.w = _gls_weights_from_cov(self.cov, shrink=0.25)

        self.t += 1
        snap = {
            "t": self.t,
            "pi": self.pi(),
            "n_eff": self.n_eff(),
            "cum_loss": {e: float(self.cum_loss[i]) for i, e in enumerate(self.experts)},
        }
        self.history.append(snap)
        return self.pi()


def honest_step(
    stacker: OnlineStacker,
    votes: Mapping[str, float],
    y: float,
    expert_losses: Mapping[str, float] | None = None,
) -> dict:
    """Prequential step: score with frozen π, then update (Dawid 1984)."""
    scored = stacker.score(votes, y)
    stacker.update(votes, expert_losses, y)
    return scored


def leaky_step(
    stacker: OnlineStacker,
    votes: Mapping[str, float],
    y: float,
    expert_losses: Mapping[str, float] | None = None,
) -> dict:
    """Leaky stacking: update first, then score with the *new* π (in-sample)."""
    stacker.update(votes, expert_losses, y)
    return stacker.score(votes, y)


def oracle_convex_combo(
    votes: np.ndarray,
    y: np.ndarray,
    *,
    n_grid: int = 17,
) -> dict:
    """Hindsight best *fixed* convex combo on a simplex grid (K small)."""
    votes = np.asarray(votes, float)
    y = np.asarray(y, float).ravel()
    t, k = votes.shape
    if k == 1:
        w = np.array([1.0])
        pred = votes[:, 0]
        return {"pi": w, "loss": float(np.mean((pred - y) ** 2))}
    if k == 2:
        grid = np.linspace(0, 1, n_grid)
        best, bw = np.inf, np.array([0.5, 0.5])
        for a in grid:
            w = np.array([a, 1.0 - a])
            pred = votes @ w
            L = float(np.mean((pred - y) ** 2))
            if L < best:
                best, bw = L, w
        return {"pi": bw, "loss": best}
    # k=3: barycentric grid
    best, bw = np.inf, np.full(k, 1.0 / k)
    xs = np.linspace(0, 1, n_grid)
    for a in xs:
        for b in xs:
            c = 1.0 - a - b
            if c < -1e-9:
                continue
            w = np.array([a, b, max(c, 0.0)])
            w = w / w.sum()
            pred = votes @ w
            L = float(np.mean((pred - y) ** 2))
            if L < best:
                best, bw = L, w
    return {"pi": bw, "loss": best}


def best_expert_loss(votes: np.ndarray, y: np.ndarray) -> dict:
    se = (votes - y.reshape(-1, 1)) ** 2
    mean = se.mean(axis=0)
    k = int(np.argmin(mean))
    return {"index": k, "loss": float(mean[k]), "per_expert": mean.tolist()}


def prequential_regret(preq_losses: Sequence[float], oracle_loss: float) -> float:
    if not preq_losses:
        return float("nan")
    return float(np.mean(preq_losses) - oracle_loss)


def pi_to_lr(
    pi: Mapping[str, float],
    experts: Sequence[str],
    *,
    beta: float = 0.10,
) -> dict[str, float]:
    """Actuator: stacking masses → next-stage LR multipliers (FWD on)."""
    from .lr_controller import alpha_to_lr

    return alpha_to_lr(pi, list(experts), beta=beta)
