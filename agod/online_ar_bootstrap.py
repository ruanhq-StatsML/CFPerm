"""Palm & Nagler online AR-bootstrap, null-adjusted on a frozen reference.

Protocol (same as the concept-drift batch-size board)::

    fit a probe on D_ref
    μ_ref = mean score of that frozen probe on D_ref mini-batches
    s_t   = score of the same frozen probe on trail batch t
    Δ_t   = s_t − μ_ref
    update OnlineARBootstrap(Δ_t)
    fire iff CI_lo > 0

This is *not* the last-two OnlineRFPerm gate. The probe here is frozen on
``D_ref``; last-two refits every hop and fires on ``e_now / e_prev ≥ γ``.
Mixing parameter ``β = √2 − 1`` is Palm & Nagler.

The class itself is scoring-agnostic: feed it ``Δ_t`` (or any scalar
stream). Classification Brier / 0-1 / PO-risk are all valid ``s_t``.
"""
from __future__ import annotations

import numpy as np

BETA = float(2.0 ** 0.5 - 1.0)


class OnlineARBootstrap:
    """Recursive AR(1) bootstrap of a running mean (Palm & Nagler)."""

    def __init__(self, n_boot: int = 200, seed: int = 0):
        self.n_boot = int(n_boot)
        self._rng = np.random.default_rng(seed)
        self.t = 0
        self.V = np.zeros(self.n_boot)
        self.Vbar = np.zeros(self.n_boot)
        self.Xstar = np.zeros(self.n_boot)
        self.mean = np.nan

    def update(self, x: float) -> None:
        t = self.t + 1
        rho = float(np.clip(1.0 - t ** (-BETA), 0.0, 1.0 - 1e-12))
        zeta = self._rng.normal(size=self.n_boot)
        self.V = 1.0 + rho * (self.V - 1.0) + np.sqrt(max(0.0, 1.0 - rho * rho)) * zeta
        if t == 1:
            self.Xstar = np.full(self.n_boot, float(x))
            self.Vbar = self.V.copy()
            self.mean = float(x)
        else:
            num = (t - 1) * self.Vbar * self.Xstar + float(x) * self.V
            den = (t - 1) * self.Vbar + self.V
            ok = np.abs(den) > 1e-15
            nxt = np.full(self.n_boot, float(x))
            nxt[ok] = num[ok] / den[ok]
            self.Xstar = nxt
            self.Vbar = (1.0 - 1.0 / t) * self.Vbar + self.V / t
            self.mean = ((t - 1) * float(self.mean) + float(x)) / t
        self.t = t

    def ci(self, alpha: float = 0.05) -> tuple[float, float]:
        if self.t < 2:
            return float("nan"), float("nan")
        lo, hi = np.quantile(self.Xstar, [alpha / 2.0, 1.0 - alpha / 2.0])
        return float(lo), float(hi)

    def rho(self) -> float:
        """Current mixing weight ``1 − t^(−β)``. ``t=0`` → 0."""
        if self.t <= 0:
            return 0.0
        return float(np.clip(1.0 - self.t ** (-BETA), 0.0, 1.0 - 1e-12))


def run_delta_bootstrap(
    delta,
    *,
    n_boot: int = 200,
    seed: int = 0,
    alpha: float = 0.05,
):
    """Walk ``Δ_t`` through the AR-bootstrap. ``fire`` = ``CI_lo > 0``."""
    delta = np.asarray(delta, dtype=float).ravel()
    boot = OnlineARBootstrap(n_boot=n_boot, seed=seed)
    rows = []
    for t, d in enumerate(delta, start=1):
        boot.update(float(d))
        lo, hi = boot.ci(alpha)
        fire = bool(np.isfinite(lo) and lo > 0.0)
        width = float(hi - lo) if np.isfinite(lo) and np.isfinite(hi) else float("nan")
        rows.append(
            {
                "t": int(t),
                "s": float(d),
                "mean": float(boot.mean),
                "lo": float(lo),
                "hi": float(hi),
                "width": width,
                "fire": fire,
            }
        )
    return rows, boot


def first_significant(rows, *, after_t: int = 0):
    """First row with ``t > after_t`` and ``CI_lo > 0``."""
    for r in rows:
        if int(r["t"]) > int(after_t) and r.get("fire"):
            return r
    return None
