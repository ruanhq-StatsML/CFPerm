"""Latency clocks for gradual / continuous-time concept drift.

A slow walk has no hop. Adjacent-window ratios go to 1, so last-two
hop_fires is the wrong WHEN. Latency is not one delay; it is several
clocks plus the excess loss paid while waiting.

No online-bootstrap. T is a batch label, not a treatment.
"""
from __future__ import annotations

from typing import Iterable

import numpy as np


def first_index(flags: Iterable, *, start: int = 0) -> int | None:
    """First True at or after start. None = never (infinite delay)."""
    for i, f in enumerate(flags):
        if i >= int(start) and f:
            return int(i)
    return None


def delay_batches(hat: int | None, t0: int) -> float | None:
    if hat is None:
        return None
    return float(hat - int(t0))


def delay_obs(hat: int | None, t0: int, n_new: int) -> float | None:
    d = delay_batches(hat, t0)
    if d is None:
        return None
    return float(d) * float(n_new)


def alpha_at(alpha_batch, hat: int | None) -> float | None:
    if hat is None:
        return None
    a = np.asarray(alpha_batch, dtype=float).ravel()
    if hat < 0 or hat >= len(a):
        return None
    return float(a[hat])


def excess_area(loss, oracle, t0: int, hat: int | None) -> float:
    """Σ (L − L_oracle) from labeled onset to detection (or to the end).

    This is the evaluable latency: extra serving loss paid while waiting.
    Point-delay is ill-posed when there is no hop.
    """
    a = np.asarray(loss, dtype=float).ravel()
    b = np.asarray(oracle, dtype=float).ravel()
    t0 = int(t0)
    hi = len(a) if hat is None else int(hat)
    hi = max(t0, min(hi, len(a)))
    if hi <= t0:
        return 0.0
    return float(np.sum(a[t0:hi] - b[t0:hi]))


def consecutive(flags, k: int = 2):
    """True once the last k flags are True (causal)."""
    f = np.asarray(list(flags), dtype=bool)
    out = np.zeros(len(f), dtype=bool)
    run = 0
    for i, v in enumerate(f):
        run = run + 1 if v else 0
        out[i] = run >= int(k)
    return out


def pre_onset_false_alarms(flags, t0: int) -> int:
    f = np.asarray(list(flags), dtype=bool)
    t0 = int(t0)
    return int(np.sum(f[:t0])) if t0 > 0 else 0


def summarize_detector(name, flags, *, t0, n_new, alpha_batch, loss, oracle) -> dict:
    hat = first_index(flags, start=t0)
    return {
        "detector": name,
        "hat_batch": hat,
        "delay_batches": delay_batches(hat, t0),
        "delay_obs": delay_obs(hat, t0, n_new),
        "alpha_at_hat": alpha_at(alpha_batch, hat),
        "false_alarms_pre": pre_onset_false_alarms(flags, t0),
        "never": hat is None,
        "area_until_hat": excess_area(loss, oracle, t0, hat),
        "area_full": excess_area(loss, oracle, t0, None),
    }


# Hop and continuous time read the same score S(τ).
# Hop is S(τ)/S(τ−W) ≥ g. Level is S(τ) vs S_ref.
# The remaining coupling is serving error vs share error.

SERVE_QUIET_SHARE_QUIET = "quiet"
SERVE_QUIET_SHARE_LOUD = "share_only"
SERVE_LOUD_SHARE_QUIET = "serve_only"
SERVE_LOUD_SHARE_LOUD = "both"


def disagreement_cell(serve_loud: bool, share_loud: bool) -> str:
    """Serving rent vs localization pointer. Four cells, not a unique decomp."""
    if serve_loud and share_loud:
        return SERVE_LOUD_SHARE_LOUD
    if share_loud:
        return SERVE_QUIET_SHARE_LOUD
    if serve_loud:
        return SERVE_LOUD_SHARE_QUIET
    return SERVE_QUIET_SHARE_QUIET


def share_error(pi_true: float, *, after_onset: bool) -> float:
    """How wrong the localization pointer is.

    After onset, true group should take the mass: 1 − π.
    Before onset, mass on that group is a false pointer: π.
    This is not Shapley error. It is pointer noise.
    """
    p = float(np.clip(pi_true, 0.0, 1.0))
    return float(1.0 - p) if after_onset else p


def lead_lag(hat_share: int | None, hat_serve: int | None) -> dict:
    """Share clock minus serving clock. None if either never fires."""
    if hat_share is None or hat_serve is None:
        return {
            "lag_batches": None,
            "share_first": None,
            "never_share": hat_share is None,
            "never_serve": hat_serve is None,
        }
    d = int(hat_share) - int(hat_serve)
    return {
        "lag_batches": d,
        "share_first": d < 0,
        "never_share": False,
        "never_serve": False,
    }
