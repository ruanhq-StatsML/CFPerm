"""Unit tests for manuscript first_k form."""
from __future__ import annotations

import numpy as np

from agod.first_k_metrics import (
    aggregate_day_forms,
    alarm_time_quantiles,
    first_k_consecutive_rej,
    first_k_form,
)


def test_first_k_oob_counting():
    # reject at indices 0,1,2,5 → first1=1, first2=2, first3=3 (1-based end)
    rej = [1, 1, 1, 0, 0, 1, 0]
    assert first_k_consecutive_rej(rej, 1) == 1.0
    assert first_k_consecutive_rej(rej, 2) == 2.0
    assert first_k_consecutive_rej(rej, 3) == 3.0
    assert np.isnan(first_k_consecutive_rej(rej, 4))


def test_first_k_never():
    assert np.isnan(first_k_consecutive_rej([0, 0, 0], 1))


def test_alarm_quantiles():
    rej = [0, 1, 0, 1, 1, 0]  # times 1,3,4
    q = alarm_time_quantiles(rej)
    assert q["n_alarm"] == 3
    assert q["median"] == 3.0
    assert q["p25"] <= q["median"] <= q["p75"]


def test_day_aggregate_detrate():
    forms = [
        {"SUM": 2, "first1": 1.0, "first2": 2.0, "first3": float("nan")},
        {"SUM": 0, "first1": float("nan"), "first2": float("nan"), "first3": float("nan")},
        {"SUM": 1, "first1": 3.0, "first2": float("nan"), "first3": float("nan")},
    ]
    agg = aggregate_day_forms(forms)
    assert agg["n_days"] == 3
    assert agg["n_alarm_days"] == 2
    assert abs(agg["det_rate"] - 2 / 3) < 1e-9


def test_first_k_form_keys():
    f = first_k_form([0, 1, 1, 0])
    for k in ("SUM", "first1", "first2", "first3", "alarm_p25", "alarm_median", "alarm_p75"):
        assert k in f
