"""Tests for soft IPTW burn policy + FLOPs ledger."""
from __future__ import annotations

import json
from pathlib import Path

from agod.soft_burn import (
    burn_soft_weights,
    soft_flops_ledger,
    summarize_burn_decisions,
    weighting_flops_proxy,
)

ROOT = Path(__file__).resolve().parents[1]


def test_weighting_negligible_vs_adapt():
    led = soft_flops_ledger(adapt_flops=1e6, n_rows=100)
    assert led["weighting_flops"] == 100
    assert led["weighting_share"] < 1e-3


def test_burn_when_beats_uniform():
    d = burn_soft_weights(
        rel_gated_sqrt=0.93, rel_gated_cbrt=0.95, soft_win_cbrt_le_sqrt=False
    )
    assert d["burn"] is True
    assert d["decision"] == "BURN_SQRT"


def test_soften_only_when_hurts():
    d = burn_soft_weights(
        rel_gated_sqrt=1.1,
        rel_gated_cbrt=1.05,
        soft_win_cbrt_le_sqrt=True,
        gate_already_on=True,
    )
    assert d["burn"] is False
    assert d["decision"] == "SOFTEN_ONLY"
    assert abs(d["alpha"] - 1.0 / 3.0) < 1e-9


def test_keep_uniform_when_gate_off():
    d = burn_soft_weights(
        rel_gated_sqrt=1.1,
        rel_gated_cbrt=1.05,
        gate_already_on=False,
    )
    assert d["decision"] == "KEEP_UNIFORM"
    assert d["burn"] is False


def test_repo_summary_smoke():
    path = ROOT / "results/agod_po_cbrt/summary.json"
    if not path.exists():
        return
    from agod.po_power_eff import power_scorecard_from_summary

    power = power_scorecard_from_summary(json.loads(path.read_text()))
    rep = summarize_burn_decisions(power["cards"])
    assert rep["n_packs"] >= 1
    assert "decision_counts" in rep
    assert weighting_flops_proxy(10) == 10.0
