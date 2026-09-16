"""Tests for bulletin landing attribution board."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from scripts.agod.bulletin_landing_attr import (
    STYLE_NAMES,
    _family_mass,
    _metric_shift,
    _path_metrics,
    _synth_halu,
    _synth_hh,
    main,
    run_agent_path_landing,
    run_audit_landing,
    run_style_landing,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "agod" / "bulletin_landing"


def test_path_metrics_hop():
    batch = np.array([0, 0, 0, 1, 1, 1])
    rng = np.random.default_rng(0)
    p = _path_metrics(6, batch, cut=1, rng=rng)
    assert p.shape == (6, 4)
    assert p[batch >= 1, 0].mean() > p[batch < 1, 0].mean()


def test_metric_shift_orders_by_abs_delta():
    pre = np.zeros((10, 10))
    post = np.zeros((10, 10))
    post[:, 6] = 1.0
    rows = _metric_shift(pre, post, STYLE_NAMES)
    assert rows[0]["metric"] == "formal"
    assert rows[0]["abs_delta"] == 1.0


def test_family_mass_normalized():
    vimp = np.array([0.5, 0.3, 0.2])
    mass = _family_mass(vimp, ["a", "b"], {"a": [0, 1], "b": [2]})
    assert abs(sum(mass.values()) - 1.0) < 1e-9
    assert mass["a"] > mass["b"]


def test_agent_audit_style_synth_board():
    hh, hu = _synth_hh(280), _synth_halu(480)
    agent = run_agent_path_landing(hu, gate=1.2, seed=0)
    audit = run_audit_landing(hh, gate=1.2, seed=0)
    style = run_style_landing(hh, seed=0)

    assert agent["detection"]["first_fire_t"] is not None
    assert agent["attribution"]["top_family"] in ("text_hash", "rag", "style", "path")
    assert abs(sum(agent["attribution"]["logo_share"].values()) - 1.0) < 1e-6
    assert "path_metric_shifts" in agent["attribution"]

    assert audit["detection"]["judge_err_ratio"] > 1.0
    assert audit["detection"]["style_domain_auc"] > 0.7
    assert "preference_vs_style" in audit["attribution"]

    assert style["detection"]["style_domain_auc"] == style["detection"]["style_domain_auc"]
    assert style["attribution"]["top_metric"] in STYLE_NAMES
    assert "另账" in style["action"]["ledger"] or "客服" in style["action"]["ledger"]


def test_main_synth_writes_artifacts(tmp_path, monkeypatch):
    # Redirect OUT by monkeypatching module constant via running with --synth
    # into the real OUT (idempotent) — assert files exist after run.
    rc = main(["--synth", "--gate", "1.2", "--seed", "0"])
    assert rc == 0
    assert (OUT / "summary.json").exists()
    assert (OUT / "REPORT.md").exists()
    summary = json.loads((OUT / "summary.json").read_text())
    assert "bulletin" in summary
    assert set(summary.keys()) >= {
        "agent_path",
        "data_audit",
        "style_drift",
        "landing_map",
    }
    assert summary["agent_path"]["attribution"]["top_family"]
    assert summary["data_audit"]["action"]["primary"]
    assert summary["style_drift"]["attribution"]["text_metric_shifts"]
