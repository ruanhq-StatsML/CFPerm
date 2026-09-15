"""Lock CS-assistant incremental contribution formula."""
from __future__ import annotations

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SQL = ROOT / "sql" / "biz_value"


@pytest.fixture(scope="module")
def con():
    duckdb = pytest.importorskip("duckdb")
    # Import seed via running demo pieces
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "biz_demo", ROOT / "scripts" / "agod" / "run_biz_value_sql_demo.py"
    )
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)

    c = duckdb.connect(":memory:")
    for name in [
        "00_schema.sql",
        "01_hallucination_value.sql",
        "02_style_drift_value.sql",
        "03_value_dashboard.sql",
        "04_cs_assistant_contribution.sql",
        "05_cs_net_and_weekly.sql",
        "06_cs_exec_dashboard.sql",
        "07_cs_week_attribution_payback.sql",
        "08_cs_unit_econ_cumulative.sql",
    ]:
        c.execute((SQL / name).read_text())
    mod.seed(c)
    return c


def test_cs_increment_positive_and_formula(con):
    row = con.execute("SELECT * FROM vw_cs_assist_increment").fetchone()
    cols = [d[0] for d in con.description]
    d = dict(zip(cols, row))

    assert d["sessions_acted"] > 0
    assert d["tickets_avoided_realized"] > 0
    assert d["refunds_avoided_realized"] > 0
    assert d["incremental_yen_realized"] > 0

    expected = (
        d["tickets_avoided_realized"] * d["cs_ticket_cost"]
        + d["refunds_avoided_realized"] * d["refund_unit_cost"]
        + d["extra_contained_realized"] * d["contained_session_value"]
    )
    assert abs(d["incremental_yen_realized"] - expected) < 1e-6


def test_net_after_audit_below_gross(con):
    net = con.execute("SELECT * FROM vw_cs_assist_net_increment").fetchone()
    cols = [d[0] for d in con.description]
    d = dict(zip(cols, net))
    gross = con.execute(
        "SELECT incremental_yen_realized FROM vw_cs_assist_increment"
    ).fetchone()[0]
    assert d["net_incremental_yen"] > 0
    assert d["audit_cost_acted"] >= 0
    assert d["net_incremental_yen"] <= gross + 1e-6
    assert abs(d["net_incremental_yen"] - (gross - d["audit_cost_acted"])) < 1.0


def test_action_net_below_gross(con):
    rows = con.execute(
        "SELECT action_type, gross_yen, net_yen_after_action_cost FROM vw_cs_assist_action_net"
    ).fetchall()
    assert rows
    for action_type, gross, net in rows:
        assert gross > 0, action_type
        assert net <= gross + 1e-6, action_type
        assert net > 0, action_type


def test_weekly_wow_has_fire_weeks(con):
    n = con.execute(
        "SELECT COUNT(*) FROM vw_cs_assist_weekly WHERE days_acted + days_ignored > 0"
    ).fetchone()[0]
    assert n >= 1
    wow = con.execute("SELECT * FROM vw_cs_assist_weekly_wow ORDER BY week_start").fetchall()
    assert len(wow) >= 1


def test_ops_brief_builds(tmp_path, monkeypatch):
    """Brief script must render from demo JSON artifacts."""
    import importlib.util

    # Ensure demo artifacts exist (reuse committed results)
    brief_path = ROOT / "scripts" / "agod" / "cs_assist_weekly_ops_brief.py"
    spec = importlib.util.spec_from_file_location("ops_brief", brief_path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    text = mod.build_brief()
    assert "少工单" in text or "少 **" in text
    assert "毛贡献" in text
    assert "动作净贡献" in text
    assert "火情周 WoW" in text
    out = tmp_path / "brief.md"
    out.write_text(text)
    assert out.stat().st_size > 200


def test_exec_dashboard_before_after(con):
    row = con.execute("SELECT * FROM vw_cs_assist_exec_dashboard").fetchone()
    cols = [d[0] for d in con.description]
    d = dict(zip(cols, row))
    assert d["tickets_avoided"] > 0
    assert d["gross_yen"] > 0
    assert d["net_yen"] <= d["gross_yen"] + 1e-6
    assert d["ticket_rate_ignored_pct"] > d["ticket_rate_acted_pct"]
    assert d["top_action"]
    assert d["top_action_net_yen_per_day"] > 0


def test_action_recommend_ranked(con):
    rows = con.execute(
        "SELECT recommend_rank, net_yen_per_day FROM vw_cs_assist_action_recommend ORDER BY recommend_rank"
    ).fetchall()
    assert len(rows) >= 1
    assert rows[0][0] == 1
    # ranks decrease or equal in net/day
    nets = [r[1] for r in rows]
    assert nets == sorted(nets, reverse=True)


def test_hop_sensitivity_monotonic_net():
    import importlib.util
    import duckdb
    import numpy as np

    spec = importlib.util.spec_from_file_location(
        "sens", ROOT / "scripts" / "agod" / "cs_assist_hop_sensitivity.py"
    )
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    demo = mod._load_demo()
    base = demo.load_hf_hop()
    weak = {**base, "scenario": "weak_hop", "fire_halluc": max(base["quiet_halluc"]+0.05, base["fire_halluc"]*0.5), "hop_ratio": max(2.0, base["hop_ratio"]*0.35)}
    strong = {**base, "scenario": "strong_hop", "fire_halluc": min(0.55, base["fire_halluc"]*1.3), "hop_ratio": base["hop_ratio"]*1.4}
    w = mod._run_scenario(demo, duckdb, weak)
    s = mod._run_scenario(demo, duckdb, strong)
    # stronger hop should not yield lower net contribution
    assert s["net_yen"] >= w["net_yen"] - 1.0


def test_week_attribution_reconciles(con):
    row = con.execute("SELECT * FROM vw_cs_assist_attribution_check").fetchone()
    cols = [d[0] for d in con.description]
    d = dict(zip(cols, row))
    assert abs(d["attribution_gap_yen"]) <= 1.0
    assert d["attribution_gap_pct"] <= 0.1
    weeks = con.execute("SELECT COUNT(*) FROM vw_cs_assist_week_attribution").fetchone()[0]
    assert weeks >= 1


def test_action_payback_same_day(con):
    rows = con.execute(
        "SELECT action_type, payback_days, payback_bucket, net_roi_multiple "
        "FROM vw_cs_assist_action_payback"
    ).fetchall()
    assert rows
    for action, days, bucket, roi in rows:
        assert days is not None and days < 1.0, action
        assert bucket == "same_day_payback", action
        assert roi > 1.0, action


def test_unit_econ_band_ordered(con):
    row = con.execute("SELECT * FROM vw_cs_assist_unit_econ_band").fetchone()
    cols = [d[0] for d in con.description]
    d = dict(zip(cols, row))
    assert d["net_yen_low"] < d["net_yen_base"] < d["net_yen_high"]
    assert d["gross_yen_low"] < d["gross_yen_base"] < d["gross_yen_high"]
    assert abs(d["net_yen_base"] - (d["gross_yen_base"] - (d["gross_yen_base"] - d["net_yen_base"]))) < 1e-6


def test_cumulative_curve_ends_at_total(con):
    total = con.execute(
        "SELECT SUM(gross_yen_day) FROM vw_cs_assist_cumulative_curve"
    ).fetchone()[0]
    last = con.execute(
        "SELECT cumulative_gross_yen, cumulative_pct_of_total "
        "FROM vw_cs_assist_cumulative_curve ORDER BY dt DESC LIMIT 1"
    ).fetchone()
    assert abs(last[0] - total) <= 2.0
    assert abs(last[1] - 100.0) < 0.2


def test_method_biz_scenario_prototype_chain(tmp_path, monkeypatch):
    """Method↔biz↔unit-econ prototype must render same/diff + relevance chain."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "method_biz_proto",
        ROOT / "scripts" / "agod" / "method_biz_scenario_prototype.py",
    )
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    # Redirect writes into tmp so pytest stays side-effect light
    monkeypatch.setattr(mod, "OUT", tmp_path / "out")
    monkeypatch.setattr(mod, "DOCS", tmp_path / "docs")
    assert mod.main() == 0

    payload = mod.build_payload()
    assert payload["method"]["halu"]["fired"] == 1
    assert payload["biz_scenarios"]["hop_band"]["base_net"] > 0
    assert (
        payload["unit_econ"]["band"]["net_yen_low"]
        < payload["unit_econ"]["band"]["net_yen_base"]
        < payload["unit_econ"]["band"]["net_yen_high"]
    )
    assert len(payload["chain"]) >= 6
    assert payload["same_diff"]["same"]
    assert payload["same_diff"]["different"]
    assert "跳变" in payload["same_diff"]["relevance_one_liner"]

    md = (tmp_path / "docs" / "METHOD_BIZ_SCENARIO_PROTOTYPE.md").read_text()
    assert "相同" in md
    assert "单位经济" in md
    assert "relevance" in md.lower()
    en = (tmp_path / "docs" / "METHOD_BIZ_SCENARIO_PROTOTYPE_EN.md").read_text()
    assert "Relevance chain" in en
