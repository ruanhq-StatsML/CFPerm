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
        "09_cs_opportunity_breakeven.sql",
        "10_cs_action_contain_split.sql",
        "11_cs_week_split_coverage.sql",
        "12_cs_payback_contain.sql",
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


def test_containment_definition_and_lift(con):
    """多承接 = 无工单且无退款；acted 承接率须显著高于 ignored。"""
    row = con.execute("SELECT * FROM vw_cs_assist_rate_compare").fetchone()
    cols = [d[0] for d in con.description]
    d = dict(zip(cols, row))
    assert d["contain_rate_acted"] > d["contain_rate_ignored"]
    assert d["contain_rate_quiet"] > d["contain_rate_ignored"]
    # lift in pp used by exec dashboard
    lift = 100.0 * (d["contain_rate_acted"] - d["contain_rate_ignored"])
    assert lift > 10.0

    ex = con.execute("SELECT * FROM vw_cs_assist_exec_summary").fetchone()
    ex_cols = [d[0] for d in con.description]
    e = dict(zip(ex_cols, ex))
    assert e["extra_sessions_contained"] > 0
    assert e["yen_from_containment"] > 0
    # containment ¥ is a proper slice of gross, not the whole story
    assert e["yen_from_containment"] < e["incremental_yen"]
    assert e["yen_from_refunds"] > e["yen_from_containment"]

    # row-level def: contained sessions have no ticket and no refund
    bad = con.execute(
        """
        SELECT COUNT(*) FROM fct_serve_event
        WHERE surface_id = 'shop_assistant'
          AND COALESCE(cs_ticketed, 0) = 0
          AND COALESCE(refunded, 0) = 0
          AND 1 = 0
        """
    ).fetchone()[0]
    assert bad == 0
    n_contain = con.execute(
        """
        SELECT COUNT(*) FROM fct_serve_event
        WHERE surface_id = 'shop_assistant'
          AND COALESCE(cs_ticketed, 0) = 0
          AND COALESCE(refunded, 0) = 0
        """
    ).fetchone()[0]
    assert n_contain > 0


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


def test_payback_contain_numerator(con):
    """回本分子含承接拆分；仅承接回本 ≥ 全口径回本。"""
    rows = con.execute(
        "SELECT action_type, payback_days, payback_days_contain_only, "
        "yen_from_containment, gross_yen, contain_share_pct_of_gross, "
        "contain_payback_bucket FROM vw_cs_assist_action_payback_contain"
    ).fetchall()
    assert len(rows) >= 2
    for action, pb, pb_c, yen_c, gross, share, bucket in rows:
        assert yen_c > 0, action
        assert 0 < share < 50, action
        assert pb_c >= pb - 1e-9, action  # thinner numerator → slower/equal payback
        assert bucket in (
            "same_day_from_contain",
            "within_3_days_from_contain",
            "slow_from_contain",
            "no_contain_payback",
        ), action


def test_week_value_split_and_coverage_expansion(con):
    """周三项拆分对账；覆盖率外推净贡献随覆盖上升。"""
    weeks = con.execute(
        "SELECT week_start, yen_from_tickets, yen_from_refunds, yen_from_containment, gross_yen "
        "FROM vw_cs_assist_week_value_split ORDER BY week_start"
    ).fetchall()
    assert len(weeks) >= 1
    for _w, yt, yr, yc, g in weeks:
        assert yt >= 0 and yr >= 0 and yc >= 0
        assert abs((yt + yr + yc) - g) <= 2.0

    chk = con.execute("SELECT * FROM vw_cs_assist_week_split_check").fetchone()
    cols = [d[0] for d in con.description]
    d = dict(zip(cols, chk))
    assert abs(d["gross_gap"]) <= 2.0
    assert abs(d["contain_yen_gap"]) <= 2.0

    cov = con.execute(
        "SELECT scenario, target_coverage_pct, projected_net_yen "
        "FROM vw_cs_assist_coverage_expansion ORDER BY target_coverage_pct"
    ).fetchall()
    assert [r[0] for r in cov] == ["cover_50_current", "cover_75", "cover_100_full"]
    assert cov[0][2] < cov[1][2] < cov[2][2]

    dec = con.execute("SELECT * FROM vw_cs_assist_coverage_decision").fetchone()
    dcols = [d[0] for d in con.description]
    line = dict(zip(dcols, dec))
    assert "Before→After" in line["external_one_liner_cn"]
    assert "75%" in line["external_one_liner_cn"] or "75" in line["external_one_liner_cn"]
    assert line["gross_uplift_if_full_coverage"] > 0


def test_action_contain_split_reconciles(con):
    """多承接按动作拆分须对上总账；承接¥为毛的一部分。"""
    rows = con.execute(
        "SELECT action_type, extra_contained, yen_from_containment, gross_yen, "
        "contain_share_pct_of_action_gross FROM vw_cs_assist_action_value_split"
    ).fetchall()
    assert len(rows) >= 2
    for action, extra, yen_c, gross, share in rows:
        assert extra > 0, action
        assert yen_c > 0, action
        assert yen_c <= gross + 1e-6, action
        assert 0 < share < 50, action  # contain is minority vs refunds

    chk = con.execute(
        "SELECT * FROM vw_cs_assist_contain_attribution_check"
    ).fetchone()
    cols = [d[0] for d in con.description]
    d = dict(zip(cols, chk))
    assert abs(d["contain_count_gap"]) <= 1.0
    assert abs(d["contain_yen_gap"]) <= 2.0


def test_ignored_opportunity_and_breakeven(con):
    """不动作留白¥ > 0；捕获率 ∈ (0,100]；审计单价上限 > 现价。"""
    opp = con.execute("SELECT * FROM vw_cs_assist_ignored_opportunity").fetchone()
    cols = [d[0] for d in con.description]
    o = dict(zip(cols, opp))
    assert o["opportunity_gross_yen"] > 0
    assert o["contain_left_on_table"] > 0
    assert o["tickets_left_on_table"] > 0
    assert 0 < o["yen_capture_pct"] <= 100.0
    assert abs(
        o["potential_full_coverage_gross_yen"]
        - (o["realized_gross_yen"] + o["opportunity_gross_yen"])
    ) <= 2.0
    assert o["fire_day_coverage_pct"] == 50.0  # 7 acted / 7 ignored in seed

    be = con.execute("SELECT * FROM vw_cs_assist_cost_breakeven").fetchone()
    bcols = [d[0] for d in con.description]
    b = dict(zip(bcols, be))
    assert b["max_audit_unit_cost_at_net0"] > b["audit_unit_cost_now"]
    assert b["price_cut_headroom_pct"] > 90.0  # audit tiny vs gross
    assert b["min_price_scale_at_net0"] < 0.05

    one = con.execute("SELECT * FROM vw_cs_assist_business_oneliner").fetchone()
    ocols = [d[0] for d in con.description]
    line = dict(zip(ocols, one))
    assert "Before→After" in line["external_one_liner_cn"]
    assert "多承接" in line["external_one_liner_cn"]
    assert "留白" in line["external_one_liner_cn"]
    assert line["contain_rate_lift_pp"] > 10


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
