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
