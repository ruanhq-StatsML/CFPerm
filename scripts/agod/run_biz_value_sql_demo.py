#!/usr/bin/env python3
"""Seed + run business-value SQL (hallucination + style drift).

Requires: pip install duckdb

Writes results/agod/biz_value_sql/{roi_*.json,dashboard.json,REPORT.md}
"""
from __future__ import annotations

import json
from datetime import date, datetime, timedelta
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SQL_DIR = ROOT / "sql" / "biz_value"
OUT = ROOT / "results" / "agod" / "biz_value_sql"
HF_HOP_PATH = ROOT / "results" / "agod" / "hf_landing" / "halu_regime_rag.json"


def _require_duckdb():
    try:
        import duckdb  # noqa: F401
    except ImportError as e:
        raise SystemExit(
            "duckdb is required: pip install duckdb\n" + str(e)
        ) from e
    import duckdb

    return duckdb


def load_hf_hop(path: Path = HF_HOP_PATH) -> dict:
    """Map HaluEval hop JSON → seed fire/rag intensity knobs."""
    defaults = {
        "source": "defaults",
        "quiet_halluc": 0.07,
        "fire_halluc": 0.41,
        "acted_halluc_scale": 0.40,
        "hop_ratio": 3.0,
        "rag_low": 0.22,
        "rag_ok": 0.60,
        "rag_threshold": 0.35,
        "precision_at_10": 0.7,
        "fired_at_cut": 0,
    }
    if not path.exists():
        return defaults
    raw = json.loads(path.read_text())
    hop = raw.get("hop_at_cut") or raw.get("hop_at_cut") or {}
    quiet = raw.get("hop_quiet") or raw.get("hop_quiet") or {}
    hop_rank = hop.get("ranking") or {}
    quiet_rank = quiet.get("ranking") or {}
    quiet_h = float(quiet_rank.get("halluc_rate", defaults["quiet_halluc"]))
    fire_h = float(hop_rank.get("halluc_rate", 0.82))
    # Business surface is milder than raw label hop; keep proportional.
    fire_h_biz = float(np.clip(fire_h * 0.5, quiet_h + 0.05, 0.55))
    ratio = float(
        hop.get("ratio")
        or (raw.get("rate_after", 3.0) / max(raw.get("rate_before", 0.08), 1e-6))
    )
    rag_top = float(
        hop_rank.get("mean_rag_hit_top10")
        or hop_rank.get("mean_rag_hit_top10")
        or 1.0
    )
    # High HF top-10 rag ⇒ generation hop dominates; still keep a low-rag
    # retrieval-refresh arm on every 3rd fire day for action split.
    rag_ok = float(np.clip(0.35 + 0.4 * rag_top, 0.45, 0.75))
    rag_low = float(np.clip(rag_ok - 0.40, 0.12, 0.30))
    return {
        "source": str(path.relative_to(ROOT)),
        "quiet_halluc": quiet_h,
        "fire_halluc": fire_h_biz,
        "acted_halluc_scale": float(np.clip(1.0 / max(ratio, 1.0), 0.25, 0.55)),
        "hop_ratio": ratio,
        "rag_low": rag_low,
        "rag_ok": rag_ok,
        "rag_threshold": 0.35,
        "precision_at_10": float(hop_rank.get("precision_at_10", 0.7)),
        "fired_at_cut": int(hop.get("fired", 0)),
        "raw_fire_halluc": fire_h,
        "dataset": raw.get("dataset") or raw.get("scenario"),
        "rate_before": float(raw.get("rate_before", quiet_h)),
        "rate_after": float(raw.get("rate_after", fire_h)),
    }


def seed(con, hop: dict | None = None) -> None:
    rng = np.random.default_rng(7)
    start = date(2026, 9, 1)
    days = 28
    hop = hop or load_hf_hop()

    con.execute(
        """
        INSERT INTO dim_surface VALUES
          ('shop_assistant','llm_assist','llm-quality'),
          ('feed_caption','content','content-growth'),
          ('ad_creative','ads','ads-creative')
        """
    )
    con.execute(
        """
        INSERT INTO dim_value_assumption VALUES
          (?, 'shop_assistant', 'cs_ticket_cost', 25.0, 'CNY', 'ops FP 2026Q3'),
          (?, 'shop_assistant', 'refund_unit_cost', 80.0, 'CNY', 'finance 2026Q3'),
          (?, 'shop_assistant', 'contained_session_value', 3.5, 'CNY', 'ops: avoided human handle'),
          (?, 'shop_assistant', 'audit_unit_cost', 6.0, 'CNY', 'ops: human audit per sample'),
          (?, 'shop_assistant', 'action_cost_retrieval_refresh_day', 40.0, 'CNY', 'ops: index/RAG refresh day'),
          (?, 'shop_assistant', 'action_cost_model_rollback_day', 80.0, 'CNY', 'ops: rollback+canary day'),
          (?, 'shop_assistant', 'action_cost_audit_topk_day', 20.0, 'CNY', 'ops: queue setup (excl per-sample)'),
          (?, 'feed_caption', 'value_per_ctr_point', 1200.0, 'CNY', 'growth proxy'),
          (?, 'feed_caption', 'brand_incident_cost', 50.0, 'CNY', 'brand ops'),
          (?, 'ad_creative', 'value_per_ctr_point', 2500.0, 'CNY', 'ads proxy'),
          (?, 'ad_creative', 'brand_incident_cost', 120.0, 'CNY', 'brand ops')
        """,
        [start] * 11,
    )

    # creatives / style clusters
    for i, cluster in enumerate(["clean_min", "loud_promo", "lifestyle", "ugc_raw"]):
        con.execute(
            "INSERT INTO dim_creative VALUES (?, ?, ?, ?, ?)",
            [f"c{i}", f"camp{i%2}", cluster, cluster, start],
        )

    serve_rows = []
    signal_rows = []
    audit_rows = []

    q_h = float(hop["quiet_halluc"])
    f_h = float(hop["fire_halluc"])
    acted_scale = float(hop["acted_halluc_scale"])
    rag_low = float(hop["rag_low"])
    rag_ok = float(hop["rag_ok"])
    rag_thr = float(hop["rag_threshold"])
    # Stronger hop ⇒ larger ticket/refund blow-up when ignored.
    ignore_ticket_bump = float(np.clip(0.03 + 0.01 * hop["hop_ratio"], 0.04, 0.12))
    ignore_refund_bump = float(np.clip(0.02 + 0.005 * hop["hop_ratio"], 0.025, 0.08))

    for d in range(days):
        dt = start + timedelta(days=d)
        # Hallucination surface: concept hop after day 14 (HF cut analogue)
        fired = 1 if d >= 14 else 0
        acted = 1 if (fired and d % 2 == 0) else 0
        # Route action by retrieval health (HF hop RAG split)
        rag = (q_h + 0.48) if d < 14 else (rag_low if d % 3 == 0 else rag_ok)
        if not fired:
            notes = None
            action = None
            act_ticket_cut = 0.0
            act_refund_cut = 0.0
            act_halluc = q_h
        elif not acted:
            notes = "ignored"
            action = None
            act_ticket_cut = 0.0
            act_refund_cut = 0.0
            act_halluc = f_h
        elif rag < rag_thr:
            notes = "acted:retrieval_refresh"
            action = "retrieval_refresh"
            act_ticket_cut = 0.028
            act_refund_cut = 0.014
            act_halluc = f_h * acted_scale * 0.85
        elif d % 4 == 0:
            notes = "acted:audit_topk"
            action = "audit_topk"
            act_ticket_cut = 0.018
            act_refund_cut = 0.009
            act_halluc = f_h * acted_scale * 1.15
        else:
            notes = "acted:model_rollback"
            action = "model_rollback"
            act_ticket_cut = 0.022
            act_refund_cut = 0.011
            act_halluc = f_h * acted_scale
        n = 200
        for i in range(n):
            # HF-driven rates: quiet / fire_ignored / fire_acted
            halluc = int(rng.random() < (q_h * 0.6 if d < 14 else act_halluc))
            ticket = int(
                rng.random()
                < (
                    0.012
                    + 0.10 * halluc
                    + (ignore_ticket_bump if (fired and not acted) else 0.0)
                    - act_ticket_cut
                )
            )
            refund = int(
                rng.random()
                < (
                    0.004
                    + 0.07 * halluc
                    + (ignore_refund_bump if (fired and not acted) else 0.0)
                    - act_refund_cut
                )
            )
            serve_rows.append(
                (
                    f"h-{d}-{i}",
                    dt,
                    datetime.combine(dt, datetime.min.time()),
                    "shop_assistant",
                    f"s-{d}-{i%40}",
                    f"u-{i%80}",
                    "m_v2" if d >= 14 else "m_v1",
                    None,
                    "text",
                    "qa",
                    halluc,
                    None,
                    float(rng.normal(80, 20)),
                    float(rng.uniform(0.2, 0.8)),
                    float(rng.uniform(0, 0.2)),
                    "zh",
                    float(np.clip(rag + rng.normal(0, 0.05), 0, 1)),
                    int(rng.random() < 0.1),
                    int(rng.random() < 0.03),
                    float(rng.uniform(0, 30)),
                    refund,
                    ticket,
                    0,
                )
            )
            if halluc and fired and i < 15:
                audit_rows.append(
                    (
                        f"a-{d}-{i}",
                        dt,
                        "shop_assistant",
                        f"h-{d}-{i}",
                        f"b-{d}",
                        float(0.5 + rng.random()),
                        i + 1,
                        "confirmed_bad" if i < 10 else "false_alarm",
                        "auditor",
                        datetime.combine(dt, datetime.min.time()),
                    )
                )
        signal_rows.append(
            (
                f"sig-h-fire-{d}",
                dt,
                datetime.combine(dt, datetime.min.time()),
                "shop_assistant",
                f"b-{d}",
                "rfperm_fire",
                "concept",
                fired,
                float(1.1 + 0.4 * fired),
                n,
                "m_v2" if d >= 14 else "m_v1",
                notes,
            )
        )
        signal_rows.append(
            (
                f"sig-h-po-{d}",
                dt,
                datetime.combine(dt, datetime.min.time()),
                "shop_assistant",
                f"b-{d}",
                "po_risk0",
                "concept",
                0,
                float(0.2 + 0.5 * fired),
                n,
                "m_v2" if d >= 14 else "m_v1",
                None,
            )
        )

        # Style surface: portrait shift after day 10, usually WITHOUT concept fire
        style_auc = 0.52 if d < 10 else 0.82
        covar_fired = 1 if d >= 10 else 0
        concept_fired_style = 1 if d >= 22 else 0  # late joint days
        for i in range(150):
            cluster_idx = (i + (3 if d >= 10 else 0)) % 4
            creative_id = f"c{cluster_idx}"
            formal = 0.3 if d < 10 else 0.75
            brand = int(rng.random() < (0.005 + 0.02 * (d >= 10)))
            ctr_p = 0.08 if d < 10 else (0.055 if cluster_idx >= 2 else 0.07)
            serve_rows.append(
                (
                    f"s-{d}-{i}",
                    dt,
                    datetime.combine(dt, datetime.min.time()),
                    "feed_caption",
                    f"fs-{d}-{i%30}",
                    f"fu-{i%60}",
                    "cap_v1",
                    creative_id,
                    "text",
                    "caption",
                    None,
                    None,
                    float(rng.normal(40, 10)),
                    float(np.clip(formal + rng.normal(0, 0.05), 0, 1)),
                    float(rng.uniform(0, 0.3)),
                    "zh",
                    None,
                    int(rng.random() < ctr_p),
                    int(rng.random() < 0.01),
                    float(rng.uniform(0, 5)),
                    0,
                    0,
                    brand,
                )
            )
        signal_rows.append(
            (
                f"sig-st-auc-{d}",
                dt,
                datetime.combine(dt, datetime.min.time()),
                "feed_caption",
                f"sb-{d}",
                "style_domain_auc",
                "covariate",
                0,
                float(style_auc),
                150,
                "cap_v1",
                "acted" if covar_fired and d % 2 == 0 else None,
            )
        )
        signal_rows.append(
            (
                f"sig-st-cov-{d}",
                dt,
                datetime.combine(dt, datetime.min.time()),
                "feed_caption",
                f"sb-{d}",
                "rfperm_fire",
                "covariate",
                covar_fired,
                float(1.0 + 0.3 * covar_fired),
                150,
                "cap_v1",
                None,
            )
        )
        signal_rows.append(
            (
                f"sig-st-con-{d}",
                dt,
                datetime.combine(dt, datetime.min.time()),
                "feed_caption",
                f"sb-{d}",
                "rfperm_fire",
                "concept",
                concept_fired_style,
                float(1.0 + 0.5 * concept_fired_style),
                150,
                "cap_v1",
                None,
            )
        )
        if concept_fired_style:
            signal_rows.append(
                (
                    f"sig-st-judge-{d}",
                    dt,
                    datetime.combine(dt, datetime.min.time()),
                    "feed_caption",
                    f"sb-{d}",
                    "judge_err_ratio",
                    "concept",
                    0,
                    1.7,
                    150,
                    "cap_v1",
                    None,
                )
            )

    con.executemany(
        """
        INSERT INTO fct_serve_event VALUES
        (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        serve_rows,
    )
    con.executemany(
        """
        INSERT INTO fct_shift_signal VALUES
        (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        signal_rows,
    )
    if audit_rows:
        con.executemany(
            """
            INSERT INTO fct_audit_queue VALUES
            (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            audit_rows,
        )


def main() -> int:
    duckdb = _require_duckdb()
    OUT.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(OUT / "biz_value.duckdb"))

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
        sql = (SQL_DIR / name).read_text()
        con.execute(sql)

    # reset seed tables
    for t in [
        "fct_audit_queue",
        "fct_shift_signal",
        "fct_serve_event",
        "dim_value_assumption",
        "dim_creative",
        "dim_surface",
    ]:
        con.execute(f"DELETE FROM {t}")

    hop = load_hf_hop()
    (OUT / "hf_hop_knobs.json").write_text(
        json.dumps(hop, ensure_ascii=False, indent=2)
    )
    seed(con, hop=hop)

    halluc = con.execute("SELECT * FROM vw_halluc_roi_rollup").fetchdf()
    style = con.execute("SELECT * FROM vw_style_roi_rollup").fetchdf()
    dash = con.execute("SELECT * FROM vw_biz_value_dashboard").fetchdf()
    gate = con.execute(
        "SELECT * FROM vw_alignment_merge_gate ORDER BY dt DESC LIMIT 10"
    ).fetchdf()
    oncall = con.execute(
        "SELECT * FROM vw_oncall_queue ORDER BY severity_proxy DESC LIMIT 15"
    ).fetchdf()

    def dump(df, path: Path):
        path.write_text(
            df.to_json(orient="records", force_ascii=False, indent=2, date_format="iso")
        )

    cs = con.execute("SELECT * FROM vw_cs_assist_exec_summary").fetchdf()
    cs_arms = con.execute(
        "SELECT * FROM vw_cs_assist_arm_stats ORDER BY surface_id, arm"
    ).fetchdf()
    cs_rates = con.execute("SELECT * FROM vw_cs_assist_rate_compare").fetchdf()
    cs_actions = con.execute(
        "SELECT * FROM vw_cs_assist_action_increment ORDER BY incremental_yen DESC"
    ).fetchdf()
    cs_action_net = con.execute(
        "SELECT * FROM vw_cs_assist_action_net ORDER BY net_yen_after_action_cost DESC"
    ).fetchdf()
    cs_traffic = con.execute(
        "SELECT * FROM vw_cs_assist_traffic_scenarios ORDER BY daily_sessions"
    ).fetchdf()
    cs_net = con.execute("SELECT * FROM vw_cs_assist_net_increment").fetchdf()
    cs_week = con.execute(
        "SELECT * FROM vw_cs_assist_weekly ORDER BY week_start"
    ).fetchdf()
    cs_wow = con.execute(
        "SELECT * FROM vw_cs_assist_weekly_wow ORDER BY week_start"
    ).fetchdf()
    cs_ledger = con.execute(
        """
        SELECT arm,
               COUNT(*) AS days,
               SUM(n_sessions) AS sessions,
               SUM(n_tickets) AS tickets,
               SUM(n_refunds) AS refunds,
               SUM(n_contained) AS contained,
               ROUND(AVG(ticket_rate), 4) AS avg_ticket_rate,
               ROUND(AVG(refund_rate), 4) AS avg_refund_rate,
               ROUND(AVG(containment_rate), 4) AS avg_containment_rate,
               ROUND(SUM(quality_cost_yen), 0) AS quality_cost_yen,
               ROUND(SUM(net_ops_value_yen), 0) AS net_ops_value_yen
        FROM vw_cs_assist_daily_ledger
        GROUP BY arm
        ORDER BY arm
        """
    ).fetchdf()

    dump(halluc, OUT / "roi_halluc.json")
    dump(style, OUT / "roi_style.json")
    dump(dash, OUT / "dashboard.json")
    dump(gate, OUT / "merge_gate_sample.json")
    dump(oncall, OUT / "oncall_sample.json")
    dump(cs, OUT / "cs_assist_increment.json")
    dump(cs_arms, OUT / "cs_assist_arms.json")
    dump(cs_rates, OUT / "cs_assist_rates.json")
    dump(cs_actions, OUT / "cs_assist_actions.json")
    dump(cs_action_net, OUT / "cs_assist_action_net.json")
    dump(cs_traffic, OUT / "cs_assist_traffic.json")
    dump(cs_net, OUT / "cs_assist_net.json")
    dump(cs_week, OUT / "cs_assist_weekly.json")
    dump(cs_wow, OUT / "cs_assist_weekly_wow.json")
    cs_exec = con.execute("SELECT * FROM vw_cs_assist_exec_dashboard").fetchdf()
    cs_rec = con.execute("SELECT * FROM vw_cs_assist_action_recommend ORDER BY recommend_rank").fetchdf()
    dump(cs_exec, OUT / "cs_assist_exec_dashboard.json")
    dump(cs_rec, OUT / "cs_assist_action_recommend.json")
    cs_week_attr = con.execute(
        "SELECT * FROM vw_cs_assist_week_attribution ORDER BY week_start"
    ).fetchdf()
    cs_payback = con.execute(
        "SELECT * FROM vw_cs_assist_action_payback ORDER BY payback_days"
    ).fetchdf()
    cs_attr_check = con.execute(
        "SELECT * FROM vw_cs_assist_attribution_check"
    ).fetchdf()
    dump(cs_week_attr, OUT / "cs_assist_week_attribution.json")
    dump(cs_payback, OUT / "cs_assist_action_payback.json")
    dump(cs_attr_check, OUT / "cs_assist_attribution_check.json")
    cs_curve = con.execute(
        "SELECT * FROM vw_cs_assist_cumulative_curve ORDER BY dt"
    ).fetchdf()
    cs_unit = con.execute(
        "SELECT * FROM vw_cs_assist_unit_econ_sensitivity ORDER BY scenario"
    ).fetchdf()
    cs_band = con.execute("SELECT * FROM vw_cs_assist_unit_econ_band").fetchdf()
    dump(cs_curve, OUT / "cs_assist_cumulative_curve.json")
    dump(cs_unit, OUT / "cs_assist_unit_econ_sensitivity.json")
    dump(cs_band, OUT / "cs_assist_unit_econ_band.json")
    # Finance CSV: day-level contribution for ledger import
    cs_curve.to_csv(OUT / "cs_assist_finance_daily.csv", index=False)
    dump(cs_ledger, OUT / "cs_assist_ledger.json")

    h = halluc.iloc[0].to_dict() if len(halluc) else {}
    s = style.iloc[0].to_dict() if len(style) else {}
    c = cs.iloc[0].to_dict() if len(cs) else {}
    r = cs_rates.iloc[0].to_dict() if len(cs_rates) else {}
    n = cs_net.iloc[0].to_dict() if len(cs_net) else {}

    def pct(x):
        try:
            return f"{100.0 * float(x):.2f}%"
        except (TypeError, ValueError):
            return str(x)

    action_rows = []
    for _, row in cs_actions.iterrows():
        action_rows.append(
            f"| {row['action_type']} | {int(row.get('n_days', row.get('days', 0)))} | {int(row['sessions'])} | "
            f"{row['tickets_avoided']} | {row['refunds_avoided']} | "
            f"{row['extra_contained']} | **¥{int(row['incremental_yen'])}** | "
            f"{row['avg_rag_hit']:.2f} |"
        )
    action_table = "\n".join(action_rows) if action_rows else "| (none) ||||||"

    action_net_rows = []
    for _, row in cs_action_net.iterrows():
        action_net_rows.append(
            f"| {row['action_type']} | {int(row['n_days'])} | "
            f"¥{int(row['gross_yen'])} | ¥{int(row['action_day_cost_yen'])} | "
            f"**¥{int(row['net_yen_after_action_cost'])}** |"
        )
    action_net_table = "\n".join(action_net_rows) if action_net_rows else "| (none) ||||"

    week_attr_rows = []
    for _, row in cs_week_attr.iterrows():
        week_attr_rows.append(
            f"| {str(row['week_start'])[:10]} | {int(row['days_acted'])} | "
            f"{int(row['sessions_acted'])} | {row['tickets_avoided']} | "
            f"{row['refunds_avoided']} | {row['extra_contained']} | "
            f"**¥{int(row['gross_yen'])}** | {row['pct_of_total_gross']}% |"
        )
    week_attr_table = "\n".join(week_attr_rows) if week_attr_rows else "| (none) |||||||"

    payback_rows = []
    for _, row in cs_payback.iterrows():
        payback_rows.append(
            f"| {row['action_type']} | {int(row['n_days'])} | "
            f"¥{int(row['gross_yen_per_day'])} | ¥{int(row['cost_yen_per_day'])} | "
            f"**{row['payback_days']} 天** | {row['payback_bucket']} | "
            f"{row['net_roi_multiple']}x |"
        )
    payback_table = "\n".join(payback_rows) if payback_rows else "| (none) ||||||"

    attr_gap = float(cs_attr_check.iloc[0]["attribution_gap_yen"]) if len(cs_attr_check) else 0
    attr_gap_pct = float(cs_attr_check.iloc[0]["attribution_gap_pct"]) if len(cs_attr_check) else 0

    unit_rows = []
    for _, row in cs_unit.iterrows():
        unit_rows.append(
            f"| {row['scenario']} | ¥{row['ticket_cost']} | ¥{row['refund_cost']} | "
            f"¥{row['contain_value']} | **¥{int(row['gross_yen'])}** | "
            f"**¥{int(row['net_yen'])}** |"
        )
    unit_table = "\n".join(unit_rows) if unit_rows else "| (none) |||||"
    band = cs_band.iloc[0].to_dict() if len(cs_band) else {}

    curve_tail = cs_curve.tail(3) if len(cs_curve) else cs_curve
    curve_rows = []
    for _, row in curve_tail.iterrows():
        curve_rows.append(
            f"| {str(row['dt'])[:10]} | ¥{int(row['gross_yen_day'])} | "
            f"**¥{int(row['cumulative_gross_yen'])}** | {row['cumulative_pct_of_total']}% |"
        )
    curve_table = "\n".join(curve_rows) if curve_rows else "| (none) |||"

    traffic_rows = []
    for _, row in cs_traffic.iterrows():
        traffic_rows.append(
            f"| {row['scenario']} | {int(row['daily_sessions']):,} | "
            f"**¥{int(row['monthly_yen']):,}** | {row['monthly_tickets']} | "
            f"{row['monthly_refunds']} | {row['monthly_extra_contained']} |"
        )
    traffic_table = "\n".join(traffic_rows)

    week_rows = []
    for _, row in cs_week.iterrows():
        week_rows.append(
            f"| {row['week_start']} | {int(row['sessions'])} | {int(row['tickets'])} | "
            f"{int(row['refunds'])} | {pct(row['ticket_rate'])} | "
            f"{pct(row['containment_rate'])} | ¥{int(row['net_contrib_yen'])} | "
            f"{int(row['days_acted'])}/{int(row['days_ignored'])} |"
        )
    week_table = "\n".join(week_rows)

    # Pull HF hop evidence if present (justify which action bucket)
    hf_note = ""
    hf_path = ROOT / "results" / "agod" / "hf_landing" / "halu_regime_rag.json"
    if hf_path.exists():
        hf = json.loads(hf_path.read_text())
        hop_cut = hf.get("hop_at_cut") or hf.get("hop_at_cut") or {}
        ranking = hop_cut.get("ranking") or {}
        hf_note = f"""
## 与 HF 幻觉子集的衔接（证据链）

HaluEval 子集原型（`results/agod/hf_landing/halu_regime_rag.json`）：

- cut 处 `fired={hop_cut.get('fired')}`，ratio≈`{hop_cut.get('ratio')}`
- 跳变后幻觉率 ≈ `{ranking.get('halluc_rate')}`；`po_risk0` P@10=`{ranking.get('precision_at_10')}`
- Top-10 平均 `rag_hit`=`{ranking.get('mean_rag_hit_top10')}` → 路由：检索缺口 vs 生成制度

本账动作拆分与之对齐：`avg_rag_hit` 低 → `retrieval_refresh`；否则 → `model_rollback` / `audit_topk`。
"""

    hop_source = hop.get("source")
    hop_quiet = hop.get("quiet_halluc")
    hop_fire = hop.get("fire_halluc")
    hop_ratio = hop.get("hop_ratio")
    hop_acted_scale = hop.get("acted_halluc_scale")
    hop_rag_low = hop.get("rag_low")
    hop_rag_ok = hop.get("rag_ok")
    hop_rag_thr = hop.get("rag_threshold")
    cs_report = f"""# 客服助手增量贡献账（可复现）

产品面：`shop_assistant`（客服助手）。  
对照：同一幻觉制度跳变期内，**动作落地** vs **同条件不动作**。  
单位经济：工单 ¥25 / 退款 ¥80 / 机器人成功承接会话 ¥3.5 / 审计件 ¥6。

## 落地产出（业务 KPI）

| 指标 | 数值 |
|------|------|
| 少产生工单 | **{c.get('tickets_avoided')}** 单 |
| 少产生退款 | **{c.get('refunds_avoided')}** 单 |
| 多机器人承接会话 | **{c.get('extra_sessions_contained')}** 次 |
| 区间增量贡献（毛） | **¥{c.get('incremental_yen')}** |
| 其中：少工单 | ¥{c.get('yen_from_tickets')} |
| 其中：少退款 | ¥{c.get('yen_from_refunds')} |
| 其中：多承接 | ¥{c.get('yen_from_containment')} |
| **扣审计后净增量** | **¥{n.get('net_incremental_yen')}** |
| 审计件数 / 成本（acted） | {n.get('audits_acted')} / ¥{n.get('audit_cost_acted')} |
| 审计 precision 代理 | {n.get('avg_audit_precision')} |
| 每千会话毛增量 | **¥{c.get('incremental_yen_per_1k_sessions')}** |
| 每千会话净增量 | **¥{n.get('net_yen_per_1k_sessions')}** |
| 30 天跑率（毛，acted 日均外推） | **¥{c.get('monthly_runrate_yen')}** / 月 |
| 承接率提升 | {c.get('containment_rate_lift_pp')} pp |
| 对照天数 acted / ignored | {c.get('days_acted')} / {c.get('days_ignored')} |
| acted 会话量 | {c.get('sessions_acted')} |

## Before → After 费率（三臂）

| 费率 | quiet（平稳） | fire+不动作 | fire+动作落地 |
|------|---------------|-------------|---------------|
| 工单率 | {pct(r.get('ticket_rate_quiet'))} | {pct(r.get('ticket_rate_ignored'))} | {pct(r.get('ticket_rate_acted'))} |
| 退款率 | {pct(r.get('refund_rate_quiet'))} | {pct(r.get('refund_rate_ignored'))} | {pct(r.get('refund_rate_acted'))} |
| 机器人承接率 | {pct(r.get('contain_rate_quiet'))} | {pct(r.get('contain_rate_ignored'))} | {pct(r.get('contain_rate_acted'))} |

读法：制度跳变后若不动作，工单/退款率显著恶化；动作落地后费率回到接近平稳期，这就是增量贡献的来源。

## 按动作类型拆贡献（增量来源）

| 动作 | 天数 | 会话 | 少工单 | 少退款 | 多承接 | 增量¥ | avg_rag_hit |
|------|------|------|--------|--------|--------|-------|-------------|
{action_table}

## 扣动作成本后净贡献

| 动作 | 天数 | 毛增量¥ | 动作日成本¥ | 净贡献¥ |
|------|------|---------|-------------|---------|
{action_net_table}

## 流量情景（月贡献外推 · 毛）

| 情景 | 日会话 | 月增量¥ | 月少工单 | 月少退款 | 月多承接 |
|------|--------|---------|----------|----------|----------|
{traffic_table}

## 周经营看板

| 周起始 | 会话 | 工单 | 退款 | 工单率 | 承接率 | 净贡献¥ | acted/ignored 天 |
|--------|------|------|------|--------|--------|---------|------------------|
{week_table}
{hf_note}
## 口径

```
毛增量¥ =
  (ignored工单率 - acted工单率) × acted会话数 × 工单单价
+ (ignored退款率 - acted退款率) × acted会话数 × 退款单价
+ (acted承接率 - ignored承接率) × acted会话数 × 承接会话价值

净增量¥ = 毛增量¥ − acted 侧审计人力成本
         （审计单价 × acted 日审计件数；不把 ignored 多烧的审计算进贡献）
```

## HF hop 驱动参数

- source: `{hop_source}`
- quiet_halluc / fire_halluc(biz): `{hop_quiet}` / `{hop_fire}`
- hop_ratio: `{hop_ratio}`；acted_halluc_scale: `{hop_acted_scale}`
- rag_low / rag_ok / thr: `{hop_rag_low}` / `{hop_rag_ok}` / `{hop_rag_thr}`

## 周归因贡献（财务对账）

| 周起始 | acted天 | 会话 | 少工单 | 少退款 | 多承接 | 毛¥ | 占总毛% |
|--------|---------|------|--------|--------|--------|-----|---------|
{week_attr_table}

周归因合计 vs 总账缺口：¥{attr_gap}（{attr_gap_pct}%）——应为 ~0。

## 动作回本天数

| 动作 | 天数 | 日均毛¥ | 日均成本¥ | 回本 | 分档 | 净ROI |
|------|------|---------|-----------|------|------|-------|
{payback_table}

读法：回本天数 < 1 = 当天回本；净ROI = 动作净贡献 / 动作日成本合计。

## 单位经济敏感度（量固定，扫单价）

| 情景 | 工单单价 | 退款单价 | 承接单价 | 毛¥ | 净¥ |
|------|----------|----------|----------|-----|-----|
{unit_table}

单价全 ±20% 净贡献带：**¥{int(band.get('net_yen_low') or 0):,} → ¥{int(band.get('net_yen_base') or 0):,} → ¥{int(band.get('net_yen_high') or 0):,}**（带宽相对 base ≈ {band.get('net_band_width_vs_base')}）。

## 累计贡献曲线（末三日）

| 日期 | 当日毛¥ | 累计毛¥ | 累计占比 |
|------|---------|---------|----------|
{curve_table}

财务日账 CSV：`results/agod/biz_value_sql/cs_assist_finance_daily.csv`

## 一句对外

客服助手在幻觉制度跳变的 {c.get('days_acted')} 个动作日里，相对同条件不动作：少了
{c.get('tickets_avoided')} 单工单、{c.get('refunds_avoided')} 单退款，
多承接 {c.get('extra_sessions_contained')} 次会话；毛贡献约 ¥{c.get('incremental_yen')}，
扣审计人力后净贡献约 ¥{n.get('net_incremental_yen')}。
生产中等流量（日 1 万会话）见上表 `prod_mid`。
"""
    (OUT / "CS_ASSISTANT_CONTRIBUTION.md").write_text(cs_report)
    docs_biz = ROOT / "docs" / "biz" / "CS_ASSISTANT_CONTRIBUTION.md"
    docs_biz.write_text(cs_report)

    # Weekly ops brief (paste-ready for biz review)
    try:
        import importlib.util

        brief_path = ROOT / "scripts" / "agod" / "cs_assist_weekly_ops_brief.py"
        spec = importlib.util.spec_from_file_location("cs_ops_brief", brief_path)
        brief_mod = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(brief_mod)
        brief = brief_mod.build_brief()
        (OUT / "CS_ASSISTANT_WEEKLY_OPS_BRIEF.md").write_text(brief)
        (ROOT / "docs" / "biz" / "CS_ASSISTANT_WEEKLY_OPS_BRIEF.md").write_text(brief)
        # Append pointer into contribution doc
        pointer = (
            "\n\n---\n\n周经营简报（WoW / 动作净贡献）："
            "`docs/biz/CS_ASSISTANT_WEEKLY_OPS_BRIEF.md`\n"
            "经营总看板 / hop 情景：`docs/biz/CS_ASSISTANT_EXEC_DASHBOARD.md`\n"
            "周归因 / 回本：见贡献账内「周归因贡献」「动作回本天数」\n"
            "单位经济带 / 累计曲线 / 财务CSV：见贡献账对应章节\n"
            "方法论 × 业务情景 × 单位经济对照原型："
            "`docs/biz/METHOD_BIZ_SCENARIO_PROTOTYPE.md`\n"
            "方法论落地 Roadmap（推理/对齐 + 回本/周归因详解）："
            "`docs/biz/METHOD_LANDING_ROADMAP.md`\n"
            "落地场景 · 多承接 · 方法异同 · 迭代更新："
            "`docs/biz/LANDING_CONTAIN_METHOD_ITER.md`\n"
        )
        docs_biz.write_text(docs_biz.read_text() + pointer)
        (OUT / "CS_ASSISTANT_CONTRIBUTION.md").write_text(
            (OUT / "CS_ASSISTANT_CONTRIBUTION.md").read_text() + pointer
        )
    except Exception as e:  # noqa: BLE001 — brief is additive; don't fail demo
        print(f"weekly ops brief skipped: {e}")

    report = f"""# Business-value SQL demo report

## 客服助手增量贡献（主结果）

- tickets avoided: `{c.get('tickets_avoided')}`
- refunds avoided: `{c.get('refunds_avoided')}`
- extra sessions contained: `{c.get('extra_sessions_contained')}`
- **gross incremental ¥: `{c.get('incremental_yen')}`**
- **net incremental ¥ (after audit): `{n.get('net_incremental_yen')}`**
- ¥ breakdown tickets/refunds/contain: `{c.get('yen_from_tickets')}` / `{c.get('yen_from_refunds')}` / `{c.get('yen_from_containment')}`
- monthly run-rate ¥ (gross): `{c.get('monthly_runrate_yen')}`
- action split: `cs_assist_actions.json`
- traffic / weekly / net: `cs_assist_traffic.json` / `cs_assist_weekly.json` / `cs_assist_net.json`
- detail: `CS_ASSISTANT_CONTRIBUTION.md` / `docs/biz/CS_ASSISTANT_CONTRIBUTION.md`

## Hallucination cost rollup (`shop_assistant`)

- est daily save act vs ignore: `{h.get('est_daily_save_act_vs_ignore')}` CNY

## Style / creative (`feed_caption`)

- style-only days: `{s.get('n_style_only_days')}`
- avg CTR style-only / quiet: `{s.get('avg_ctr_style_only')}` / `{s.get('avg_ctr_quiet')}`

## Artifacts

- `sql/biz_value/04_cs_assistant_contribution.sql`
- `sql/biz_value/05_cs_net_and_weekly.sql`
- `results/agod/biz_value_sql/CS_ASSISTANT_CONTRIBUTION.md`
"""
    (OUT / "REPORT.md").write_text(report)
    print(cs_report)
    print(report)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
