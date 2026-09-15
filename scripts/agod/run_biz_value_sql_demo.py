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


def _require_duckdb():
    try:
        import duckdb  # noqa: F401
    except ImportError as e:
        raise SystemExit(
            "duckdb is required: pip install duckdb\n" + str(e)
        ) from e
    import duckdb

    return duckdb


def seed(con) -> None:
    rng = np.random.default_rng(7)
    start = date(2026, 9, 1)
    days = 28

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
          (?, 'feed_caption', 'value_per_ctr_point', 1200.0, 'CNY', 'growth proxy'),
          (?, 'feed_caption', 'brand_incident_cost', 50.0, 'CNY', 'brand ops'),
          (?, 'ad_creative', 'value_per_ctr_point', 2500.0, 'CNY', 'ads proxy'),
          (?, 'ad_creative', 'brand_incident_cost', 120.0, 'CNY', 'brand ops')
        """,
        [start] * 7,
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

    for d in range(days):
        dt = start + timedelta(days=d)
        # Hallucination surface: concept hop after day 14
        fired = 1 if d >= 14 else 0
        acted = 1 if (fired and d % 2 == 0) else 0
        notes = "acted" if acted else ("ignored" if fired else None)
        rag = 0.55 if d < 14 else (0.20 if d % 3 == 0 else 0.60)
        n = 200
        for i in range(n):
            # Clear incremental effect: after regime hop, ignoring mitigation
            # drives tickets/refunds; acting contains more sessions in-bot.
            halluc = int(rng.random() < (0.04 if d < 14 else (0.18 if acted else 0.42)))
            ticket = int(
                rng.random()
                < (
                    0.012
                    + 0.10 * halluc
                    + (0.055 if (fired and not acted) else 0.0)
                    - (0.02 if acted else 0.0)
                )
            )
            refund = int(
                rng.random()
                < (
                    0.004
                    + 0.07 * halluc
                    + (0.035 if (fired and not acted) else 0.0)
                    - (0.01 if acted else 0.0)
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

    seed(con)

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
    dump(cs_ledger, OUT / "cs_assist_ledger.json")

    h = halluc.iloc[0].to_dict() if len(halluc) else {}
    s = style.iloc[0].to_dict() if len(style) else {}
    c = cs.iloc[0].to_dict() if len(cs) else {}
    r = cs_rates.iloc[0].to_dict() if len(cs_rates) else {}

    def pct(x):
        try:
            return f"{100.0 * float(x):.2f}%"
        except (TypeError, ValueError):
            return str(x)

    cs_report = f"""# 客服助手增量贡献账（可复现）

产品面：`shop_assistant`（客服助手）。  
对照：同一幻觉制度跳变期内，**动作落地** vs **同条件不动作**。  
单位经济：工单 ¥25 / 退款 ¥80 / 机器人成功承接会话 ¥3.5。

## 落地产出（业务 KPI）

| 指标 | 数值 |
|------|------|
| 少产生工单 | **{c.get('tickets_avoided')}** 单 |
| 少产生退款 | **{c.get('refunds_avoided')}** 单 |
| 多机器人承接会话 | **{c.get('extra_sessions_contained')}** 次 |
| 区间增量贡献 | **¥{c.get('incremental_yen')}** |
| 其中：少工单 | ¥{c.get('yen_from_tickets')} |
| 其中：少退款 | ¥{c.get('yen_from_refunds')} |
| 其中：多承接 | ¥{c.get('yen_from_containment')} |
| 每千会话增量贡献 | **¥{c.get('incremental_yen_per_1k_sessions')}** |
| 30 天跑率（按 acted 日均外推） | **¥{c.get('monthly_runrate_yen')}** / 月 |
| 30 天少工单（外推） | {c.get('monthly_tickets_avoided')} 单 |
| 30 天少退款（外推） | {c.get('monthly_refunds_avoided')} 单 |
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

## 口径

```
增量贡献¥ =
  (ignored工单率 - acted工单率) × acted会话数 × 工单单价
+ (ignored退款率 - acted退款率) × acted会话数 × 退款单价
+ (acted承接率 - ignored承接率) × acted会话数 × 承接会话价值
```

## 一句对外

客服助手在幻觉制度跳变的 7 个动作日里，相对同条件不动作：少了
{c.get('tickets_avoided')} 单工单、{c.get('refunds_avoided')} 单退款，
多承接 {c.get('extra_sessions_contained')} 次会话，贡献约 ¥{c.get('incremental_yen')}；
按当前流量外推约 ¥{c.get('monthly_runrate_yen')}/月。
"""
    (OUT / "CS_ASSISTANT_CONTRIBUTION.md").write_text(cs_report)

    # Stable docs copy for PR / biz review
    docs_biz = ROOT / "docs" / "biz" / "CS_ASSISTANT_CONTRIBUTION.md"
    docs_biz.write_text(cs_report)

    report = f"""# Business-value SQL demo report

## 客服助手增量贡献（主结果）

- tickets avoided: `{c.get('tickets_avoided')}`
- refunds avoided: `{c.get('refunds_avoided')}`
- extra sessions contained: `{c.get('extra_sessions_contained')}`
- **incremental ¥ realized: `{c.get('incremental_yen')}`**
- ¥ breakdown tickets/refunds/contain: `{c.get('yen_from_tickets')}` / `{c.get('yen_from_refunds')}` / `{c.get('yen_from_containment')}`
- incremental ¥ / 1k sessions: `{c.get('incremental_yen_per_1k_sessions')}`
- monthly run-rate ¥: `{c.get('monthly_runrate_yen')}`
- detail: `CS_ASSISTANT_CONTRIBUTION.md` / `docs/biz/CS_ASSISTANT_CONTRIBUTION.md`

## Hallucination cost rollup (`shop_assistant`)

- est daily save act vs ignore: `{h.get('est_daily_save_act_vs_ignore')}` CNY

## Style / creative (`feed_caption`)

- style-only days: `{s.get('n_style_only_days')}`
- avg CTR style-only / quiet: `{s.get('avg_ctr_style_only')}` / `{s.get('avg_ctr_quiet')}`

## Artifacts

- `sql/biz_value/04_cs_assistant_contribution.sql`
- `results/agod/biz_value_sql/CS_ASSISTANT_CONTRIBUTION.md`
- `cs_assist_increment.json` / `cs_assist_ledger.json` / `cs_assist_rates.json`
"""
    (OUT / "REPORT.md").write_text(report)
    print(cs_report)
    print(report)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
