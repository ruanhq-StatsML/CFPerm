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
          (?, 'feed_caption', 'value_per_ctr_point', 1200.0, 'CNY', 'growth proxy'),
          (?, 'feed_caption', 'brand_incident_cost', 50.0, 'CNY', 'brand ops'),
          (?, 'ad_creative', 'value_per_ctr_point', 2500.0, 'CNY', 'ads proxy'),
          (?, 'ad_creative', 'brand_incident_cost', 120.0, 'CNY', 'brand ops')
        """,
        [start] * 6,
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
            halluc = int(rng.random() < (0.05 if d < 14 else 0.35))
            ticket = int(rng.random() < (0.01 + 0.08 * halluc + 0.04 * fired * (1 - acted)))
            refund = int(rng.random() < (0.005 + 0.06 * halluc + 0.03 * fired * (1 - acted)))
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
        path.write_text(df.to_json(orient="records", force_ascii=False, indent=2))

    dump(halluc, OUT / "roi_halluc.json")
    dump(style, OUT / "roi_style.json")
    dump(dash, OUT / "dashboard.json")
    dump(gate, OUT / "merge_gate_sample.json")
    dump(oncall, OUT / "oncall_sample.json")

    h = halluc.iloc[0].to_dict() if len(halluc) else {}
    s = style.iloc[0].to_dict() if len(style) else {}
    report = f"""# Business-value SQL demo report

Seeded DuckDB run for hallucination + style/creative surfaces.

## Hallucination ROI (`shop_assistant`)

- fire acted days: `{h.get('n_fire_acted_days')}`
- fire ignored days: `{h.get('n_fire_ignored_days')}`
- avg quality cost fire+acted: `{h.get('avg_cost_fire_acted')}`
- avg quality cost fire+ignored: `{h.get('avg_cost_fire_ignored')}`
- **est daily save act vs ignore: `{h.get('est_daily_save_act_vs_ignore')}` CNY**

Narrative: acting on concept-fire (√PO / rollback / Top-k audit) vs ignoring
reduces CS+refund cost on the same fire regime — this is the value sentence.

## Style / creative ROI (`feed_caption`)

- style-only days: `{s.get('n_style_only_days')}`
- concept-only days: `{s.get('n_concept_only_days')}`
- joint days: `{s.get('n_joint_days')}`
- avg CTR style-only / quiet: `{s.get('avg_ctr_style_only')}` / `{s.get('avg_ctr_quiet')}`
- avg brand cost style-only / quiet: `{s.get('avg_brand_cost_style_only')}` / `{s.get('avg_brand_cost_quiet')}`

Narrative: high `style_domain_auc` without concept-fire → open **creative mix**
ticket only; do not retrain preference head.

## Artifacts

- DuckDB: `results/agod/biz_value_sql/biz_value.duckdb`
- JSON: `roi_halluc.json`, `roi_style.json`, `dashboard.json`, `merge_gate_sample.json`
- SQL package: `sql/biz_value/`
- Mapping doc: `docs/biz/HALLUC_STYLE_BIZ_SQL_MAP.md`
"""
    (OUT / "REPORT.md").write_text(report)
    print(report)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
