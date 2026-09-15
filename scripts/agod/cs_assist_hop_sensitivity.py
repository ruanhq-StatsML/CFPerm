#!/usr/bin/env python3
"""Hop-intensity sensitivity for 客服助手 contribution.

Re-seeds in-memory with weak / base / strong HF-hop knobs and reports
before→after ¥ so biz sees a contribution band, not a single point.

Usage:
  python3 scripts/agod/cs_assist_hop_sensitivity.py
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SQL = ROOT / "sql" / "biz_value"
OUT = ROOT / "results" / "agod" / "biz_value_sql"
DOCS = ROOT / "docs" / "biz"
SQL_FILES = [
    "00_schema.sql",
    "01_hallucination_value.sql",
    "02_style_drift_value.sql",
    "03_value_dashboard.sql",
    "04_cs_assistant_contribution.sql",
    "05_cs_net_and_weekly.sql",
    "06_cs_exec_dashboard.sql",
    "07_cs_week_attribution_payback.sql",
    "08_cs_unit_econ_cumulative.sql",
]


def _load_demo():
    spec = importlib.util.spec_from_file_location(
        "biz_demo", ROOT / "scripts" / "agod" / "run_biz_value_sql_demo.py"
    )
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def _run_scenario(mod, duckdb, hop: dict) -> dict:
    con = duckdb.connect(":memory:")
    for name in SQL_FILES:
        con.execute((SQL / name).read_text())
    mod.seed(con, hop=hop)
    row = con.execute("SELECT * FROM vw_cs_assist_exec_dashboard").fetchone()
    cols = [d[0] for d in con.description]
    d = dict(zip(cols, row))
    d["scenario"] = hop["scenario"]
    d["hop_ratio"] = hop["hop_ratio"]
    d["fire_halluc"] = hop["fire_halluc"]
    d["quiet_halluc"] = hop["quiet_halluc"]
    con.close()
    return d


def main() -> int:
    duckdb = __import__("duckdb")
    mod = _load_demo()
    base = mod.load_hf_hop()

    scenarios = [
        {
            **base,
            "scenario": "weak_hop",
            "fire_halluc": float(
                np.clip(base["fire_halluc"] * 0.55, base["quiet_halluc"] + 0.05, 0.35)
            ),
            "hop_ratio": max(2.0, base["hop_ratio"] * 0.35),
            "acted_halluc_scale": float(
                np.clip(base["acted_halluc_scale"] * 1.2, 0.3, 0.6)
            ),
        },
        {**base, "scenario": "base_hop"},
        {
            **base,
            "scenario": "strong_hop",
            "fire_halluc": float(np.clip(base["fire_halluc"] * 1.25, 0.35, 0.55)),
            "hop_ratio": base["hop_ratio"] * 1.35,
            "acted_halluc_scale": float(
                np.clip(base["acted_halluc_scale"] * 0.85, 0.2, 0.45)
            ),
        },
    ]

    rows = [_run_scenario(mod, duckdb, h) for h in scenarios]
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "cs_assist_hop_sensitivity.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2, default=str)
    )

    def yen(x):
        try:
            return f"¥{int(float(x)):,}"
        except (TypeError, ValueError):
            return str(x)

    table = "\n".join(
        f"| {r['scenario']} | {r['fire_halluc']:.2f} | {r['hop_ratio']:.1f} | "
        f"{r.get('tickets_avoided')} | {r.get('refunds_avoided')} | "
        f"{r.get('extra_sessions_contained')} | **{yen(r.get('gross_yen'))}** | "
        f"**{yen(r.get('net_yen'))}** | {yen(r.get('prod_mid_monthly_net_yen'))} |"
        for r in rows
    )
    base_row = next(r for r in rows if r["scenario"] == "base_hop")
    weak = next(r for r in rows if r["scenario"] == "weak_hop")
    strong = next(r for r in rows if r["scenario"] == "strong_hop")

    md = f"""# 客服助手经营总看板（大迭代）

## 对外一句

幻觉制度跳变期内，客服助手相对不动作：少 **{base_row.get('tickets_avoided')}** 单工单、
少 **{base_row.get('refunds_avoided')}** 单退款、多承接 **{base_row.get('extra_sessions_contained')}** 次；
毛贡献 **{yen(base_row.get('gross_yen'))}**，扣审计后净贡献 **{yen(base_row.get('net_yen'))}**；
优先动作 **`{base_row.get('top_action')}`**（日均净 ~{yen(base_row.get('top_action_net_yen_per_day'))}）。
生产中等流量（日 1 万会话）月净约 **{yen(base_row.get('prod_mid_monthly_net_yen'))}**。

## Before → After（base hop）

| 费率 | quiet | fire+不动作 | fire+动作 | 节省/提升 |
|------|-------|-------------|-----------|-----------|
| 工单率 | {base_row.get('ticket_rate_quiet_pct')}% | {base_row.get('ticket_rate_ignored_pct')}% | {base_row.get('ticket_rate_acted_pct')}% | **{base_row.get('ticket_rate_saved_pp')} pp** |
| 退款率节省 | — | — | — | **{base_row.get('refund_rate_saved_pp')} pp** |
| 承接率提升 | — | — | — | **{base_row.get('contain_rate_lift_pp')} pp** |

## HF hop 强度情景（贡献带）

| 情景 | fire_halluc | hop_ratio | 少工单 | 少退款 | 多承接 | 毛¥ | 净¥ | 日1万会话月净 |
|------|-------------|-----------|--------|--------|--------|-----|-----|----------------|
{table}

读法：hop 越强，不动作代价越大，动作落地的增量贡献越高。
弱→强净贡献带：**{yen(weak.get('net_yen'))} → {yen(base_row.get('net_yen'))} → {yen(strong.get('net_yen'))}**。

## 优先动作

- `{base_row.get('top_action')}`：{base_row.get('top_action_playbook')}
- 日均净贡献约 {yen(base_row.get('top_action_net_yen_per_day'))}

## 口径

```
毛¥ = 费率差 × acted会话 × 单价
净¥ = 毛¥ − acted 审计人力
动作净¥/日 = (动作毛¥ − 动作日成本 − 审计分摊) / 天数
```

产物：`vw_cs_assist_exec_dashboard` / `cs_assist_hop_sensitivity.json`
"""
    (OUT / "CS_ASSISTANT_EXEC_DASHBOARD.md").write_text(md)
    (DOCS / "CS_ASSISTANT_EXEC_DASHBOARD.md").write_text(md)
    print(md)
    print(f"wrote {OUT / 'CS_ASSISTANT_EXEC_DASHBOARD.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
