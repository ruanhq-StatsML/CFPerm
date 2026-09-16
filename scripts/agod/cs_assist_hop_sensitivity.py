#!/usr/bin/env python3
"""HF hop / rag knobs → re-seed → 可对账贡献带（毛 / 净 / 全成本净 / 动作混比）.

Business question: if the HF hop is weaker/stronger, or rag_hit shifts,
how much ¥ moves — and does action mix flip to retrieval_refresh?

Usage:
  PYTHONPATH=. python3 scripts/agod/cs_assist_hop_sensitivity.py
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
    "09_cs_opportunity_breakeven.sql",
    "10_cs_action_contain_split.sql",
    "11_cs_week_split_coverage.sql",
    "12_cs_payback_contain.sql",
    "13_cs_hf_knob_bridge.sql",
    "14_cs_marginal_day.sql",
    "15_cs_payback_contain_price_stress.sql",
    "16_cs_detection_delay_profit.sql",
    "17_cs_fully_loaded_capture.sql",
    "18_cs_ops_onepager.sql",
    "19_cs_weekly_fully_loaded.sql",
    "20_cs_hf_hop_yen_band.sql",
]


def _load_demo():
    spec = importlib.util.spec_from_file_location(
        "biz_demo", ROOT / "scripts" / "agod" / "run_biz_value_sql_demo.py"
    )
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def _action_mix(con) -> dict:
    rows = con.execute(
        """
        SELECT action_type, n_days, sessions, incremental_yen, avg_rag_hit
        FROM vw_cs_assist_action_increment
        ORDER BY action_type
        """
    ).fetchall()
    out = {}
    total_days = sum(r[1] for r in rows) or 1
    for action, n_days, sessions, yen, rag in rows:
        out[action] = {
            "n_days": int(n_days),
            "sessions": int(sessions),
            "gross_yen": float(yen),
            "avg_rag_hit": float(rag) if rag is not None else None,
            "day_share": round(float(n_days) / total_days, 3),
        }
    return out


def run_scenario(mod, duckdb, hop: dict) -> dict:
    con = duckdb.connect(":memory:")
    for name in SQL_FILES:
        con.execute((SQL / name).read_text())
    mod.seed(con, hop=hop)

    dash = con.execute("SELECT * FROM vw_cs_assist_exec_dashboard").fetchone()
    dcols = [d[0] for d in con.description]
    d = dict(zip(dcols, dash))

    full = con.execute("SELECT * FROM vw_cs_assist_fully_loaded_capture").fetchone()
    fcols = [c[0] for c in con.description]
    f = dict(zip(fcols, full))

    rates = con.execute("SELECT * FROM vw_cs_assist_rate_compare").fetchone()
    rcols = [c[0] for c in con.description]
    r = dict(zip(rcols, rates))

    mix = _action_mix(con)
    top_retrieval_share = float((mix.get("retrieval_refresh") or {}).get("day_share") or 0.0)

    con.close()
    return {
        "scenario": hop["scenario"],
        "hop_ratio": float(hop["hop_ratio"]),
        "fire_halluc": float(hop["fire_halluc"]),
        "quiet_halluc": float(hop["quiet_halluc"]),
        "rag_low": float(hop["rag_low"]),
        "rag_ok": float(hop["rag_ok"]),
        "rag_threshold": float(hop["rag_threshold"]),
        "tickets_avoided": d.get("tickets_avoided"),
        "refunds_avoided": d.get("refunds_avoided"),
        "extra_sessions_contained": d.get("extra_sessions_contained"),
        "gross_yen": float(d.get("gross_yen") or 0),
        "net_yen": float(d.get("net_yen") or 0),
        "fully_loaded_net_yen": float(f.get("fully_loaded_net_yen") or 0),
        "audit_cost_yen": float(f.get("audit_cost_yen") or 0),
        "action_day_cost_yen": float(f.get("action_day_cost_yen") or 0),
        "ticket_rate_ignored_pct": float(r.get("ticket_rate_ignored") or 0) * 100.0,
        "ticket_rate_acted_pct": float(r.get("ticket_rate_acted") or 0) * 100.0,
        "contain_lift_pp": float(
            (r.get("contain_rate_acted") or 0) - (r.get("contain_rate_ignored") or 0)
        )
        * 100.0,
        "top_action": d.get("top_action"),
        "retrieval_day_share": top_retrieval_share,
        "action_mix": mix,
        "prod_mid_monthly_net_yen": d.get("prod_mid_monthly_net_yen"),
    }


def build_scenarios(base: dict) -> list[dict]:
    return [
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
        # rag worse → more days below thr → retrieval_refresh share ↑
        {
            **base,
            "scenario": "rag_worse",
            "rag_low": float(np.clip(base["rag_low"] * 0.4, 0.05, 0.20)),
            "rag_ok": float(
                np.clip(min(base["rag_ok"] * 0.45, base["rag_threshold"] - 0.05), 0.10, 0.40)
            ),
            "rag_threshold": float(np.clip(base["rag_threshold"] + 0.20, 0.45, 0.70)),
        },
        # rag better → fewer retrieval days
        {
            **base,
            "scenario": "rag_better",
            "rag_low": float(np.clip(base["rag_low"] + 0.25, 0.35, 0.70)),
            "rag_ok": float(np.clip(base["rag_ok"] + 0.15, 0.55, 0.90)),
            "rag_threshold": float(np.clip(base["rag_threshold"] - 0.10, 0.20, 0.35)),
        },
    ]


def attach_deltas(rows: list[dict]) -> list[dict]:
    base = next(r for r in rows if r["scenario"] == "base_hop")
    out = []
    for r in rows:
        x = dict(r)
        x["delta_gross_vs_base"] = round(r["gross_yen"] - base["gross_yen"], 0)
        x["delta_net_vs_base"] = round(r["net_yen"] - base["net_yen"], 0)
        x["delta_fully_loaded_vs_base"] = round(
            r["fully_loaded_net_yen"] - base["fully_loaded_net_yen"], 0
        )
        x["delta_retrieval_share_vs_base"] = round(
            r["retrieval_day_share"] - base["retrieval_day_share"], 3
        )
        out.append(x)
    return out


def run_band() -> dict:
    duckdb = __import__("duckdb")
    mod = _load_demo()
    base = mod.load_hf_hop()
    rows = attach_deltas(
        [run_scenario(mod, duckdb, h) for h in build_scenarios(base)]
    )
    by = {r["scenario"]: r for r in rows}
    weak, base_r, strong = by["weak_hop"], by["base_hop"], by["strong_hop"]
    rag_w, rag_b = by["rag_worse"], by["rag_better"]
    summary = {
        "base_gross_yen": base_r["gross_yen"],
        "base_net_yen": base_r["net_yen"],
        "base_fully_loaded_net_yen": base_r["fully_loaded_net_yen"],
        "weak_to_strong_net_band": [
            weak["net_yen"],
            base_r["net_yen"],
            strong["net_yen"],
        ],
        "weak_to_strong_fully_loaded_band": [
            weak["fully_loaded_net_yen"],
            base_r["fully_loaded_net_yen"],
            strong["fully_loaded_net_yen"],
        ],
        "rag_worse_retrieval_share": rag_w["retrieval_day_share"],
        "rag_better_retrieval_share": rag_b["retrieval_day_share"],
        "base_retrieval_share": base_r["retrieval_day_share"],
        "external_one_liner_cn": (
            f"HF hop 弱→强：全成本净¥{int(weak['fully_loaded_net_yen'])}"
            f"→{int(base_r['fully_loaded_net_yen'])}"
            f"→{int(strong['fully_loaded_net_yen'])}；"
            f"rag 变差时 retrieval 天占比 "
            f"{base_r['retrieval_day_share']:.0%}→{rag_w['retrieval_day_share']:.0%}，"
            f"变好时→{rag_b['retrieval_day_share']:.0%}。"
        ),
    }
    return {"rows": rows, "summary": summary}


def _yen(x) -> str:
    try:
        return f"¥{int(float(x)):,}"
    except (TypeError, ValueError):
        return str(x)


def render_exec_dashboard(band: dict) -> str:
    rows = band["rows"]
    by = {r["scenario"]: r for r in rows}
    base_row = by["base_hop"]
    weak, strong = by["weak_hop"], by["strong_hop"]
    hop_table = "\n".join(
        f"| {r['scenario']} | {r['fire_halluc']:.2f} | {r['hop_ratio']:.1f} | "
        f"{r.get('tickets_avoided')} | {r.get('refunds_avoided')} | "
        f"{r.get('extra_sessions_contained')} | **{_yen(r.get('gross_yen'))}** | "
        f"**{_yen(r.get('net_yen'))}** | **{_yen(r.get('fully_loaded_net_yen'))}** | "
        f"{_yen(r.get('delta_fully_loaded_vs_base'))} |"
        for r in rows
        if r["scenario"] in ("weak_hop", "base_hop", "strong_hop")
    )
    rag_table = "\n".join(
        f"| {r['scenario']} | {r['rag_low']:.2f}/{r['rag_ok']:.2f} | "
        f"{r['rag_threshold']:.2f} | {r['retrieval_day_share']:.0%} | "
        f"{r.get('top_action')} | **{_yen(r.get('fully_loaded_net_yen'))}** | "
        f"{r['delta_retrieval_share_vs_base']:+.0%} |"
        for r in rows
        if r["scenario"] in ("rag_worse", "base_hop", "rag_better")
    )
    return f"""# 客服助手经营总看板（大迭代）

## 对外一句

幻觉制度跳变期内，客服助手相对不动作：少 **{base_row.get('tickets_avoided')}** 单工单、
少 **{base_row.get('refunds_avoided')}** 单退款、多承接 **{base_row.get('extra_sessions_contained')}** 次；
毛贡献 **{_yen(base_row.get('gross_yen'))}**，扣审计后净贡献 **{_yen(base_row.get('net_yen'))}**，
全成本净 **{_yen(base_row.get('fully_loaded_net_yen'))}**；
优先动作 **`{base_row.get('top_action')}`**。
{band['summary']['external_one_liner_cn']}

## Before → After（base hop）

| 费率 | quiet→ignored→acted（工单） | 承接提升 |
|------|------------------------------|----------|
| base | {base_row.get('ticket_rate_ignored_pct'):.2f}%→{base_row.get('ticket_rate_acted_pct'):.2f}% | **{base_row.get('contain_lift_pp'):.2f} pp** |

## HF hop 强度 → 全成本净带

| 情景 | fire_halluc | hop_ratio | 少工单 | 少退款 | 多承接 | 毛¥ | 净¥ | 全成本净¥ | Δ全成本 vs base |
|------|-------------|-----------|--------|--------|--------|-----|-----|-----------|-----------------|
{hop_table}

弱→强全成本净带：**{_yen(weak.get('fully_loaded_net_yen'))} → {_yen(base_row.get('fully_loaded_net_yen'))} → {_yen(strong.get('fully_loaded_net_yen'))}**。

## rag_hit knobs → 动作混比

| 情景 | rag_low/ok | thr | retrieval天占比 | 优先动作 | 全成本净¥ | Δ retrieval占比 |
|------|------------|-----|-----------------|----------|-----------|-----------------|
{rag_table}

读法：HF `rag_*` 直接进 seed；thr 抬高 / rag 变差 → 更多天走 `retrieval_refresh`。

## 口径

```
seed: hop_ratio → ignore 工单/退款 bump；fire_halluc → 火情标签率；
      rag_low/ok/thr → 动作路由（retrieval vs rollback/audit）
毛¥ = 费率差 × acted会话 × 单价
净¥ = 毛¥ − acted 审计人力
全成本净¥ = 毛¥ − 审计 − 动作日成本
```

产物：`cs_assist_hop_sensitivity.json` / `cs_assist_hop_yen_band_summary.json`
"""


def contribution_section(band: dict) -> str:
    """Markdown section injected into CS_ASSISTANT_CONTRIBUTION.md."""
    rows = band["rows"]
    by = {r["scenario"]: r for r in rows}
    hop_lines = []
    for key in ("weak_hop", "base_hop", "strong_hop"):
        r = by[key]
        hop_lines.append(
            f"| {r['scenario']} | {r['fire_halluc']:.2f} | {r['hop_ratio']:.1f} | "
            f"{r['tickets_avoided']} | {r['refunds_avoided']} | "
            f"{r['extra_sessions_contained']} | **¥{int(r['gross_yen'])}** | "
            f"**¥{int(r['net_yen'])}** | **¥{int(r['fully_loaded_net_yen'])}** | "
            f"¥{int(r['delta_fully_loaded_vs_base']):+d} |"
        )
    rag_lines = []
    for key in ("rag_worse", "base_hop", "rag_better"):
        r = by[key]
        rag_lines.append(
            f"| {r['scenario']} | {r['rag_low']:.2f}/{r['rag_ok']:.2f} | "
            f"{r['rag_threshold']:.2f} | {r['retrieval_day_share']:.0%} | "
            f"{r['top_action']} | **¥{int(r['fully_loaded_net_yen'])}** | "
            f"{r['delta_retrieval_share_vs_base']:+.0%} |"
        )
    return f"""## HF hop / rag 驱动 seed → 贡献带（可对账）

口径：同一套 SQL，只改 HF knobs 重 seed；Δ 相对 `base_hop`。

| 情景 | fire_halluc | hop_ratio | 少工单 | 少退款 | 多承接 | 毛¥ | 净¥ | 全成本净¥ | Δ全成本 |
|------|-------------|-----------|--------|--------|--------|-----|-----|-----------|---------|
{chr(10).join(hop_lines)}

| 情景 | rag_low/ok | thr | retrieval天占比 | 优先动作 | 全成本净¥ | Δ retrieval占比 |
|------|------------|-----|-----------------|----------|-----------|-----------------|
{chr(10).join(rag_lines)}

对外一句：{band['summary']['external_one_liner_cn']}
"""


def main() -> int:
    band = run_band()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "cs_assist_hop_sensitivity.json").write_text(
        json.dumps(band["rows"], ensure_ascii=False, indent=2, default=str)
    )
    (OUT / "cs_assist_hop_yen_band_summary.json").write_text(
        json.dumps(band["summary"], ensure_ascii=False, indent=2, default=str)
    )
    md = render_exec_dashboard(band)
    (OUT / "CS_ASSISTANT_EXEC_DASHBOARD.md").write_text(md)
    (DOCS / "CS_ASSISTANT_EXEC_DASHBOARD.md").write_text(md)
    print(md)
    print(f"wrote {OUT / 'cs_assist_hop_yen_band_summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
