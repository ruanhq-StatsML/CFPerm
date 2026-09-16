#!/usr/bin/env python3
"""Generate a weekly ops brief for 客服助手 incremental contribution.

Reads results from the biz-value SQL demo (JSON) and writes a CN brief that
ops/product can paste into the weekly review — tickets, refunds, ¥, WoW.

Usage:
  python3 scripts/agod/cs_assist_weekly_ops_brief.py
  python3 scripts/agod/cs_assist_weekly_ops_brief.py --out docs/biz/CS_ASSISTANT_WEEKLY_OPS_BRIEF.md
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
IN = ROOT / "results" / "agod" / "biz_value_sql"
DEFAULT_OUT = IN / "CS_ASSISTANT_WEEKLY_OPS_BRIEF.md"
DOCS_OUT = ROOT / "docs" / "biz" / "CS_ASSISTANT_WEEKLY_OPS_BRIEF.md"


def _load(name: str) -> list | dict:
    path = IN / name
    if not path.exists():
        raise SystemExit(f"missing {path}; run scripts/agod/run_biz_value_sql_demo.py first")
    return json.loads(path.read_text())


def _pct(x) -> str:
    try:
        return f"{100.0 * float(x):.2f}%"
    except (TypeError, ValueError):
        return str(x)


def _week_label(ws) -> str:
    s = str(ws)[:10]
    return s


def build_brief() -> str:
    weekly = _load("cs_assist_weekly.json")
    inc = _load("cs_assist_increment.json")
    net = _load("cs_assist_net.json")
    actions = _load("cs_assist_action_net.json")
    hop = _load("hf_hop_knobs.json") if (IN / "hf_hop_knobs.json").exists() else {}
    traffic = _load("cs_assist_traffic.json")

    c = inc[0] if isinstance(inc, list) and inc else (inc or {})
    n = net[0] if isinstance(net, list) and net else (net or {})

    # WoW on fire weeks (days_acted + days_ignored > 0)
    fire_weeks = [
        w for w in weekly
        if float(w.get("days_acted") or 0) + float(w.get("days_ignored") or 0) > 0
    ]
    wow_rows = []
    for i, w in enumerate(fire_weeks):
        prev = fire_weeks[i - 1] if i else None
        d_ticket = (
            float(w["ticket_rate"]) - float(prev["ticket_rate"]) if prev else None
        )
        d_net = (
            float(w["net_contrib_yen"]) - float(prev["net_contrib_yen"]) if prev else None
        )
        wow_rows.append(
            f"| {_week_label(w['week_start'])} | {int(w['sessions'])} | "
            f"{_pct(w['ticket_rate'])} | {_pct(w['containment_rate'])} | "
            f"{int(w['days_acted'])}/{int(w['days_ignored'])} | "
            f"¥{int(w['net_contrib_yen'])} | "
            f"{('+' if d_ticket and d_ticket >= 0 else '') + f'{100*d_ticket:.2f}pp' if d_ticket is not None else '—'} | "
            f"{('+' if d_net and d_net >= 0 else '') + f'¥{int(d_net)}' if d_net is not None else '—'} |"
        )
    wow_table = "\n".join(wow_rows) if wow_rows else "| (no fire weeks) |||||||"

    action_rows = []
    for a in actions:
        action_rows.append(
            f"| {a['action_type']} | {int(a['n_days'])} | "
            f"¥{int(a['gross_yen'])} | ¥{int(a['action_day_cost_yen'])} | "
            f"**¥{int(a['net_yen_after_action_cost'])}** |"
        )
    action_table = "\n".join(action_rows) if action_rows else "| (none) ||||"

    prod = next(
        (
            t
            for t in traffic
            if t.get("scenario") in ("prod_mid", "prod_mid", "prod_mid")
            or int(t.get("daily_sessions") or 0) == 10000
        ),
        None,
    )
    prod_line = (
        f"生产中等流量（日 {int(prod['daily_sessions']):,} 会话）月毛贡献约 **¥{int(prod['monthly_yen']):,}**。"
        if prod
        else ""
    )

    # Ops call: which fire week hurt most / which action to prioritize
    worst = min(fire_weeks, key=lambda w: float(w["net_contrib_yen"])) if fire_weeks else None
    best_action = (
        max(actions, key=lambda a: float(a["net_yen_after_action_cost"])) if actions else None
    )

    callouts = []
    if worst is not None:
        callouts.append(
            f"- 最差火情周 `{_week_label(worst['week_start'])}`：净贡献 ¥{int(worst['net_contrib_yen'])}，"
            f"工单率 {_pct(worst['ticket_rate'])}，acted/ignored={int(worst['days_acted'])}/{int(worst['days_ignored'])}。"
        )
    if best_action is not None:
        callouts.append(
            f"- 优先动作 `{best_action['action_type']}`：扣动作成本后净贡献 ¥{int(best_action['net_yen_after_action_cost'])}"
            f"（毛 ¥{int(best_action['gross_yen'])} − 日成本 ¥{int(best_action['action_day_cost_yen'])}）。"
        )
    callouts.append(
        f"- 区间对照：少工单 {c.get('tickets_avoided')} / 少退款 {c.get('refunds_avoided')} / "
        f"多承接 {c.get('extra_sessions_contained')}；毛 ¥{c.get('incremental_yen')} → 净 ¥{n.get('net_incremental_yen')}。"
    )
    if hop:
        callouts.append(
            f"- HF hop 驱动：quiet={hop.get('quiet_halluc')} / fire(biz)={hop.get('fire_halluc')} / "
            f"ratio={hop.get('hop_ratio')}（`{hop.get('source')}`）。"
        )

    return f"""# 客服助手周经营简报（可复现）

生成自 `results/agod/biz_value_sql/*.json`。只讲工单 / 退款 / 承接 / ¥。

## 本周会一句

在幻觉制度跳变对照下，客服助手相对不动作：少 **{c.get('tickets_avoided')}** 单工单、
少 **{c.get('refunds_avoided')}** 单退款、多承接 **{c.get('extra_sessions_contained')}** 次；
毛贡献 **¥{c.get('incremental_yen')}**，扣审计后净贡献 **¥{n.get('net_incremental_yen')}**。
{prod_line}

## 运营点拨

{chr(10).join(callouts)}

## 火情周 WoW

| 周起始 | 会话 | 工单率 | 承接率 | acted/ignored | 净贡献¥ | 工单率Δ | 净贡献Δ |
|--------|------|--------|--------|----------------|---------|---------|---------|
{wow_table}

## 动作净贡献（扣动作日成本）

| 动作 | 天数 | 毛¥ | 动作日成本¥ | 净¥ |
|------|------|-----|-------------|-----|
{action_table}

## 口径（钉死）

```
毛增量¥ = 费率差 × acted会话 × 单价（工单/退款/承接）
净增量¥ = 毛增量¥ − acted 审计人力
动作净¥ = 该动作毛增量¥ − 动作日成本 −（audit_topk 的审计件成本分摊）
```

数据：`cs_assist_weekly.json` / `cs_assist_increment.json` / `cs_assist_action_net.json`
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--also-docs", action="store_true", default=True)
    args = ap.parse_args()
    brief = build_brief()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(brief)
    if args.also_docs:
        DOCS_OUT.write_text(brief)
        print(f"wrote {DOCS_OUT}")
    print(brief)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
