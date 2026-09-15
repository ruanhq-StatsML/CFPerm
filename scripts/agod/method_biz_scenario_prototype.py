#!/usr/bin/env python3
"""Method ↔ Biz-scenario ↔ Unit-econ relevance prototype.

Pulls HF landing signals (OnlineRFPerm / PO-risk / style AUC) and the
客服助手 scenario suite (hop intensity + unit-econ band) into one
readable map: what is the SAME, what DIFFERS, and why each layer is
relevant to the next.

Does NOT claim AUROC = business value. Value starts only after
act-vs-ignore volumes × unit prices.

Usage:
  python3 scripts/agod/method_biz_scenario_prototype.py
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HF = ROOT / "results" / "agod" / "hf_landing"
BIZ = ROOT / "results" / "agod" / "biz_value_sql"
OUT = ROOT / "results" / "agod" / "method_biz_proto"
DOCS = ROOT / "docs" / "biz"


def _load(path: Path, default=None):
    if not path.exists():
        return default
    return json.loads(path.read_text())


def _yen(x) -> str:
    try:
        return f"¥{int(float(x)):,}"
    except (TypeError, ValueError):
        return str(x)


def _pct(x, digits=1) -> str:
    try:
        return f"{100.0 * float(x):.{digits}f}%"
    except (TypeError, ValueError):
        return str(x)


def build_payload() -> dict:
    halu = _load(HF / "halu_regime_rag.json", {}) or {}
    hh = _load(HF / "hh_judge_style.json", {}) or {}
    hop_sens = _load(BIZ / "cs_assist_hop_sensitivity.json", []) or []
    unit_band = (_load(BIZ / "cs_assist_unit_econ_band.json", []) or [{}])[0]
    unit_sens = _load(BIZ / "cs_assist_unit_econ_sensitivity.json", []) or []
    exec_dash = (_load(BIZ / "cs_assist_exec_dashboard.json", []) or [{}])[0]
    payback = _load(BIZ / "cs_assist_action_payback.json", []) or []
    traffic = _load(BIZ / "cs_assist_traffic.json", []) or []
    week_attr = _load(BIZ / "cs_assist_week_attribution.json", []) or []
    knobs = _load(BIZ / "hf_hop_knobs.json", {}) or {}

    hop = halu.get("hop_at_cut") or {}
    hop_rank = hop.get("ranking") or {}
    quiet = (halu.get("hop_quiet") or {}).get("ranking") or {}

    by_scen = {r.get("scenario"): r for r in hop_sens}
    base = by_scen.get("base_hop") or exec_dash or {}
    weak = by_scen.get("weak_hop") or {}
    strong = by_scen.get("strong_hop") or {}

    prod = next((t for t in traffic if t.get("scenario") == "prod_mid"), {})

    # Relevance chain rows
    chain = [
        {
            "layer": "方法论 · OnlineRFPerm",
            "object": "P(Y|X) / P(X) 是否进入新制度",
            "signal": f"fired={hop.get('fired')}, ratio≈{float(hop.get('ratio') or 0):.2f}",
            "same_as_biz": "给「要不要动作」打时间戳；不直接给钱",
            "diff_from_biz": "不产出工单/退款/¥；AUROC≠贡献",
            "relevance": "fire → 打开 act-vs-ignore 对照窗（制度跳变期）",
        },
        {
            "layer": "方法论 · po_risk0",
            "object": "跳变下谁该进审计队列",
            "signal": f"P@10={hop_rank.get('precision_at_10')}, auroc={hop_rank.get('auroc_po_risk0')}",
            "same_as_biz": "决定审计人力花在哪",
            "diff_from_biz": "排序器不是事实核查器；P@10≠少退款数",
            "relevance": "驱动 audit_topk 动作臂 + 审计件成本（净贡献扣减）",
        },
        {
            "layer": "方法论 · RAG / router",
            "object": "检索缺口 vs 生成制度 vs 路由漂移",
            "signal": (
                f"rag_top10={hop_rank.get('mean_rag_hit_top10')}, "
                f"domain_auc≈{float((halu.get('router_task_shift') or {}).get('domain_auc_pre_post') or 0):.2f}"
            ),
            "same_as_biz": "决定「刷检索」还是「回滚模型」",
            "diff_from_biz": "不解释 GMV 因果",
            "relevance": "路由到 retrieval_refresh / model_rollback 动作拆分",
        },
        {
            "layer": "方法论 · style_domain_auc（对照面）",
            "object": "P(X) 文风/画风画像",
            "signal": f"style_auc≈{float(hh.get('style_domain_auc') or 0):.3f}, judge_err_ratio≈{float(hh.get('judge_err_ratio') or 0):.2f}",
            "same_as_biz": "也是制度/画像监控",
            "diff_from_biz": "客服助手主账不走 CTR/品牌路径；勿与幻觉 hop 混账",
            "relevance": "证明要拆轴：风格工单 ≠ 客服工单；防误伤偏好头",
        },
        {
            "layer": "业务情景 · hop 强度",
            "object": "同一 fire 下跳变强弱 → 量差",
            "signal": (
                f"weak/base/strong net = "
                f"{_yen(weak.get('net_yen'))} / {_yen(base.get('net_yen'))} / {_yen(strong.get('net_yen'))}"
            ),
            "same_as_biz": "仍用 act-vs-ignore 对照",
            "diff_from_method": "把 ratio/幻觉率翻译成工单·退款·承接量",
            "relevance": "方法论 hop 强度的业务可感知带（量）",
        },
        {
            "layer": "业务情景 · 流量缩放",
            "object": "每千会话单价 × 日会话",
            "signal": f"prod_mid 日{int(prod.get('daily_sessions') or 0):,} → 月毛 {_yen(prod.get('monthly_yen'))}",
            "same_as_biz": "线性外推，不改费率",
            "diff_from_method": "方法层没有「会话量」概念",
            "relevance": "把子集原型外推到生产量级",
        },
        {
            "layer": "单位经济 · 单价敏感度",
            "object": "量固定，扫工单/退款/承接单价",
            "signal": (
                f"净带 {_yen(unit_band.get('net_yen_low'))} → "
                f"{_yen(unit_band.get('net_yen_base'))} → "
                f"{_yen(unit_band.get('net_yen_high'))}"
            ),
            "same_as_biz": "仍基于同一套少工单/退款/承接量",
            "diff_from_method": "与 fire/AUROC 无关；财务假设层",
            "relevance": "回答「单价变了贡献还站不站得住」——方法层答不了",
        },
    ]

    # Compact matrix for JSON
    same_diff = {
        "same": [
            "都对「制度/画像是否变了」敏感",
            "都需要对照臂（quiet / fire_ignored / fire_acted），不能只看绝对值",
            "动作路由都依赖轴拆分（concept vs covariate / RAG vs generation）",
        ],
        "different": [
            "方法论输出：fired / ratio / AUROC / P@10 / style_auc（诊断）",
            "业务情景输出：少工单 / 少退款 / 多承接 / 动作净¥（运营）",
            "单位经济输出：同一量下的 ¥ 带（财务假设）",
            "方法论不保证因果归因；业务账是对照代理，不是 RCT",
        ],
        "relevance_one_liner": (
            "OnlineRFPerm 决定「哪段日子算跳变对照窗」；"
            "hop 情景把跳变强度变成量；"
            "单位经济把量变成可对账的 ¥ 带——三层叠乘，缺一不可。"
        ),
    }

    return {
        "method": {
            "halu": {
                "dataset": halu.get("dataset"),
                "fired": hop.get("fired"),
                "ratio": hop.get("ratio"),
                "halluc_quiet": quiet.get("halluc_rate"),
                "halluc_fire": hop_rank.get("halluc_rate"),
                "po_p10": hop_rank.get("precision_at_10"),
                "po_auroc": hop_rank.get("auroc_po_risk0"),
                "rag_top10": hop_rank.get("mean_rag_hit_top10"),
                "reading": halu.get("reading"),
            },
            "hh_style": {
                "dataset": hh.get("dataset"),
                "style_domain_auc": hh.get("style_domain_auc"),
                "judge_err_ratio": hh.get("judge_err_ratio"),
                "fired": (hh.get("hop_at_cut") or {}).get("fired"),
                "reading": hh.get("reading"),
            },
            "knobs_into_seed": knobs,
        },
        "biz_scenarios": {
            "hop_band": {
                "weak_net": weak.get("net_yen"),
                "base_net": base.get("net_yen"),
                "strong_net": strong.get("net_yen"),
                "base_tickets": base.get("tickets_avoided"),
                "base_refunds": base.get("refunds_avoided"),
                "base_contained": base.get("extra_sessions_contained"),
                "base_gross": base.get("gross_yen"),
                "top_action": base.get("top_action"),
            },
            "before_after": {
                "ticket_quiet": base.get("ticket_rate_quiet_pct"),
                "ticket_ignored": base.get("ticket_rate_ignored_pct"),
                "ticket_acted": base.get("ticket_rate_acted_pct"),
                "contain_lift_pp": base.get("contain_rate_lift_pp"),
            },
            "traffic_prod_mid_monthly_yen": prod.get("monthly_yen"),
            "week_attribution": week_attr,
            "action_payback": payback,
        },
        "unit_econ": {
            "band": unit_band,
            "scenarios": unit_sens,
        },
        "chain": chain,
        "same_diff": same_diff,
        "external_one_liner_cn": (
            f"方法论确认幻觉制度跳变（ratio≈{float(hop.get('ratio') or 0):.1f}）后，"
            f"客服助手相对不动作少 {base.get('tickets_avoided')} 工单 / "
            f"{base.get('refunds_avoided')} 退款 / 多承接 {base.get('extra_sessions_contained')}；"
            f"净贡献 {_yen(base.get('net_yen'))}，"
            f"单价±20% 净带 {_yen(unit_band.get('net_yen_low'))}–{_yen(unit_band.get('net_yen_high'))}；"
            f"优先动作 {base.get('top_action')}。"
        ),
    }


def render_cn(p: dict) -> str:
    m = p["method"]["halu"]
    hh = p["method"]["hh_style"]
    b = p["biz_scenarios"]["hop_band"]
    ba = p["biz_scenarios"]["before_after"]
    ue = p["unit_econ"]["band"]
    sd = p["same_diff"]

    hop_table = (
        f"| weak_hop | {_yen(b.get('weak_net'))} |\n"
        f"| base_hop | **{_yen(b.get('base_net'))}** |\n"
        f"| strong_hop | {_yen(b.get('strong_net'))} |"
    )

    chain_rows = []
    for row in p["chain"]:
        chain_rows.append(
            f"| {row['layer']} | {row['object']} | `{row['signal']}` | "
            f"{row.get('same_as_biz') or row.get('same_as_method', '')} | "
            f"{row.get('diff_from_biz') or row.get('diff_from_method', '')} | "
            f"{row['relevance']} |"
        )
    chain_table = "\n".join(chain_rows)

    pay_rows = []
    for a in p["biz_scenarios"].get("action_payback") or []:
        pay_rows.append(
            f"| {a.get('action_type')} | {a.get('payback_days')} 天 | "
            f"{a.get('payback_bucket')} | {a.get('net_roi_multiple')}x |"
        )
    pay_table = "\n".join(pay_rows) if pay_rows else "| (run biz demo first) |||"

    week_rows = []
    for w in p["biz_scenarios"].get("week_attribution") or []:
        week_rows.append(
            f"| {str(w.get('week_start'))[:10]} | {w.get('tickets_avoided')} | "
            f"{w.get('refunds_avoided')} | **{_yen(w.get('gross_yen'))}** | "
            f"{w.get('pct_of_total_gross')}% |"
        )
    week_table = "\n".join(week_rows) if week_rows else "| (run biz demo first) ||||"

    unit_rows = []
    for s in p["unit_econ"].get("scenarios") or []:
        if s.get("scenario") in ("base", "all_minus20", "all_plus20", "ticket_plus20", "refund_plus20"):
            unit_rows.append(
                f"| {s.get('scenario')} | ¥{s.get('ticket_cost')} | ¥{s.get('refund_cost')} | "
                f"¥{s.get('contain_value')} | {_yen(s.get('gross_yen'))} | {_yen(s.get('net_yen'))} |"
            )
    unit_table = "\n".join(unit_rows)

    return f"""# 方法论 × 业务情景 × 单位经济 — 全面对照原型

> 目标：把 OnlineRFPerm / PO-risk / 风格轴 与「客服助手全面情景 + 单位经济敏感度」放在同一张图上，
> 讲清**异同**与**relevance**。价值不来自 AUROC，来自 **act-vs-ignore 量 × 单价**。

## 对外一句

{p['external_one_liner_cn']}

## 1. 三层叠乘（relevance 主链）

```text
[方法论] OnlineRFPerm fire / po_risk0 / RAG轴
    → 打开「制度跳变对照窗」+ 动作路由
[业务情景] hop弱/基/强 × 流量 × 动作拆分
    → 少工单 / 少退款 / 多承接（量）
[单位经济] 工单¥ / 退款¥ / 承接¥ （±20%扫）
    → 毛/净贡献¥带（钱）
```

**一句话 relevance：**  
{sd['relevance_one_liner']}

## 2. 异同矩阵

### 相同（same）

{chr(10).join(f'- {x}' for x in sd['same'])}

### 不同（different）

{chr(10).join(f'- {x}' for x in sd['different'])}

### 对照表（层 × 对象 × 信号 × 同/异 × relevance）

| 层 | 对象 | 信号（本原型） | 与业务的同 | 与业务的异 | relevance |
|----|------|----------------|------------|------------|-----------|
{chain_table}

## 3. 方法论侧（HF 子集落地）

### 3.1 幻觉制度（HaluEval）— 客服助手主账来源

| 项 | 值 |
|----|-----|
| dataset | `{m.get('dataset')}` |
| cut fire | `{m.get('fired')}`，ratio≈`{float(m.get('ratio') or 0):.2f}` |
| 幻觉率 quiet → fire | `{m.get('halluc_quiet')}` → `{m.get('halluc_fire')}` |
| po_risk0 P@10 / AUROC | `{m.get('po_p10')}` / `{m.get('po_auroc')}` |
| Top-10 rag_hit | `{m.get('rag_top10')}` |

读法（方法层）：{m.get('reading')}

**与业务的衔接：** `fired+ratio` → seed 的 `fire_halluc/hop_ratio`；`rag_hit` → `retrieval_refresh` vs `model_rollback`；`P@10` → 审计成本进入净贡献。

### 3.2 文风/偏好轴（HH-RLHF）— 对照面（不要混进客服账）

| 项 | 值 |
|----|-----|
| style_domain_auc | `{hh.get('style_domain_auc')}` |
| judge_err_ratio | `{hh.get('judge_err_ratio')}` |
| concept fire | `{(hh.get('fired'))}` |

读法：{hh.get('reading')}

**relevance：** 高 style_auc + 有/无 concept-fire 决定「只调素材」还是「动偏好」——**钱走 CTR/品牌路径，不是客服工单路径**。本原型把它放在对照列，防止方法论信号被当成同一本账。

## 4. 业务情景侧（全面场景）

### 4.1 Before → After（base hop）

| 费率 | quiet | fire+不动作 | fire+动作 |
|------|-------|-------------|-----------|
| 工单率 | {ba.get('ticket_quiet')}% | {ba.get('ticket_ignored')}% | {ba.get('ticket_acted')}% |
| 承接率提升 | — | — | **{ba.get('contain_lift_pp')} pp** |

### 4.2 Hop 强度情景（方法论 ratio 的业务翻译）

| 情景 | 净贡献¥ |
|------|---------|
{hop_table}

量（base）：少工单 **{b.get('base_tickets')}** / 少退款 **{b.get('base_refunds')}** / 多承接 **{b.get('base_contained')}** → 毛 **{_yen(b.get('base_gross'))}**，优先动作 **`{b.get('top_action')}`**。

生产中等流量月毛外推：**{_yen(p['biz_scenarios'].get('traffic_prod_mid_monthly_yen'))}**。

### 4.3 周归因（对账）

| 周 | 少工单 | 少退款 | 毛¥ | 占比 |
|----|--------|--------|-----|------|
{week_table}

### 4.4 动作回本

| 动作 | 回本 | 分档 | 净ROI |
|------|------|------|-------|
{pay_table}

## 5. 单位经济敏感度（财务层）

量固定（少工单/退款/承接不变），只扫单价：

| 情景 | 工单¥ | 退款¥ | 承接¥ | 毛¥ | 净¥ |
|------|-------|-------|-------|-----|-----|
{unit_table}

**净贡献带（全单价 ±20%）：** {_yen(ue.get('net_yen_low'))} → **{_yen(ue.get('net_yen_base'))}** → {_yen(ue.get('net_yen_high'))}  
（带宽 / base ≈ {ue.get('net_band_width_vs_base')}）

### 与方法论的异同（专指单位经济）

| | 内容 |
|--|------|
| **同** | 仍站在同一套 act-vs-ignore 量上；不另编故事 |
| **异** | 完全不依赖 fire/AUROC；回答的是财务假设，不是检测力 |
| **relevance** | 方法层证明「该动作」；情景层证明「有多少量」；单位经济证明「换成你们财务口径还值不值」 |

## 6. 一张总图：谁回答什么问题

| 问题 | 该问哪一层 | 本原型答案 |
|------|------------|------------|
| 模型/标签制度跳了吗？ | 方法论 OnlineRFPerm | fire=1，ratio≈{float(m.get('ratio') or 0):.1f} |
| 该审哪些样本？ | 方法论 po_risk0 | P@10={m.get('po_p10')} |
| 刷检索还是回滚？ | 方法论 RAG + 业务动作拆分 | 优先 `{b.get('top_action')}` |
| 动作相对不动作少多少单？ | 业务情景 | {b.get('base_tickets')} 工单 / {b.get('base_refunds')} 退款 |
| 值多少钱？ | 单位经济 | 净 {_yen(b.get('base_net'))}；带 {_yen(ue.get('net_yen_low'))}–{_yen(ue.get('net_yen_high'))} |
| 文风漂了要不要动客服账？ | 方法论风格轴（对照） | **不要混账**；走素材/CTR 路径 |

## 7. 怎么跑

```bash
# 1) 业务账 + 单位经济（若尚未跑）
PYTHONPATH=. python3 scripts/agod/run_biz_value_sql_demo.py
python3 scripts/agod/cs_assist_hop_sensitivity.py

# 2) 本对照原型
python3 scripts/agod/method_biz_scenario_prototype.py
# → docs/biz/METHOD_BIZ_SCENARIO_PROTOTYPE.md
# → results/agod/method_biz_proto/summary.json
```

## 8. 边界（必须写死）

1. OnlineRFPerm / PO-risk **不是**事实核查器，也**不是**因果归因。
2. 业务增量是 **fire 期内 acted vs ignored 的对照代理**，不是 RCT。
3. 单位经济带只扫单价；流量情景只缩放会话量——两者正交，不要合成一个「玄学 AUROC¥」。
"""


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    payload = build_payload()
    (OUT / "summary.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=str)
    )
    md = render_cn(payload)
    (OUT / "METHOD_BIZ_SCENARIO_PROTOTYPE.md").write_text(md)
    (DOCS / "METHOD_BIZ_SCENARIO_PROTOTYPE.md").write_text(md)

    # Short EN mirror for toolkit readers
    en = f"""# Method × Biz-scenario × Unit-econ Prototype

## One-liner

{payload['external_one_liner_cn']}

## Relevance chain

1. **Method (OnlineRFPerm / PO / RAG / style AUC)** — opens the regime window and routes actions.
2. **Biz scenarios (hop weak/base/strong × traffic × actions)** — turns the window into ticket/refund/containment volumes.
3. **Unit economics (±20% prices)** — turns volumes into a ¥ band under finance assumptions.

## Same

{chr(10).join('- ' + x for x in payload['same_diff']['same'])}

## Different

{chr(10).join('- ' + x for x in payload['same_diff']['different'])}

## Relevance

{payload['same_diff']['relevance_one_liner']}

See CN full report: `docs/biz/METHOD_BIZ_SCENARIO_PROTOTYPE.md`
"""
    (OUT / "METHOD_BIZ_SCENARIO_PROTOTYPE_EN.md").write_text(en)
    (DOCS / "METHOD_BIZ_SCENARIO_PROTOTYPE_EN.md").write_text(en)

    print(md)
    print(f"wrote {DOCS / 'METHOD_BIZ_SCENARIO_PROTOTYPE.md'}")
    print(f"wrote {OUT / 'summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
