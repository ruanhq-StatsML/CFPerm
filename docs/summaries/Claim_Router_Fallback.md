# Claim router：防 overclaim 的全局检查与 fallback

> 问题：scorecard / RSI 文案容易 **overclaim**（skill、drivers、win、promote、ROI）。  
> 办法：所有 interpretive reading 先过 **`agod/claim_router.py`** ——  
> **分级检查 → 不够证据就降级 → L4（因果/上线/ROI）永远 fallback。**

相关：[`Coverage_Confidence_Gray_Rollback.md`](./Coverage_Confidence_Gray_Rollback.md)（诚实分层）·  
[`Sample_Chunk_Board_Logic.md`](./Sample_Chunk_Board_Logic.md)

---

## 1. 等级（只升不造）

| Level | 能说什么 | 解锁条件（例） |
|---|---|---|
| **L0** descriptive | 数字本身 | 总有 |
| **L1** diagnostic | 比较 / 曲线形状 / Δexcess | 有限指标 |
| **L2** evidential | excess / null / burn 码 | excess≥0.05 且（无 CI 或 CI_lo>0）；若有 partial 则 partial≥0.03 |
| **L3** action hint | SOFTEN / KEEP / 可调查 | `allow_action_hint` + 合法 burn_decision |
| **L4** forbidden | 因果驱动、promote 生产、ROI/省人分钟 | **永不**从诊断板放出 |

---

## 2. 路由

```text
proposed text
    │
    ├─ scan_overclaim_text  → flags (skill/beats/stable/alarm/causal/…)
    ├─ max_supported_level(evidence)
    │
    ├─ need ≤ supported 且非 L4 → 原文放行
    └─ else → fallback_text（降级话术）
```

API：

- `route_claim(proposed, Evidence(...))` → `{ok, level, text, flags, fallback_applied}`
- `harden_reading(...)` → 只要安全 `text`
- `route_pack_interpretation(...)` → 相邻窗板 reading（去掉 “persistent drivers”）
- `route_ml_suite_claims(suite)` → ML ops 旁路 headline

---

## 3. 已挂载出口

| 出口 | 行为 |
|---|---|
| `interpret_pack`（adjacent board） | 经 `route_pack_interpretation`；correlates ≠ drivers |
| `ml_method_ops` readings | curve / PH / suite `claims` 经 router |
| `summarize_burn_decisions` headline | L3 hint + evidence；禁 ship/ROI 语 |

---

## 4. 原则

> **能写进 JSON 的主张 ≤ evidence 解锁的最高级。**  
> Overclaim 不是文风问题，是路由失败——fallback 比漂亮 headline 优先。
