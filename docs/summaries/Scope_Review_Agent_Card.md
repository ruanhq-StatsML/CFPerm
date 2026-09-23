# Scope 收窄：先走「审出加速包」一条

> 决策：OOD 账本先冻结扩面；**只走一条**可现网变现的支路。  
> 选定：**审出加速包**（审核 agent 吃线索卡 → 降人分钟）。  
> ETA / 规则增敏 / B1–B15 → **Backlog**，本 scope 不做。

---

## 1. 为什么先走这个

| 候选 | 现网依赖 | 我们已有弹药 | 2 周可点？ |
|---|---|---|---|
| **审出加速（选中）** | 审单/工单上下文 | `summary.json` + `direction` + tip | **是**（导出卡即可） |
| 时效稳健 ETA | ETA 服务改输出 | 波次 \(R_t\) 概念有，Tencent 上 harness 未挂 | 偏后 |
| 规则增敏 | 规则引擎 API | 需对方排期 | 偏后 |
| 支付失败 / 退货边… | 新数据源 | 无 | 杀掉本迭代 |

商业 KPI 只盯一个：**件均审核耗时 ↓**（辅：人审「有用」率）。

---

## 2. 本 scope 边界（In / Out）

**In（必须做完才算走通）：**
1. 从 localize `summary.json` 生成 **审核 agent 上下文卡**（JSON + 短中文）  
2. 卡内固定字段：`K*` 摘要、`direction`、`tips`、推荐审出话术桶  
3. CLI 一键：`summary.json` → `review_agent_card.json` / `.md` / `ticket_custom_fields.json`  
4. 任务卡验收句：打开卡文件即能当审单上下文粘贴  

**Out（明确不做）：**
- 改 MMD/PO/cmean / Drill 门  
- ETA 门控、规则插队、名单 API  
- math-Eval、AUC、GT hit 报表  
- 新 OOD 数据源（支付/退货）  

---

## 3. MVP 走法（本周）

```text
已有 run_*_localize → summary.json(+direction)
        ↓
scripts/tencent_gr/export_review_agent_card.py
        ↓
review_agent_card.json  +  review_agent_card.md
        ↓
审核 agent / 人审 粘贴进工单（人工也行）
        ↓
回流：有用/没用 记一笔（下个 scope）
```

成功标准：
- [ ] CLI 对一份真实 `summary.json` 产出卡  
- [ ] 卡含 `sign_Dy` + `tip_signs` + 支撑摘要  
- [ ] README/任务卡写明「禁止当定罪结论」  

---

## 4. Backlog（冻住，不并行）

- 时效稳健包（ETA）  
- 规则增敏包  
- OOD B1–B15（支付失败、退货边、运力图…）  

解冻条件：审出加速卡 **被真人审至少用过 1 次**（有用/没用打标即可）。

---

## 5. Takeaway

点子已经够多；scope = **只做审出加速导出卡**。  
走通这一条再谈 ETA / 规则。
