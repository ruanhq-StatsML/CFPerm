# Scope：审出加速包（业务场景 + 边界）

> **一句话**：图谱分布变动波次 → Drill 出 \(K^\star\) → tip+正负向线索卡 → **贴进审核工单**，降低件均人分钟。  
> 不是定罪模型，不是自动封禁；是 **L1 盯梢级** 的可读分流上下文。  
> 底座算法冻结：Drill = rank-avg(cmean, MMD, PO)；本 scope **不改门**。

相关：[`Wave_Drill_Clue_Pipeline.md`](Wave_Drill_Clue_Pipeline.md) · [`AntiFraud_Tip_Sign_Actions.md`](AntiFraud_Tip_Sign_Actions.md) · [`TenMin_OOD_Iteration_Opportunities.md`](TenMin_OOD_Iteration_Opportunities.md) · [`Business_Scenarios_Brush_vs_Inject.md`](Business_Scenarios_Brush_vs_Inject.md)（**业务场景主文：刷量 vs 灌入 vs 漂移**） · [`Feature_Brush_LLM_SelfEvolution.md`](Feature_Brush_LLM_SelfEvolution.md) · [`Direction_UR_Gray_Unfreeze.md`](Direction_UR_Gray_Unfreeze.md)（Direction→UR→灰度→解冻ETA/规则） · [`Review_Accel_Personas_Feedback_Loop.md`](Review_Accel_Personas_Feedback_Loop.md)

---

## 0. 业务舞台（先看清再定 scope）

### 0.1 数据与实体（TencentGR 子集）
| 层 | 是什么 | 在卡上怎么出现 |
|---|---|---|
| 交互边 | \((u,i,t)\) 曝光/点击流 | 支撑边数 `n_localized_edges` |
| Item / 商户侧 | 路径份额、归因、共现邻域 | tip：`i_share_*` / `i_credit_*` / `i_n_covisit_*` |
| User 侧 | 活跃跨度、事件规模 | tip：`u_span_sec` 等 |
| 波次 | W1↔W2（或连续 \(R_t=1\)） | `gap_days` → SLA；`localize_k` → 支撑规模 |
| 结局 | 块上 click/convert 均值差 | `sign_Dy` ∈ {pos, neg, flat} |

当前样例跑（`tencent_gr_w1w2_mmd_po_fsds`）：`localize_k=200`，`gap_days=30`，顶 tip 以 **路径线性份额/归因 + 共现邻域** 为主；`sign_Dy` 缺省按 **flat** 读（慎升强动作）。

### 0.2 角色与痛点
| 角色 | 今天怎么干活 | 痛点 | 本包给什么 |
|---|---|---|---|
| 监测/算法 | 出 `summary.json` | 数字堆在结果目录，审单看不见 | 一键导出粘贴卡 |
| 审核同学 | 人工翻图、翻规则、翻历史工单 | 件均耗时长、不知道先看哪类案 | **队列桶 + 读法 + tip 话术** |
| 审核 Agent | 吃工单上下文 | 缺结构化图谱线索 | `paste_for_agent` + `ticket_custom_fields` |
| 运营/风控 | 定 L1/L2 灰度 | 怕误伤、难回滚 | kill_switch / 按队列 sample_rate |
| 调查（L3） | 开案 | 本包 **不直接开案** | 只提供时间线与特征卡（Out） |

角色扩表、痛点→下一步、真人反馈闭环详见 [`Review_Accel_Personas_Feedback_Loop.md`](Review_Accel_Personas_Feedback_Loop.md)（R1–R11 + §3 闭环状态机）。

### 0.3 动作级（与全站一致，本包默认停在 L1）
| 级 | 含义 | 本 scope |
|---|---|---|
| L0 | 观察 / 灰度未命中 | `GRAY_DISABLED` 时强制 |
| **L1** | 盯梢、加监、提审单优先级 | **默认交付** |
| L2 | 限频/验证码等可回滚收紧 | Out（解冻后另包） |
| L3 | 开案冻结 | Out |

---

## 1. 为什么只走「审出加速」

| 候选包 | 现网依赖 | 我们已有 | 2 周可点？ | 本迭代 |
|---|---|---|---|---|
| **审出加速（选中）** | 工单/审核上下文可粘贴 | `summary`+tip+导出 CLI 全家桶 | **是** | **In** |
| 时效稳健 ETA | ETA 服务改输出 | 仅卡上 `eta_soft_hint` 只读钩子 | 偏后 | Backlog |
| 规则增敏 | 规则引擎 API | 概念有，无对接 | 偏后 | Backlog |
| 资损/支付失败边 | 新数据源 | 无 | 否 | Kill |

**商业 KPI（唯一主指标）**：件均审核耗时 ↓  
**辅指标**：人审「有用」率（`useful` / (`useful`+`not_useful`)）按 `queue_bucket` 切片。

---

## 2. 业务场景（先捋清；详文另册）

> **主文**：[`Business_Scenarios_Brush_vs_Inject.md`](Business_Scenarios_Brush_vs_Inject.md)  
> 读卡顺序：`sign_Dy`（极性）→ 主 tip 桶 → 人可改队。其它（角色/闭环）都挂这张图下。

### 2.1 三族（案由轴）
| 族 | `sign_Dy` | 一句话 | 人审先看 |
|---|---|---|---|
| **S1 刷量族** | **pos** | 块上成功变多：互点/养号/末跳，或真爆款 | 成功是否压在末跳；别当爆款误伤 |
| **S2 灌入族** | **neg** | 块上成功变少：劣质短会话 / 劫持残留 | 差流是否灌进 \(K^\star\)；别甩锅商户 |
| **S3 漂移族** | **flat** | 结构漂了、成功没动（**当前样例**） | 活动/推荐/库存；**慎升强动作** |

刷量 ≠ 灌入：极性相反，查法相反；**不能「有 tip 就叫刷量」**。

### 2.2 修饰（不单开案由）
| ID | 作用 |
|---|---|
| S4 共点 | covisit 顶 → 抽检↑，不对外认定团伙 |
| S5 错配 | 帮刷量族拆真爆款 vs 刷热 |
| S6 内容复用 | 换词典，三族轴不变 |
| S7 SLA | 同队先清老波次 |

### 2.3 字段速查
| 场景 | 关键字段 | 默认级 |
|---|---|---|
| S1 | `sign_Dy=pos` + last/share 或短刷 | L1 |
| S2 | `sign_Dy=neg` + 短 span / 差流 tip | L1 |
| S3 | `sign_Dy=flat` + linear/规模 | L1 慎升 |
| S4–S7 | covisit / mismatch / overlay / gap_days | 修饰或调度 |

---

## 3. Scope 边界（In / Out / 解冻）

### 3.1 In（本包必须能演示）
1. `summary.json` → `review_agent_card.{json,md}` + `ticket_custom_fields.json`  
2. 卡含：支撑摘要、`direction`（缺省 flat）、tips、队列桶、读法、L0/L1、免责声明  
3. 工单字段可 dry-run POST；灰度 `kill_switch` + 按队列 `sample_rate`  
4. 有用/没用回流 + 按队列 useful_rate + sample_rate **建议**（人工 cp，禁止自动覆盖）  
5. 波次批处理出卡、波次差分、e2e smoke、内容 overlay（已有）  
6. human-use gate：真人打标前 **feature freeze**

### 3.2 Out（明确不做）
- 改 MMD / PO / cmean / Drill 门公式  
- 自动 L2/L3、自动封禁、对外「团伙认定」通报  
- ETA 服务改输出、规则引擎插队、名单 API 写库（卡上 `eta_soft_hint` 只读钩子除外）  
- math-Eval / AUC ablation / GT hit 当商业验收  
- 新数据源（支付失败、退货边、运力图）— 解冻后另开包  
- 审出卡 Web UI、tip 排序小模型、自动覆盖现网 flags  

### 3.3 解冻条件（唯一）
```text
check_review_accel_human_gate.py → READY
= 至少 1 条 label∈{useful,not_useful} 且 reviewer 非 smoke/ood_timer/agent/ci/bot/test
```
READY 后才议：ETA 包、规则增敏、B1–B15、新行业 overlay。

---

## 4. 交付链路（现网可点路径）

```text
波次 R_t=1 / 两窗 localize
        ↓
summary.json (+ direction 理想必有)
        ↓
export_review_agent_card.py  [--tip-overlay] [--flags]
        ↓
卡 JSON/MD + ticket_custom_fields
        ↓
人审 / 审核 Agent 粘贴工单  （L1）
        ↓
log_review_card_feedback.py --reviewer <真名>
        ↓
summarize → useful_rate → suggest sample_rate（人工确认）
        ↓
human-gate READY？ → 解冻其它包 : 继续只催打标
```

一键演示：`smoke_review_accel_pack.py`  
真人交接：`artifacts/review_accel_human_handoff.md`  
等待看板：`artifacts/review_accel_waiting_board.md`

---

## 5. 验收（怎么算「效果」）

| 层 | 标准 | 现状 |
|---|---|---|
| 工程 | smoke 五步绿；CI 绿；卡可读非定罪 | **已满足** |
| 接通 | 工单字段 dry-run；灰度可回滚 | **已满足** |
| 业务 | ≥1 真人 useful/not_useful | **NOT_READY** |
| KPI | 件均耗时有前后对比 或 有用率按队列可汇报 | 待真人样本 |

**没有真人打标，不得宣称「人分钟下降」已验证。**

---

## 6. 决策摘要

| 项 | 决定 |
|---|---|
| 走哪条 | 只走审出加速包（S1–S7 场景词典 + L1 卡） |
| 主 KPI | 件均审核耗时 ↓（辅：有用率） |
| 默认动作 | L1_watch；flat → 慎升强动作 |
| 冻结 | ETA / 规则 / B1–B15 / 新 feature（至 READY） |
| 节奏 | 10min OOD 已停；6h 只巡检 human-gate |

Takeaway：业务场景已经够具体；**scope = 把 S1–S7 的读法稳稳送进审单**，而不是再扩算法面。下一步是 **1 位真人打标**，不是 B35。
