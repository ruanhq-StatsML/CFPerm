# Scope：审出加速包（业务场景 + 边界）

> **一句话**：图谱分布变动波次 → Drill 出 \(K^\star\) → tip+正负向线索卡 → **贴进审核工单**，降低件均人分钟。  
> 不是定罪模型，不是自动封禁；是 **L1 盯梢级** 的可读分流上下文。  
> 底座算法冻结：Drill = rank-avg(cmean, MMD, PO)；本 scope **不改门**。

相关：[`Wave_Drill_Clue_Pipeline.md`](Wave_Drill_Clue_Pipeline.md) · [`AntiFraud_Tip_Sign_Actions.md`](AntiFraud_Tip_Sign_Actions.md) · [`TenMin_OOD_Iteration_Opportunities.md`](TenMin_OOD_Iteration_Opportunities.md)

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

## 2. 具体业务场景（审核怎么用这张卡）

每张卡固定三件套：**`sign_Dy` 读法** + **主 tip 队列桶** + **建议动作级（默认 L1）**。  
下列场景是运营词典，不是自动分类器；人审可改队。

### S1 · 刷量 / 互点 / 末跳操控（偏 pos）
| | |
|---|---|
| 触发读法 | `sign_Dy=pos` + tip 含 `i_share_last` / `i_credit_last`（+）或短跨度刷 |
| 建议队列 | 末次份额/归因偏高（末跳嫌疑）；或活跃跨度异常（短刷） |
| 审核动作 | 优先看成功路径是否压在末跳商品；抽检互点/养号；**不**直接封 |
| 禁区 | 无 tip 符号、仅 Top ID 名单 → 不得升 L2 |

### S2 · 劣质灌入 / 劫持残留（偏 neg）
| | |
|---|---|
| 触发读法 | `sign_Dy=neg` + tip 活跃跨度变短或共现/触达异常 |
| 建议队列 | 劣质灌入 / 劫持残留优先 |
| 审核动作 | 看失败/低转化会话是否灌进 \(K^\star\)；对照供给是否被打压后残留 |
| 禁区 | 把「成功率掉」直接当商户作恶 |

### S3 · 供给 / 分布漂移（flat，当前样例）
| | |
|---|---|
| 触发读法 | `sign_Dy=flat`（或 Dy 缺失）且 tip 多为路径份额/规模类 |
| 建议队列 | 路径线性份额/归因变动；触达/曝光规模；共现邻域 |
| 审核动作 | **先当供给或流量结构漂了**：核对活动、库存、推荐策略；慎升强动作 |
| 卡点 | 当前样例即此状——产品上要能「读得懂的平」，而不是空白 |

### S4 · 团伙共点（结构）
| | |
|---|---|
| 触发读法 | tip 顶 `i_n_covisit_neighbors` / `i_log1p_n_covisit` |
| 建议队列 | 共现邻域变动（团伙共点） |
| 审核动作 | 看 \(K^\star\) 内共现是否突然变密；与 S1/S2 叠加时提高抽检，仍停 L1 |
| 禁区 | 把 \(K^\star\) 对外通报为「已认定欺诈团伙」 |

### S5 · 热度-活跃错配
| | |
|---|---|
| 触发读法 | tip `ui_pop_mismatch` |
| 建议队列 | 热度-活跃错配 |
| 审核动作 | 区分「真爆款」vs「刷热度」；需结合 `sign_Dy` |

### S6 · 内容审核复用（同一 UI，换词典）
| | |
|---|---|
| 触发 | `--tip-overlay configs/review_tip_bucket_overlay_content.json` |
| 建议队列 | 互粉同发 / 种草入口 / 导流末跳 等话术 |
| 审核动作 | 内容安审队列复用同一卡组件；**不**新开算法门 |
| 边界 | 解冻前不扩金融/运力等新 overlay 包 |

### S7 · 老波次清队（SLA）
| | |
|---|---|
| 触发 | `gap_days` 大 → `sla_urgency=urgent|tight` |
| 审核动作 | 同队列内优先清老波次；灰度关闭则 SLA=`deferred` 不催审 |

**场景 × 字段速查**

| 场景 | 关键字段 | 默认级 |
|---|---|---|
| S1 刷量末跳 | `sign_Dy=pos`, tip last/share | L1 |
| S2 劣质灌入 | `sign_Dy=neg`, tip span/灌入 | L1 |
| S3 供给漂移 | `sign_Dy=flat`, 路径/规模 tip | L1（慎升） |
| S4 共点 | covisit tip | L1 |
| S5 错配 | `ui_pop_mismatch` | L1 |
| S6 内容复用 | tip overlay industry | L1 |
| S7 清队 | `graph_shift_sla_level` | 调度 |

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
