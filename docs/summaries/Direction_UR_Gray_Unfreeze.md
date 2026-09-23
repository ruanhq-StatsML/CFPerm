# Direction → 队列有用率调灰度/话术 → 解冻 ETA/规则

> 现状怎么评价；「先补 Direction、再按队列有用率调灰度或话术」**怎么评估 / 刻画 / 预测**；  
> **解冻 ETA / 规则** 的逻辑具体怎么走。  
> 场景轴仍是 pos刷量 / neg灌入 / flat漂移。

相关：[`Feature_Brush_LLM_SelfEvolution.md`](Feature_Brush_LLM_SelfEvolution.md) · [`Business_Scenarios_Brush_vs_Inject.md`](Business_Scenarios_Brush_vs_Inject.md) · [`Review_Accel_Personas_Feedback_Loop.md`](Review_Accel_Personas_Feedback_Loop.md) · [`Coverage_Confidence_Gray_Rollback.md`](Coverage_Confidence_Gray_Rollback.md)（覆盖/置信/回滚可辩护版：Cover 实测、C0≠C2、期望伤害、L-roll 爆炸半径、升慢降快）

---

## 0. 现在效果怎么样（一屏）

| 层 | 现状 | 读法 |
|---|---|---|
| 工程 | smoke / CI / 出卡 / 灰度 CLI **绿** | 可演示 |
| Direction | 当前 `summary.json` **无 `direction` 键** | 卡只能假 flat → **刷量/灌入族展不开** |
| 真人反馈 | n_human=0（仅 ood_timer smoke） | gate **NOT_READY**；有用率不可信 |
| 商业 KPI | 无人分钟前后差 | Score≈1（仅接通） |
| 下一步主链 | ①补 Direction → ②真人打标 → ③按队列调灰度/话术 → ④才议解冻 | 见下 |

**效果结论：** 管道通了，**极性与进化燃料都还没接上**；谈「调灰度/解冻」之前必须先过 Direction + 真人两关。

---

## 1. 阶段 A：先补 Direction

### 1.1 它解决什么
没有 `sign_Dy` / `tip_signs`，审出卡无法区分 S1/S2/S3，后面所有「按队列有用率调话术」都会在 **假 flat** 上噪声训练。

代码已有：`direction_report.build_direction_dict`；`run_w1w2_mmd_po_localize_fsds.py` 写入路径存在——**当前落地 summary 缺键 = 跑数/落盘缺口，不是理论缺口。**

### 1.2 怎么刻画（指标）
| 符号 | 定义 | 好长相 |
|---|---|---|
| `dir_cover` | 有 `direction.sign_Dy` 的 summary 占比 | → 1.0 |
| `dy_missing_rate` | `Dy is null` 且 y 列本应存在的比例 | ↓ |
| `tip_sign_cover` | tip 中非 `"0"` 符号占比 | 波次有真实移位时应 >0 |
| `family_entropy` | 出卡后 S1/S2/S3 分布熵 | 不应长期塌缩到仅 S3 |
| `attach_latency` | localize 结束→direction 写入 | 同一次 run 内完成 |

### 1.3 怎么评估（过线）
| 门 | 标准 | 不过线则 |
|---|---|---|
| A0 | 新跑 `summary` 必含 `direction` 对象 | 修落盘；禁止当「刷量 demo」 |
| A1 | 演示集至少 1 张 `sign_Dy=pos` + last/share tip（真 S1） | 只讲漂移，不讲刷量已验证 |
| A2 | 卡上 `dy_missing` 显式标出（若 Dy 真算不出） | 禁止静默当 flat 成功 |

### 1.4 怎么预测（下一步值不值得）
| 预测问题 | 判据 | 动作 |
|---|---|---|
| 补 Direction 能否解锁刷量讲述？ | A0∧(A1∨能重跑出 pos) | **P0 立刻做** |
| 只改文案不补 Dy？ | `dir_cover=0` | **预测无效**（kill） |
| 补完就能解冻 ETA？ | 仍无真人 | **否**；只升到「可讲三族」 |

**预测句：** Direction 是 **案由可辨识** 的前置；ROI 预测在此阶段 = 0，只预测「分流可解释性↑」。

---

## 2. 阶段 B：按队列有用率调灰度 / 话术

前置：`human-gate READY`（≥1 真人）且建议 **n≥5** 再动灰度；n=1 只动话术草稿。

### 2.1 刻画（队列切片）
对每条反馈：`label`, `queue_bucket`, `sign_Dy`, `reviewer`, `note`。

\[
\mathrm{UR}(q)=\frac{\#useful(q)}{\#useful(q)+\#not\_useful(q)}
\]

| 符号 | 含义 | 用途 |
|---|---|---|
| \(\mathrm{UR}(q)\) | 队列 \(q\) 有用率 | 调 `queue_sample_rates` / 话术 |
| \(n(q)\) | 该队样本量 | 小样本不调灰度 |
| `wrong_family_rate` | note/改队暗示走错 S1↔S3 等 | 优先改读法词典 |
| `sample_rate(q)` | 当前灰度 | 建议升降的自变量 |
| `ΔUR` 周同比 | 进化是否生效 | 后验 |

脚本：`summarize_review_feedback.py` → by_queue；`suggest_review_card_sample_rate.py` → 全局建议（可扩成 per-queue）。

### 2.2 评估规则（决策表）
| 条件 | 评估结论 | 动作（人工确认） |
|---|---|---|
| \(n(q)<5\) | 证据不足 | **只收集**；可改 1 句话术 A/B，不改灰度 |
| \(\mathrm{UR}(q)\ge 0.7\) 且 \(n\ge5\) | 队列健康 | sample_rate 维持或略升；话术固化进场景表 |
| \(0.4\le \mathrm{UR}<0.7\) | 中间带 | 先改话术/对照表，再观察一周 |
| \(\mathrm{UR}<0.4\) 且 \(n\ge5\) | 队列有害或难读 | **降** `queue_sample_rates[q]`；重写桶文案 |
| \(\mathrm{UR}<0.2\) 且多队差 | 系统级 | 考虑全局 `kill_switch` 或降默认 sample_rate |
| note 高频「应去刷量/漂移」 | 词典错映 | 改 `TIP_BUCKETS` / 场景对照，**先于**改灰度 |

**灰度 vs 话术谁先：**  
- 走错族 / 看不懂 → **先话术**  
- 话术改过仍低 UR → **再降流量**  
禁止：无 UR 证据自动 `cp` suggested flags。

### 2.3 预测（调之前估效果）
对计划「队列 \(q\) 降 sample_rate：\(r\to r'\)」：

| 预测对象 | 粗糙公式/启发 | 验收 |
|---|---|---|
| 曝光次数 | \(\propto r\) | 工单带卡数 ↓ |
| 期望有用命中 | \(\approx r\cdot\mathrm{UR}(q)\cdot N\) | 周报 useful 绝对数 |
| 误导成本 | \(\propto r\cdot(1-\mathrm{UR})\) | not_useful ↓ 则预测对 |
| 人分钟 | 仅当 useful 场景省时有粗数后才估 | 自评分钟或工单时长 |

**预测原则：** 降 \(r\) 先预测「减少有害曝光」，不预测「AUC↑」。

### 2.4 阶段 B 状态机
```text
READY?
  N → 停在催打标（现在）
  Y → 攒 n≥5 by 队列
        → UR 决策表 → 话术草稿 / sample_rate 建议
        → 人确认写入 configs
        → 下一周 ΔUR 后验
        → 稳定后才进入解冻评估（§3）
```

---

## 3. 解冻 ETA / 规则：逻辑具体谈

### 3.1 为什么要闸门
ETA 包、规则增敏包会碰 **另一条现网系统**；在审出卡还是假 flat、无人用时解冻 = Score 造假 + 误伤面扩大。  
闸门服务一个商业问题：**审出主链是否已证明「图谱变动卡有人吃、吃了不害人」？**

### 3.2 解冻层级（不要一次全开）
| 级 | 解冻什么 | 仍禁止 |
|---|---|---|
| **U0** | 无 | 一切旁路 |
| **U1 钩子可读** | 已有：卡上 `eta_soft_hint` 只读字段 | ETA 服务改输出 |
| **U2 ETA 软门控** | ETA 在 `eta_soft_hint=lower_confidence` 时降置信展示 | 改报价引擎主模型 |
| **U3 规则∩K\*** | 规则命中且实体∈支撑时插队 | 规则改写 Drill |
| **U4 名单 TTL** | \(K^\star\) 灰名单短 TTL | 永久黑名单自动化 |

当前政策：**停在 U0/U1**；U2+ 要过下面清单。

### 3.3 解冻清单（必须全满足才 U2）
| # | 门 | 标准 | 刻画 |
|---|---|---|---|
| G1 | 人用 | `human-gate READY` | n_human≥1（建议解冻前 ≥10） |
| G2 | 极性真 | `dir_cover≥0.95` 且演示集含 S1 与 S3 各≥1 | 非永假 flat |
| G3 | 不害人 | 近 2 周 overall \(\mathrm{UR}\ge0.5\) 且无「灾难队」\(\mathrm{UR}<0.3\) 且 \(n\ge5\) | summarize 周报 |
| G4 | 灰度可控 | kill_switch 演练过一次；可按队列关 | flags 演练记录 |
| G5 | 场景共识 | 值班长书面确认 S1/S2/S3 读法 | 场景主文签字/PR |
| G6 | 旁路契约 | ETA/规则 owner 接受「只读 hint / ∩K\*，不改 Drill」 | 任务卡验收人 |
| G7 | 回滚 | 旁路独立 flag，关掉后 ETA/规则回到解冻前 | 演练 |

**U3 规则** 在 G1–G7 之上再加：  
G8：至少 1 条「规则命中 ∩ 卡上支撑」的试点工单有用；  
G9：误伤抽检（非 \(K^\star\) 被插队）低于约定阈值。

### 3.4 解冻预测（开之前估）
| 旁路 | 预测收益 | 预测风险 | 先做的小实验 |
|---|---|---|---|
| U2 ETA | 异常波次日客诉↓ | 误降置信→物流投诉 | 仅 `lower_confidence` 流量 5% 灰 |
| U3 规则 | 真案更早进队 | 误插队增加人审爆量 | 单规则 × 单队列 1 周 |
| U4 名单 | 漏放↓ | 误伤商户 | TTL 极短 + 人工撤销入口 |

若 G3 的 UR 不稳，**预测 U2/U3 ROI 为负** → kill 解冻提案。

### 3.5 解冻状态机
```text
U0/U1（现在）
  │ G1 真人 + G2 Direction
  ▼
审出主链可讲三族、可调队列灰度（§2）
  │ G3–G7
  ▼
U2 ETA 软门控（小流量）
  │ G8–G9
  ▼
U3 规则∩K* →（更后）U4 名单
任一时刻 kill_switch / 旁路 flag → 退回 U1
```

### 3.6 和「调灰度/话术」的衔接
- **话术/队列 UR**：证明审出主链可读 → 给 G3、G5 供数  
- **Direction**：证明案由可分 → 给 G2；ETA 才知道何时 `lower_confidence` 该亮  
- **解冻**：主链稳定后的 **横向复制**，不是替代主链  

顺序固定：

```text
Direction → 真人 UR → 调话术/队列灰度 →（稳定）→ 解冻 ETA → 解冻规则
```

跳步解冻 = 预测失败。

---

## 4. 一页评估看板（建议每周填）

| 周 | dir_cover | n_human | UR_overall | 最差队 q | UR(q) | 灰度动作 | 话术动作 | 解冻级 |
|---|---|---|---|---|---|---|---|---|
| 今 | 0 | 0 | — | — | — | freeze | — | U0/U1 |
| … | | | | | | | | |

今周预测：只做 **A0 Direction 落盘 + 1 真人打标**；不预测解冻。

---

## 5. Takeaway

1. **效果：** 工程通、Direction 空、无人用 → 刷量逻辑与进化都还没开火。  
2. **补 Direction：** 用 dir_cover / 真 S1 样例评估；预测只解锁「可讲三族」。  
3. **按队 UR 调灰度/话术：** \(n\ge5\) 才动灰度；走错族先改话术；用 ΔUR 后验。  
4. **解冻 ETA/规则：** G1–G7 全过才 U2；规则再加 G8–G9；独立 flag 可瞬间退回 U1。  
5. 顺序不可倒：`Direction → UR → 灰度/话术 → ETA → 规则`。
