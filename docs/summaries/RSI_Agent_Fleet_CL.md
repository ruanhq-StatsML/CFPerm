# RSI → 连续学习 · 180 智能体均力赋能

> 场景换轨：不再堆 soft-burn / PO FLOPs 细账。  
> 方法同构：**发散 → justify → 落地 → 反思**，对象改成  
> **连续学习（CL）如何让 ~180 个 agent 平均发力**。  
> 样例车道：约 **10 个反欺诈**、约 **10 个 CUPED**，其余按 specialty 锁仓。

相关：[`Twenty_Agent_Commercial_Landing.md`](./Twenty_Agent_Commercial_Landing.md) ·  
[`Motivate_Agents_Commercial_Value.md`](./Motivate_Agents_Commercial_Value.md) ·  
[`AntiFraud_Eval_Char_Model.md`](./AntiFraud_Eval_Char_Model.md) ·  
代码：`agod/agent_fleet.py`

---

## 0. 为什么换场景

| 旧 RSI 焦点 | 问题 | 新焦点 |
|---|---|---|
| soft IPTW √/∛ · duty FLOPs | 对导师交付偏超纲 | **舰队均力 + CL 反馈** |
| 单模型效率 scorecard | 难直接派 180 个执行体 | **cohort 锁仓 + 商业 Score** |
| 数学 ablation | 激励扎堆刷同一指标 | **反扎堆税 + soft cap** |

一句话：

> **RSI 的「效率」在这里 = 单位 agent 产能是否均匀、是否被反馈拉到正确 specialty，而不是 α 省了多少 FLOPs。**

---

## 1. 舰队结构（默认 N≈180）

| Cohort | n | 主交付（可点） | CL 反馈 |
|---|---:|---|---|
| `antifraud` | 10 | 波次→工单 / \(K^\star\) 名单 / 线索卡 | 人审 **useful rate** |
| `cuped` | 10 | 实验预周期特征、θ、方差比看板 | **Var(Y_cuped)/Var(Y)** |
| `review_accel` | 20 | 审出卡 / 回流 UI | 有用/没用 |
| `localize_fsds` | 20 | Drill / tip+direction | Lead / 命中粗数 |
| `board_stream` | 15 | 相邻窗板 / stream packs | excess 技能（旁路） |
| `eta_gate` | 10 | \(R_t=1\) 时 ETA 行为 | 异常日客诉 |
| `gray_rollback` | 15 | L1/L2 flag + 回滚 | 客诉 / 解封 |
| `data_contract` | 15 | OpenAPI / 审计 | 调用次数 |
| `roi_ops` | 15 | 成本仪表 / 周报机器人 | ROI 粗数 |
| `flex_reserve` | 50 | 浮动运力；**禁止**整队贴到单一 specialty | 再分配源 |

**Soft cap：** 任一 specialty 并发 WIP ≤ 约 12% 舰队（`flex_reserve` 除外）。  
**激励：** Score = \(1\{\mathrm{shipped}\}·(1_{\mathrm{wired}}+1_{\mathrm{used}}+1_{\mathrm{roi}})\) — 与 [`Motivate_Agents_Commercial_Value.md`](./Motivate_Agents_Commercial_Value.md) 同构。

---

## 2. 两个样例车道怎么「赋能」

### 2.1 反欺诈 ×10（不是 180 人一起刷 AUC）

```text
监测 R_t → Drill K* → tip+direction → 工单/名单/规则
                              ↑
                     10 个 antifraud agent 均分：
                     工单适配 / 名单 TTL / 规则插队 / 线索包装 /
                     L1 加监 / L2 灰度 / 假尖端过滤 / 降噪 /
                     误伤抽检 / 回流标注
```

CL：人审点「有用」→ `useful_rate`；lift vs baseline >0 才给 `used`/`roi` 分。  
**禁止：** 10 人同时改同一 Drill 门或只交 hit@k 表。

### 2.2 CUPED ×10（实验方差，不是又一套模型）

\[
Y^{\mathrm{cuped}} = Y - \theta\big(X_{\mathrm{pre}}-\mathbb{E}[X_{\mathrm{pre}}]\big),\quad
\theta=\frac{\mathrm{Cov}(Y,X)}{\mathrm{Var}(X)}
\]

10 车均分：预周期特征契约、θ 估计服务、方差比看板、实验平台接头、  
guardrail（θ 不稳态告警）、文档/培训、回放夹具、灰度开关、周报、审计。

CL：`ratio = Var(Y_cuped)/Var(Y)`；\<1 才算 `used`，\<0.85 才算 `roi`。  
**禁止：** 把 CUPED 车道改成「再训一个预测模型刷 AUC」。

---

## 3. RSI 单轮（舰队版，仍掐表）

| 分 | 动作 | 产出 |
|---:|---|---|
| 0–1 | 抽约束卡（均力 / 反扎堆 / 禁止 math-Eval） | 本轮边界 |
| 1–4 | 发散：哪个 cohort 反馈最差？ | 1 句 idea |
| 4–7 | Justify：Gini/HHI + useful/CUPED ratio | 可检量 |
| 7–9 | 落地：改 mix / soft cap / 任务卡 ≤1 | diff |
| 9–10 | 账本 + scorecard | `results/agod_agent_fleet/` |

**反思规则（代码里 `continuous_learn_tick`）：**

1. fraud lift≥0 且 CUPED ratio≤1 且 WIP Gini 低 → **keep specialty mix**  
2. 否则从 `flex_reserve` **rehome** 到 mean_score 最低的 specialty（不是涌向最吵的车道）

---

## 4. 观测契约（对齐 ToT / board）

| 符号 | 定义 |
|---|---|
| X | cohort 反馈特征（useful、CUPED 预周期 X、WIP） |
| Y | 商业 Score（接通·真用·ROI） |
| intermediate | 舰队 Thought = {mix, gini, rehome} — **≠ Ŷ** |

与边∥流 X 的关系：反欺诈车道仍可吃 `edge∥stream` 特征进线索机；  
舰队层 RSI **编排产能**，不替代特征拼接。

---

## 5. 跑法与读数

```bash
PYTHONPATH=. python3 scripts/run_agent_fleet_scorecard.py \
  --out results/agod_agent_fleet
```

读法：

| 信号 | 期望 |
|---|---|
| even vs pile-on | even 的 WIP Gini 更低、mean Score ≥ pile |
| antifraud useful lift | specialty-10 下应为正（相对稀释的扎堆） |
| CUPED ratio | \<1 才算实验车道赋能成功 |
| rehome | 只从 flex → underserved，禁止反向抽空 CUPED/反欺诈配额 |

---

## 6. Takeaway

> **180 人平均发力 = 锁仓 specialty + 商业 Score + CL 反馈再分配。**  
> 10 个反欺诈救命审分钟；10 个 CUPED 省实验方差；其余车道同样均分。  
> RSI 在这里 polish 的是 **舰队公平与反馈闭环**，不是又一张 PO α 表。
