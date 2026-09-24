# RSI 20 分钟 · ML / 统计方法机会（不改核、不锁商业）

> **方式换轨：** 不再堆商业车道 / agent 均力叙事。  
> **底座不动：** `excess_auc` · ECE · probe_eff · 块 bootstrap · reservoir ·  
> PO `rank_eff`/`mse_eff` · soft-burn · ToT —— **只组合、不改写**。  
> **节奏：** 每 **20 分钟**一轮：发散 ML 点子 → justify（H0/效率）→ 落地 ≤1 模块 → 谈效果。

账本：[`RSI_Iteration_Log.md`](./RSI_Iteration_Log.md) ·  
协议父本：[`RSI_TenMin_Model_Stats_Loop.md`](./RSI_TenMin_Model_Stats_Loop.md) ·  
代码：`agod/ml_method_ops.py`

---

## 0. 本轮主张（为什么换方式）

| 旧做法 | 问题 | 新做法 |
|---|---|---|
| soft-burn / α FLOPs 细抠 | 容易超纲、动核 | **旁路 ML 算子**接在 excess 旁 |
| 180-agent 商业均力 | 偏交付、弱方法 | **方法机会账本**（可检统计量） |
| 只报 AUC / 接通 | 不可识别技能 | 继续用你的 **null / excess / ECE** |

一句话：

> **提升 = 在同一 H0 语言下，多三个 ML 诊断：偏相关技能、样本量曲线、技能流变点。**

---

## 1. 20 分钟 Protocol

| 分 | 动作 | 产出 |
|---:|---|---|
| 0–2 | 抽方法约束卡（下表） | 本轮不许重复 |
| 2–7 | 发散 1 个 **ML/统计** 点子（跨 pack） | idea + 假设 |
| 7–12 | Justify：H0 / 可检量 / 与 excess 关系 | 公式或伪码 |
| 12–18 | 落地 ≤1 旁路模块 + 测试；**禁止改** null/PO 核 | green tests |
| 18–20 | 写效果段 + 账本 + commit | log 行 |

**方法约束卡（轮流）：**
1. null / 置换 / 重采样  
2. 校准（ECE / Brier / isotonic）  
3. 相对 FLOPs / 样本效率  
4. 跨非广告 pack  
5. PO / IPTW / rfperm **旁路**稳健（不改核）  
6. 流式变点 / 遗忘  
7. 只 polish 读法（本轮）

---

## 2. 机会 elaborations（侧重 ML）

### M1. Partial excess（混杂偏出）— **Iter 本轮落地**

**问题：** 高 `excess_auc` 常被 volume / intensity 单调特征抬起。  
**方法：** 对探针分数做 \(s \leftarrow s - \hat\beta c\)（\(c\)=混杂），再用**同一** label-permutation null 算 `excess_partial`。

\[
\Delta_{\mathrm{excess}} = \mathrm{excess}_{raw} - \mathrm{excess}_{partial}
\]

| 读法 | 含义 |
|---|---|
| \(\Delta\approx 0\) | 技能不只是混杂线性足迹 |
| \(\Delta\gg 0\) 且 partial≈0 | 「技能」主要是混杂 — 故事降级（仍可作 ops 信号） |
| R²(c→s) 高 | 探针几乎在读 volume |

**效果期望：** 买量基线 pack 上 raw excess 虚高、partial 塌缩 → 与「AUC≠限投」同向，但是 **ML 语言**（partialling-out），不是业务口号。

---

### M2. Excess learning curve（样本效率）— **落地**

固定已训分数，在测试窗嵌套前缀 \(n\cdot f\) 上估 excess。

| 读法 | 含义 |
|---|---|
| 20%→100% excess 几乎不变 | 估计不饥渴，小窗可报 |
| 明显上升 | 稀标签 / 尾部决定技能，别早停宣称 |
| 乱抖 | 先修 null 次数 / 类别平衡再比 pack |

**效果：** 给 `probe_eff` 补一维「要多少 eval mass」——算力账的 **样本轴**，不是 α 轴。

---

### M3. Page–Hinkley on excess stream（技能变点）— **落地**

对相邻对序列 \(\{excess_t\}\) 跑 PH（检测持续下降）。

| 读法 | 含义 |
|---|---|
| 无报警 | 技能水平稳态（仍用块 bootstrap 谈均值方差） |
| 报警在 index \(t^\star\) | **技能本身**漂移 — 触发「换窗 / 重选特征 / 停用探针」的方法门，而非业务工单 |

**效果：** 把 OnlineRFPerm 的「数据漂移门」对偶成「**探针技能漂移门**」。

---

### M4–M10. Backlog（下一轮 20min 可抽）

| ID | ML 点子 | Justify | 与现栈关系 |
|---|---|---|---|
| M4 | Isotonic / Platt 后校准再报 ECE | 校准可修 vs 不可修 | 旁路 ECE，不改 enrich |
| M5 | Spearman(hardness, score) vs AUROC 分歧 | 排序族不一致 → 硬样本不稳 | 接 `hard_rank_metrics` |
| M6 | DiffusionDB token excess vs Tencent 对照 | 弱特征族 null 压力测试 | R4 backlog |
| M7 | Cross-fit PO（已有 `obs_po_cv`）效率表 | 样本复用 vs 泄漏 | 只做 scorecard |
| M8 | Conformal set size @ excess 分位 | 高 excess 是否更窄 | 新算子 |
| M9 | Forgetting-factor IPTW（λ-折扣） | 连续学习权重 | 旁路 soft_burn |
| M10 | Freeze Pareto 的 bootstrap 带 | 条件 Pareto 不确定性 | 旁路 freeze_eff |

---

## 3. 本轮效果（合成 + 可复现）

跑：

```bash
PYTHONPATH=. python3 scripts/run_ml_method_scorecard.py \
  --out results/agod_ml_method_ops
```

读数（seed=0 合成）：

| 诊断 | 结果 | 提升含义 |
|---|---|---|
| Partial excess | volume 混杂下 \(\Delta_{\mathrm{excess}}\) 大、partial↓ | 识别「假技能」 |
| Learning curve | 前缀稳定 / 上升可分 | 决定评多少行再立论 |
| Page–Hinkley | 人为注入 skill drop → alarm | 技能流可监 |

**不宣称：** 改写了 null；省了 PO FLOPs；商业接通。

---

## 4. Takeaway

> 20min RSI = 在你的 **excess / ECE / FLOPs** 方法论上，加 **偏相关 · 学习曲线 · 变点** 三个 ML 旁路。  
> 效果可谈、可测、不超核。下一轮从 M4–M10 抽一张约束卡继续。
