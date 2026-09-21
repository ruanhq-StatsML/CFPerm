# 反欺诈评估–刻画模型 elaborate

> 在已直观的三块上再立一层 **模型**：  
> **(1) 连续时间特征** · **(2) 连续图谱变动** · **(3) 监测波次**  
> → 用于 **评估（准不准）** 与 **刻画（像哪类、往哪漂）**。  
> 底层尺子仍是 MMD+PO+cmean；本层不替代 Drill 门。

相关：[`Wave_Drill_Clue_Pipeline.md`](Wave_Drill_Clue_Pipeline.md) ·  
[`AntiFraud_Tip_Sign_Actions.md`](AntiFraud_Tip_Sign_Actions.md) ·  
[`Tip_Sign_Algorithm_Elaboration.md`](Tip_Sign_Algorithm_Elaboration.md)

---

## 0. 要做什么模型（两件事）

| | **评估模型** | **刻画模型** |
|---|---|---|
| 问 | 波次 / \(K^\star\) / tip 是否对准欺诈标签？ | 这一波是什么模式、严重度、演变？ |
| 监督 | GT 团伙 / 案件 ID / 规则命中（若有） | 弱标签 + tip 词典 + 无监督分型 |
| 输出 | AUC、hit@k、Lead、波次 P/R | 类型、风险分、方向卡、时间线 |
| 与线索机关系 | 外层打分器 / 校准器 | 外层表征器 / 分桶器 |

线索机：\(R_t\to K^\star\to J^\star+\mathrm{direction}\)。  
反欺诈模型：吃线索机产物当 **特征**，再评估与刻画。

```text
连续时间 FE ──┐
图谱 shift ──┼─→  z_t / z_wave / z_K*  ─→  评估头 + 刻画头
监测波次   ──┘
```

---

## 1. 三路特征（模型输入）

### 1.1 连续时间特征 \(x^{\mathrm{CT}}\)
窗 \([t,t+\Delta t)\) 内 `feature_engineer` 的边/用户/物品统计：  
`u_span_sec`、`i_share_*`、`i_n_covisit_*`、`ui_pop_mismatch`…  
可再加：滑动均值/方差、相对 \(t-h\) 的差分 \(\Delta x\)。

### 1.2 连续图谱变动特征 \(x^{\mathrm{GS}}\)（graph shift）
在 ref vs cur（或 \(t\!-\!1\) vs \(t\)）上：

| 特征 | 定义 |
|---|---|
| \(\|\delta\|\) | 一阶均值漂移能量 |
| \(\mathrm{MMD}^2\) | 全分布差 |
| PO-gap / \(\mathrm{mean}(\hat\tau^2)\) | 风险窗差 |
| \(\mathrm{share}_j=\delta_j^2/\|\delta\|^2\) | 列能量份额 |
| \(\mathrm{sign}(\delta_j)\) 编码 | \(+1/0/-1\)（进刻画，慎进纯评估） |
| \(R(K),\alpha,\mathrm{eff}\) | 支撑尖度 / 质量 |

### 1.3 监测波次特征 \(x^{\mathrm{W}}\)
| 特征 | 定义 |
|---|---|
| \(T_t^{\mathrm{MMD}},T_t^{\mathrm{PO}},T_t^{\mathrm{cmean}}\) | 当步监测标量 |
| \(R_t\)、财富/α-investing 状态 | 拒识状态 |
| 波次年龄 \(t-t_a\)、峰值 \(T\)、累计 reject 数 | 波次形态 |
| Lead vs MSE/其它流 | 早晚 |

**拼接：** 边级 / 实体级用 \(x^{\mathrm{CT}}\|x^{\mathrm{GS}}\)；波次级用 \(x^{\mathrm{W}}\|\) 聚合后的 GS。

---

## 2. 模型结构（建议双头）

### 2.1 评估头 \(f_{\mathrm{eval}}\)（有 GT 时）
目标：预测「该波次 / 该 \(K^\star\) 是否对应欺诈案件」。

\[
\hat p = \sigma\big(f_{\mathrm{eval}}(z_{\mathrm{wave}}, z_{K^\star})\big),\quad
z=\mathrm{MLP}/\mathrm{HGB}([x^{\mathrm{W}};x^{\mathrm{GS}}; \mathrm{pool}(x^{\mathrm{CT}}\mid K^\star)]).
\]

标签来源（由易到难）：
1. 规则命中 ∩ \(K^\star\)（弱）  
2. 人审案件 ID ↔ 实体（中）  
3. 官方 GT item/order 列表（`gt_subset_evaluator`，强）

**评估指标（算法）：**
| 层级 | 指标 |
|---|---|
| 波次 | Precision/Recall of \(R_t\) vs 案件日；Lead |
| 支撑 | hit@k / P@k / R@k of \(K^\star\) vs GT（已有） |
| tip | tip 命中词典准确率；符号与案件类型一致率 |
| 校准 | ECE of \(\hat p\)；AR 与案件密度对照 |

无 GT 时：**不装评估头**，只跑刻画 + 内在一致性（见 §3.2）。

### 2.2 刻画头 \(f_{\mathrm{char}}\)（可无 GT）
目标：把波次+\(K^\star\) 映到 **模式空间**。

\[
c = f_{\mathrm{char}}(z),\quad
c\in\{\text{刷量末跳},\,\text{劣质灌入},\,\text{供给漂移},\,\ldots\}\ \text{或连续 embedding}.
\]

输入强调：
- \(\mathrm{sign}(\Delta\bar y)\)、\(\mathrm{tip\_signs}\)（direction JSON）  
- Top tip 的 \(\delta_j^2\) 份额  
- 波次形状（短尖峰 vs 长高原）

实现选项：
1. **规则词典**：tip 组合 → 类型（可解释，先上）  
2. **浅层分类**：HGB on hand features（有弱标时）  
3. **对比学习**：同波次切片靠近、跨波次推远（无标刻画）

输出卡（进 `summary.json`）：

```json
{
  "fraud_model": {
    "eval": {"p_hat": 0.73, "gt_hit@100": 0.41},
    "char": {
      "type": "末跳份额升高+成功率升",
      "severity": 0.66,
      "direction": {"sign_Dy": "pos", "tip_signs": {"i_share_last": "+"}},
      "wave_shape": "short_spike"
    }
  }
}
```

---

## 3. 进一步评估与刻画：协议

### 3.1 离线评估协议（有 GT）
1. 时间切：train 波次 / test 波次（禁止未来泄漏）  
2. 线索机只在 train 调阈值（\(\varepsilon,\tau,R^\star\)）；test 冻结  
3. \(f_{\mathrm{eval}}\) fit on train 波次标签；报 test AUC、P/R、hit@k  
4. Ablation：`CT only` / `GS only` / `Wave only` / `CT+GS+Wave`  
5. 对照：纯规则、纯全局 FSDS（无 Drill）、随机 \(K\)

### 3.2 无 GT 时的内在刻画评估
- **稳定性**：多种子 \(K^\star\) Jaccard、tip \(\pi\)  
- **可分性**：不同类型词典簇的轮廓系数  
- **方向一致**：同类型案件应同号 \(\mathrm{sign}(\Delta\bar y)\) 比例  
- **波次质量**：假尖端率（\(R(K)\approx 1\) 仍报警的比例）

### 3.3 刻画维度（建议固定报表）
| 维 | 量 |
|---|---|
| 时间 | 波次长度、Lead、是否复发 |
| 结构 | \(\|K^\star\|\)、\(\alpha\)、键=merchant/user |
| 分布 | \(\|\delta\|\)、MMD、PO |
| 方向 | \(\mathrm{sign}(\Delta\bar y)\)、tip_signs |
| 模式 | char type / embedding 近邻历史波次 |

---

## 4. 训练目标（可分可合）

**评估（有标）：**  
\(\mathcal{L}_{\mathrm{eval}}=\mathrm{BCE}(\hat p, y_{\mathrm{fraud}})\)  
或支撑级：listwise / pairwise ranking so that GT entities rank high in `score(v)`（可微调 rank-avg 权重，**默认不改锁定平均**）。

**刻画（弱标/无标）：**  
- 监督：\(\mathcal{L}_{\mathrm{char}}=\mathrm{CE}(c, c^{\mathrm{dict}})\)  
- 无标：InfoNCE on wave slices；或仅聚类 + 人工命名  

**联合（可选）：**  
\(\mathcal{L}=\mathcal{L}_{\mathrm{eval}}+\lambda\mathcal{L}_{\mathrm{char}}\)，共享 \(z\) 编码器；Drill 分数仍外部冻结，避免用标签改 MMD 门。

---

## 5. 与现有代码的落点

| 块 | 现状 | 模型层下一步 |
|---|---|---|
| CT 特征 | `feature_engineer` | 滑窗差分、波次 pool |
| GS | `rank-avg(cmean,MMD,PO)` + FSDS | 导出 \(x^{\mathrm{GS}}\) 向量进矩阵 |
| 波次 | OnlineRFPerm 骨架（Grad 文档） | TencentGR 上 \(\Delta t\) harness 产 \(x^{\mathrm{W}}\) |
| GT 评估 | `gt_subset_evaluator` hit@k | 接 \(f_{\mathrm{eval}}\) 的标签与 AUC |
| direction | `summary.json["direction"]` | 进 \(f_{\mathrm{char}}\) 必选输入 |

建议原型脚本名：`scripts/tencent_gr/run_antifraud_eval_char.py`  
读入线索卡 JSONL（每波次一行）→ 训/评双头 → 写 `fraud_model` 块。

---

## 6. 边界（讲武德）

- 模型 \(\hat p\) = **案件对齐概率 / 风险分**，不是因果作弊概率  
- 刻画类型来自 tip+方向+词典，可错分；需人审校准词典  
- 不把 \(\hat p\) 写回 OnlineRFPerm 的 \(T_t\)（避免循环门控）  
- Bench 准 ≠ 可布控；布控仍要 direction 完整

---

## 7. Takeaway

反欺诈模型 = 在 **CT 特征 + 图谱变动 + 监测波次** 三路上建  
**评估头（对准案件/GT）** + **刻画头（模式/方向/波次形态）**。  
线索机负责发现；模型负责 **打分与分型**——进一步评估与刻画都落在这一层。
