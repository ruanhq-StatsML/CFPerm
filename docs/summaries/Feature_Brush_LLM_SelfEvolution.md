# 选出的特征 · 刷量逻辑 · LLM 出点子 · 评估/预测 · Self-evolution

> 回答四件事：  
> 1）「选特征 + 刷量逻辑」完成了吗、缺什么；  
> 2）LLM/agent **怎么自己想出新点子**；  
> 3）**下一步如何评估与预测**；  
> 4）系统如何 **self-evolution**（具体闭环，不是口号）。  
> 场景轴仍以 [`Business_Scenarios_Brush_vs_Inject.md`](Business_Scenarios_Brush_vs_Inject.md) 为准。

---

## 0. 完成度（先给结论）

| 块 | 状态 | 说明 |
|---|---|---|
| 行：Drill \(K^\star\) | **算法完成** | rank-avg(cmean,MMD,PO)；门冻结 |
| 列：FSDS 选 tip \(J^\star\) | **算法完成** | 仅在 \(K^\star\) 上 Scaler→Var→SelectKBest→HGB/LR |
| 刷量/灌入/漂移读法 | **词典完成** | pos/neg/flat → S1/S2/S3 |
| 审出卡导出 + 灰度 + 回流 CLI | **工程完成** | smoke 五步绿 |
| `direction` 写入样例 summary | **有缺口** | 当前样例可无 direction → 全按 flat，刷量逻辑**展不开** |
| 真人 useful/not_useful | **未完成** | gate NOT_READY；无人进化信号 |
| 自动改 tip 词典 / 自动升 L2 | **刻意不做** | 进化只经人工确认 |

**一句话：** 选特征与刷量**算法+词典**齐了；**业务验证与 self-evolution 燃料（真人反馈 + 可靠 Dy）还缺。**

---

## 1. 选出的特征：算法在干什么（固定流水线）

### 1.1 硬顺序（不许反）
```text
(可选) 波次 R_t=1
    → Drill 定 K*          【行：谁在共变】
    → 仅在 K* 边上 FSDS    【列：哪些特征 tip】
    → cmean 给 tip 符号     【升/降】
    → Δȳ → sign_Dy          【刷量/灌入/漂移极性】
```

列尖 **从不**改行门；行尖 **推迟** FSDS 到停钻之后。

### 1.2 FSDS 选 tip（代码锚点）
`run_w1w2_mmd_po_localize_fsds.py`：

```text
StandardScaler → VarianceThreshold → SelectKBest(f_classif, k)
  → HistGradientBoosting / LogReg
  → ranking → top_features = J*
```

- 监督标签：\(K^\star\) 内 click/convert（稀时可用代理）  
- 输出进 `summary.json` 的 `fsds_* .top_features`  
- 审出卡：`export_review_agent_card` 读 top_features → tip 列表 → `TIP_BUCKETS` 话术桶  

### 1.3 符号从哪来
| 符号 | 定义 | 业务用处 |
|---|---|---|
| \(\mathrm{sign}(\delta_j)\) | tip 在 \(K^\star\) 上均值 cur−ref | tip「升高/降低」 |
| \(\mathrm{sign}(\Delta\bar y)=\) `sign_Dy` | 块上成功率 cur−ref | **S1/S2/S3 族** |
| tip 能量 \(\delta_j^2\) | 排序用 | 谁上卡，不单独定案由 |

`direction_report.build_direction_dict` 应写入 `summary.direction`；**缺则卡只能 flat 读。**

### 1.4 当前样例实际选出的特征
顶 tip（W1/W2 一致倾向）：

`i_share_linear`, `i_credit_linear`, `i_share_first`, `i_credit_first`,  
`i_log1p_n_users`, `i_log1p_n_covisit`, `i_n_covisit_neighbors`, …

→ 路径份额/归因 + 规模/共现；**没有**可靠 `i_share_last=+` 且 `sign_Dy=pos` 时，**不能叫刷量案**。

---

## 2. 刷量逻辑：从特征到案由（怎么判）

### 2.1 判定式（给人审 / Agent 的规则，非自动封禁）
```text
若 sign_Dy = pos:
    主案由候选 = 刷量族(S1)
    若 tip 含 i_share_last / i_credit_last 且符号 +:
        子队 = 末跳操控
    elif tip 含 u_span_sec 且偏短刷:
        子队 = 短窗狂点
    else:
        子队 = 弱刷量或真爆款 → 必须人工辨爆款
若 sign_Dy = neg → 灌入族(S2) …
若 sign_Dy = flat / 缺失 → 漂移族(S3)，禁止自称刷量
```

### 2.2 「刷量」在算法里到底证明了什么
| 证明了 | **没**证明 |
|---|---|
| \(K^\star\) 上成功更多 + 某些图特征同向漂 | 法律意义上的刷单团伙 |
| 该进哪条 L1 审核抽屉 | 该封号 / 该 L2 |

### 2.3 未完成项（刷量逻辑要「能用」还差）
1. 跑 localize 后 **强制 attach direction**（否则永假 flat）  
2. 至少 1 张 **真 pos + last tip** 样例卡给人审对照 S1  
3. 真人打标：S1 卡 useful？走错成 S3？→ 回写词典  

---

## 3. LLM 如何自己 come-up-with 新 idea

这里的「LLM」= 本仓库里的 **OOD / 审出 agent 流程**（不是改 Drill 的模型）。

### 3.1 出点子的约束机（强制结构，防胡想）
每轮（曾用 10min timer；NOT_READY 时改为 6h 只巡检）强制：

| 步 | 规则 | 产出 |
|---|---|---|
| 1 | 抽 1 张约束卡（只谈接口 / 跨系统 / 人分钟 / 可回滚 / 复用 summary…） | 边界 |
| 2 | 写 **1 个**现栈没有的点子（禁 AUC/公式证明） | idea 一句 |
| 3 | 串到波次·Drill·线索·审核·ETA·名单 **≥2 节点** | links |
| 4 | 商业过滤：2 周内现网可点？ | keep / kill |
| 5 | 写 `ood_opportunity_log.jsonl` | 账本 |
| 6 | keep 才改表/代码；且 **scope_gate=审出加速** | 防扩面 |
| 7 | human-gate NOT_READY → **kill 新 feature**，只催打标 | 防空转 |

点子不是自由生成，是 **约束卡 × 骨架边 × 商业过滤 × scope 闸门** 的搜索。

### 3.2 点子从哪「长」出来（启发式，可复述）
LLM 被要求在下列空位填边，而不是空想模型：

```text
[监测波次] ──?──► [Drill K*] ──?──► [线索卡]
                      │                  │
                      ▼                  ▼
                   名单/规则          审核agent / 工单
                                         │
                                         ▼
                                    有用回流 → 灰度
```

常见 keep 模式（已落地 B16–B34）：  
缺工单字段 → 补 ticket fields；缺回流 → feedback CLI；缺回滚 → kill_switch；  
缺批量 → batch export；缺双 agent 钩子 → eta_soft_hint（只读）；  
缺真人 → handoff / gate / **停 10min 堆 feature**。

### 3.3 什么叫「好 idea」（预测它值不值得做）
在动手前用一张卡预测（写进 log 的 keep 字段）：

| 问题 | 过线 |
|---|---|
| 2 周内谁能点？ | 有工单/CLI/开关 |
| 串 ≥2 节点？ | 是 |
| 改 Drill 门？ | 否 → 否则 kill |
| 无真人时还加 feature？ | gate NOT_READY → kill |
| 服务人分钟或客诉？ | 至少一条 |
| 可灰度回滚？ | 优先 keep |

**预测分数（与 Motivate 一致）：**  
\(\mathrm{Score}=1_{\text{可点}}(1_{\text{接通}}+1_{\text{有人用}}+1_{\text{ROI粗数}})\)  
出点子阶段只能估「接通」；「有人用/ROI」必须靠真人实验后验。

---

## 4. 下一步如何评估与预测（具体）

### 4.1 三层评估（别混）
| 层 | 问什么 | 指标 | 谁做 |
|---|---|---|---|
| **A 算法层** | tip/K* 稳不稳 | 残差 \(R(K)\)、tip 稳定性跨波次 | 离线脚本 |
| **B 产品层** | 卡能不能分流 | 队列走错率、useful_rate by bucket | 真人打标 |
| **C 商业层** | 人分钟降没降 | 件均耗时前后差、有用率 | 值班长周报 |

刷量逻辑是否「对」，**B/C 说了算**，不是 AUC。

### 4.2 预测「下一波该干嘛」（决策树）
```text
human-gate?
  NO  → 唯一预测：催 1 名真人打标；kill 新 feature
  YES → direction 样例是否常空?
          YES → 预测优先：强制 attach direction（否则 S1 刷量展不开）
          NO  → useful_rate by 队列
                  低有用队列 → 降 queue_sample_rates 或改 TIP_BUCKETS 话术
                  高有用 + 要客诉 → 才预测解冻 ETA 钩子真接通
                  高有用 + 要降误伤 → 才预测规则∩K* 包
```

### 4.3 对「刷量检出」的预测实验（最小）
| 假设 | 设计 | 成功标准 |
|---|---|---|
| pos+last tip 的卡能让审核少走错队 | 造/选出 1 张真 S1 卡 vs 1 张 S3 卡，盲贴 | useful↑ 且 note 少「走错族」 |
| flat 卡自称刷量会害人 | 故意错标读法做对照 | not_useful↑ → 验证词典必要 |

---

## 5. Self-evolution：具体怎么进化（不是自动变聪明）

### 5.1 进化的是什么 / 不是什么
| 会进化（人工确认后） | 不会自动进化 |
|---|---|
| tip → 话术桶 `TIP_BUCKETS` / overlay | Drill 门公式 |
| 队列 `queue_sample_rates`、kill_switch | 自动 L2/L3 |
| 场景对照表（哪组 tip+Dy → 哪族） | 无反馈的「新行业包」乱扩 |
| OOD 账本 keep/kill 记忆 | 用 AUC 改 \(K^\star\) |

### 5.2 四条可跑的进化环
**环 ① 审出词典进化（主）**
```text
出卡 → 人审 → useful/not_useful + note
  → summarize by queue_bucket
  → 低有用：改话术 / 改队映射 / 降 sample_rate
  → 高有用：固化进 Business_Scenarios 对照表
  → 下一波出卡用新词典
```

**环 ② 灰度进化**
```text
suggest_review_card_sample_rate(useful_rate)
  → flags.suggested.json
  → 人 cp 到 configs（禁止自动）
  → 再出卡观察有用率是否回升
```

**环 ③ 点子进化（OOD）**
```text
约束卡 → idea → keep/kill 写入 jsonl
  → keep 才落地脚本
  → gate NOT_READY 时进化目标切换为「催打标」而非「新功能」
  → READY 后恢复 10min OOD，优先解冻 ETA/规则边
```

**环 ④ 特征侧弱进化（算法同学，仍人工）**
```text
跨波次 tip 集合 Jaccard 低 → 查 FE 漂移或 K* 漂
  → 不自动改 SelectKBest 的 k
  → 若人审总说「线性份额看不懂」→ 优先改桶文案，不是加模型
```

### 5.3 Self-evolution 状态机（落地版）
```text
[冷启动] 词典手工 S1–S3
    │ 工程 smoke 绿
    ▼
[可演示] Score≈1
    │ 缺真人 → freeze feature；6h 只巡检 gate
    ▼
[有人用] n_human≥1 → READY；Score→2
    │ 攒 ≥5 条 → 调队列灰度；改话术
    ▼
[可卖粗数] 有用率/耗时粗数 → Score→3
    │ 才允许开 ETA/规则 旁路进化
    ▼
[旁路进化] 仍禁止改 Drill 门
```

### 5.4 LLM 在进化里的角色边界
| LLM/agent 可做 | 必须人做 |
|---|---|
| 在约束下提案、写 CLI、改字典草稿、跑 gate/smoke | 打 useful；cp 灰度；解冻旁路；认定欺诈 |
| 根据 not_useful note 提议「改哪句读法」 | 接受或拒绝提议 |
| 预测下一优先债（direction / 真人 / 某队列） | 拍板资源 |

**Self-evolution = 反馈数据驱动词典与灰度；LLM 是提议者与执行工，不是自律法官。**

---

## 6. 你现在该盯的「完成清单」

1. □ 样例 `summary` **带 direction**（刷量逻辑可演示）  
2. □ 1 张真 **S1（pos+末跳 tip）** 对照卡  
3. □ 1 名真人 `--reviewer` 打标 → gate READY  
4. □ ≥5 条后首张队列有用率周报 → 第一次词典/灰度进化  
5. □ 才讨论 ETA/规则 self-evolution 旁路  

1–2 是算法展示债；3–4 是 self-evolution 点火；没有 3，谈进化都是空转。

---

## 7. Takeaway

1. **选特征** = \(K^\star\) 上 FSDS → \(J^\star\)；**刷量** = `sign_Dy=pos` + tip 抽屉，不是 tip 本身。  
2. 当前样例特征偏路径/规模 + Dy 空 → **漂移讲述**，刷量链路文档齐、样例未点亮。  
3. LLM 出点子 = 约束卡 × 骨架边 × 2 周可点 × scope/gate 闸门。  
4. 评估分 A/B/C 三层；下一步预测跟 gate/Dy/useful_rate 决策树走。  
5. Self-evolution 四环：词典、灰度、OOD 账本、弱特征体检——**全部人工确认**，Drill 永不自进化。
