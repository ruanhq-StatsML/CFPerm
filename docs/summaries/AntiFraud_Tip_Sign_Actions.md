# Tip 特征 / 正负向 · 反欺诈排查布控逻辑

> 两套逻辑不要混：  
> **(A) Benchmark 准确率** = 离线 hit@k / AUC 等，评「选得准不准」；  
> **(B) 排查/布控线索** = 线上交付 \(K^\star\) + tip + 正负向，驱动 **人/系统动作**。  
> 本文谈 (B)：actions 有哪些、tip 与符号怎么读、反欺诈场景怎么用。

相关：[`Graph_Shift_Method_Elaboration.md`](Graph_Shift_Method_Elaboration.md) ·
算法层（无行动项）：[`Tip_Sign_Algorithm_Elaboration.md`](Tip_Sign_Algorithm_Elaboration.md)

---

## 0. 先分清：准确率 ≠ 布控逻辑

| | Benchmark | 排查/布控线索 |
|---|---|---|
| 问什么 | 算法是否找回 GT 子集 / 预测得好不好 | 这一波 shift 要不要看、看谁、往哪查 |
| 指标 | hit@k、P/R、AUC、Lead | \(R_t\)、\(K^\star\)、tip 卡、\(\mathrm{sign}\) |
| 输出 | 表上的数 | **可执行动作**（见下） |
| 失败形态 | 分数低 | 线索不可读 / 无方向 / 只有 ID 名单 |

Bench 做过关只说明「尺子有区分度」；布控还要 **可读的 tip + 正负向**，否则运营无法下手。

---

## 1. 系统能发出的 actions（动作清单）

每次 graph shift 交付一条线索卡后，动作分四级（由轻到重）：

### L0 — 观察（默认）
- \(R_t=0\) 或 \(\|\delta\|\) 小：不派单  
- 只记 \(T_t^{\mathrm{MMD/PO/cmean}}\) 曲线

### L1 — 盯梢 / 加监（线索级）
触发：\(R_t=1\) 且有 \(K^\star\)  
动作举例：
- 把 \(K^\star\) 内 user / merchant / item 进 **关注名单**（TTL=若干 \(\Delta t\)）  
- 提高抽检率 / 人审队列优先级  
- 看板挂出 tip 卡片（见 §2）

### L2 — 策略收紧（可回滚）
触发：同一 \(K^\star\) 多日加剧（MMD/PO/cmean 持续高）+ 负向结局或业务规则命中  
动作举例：
- 对该支撑降权曝光 / 限频 / 验证码 / 延迟结算  
- 要求补充核验（不直接永久封）

### L3 — 案件升级（人工）
触发：L2 后仍尖 + tip 与已知欺诈模式同向（运营词典映射）  
动作举例：
- 开案、冻结资金链路、跨团伙并案  
- **算法到此为止**；定性归调查，不声称 ATE

**禁止当 action 的：** 仅凭 FSDS AUC 高就封号；仅凭 Top ID 无 tip/符号就布控。

---

## 2. Tip 特征逻辑（WHAT）

### 2.1 Tip 从哪来
1. 先定支撑 \(K^\star\)（shift 尖端实体块）  
2. **只在 \(K^\star\) 边**上跑一次 FSDS：  
   `Scaler → Var → SelectKBest(F) → HGB/LR`  
3. Top-\(k\) 特征 = **tip 集 \(J^\star\)**

Tip = 「在这个异常块上，最能分开 convert/非 convert（或代理标签）的图特征」——  
不是全图全局重要特征，也不是因果效应最大的特征。

### 2.2 原型里常见 tip（TencentGR）
| tip（例） | 可读含义（布控话术） |
|---|---|
| `u_span_sec` | 用户活跃时间跨度异常（短窗刷 / 长挂机） |
| `i_credit_last` / `i_share_last` | 成功路径压在末次物品上（末跳劫持感） |
| `i_n_covisit_neighbors` | 共现邻域变密/变稀（团伙共点） |
| `ui_pop_mismatch` | 用户活跃度 vs 物品热度错配 |

映射到运营词典后才能进 L2/L3；裸特征名不派策略。

### 2.3 Tip 与「准确率」的关系
- Bench：tip 选得好 → 局部 AUC/GT hit 更高（评测）  
- 布控：tip 是 **解释与分流依据**（这类异常走「末次归因」队列 vs 「活跃跨度」队列）

---

## 3. 正负向逻辑（方向，必报）

幅度 \(\|\delta\|\) / F-score **不够**；方向决定动作极性。

### 3.1 结局方向 \(\Delta\bar y\)
\[
\Delta\bar y(K^\star)=\bar y_{\mathrm{cur}}(K^\star)-\bar y_{\mathrm{ref}}(K^\star)
\quad(y=\text{click/convert})
\]

| \(\Delta\bar y\) | 标签 | 反欺诈读法（例） |
|---|---|---|
| \(>0\) | **正向** | 块上成功变多：刷量/互点嫌疑，或正常爆款（需 tip 辨别） |
| \(<0\) | **负向** | 块上成功变少：流量劫持、劣质灌入、被打压后残留 |
| \(\approx 0\) | flat | X 漂了但 click 率没动：先当分布/供给漂移，慎升 L2 |

### 3.2 Tip 方向 \(\mathrm{sign}(\delta_j)\)
\[
\mathrm{sign}(\delta_j)=\mathrm{sign}(\mu_{\mathrm{cur},j}-\mu_{\mathrm{ref},j})
\quad j\in J^\star
\]

| | 含义 |
|---|---|
| tip **+** | 该图特征在 \(K^\star\) 上 **升高**（如 span 变长、末次份额变大） |
| tip **−** | 该特征 **降低** |

能量仍按 \(\delta_j^2\) 排序，但卡片必须打印符号。  
例：`i_share_last (+)` + \(\Delta\bar y>0\) → 「末跳占成功变多且块上 click 升」→ 偏刷量/末跳操控队列；  
`u_span_sec (−)` + \(\Delta\bar y<0\) → 「活跃跨度变短且成功率掉」→ 偏劣质短会话灌入。

### 3.3 实体共向 \(s_v=\langle D_{v:},\delta\rangle\)（可选）
- \(s_v>0\)：该实体与父 shift **同向**（主嫌）  
- \(s_v<0\)：对冲/代偿（可能是对照块，勿当主嫌）

### 3.4 交付条（布控最小集）
```text
[Action-card]
  fire=R_t; level=L0|L1|L2|L3
  K*=...; n=...
  Dy=... (pos|neg|flat)
  tips={ u_span_sec:+, i_share_last:+, ... }
  next= watchlist | tighten | case
```

---

## 4. 反欺诈场景逻辑（怎么用这套线索）

### 4.1 角色
| 角色 | 用什么 |
|---|---|
| 监测 | \(R_t\)、Lead、AR（何时有事） |
| 研判 | \(K^\star\) + tip + 正负向（看谁、像哪类案） |
| 策略 | L1–L2 可回滚动作 |
| 调查 | L3 开案；算法只提供时间线与特征卡 |

### 4.2 典型故事线（连续时间）
1. **日切 \(R_t=1\)** → 出卡  
2. Drill 到 merchant/user 块 \(K^\star\)  
3. 读 \(\Delta\bar y\) 与 tip 符号 → 分到「刷量 / 灌入 / 劫持」等运营桶  
4. L1 关注；若连续多日 `score` 加剧 → L2 限频  
5. 人审确认模式 → L3；同时回写「tip 词典」（哪些 tip 组合对应哪类欺诈）

### 4.3 和规则引擎的关系
- 规则：硬条件（设备、支付、黑词）  
- Graph shift 线索：软的 **结构+分布** 异常波次  
- 合流：规则命中 ∩ \(K^\star\) → 优先；仅 shift 无线索词典 → 停在 L1

### 4.4 明确不做什么
- 不把 \(K^\star\) 当「已认定欺诈团伙」对外通报  
- 不用全局 AUC 替代单案卡片  
- 不在无符号 tip 时自动 L2  

---

## 5. Takeaway

1. Bench 准 ≠ 布控完成；布控要 **动作 + tip + 正负向**。  
2. Tip = \(K^\star\) 上 FSDS 选出的图特征；符号 = cmean 的 \(\mathrm{sign}(\delta_j)\) 与 \(\Delta\bar y\)。  
3. Actions：L0 观察 → L1 关注 → L2 可回滚收紧 → L3 人工开案。  
4. 反欺诈里本方法 = **连续图谱变动线索机**，与规则/人审互补。
