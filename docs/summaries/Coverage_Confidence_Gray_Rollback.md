# 落地效果 · 覆盖/置信 · 灰度回滚（可辩护版）

> 回答：实际落地效果如何；**覆盖（coverage）与置信（confidence）**怎么定义、怎么量；  
> **灰度回滚**为什么这样设计才站得住脚。  
> 代码锚点：`export_review_agent_card.flags_allow` · `configs/review_agent_card_flags.json` · `suggest_review_card_sample_rate.py`。

相关：[`Direction_UR_Gray_Unfreeze.md`](Direction_UR_Gray_Unfreeze.md)

---

## 0. 落地效果（诚实分层）

| 层 | 落地了什么 | 证据 | 还没落地 |
|---|---|---|---|
| **L-eng** | 出卡→字段→灰度门→建议→dry-run POST | pytest / smoke / CI | — |
| **L-gray** | `sample_rate` / `kill_switch` / 按队列覆盖；拒绝自动覆盖 | 单测 + 实测 cover≈rate | 现网 SOAR 真 POST |
| **L-polarity** | 词典 S1/S2/S3 | 文档 | 样例 summary **缺 direction** |
| **L-human** | feedback CLI + gate | 脚本在 | **n_human=0** |
| **L-biz** | KPI 人分钟 | 公式有 | **无前后实测** |

**一句话：** 灰度与回滚在**工程上可辩护且可演练**；覆盖/置信在卡注入层已实现；**业务效果尚未被真人数据证明**——这不削弱回滚设计的正当性，但限制「效果很好」的表述。

---

## 1. 覆盖（Coverage）：注入层暴露了多少

### 1.1 定义（可计算）
对一次导出，是否把卡建议写入工单侧：

\[
\mathrm{Cover}
=\Pr(\texttt{gray\_flags.allow}=\mathrm{true})
\]

实现：确定性哈希采样（非随机抖）：

```text
h = md5(source [:: queue])[:8] % 10000
sampled = (h/10000) < sample_rate
allow = enabled ∧ ¬kill_switch ∧ sampled
```

| 旋钮 | 作用 | 回滚含义 |
|---|---|---|
| `sample_rate∈[0,1]` | 全局覆盖目标 | 降 r ⇒ 覆盖近似按比例降 |
| `queue_sample_rates[q]` | 单队列覆盖覆盖全局 | 只关「有害队」 |
| `enabled=false` | 全关但非紧急 | 软关 |
| `kill_switch=true` | 紧急全关 | **硬回滚**（优先） |

### 1.2 实测（本机 1000 个 source）
| 设定 `sample_rate` | 经验 Cover |
|---|---|
| 0.0 | 0.00 |
| 0.1 | ≈0.09 |
| 0.3 | ≈0.31 |
| 0.5 | ≈0.50 |
| 1.0 | 1.00 |

→ **覆盖可控、可预期**：调 r 近似线性调暴露面；可辩护为「灰度流量旋钮」，不是装饰字段。

### 1.3 覆盖 ≠ 业务覆盖
| | 含义 |
|---|---|
| 注入覆盖 | 卡是否 allow 进工单字段 |
| Direction 覆盖 `dir_cover` | summary 是否带极性（现状样例 ≈0） |
| 人用覆盖 | 真人打标占比（现状 0） |

谈「全面上线」必须三者分开报；**只报 sample_rate=1 会高估业务覆盖**。

---

## 2. 置信（Confidence）：我们有多敢放量

这里的置信不是模型 AUC，而是 **放量置信 = 证据是否够支撑当前 Cover**。

### 2.1 三层置信（可辩护框架）
| 层 | 符号 | 证据 | 允许的 Cover |
|---|---|---|---|
| **C0 工程置信** | smoke/CI 绿 | 自动测试 | 可在预发 r=1 演练 |
| **C1 极性置信** | `dir_cover`、非假 flat | direction 落盘 | 才敢讲 S1/S2 |
| **C2 人用置信** | \(n,\mathrm{UR}\) | 真人反馈 | 才敢按 UR 调 r |
| **C3 旁路置信** | G1–G7 | 解冻清单 | 才敢 ETA/规则 |

**规则：** 置信层级决定天花板 Cover，不能用 C0 证明 C2。

\[
r_{\max}
=
\begin{cases}
0 & \text{kill\_switch}\\
r_{\mathrm{trial}}\le 0.1 & \mathrm{C0}\land\neg\mathrm{C2}\quad\text{（仅预发）}\\
f(\mathrm{UR},n) & \mathrm{C2}\\
\end{cases}
\]

### 2.2 从有用率到放量置信（与 suggest 对齐）
`suggest_review_card_sample_rate`（**只建议、不自动写现网**）：

| 证据 | 动作建议 | 置信解释 |
|---|---|---|
| \(n<\min\_n(=5)\) | hold | **证据不足 → 不调 Cover** |
| \(\mathrm{UR}<\mathrm{kill\_rate}(0.2)\) | kill_switch 建议 on + r≤0.1 | 置信崩塌 → 硬回滚建议 |
| \(\mathrm{UR}<\mathrm{low}(0.4)\) | r ← max(0.1, 0.5r) | 置信低 → 收缩暴露 |
| \(\mathrm{UR}\ge\mathrm{high}(0.7)\) | r 略升、关 kill | 置信高 → 谨慎放量 |
| 中间带 | hold | 先改话术再动 r |

**可辩护点：**  
1. 小样本不调（避免噪声驱动回滚抖动）；  
2. 降流量有下限 0.1，极端才 kill；  
3. 升流量更慢（×1.25 / +0.1），非对称——**损人优先于抢覆盖**。

### 2.3 卡级「置信展示」（给审核看的）
| 字段 | 含义 |
|---|---|
| `gray_flags.allow` | 本卡是否建议写入 |
| `gray_flags.reason` | ok / not_sampled / kill_switch / queue_not_sampled… |
| `suggested_action_level` | allow→L1_watch；否则 **L0_observe** |
| `queue_bucket=GRAY_DISABLED` | 回滚态显式可见，避免假装在审 |

→ 回滚不是静默丢卡，而是 **降级为观察**，可审计。

---

## 3. 灰度回滚逻辑（Justifiable 设计）

### 3.1 回滚阶梯（由轻到重）
```text
L-roll-1  单队列 queue_sample_rates[q]=0     【外科手术】
L-roll-2  全局 sample_rate 减半/降到 0.1   【收缩】
L-roll-3  enabled=false                      【软关】
L-roll-4  kill_switch=true                   【硬关：优先应急】
```

实测：`kill_switch` → `allow=false, reason=kill_switch, L0`；  
`queue_rate=0` → `queue_not_sampled, GRAY_DISABLED`。

### 3.2 为什么这样可辩护（对照常见翻车）
| 原则 | 我们怎么做 | 反面 |
|---|---|---|
| **可预期覆盖** | md5 稳定采样，同 source 同结果 | 真随机导致复查不一致 |
| **可局部回滚** | 按队列覆盖 | 一队差就全局杀 |
| **硬开关独立** | kill_switch 盖过 enabled/rate | 和 sample_rate 缠在一起难关 |
| **人在回路** | suggested.json 不自动 cp | 反馈噪声自动关量 |
| **状态可见** | reason + L0 + GRAY_DISABLED | 静默失败像系统坏了 |
| **不碰算法门** | 只闸注入，不改 Drill | 回滚误伤定位逻辑 |
| **非对称调参** | 降快升慢 + 小样本 hold | 有用率抖一下就全开 |

### 3.3 回滚决策表（值班用）
| 触发 | 置信判断 | 回滚动作 | 验收 |
|---|---|---|---|
| 单一队列 UR&lt;0.4，n≥5 | 局部差 | L-roll-1 该队 r=0 或 0.1 | 该队工单带卡↓，其它队不变 |
| 多队 UR 差 / 客诉「卡干扰」 | 系统差 | L-roll-2 或 L-roll-4 | 全局 allow≈0 |
| 误 POST / 字段事故 | 接通事故 | **立即 L-roll-4** | 新导出全 L0 |
| Direction 大面积缺失 | 极性不可信 | 维持低 r；禁止解冻旁路 | dir_cover 修复前不升 r |
| 误触 kill | — | kill_switch=false 恢复 | 演练记录 |

### 3.4 回滚演练（证明「能回得来」）
最低限度演练清单（解冻旁路前必做，对应解冻 G4/G7）：

1. r=1 导出 → allow=true  
2. kill_switch=true 再导出 → allow=false / L0 / reason=kill_switch  
3. kill=false，某队 r=0 → 仅该队 GRAY_DISABLED  
4. 确认 **Drill/summary 未改**（只改 flags）  

演练留一份前后 `gray_flags` JSON，即「回滚可辩护」证据。

---

## 4. 和「效果」怎么一起报（避免过度宣称）

推荐对外三句话：

1. **注入覆盖** 已由 `sample_rate` 线性可控（实测≈设定）。  
2. **回滚** 支持队列级与 kill 硬关，状态对审核可见（L0）。  
3. **业务效果** 仍待 Direction 落盘 + 真人 UR；在此之前 Cover↑ **不**等于人分钟↓。

---

## 5. Takeaway

1. 落地强项在 **可控覆盖 + 可演练回滚**；弱项在 **极性与真人置信**。  
2. Coverage = Pr(allow)；Confidence = C0–C3 证据阶梯，决定 \(r_{\max}\)。  
3. 回滚阶梯 L1 队列 → L4 kill；人确认；不改 Drill。  
4. suggest 阈值（5/0.2/0.4/0.7）把「何时收缩/硬关」写死，升慢降快，可辩护。  
5. 下一步要「效果好看」：先补 Direction 与真人 UR，再让 Cover 跟着 C2 走——不是先把 r 拉满。
