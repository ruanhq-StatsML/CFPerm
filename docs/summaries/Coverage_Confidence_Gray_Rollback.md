# 落地效果 · 覆盖/置信 · 灰度回滚（可辩护版）

> 回答：实际落地效果如何；**覆盖（coverage）与置信（confidence）**怎么定义、怎么量、凭什么这样量；  
> **灰度回滚**为什么这样设计才站得住脚（对照翻车、对照期望伤害、对照实测）。  
> 代码锚点：`export_review_agent_card.flags_allow` · `configs/review_agent_card_flags.json` · `suggest_review_card_sample_rate.py`。

相关：[`Direction_UR_Gray_Unfreeze.md`](Direction_UR_Gray_Unfreeze.md) · [`Scope_Review_Agent_Card.md`](Scope_Review_Agent_Card.md)

---

## 0. 落地效果（诚实分层，先定调）

| 层 | 落地了什么 | 证据 | 还没落地 | 敢不敢说「效果好」 |
|---|---|---|---|---|
| **L-eng** | 出卡→字段→灰度门→建议→dry-run POST | pytest / smoke / CI | — | 工程效果：**是** |
| **L-gray** | `sample_rate` / `kill_switch` / 按队列覆盖；拒绝自动覆盖 | 单测 + 本机 Cover≈rate（§1.2） | 现网 SOAR 真 POST | 灰度可控：**是** |
| **L-polarity** | 词典 S1/S2/S3 | 文档 | 样例 summary **缺 direction** | 案由可分：**否** |
| **L-human** | feedback CLI + gate | 脚本在 | **n_human=0** | 人用置信：**否** |
| **L-biz** | KPI 人分钟公式 | 公式有 | **无前后实测** | 商业效果：**否** |

**一句话：** 灰度与回滚在**工程上可辩护且可演练**；覆盖/置信在卡注入层已实现；**业务效果尚未被真人数据证明**——这不削弱回滚设计的正当性，但限制「效果很好」的表述。

**可辩护边界：** 本文论证的是「Cover 旋钮真的控得住 + 出事能退得回来 + 升量比降量更苛刻」；**不是**「已经省了人分钟」。把两件事绑在一起报，才是不可辩护。

---

## 1. 覆盖（Coverage）：注入层暴露了多少

### 1.1 定义（可计算、可审计）

对一次导出，是否把卡建议写入工单侧：

\[
\mathrm{Cover}
=\Pr(\texttt{gray\_flags.allow}=\mathrm{true})
\]

实现（确定性哈希，非随机抖）：

```text
seed   = source                    # 无队列覆盖时
       | source::queue_bucket      # 有 queue_sample_rates 覆盖时
h      = md5(seed)[:8] % 10000
sampled = (h / 10000) < sample_rate
allow   = enabled ∧ ¬kill_switch ∧ sampled
```

| 旋钮 | 作用 | 回滚含义 | 为什么独立存在才可辩护 |
|---|---|---|---|
| `sample_rate∈[0,1]` | 全局覆盖目标 | 降 r ⇒ 覆盖近似按比例降 | 日常收缩不必动紧急开关 |
| `queue_sample_rates[q]` | 单队列覆盖覆盖全局 | 只关「有害队」 | 一队差不牵连健康队（§3.2 爆炸半径） |
| `enabled=false` | 全关但非紧急 | 软关 | 与 kill 语义分开，避免「关了还以为在应急」 |
| `kill_switch=true` | 紧急全关 | **硬回滚**（优先） | 盖过 enabled/rate，一条命令灭火 |

### 1.2 实测（本机 1000 个 source，`flags_allow`）

| 设定 `sample_rate` | 经验 Cover | 偏差 Cover−r |
|---|---|---|
| 0.00 | 0.000 | 0.000 |
| 0.05 | 0.045 | −0.005 |
| 0.10 | 0.085 | −0.015 |
| 0.20 | 0.179 | −0.021 |
| 0.30 | 0.282 | −0.018 |
| 0.50 | 0.504 | +0.004 |
| 0.70 | 0.690 | −0.010 |
| 1.00 | 1.000 | 0.000 |

补充实测：

| 设定 | 结果 |
|---|---|
| `kill_switch=true`，r=1 | Cover=0，reason 全为 `kill_switch` |
| `queue_sample_rates[brush_pos]=0`，全局 r=1 | 该队 Cover=0，reason=`queue_not_sampled` |
| 同 source 两次采样 | **完全一致**（稳定可复查） |
| 全局 r=0.5 + 坏队 r=0 | 坏队 Cover=0；好队 / 无队 ≈0.50（外科手术成立） |

→ **覆盖可控、可预期**：调 r 近似线性调暴露面；可辩护为「灰度流量旋钮」，不是装饰字段。

### 1.3 为什么用 md5 而不是 `random()`（可辩护点）

| 要求 | md5 稳定采样 | 真随机 |
|---|---|---|
| 复查同一工单 | 同 source → 同 allow | 可能今天进、明天出 |
| 审计「当时为什么进灰」 | 可由 seed 复算 | 难复盘 |
| 调 r 后新旧集合关系 | 单调近似嵌套（r↑ 多进、少出） | 集合抖动 |
| A/B 与值班对账 | 可复现 | 「抽运气」不可辩护 |

灰度是**风控暴露面**，不是广告 CTR 实验；稳定性优先于无偏抖动。

### 1.4 覆盖 ≠ 业务覆盖（三层必须拆报）

| 层 | 含义 | 现状量级 | 误报风险 |
|---|---|---|---|
| **注入覆盖** Cover | 卡是否 allow 进工单字段 | 可由 r 拉到 1 | 单独报 r=1 像「全上线」 |
| **Direction 覆盖** `dir_cover` | summary 是否带极性 | 样例 ≈0 | 假 flat → 刷量/灌入讲不通 |
| **人用覆盖** | 真人打标占比 | 0 | Score≈1 只是接通，不是有用 |

谈「全面上线」必须三者分开报；**只报 sample_rate=1 会高估业务覆盖**——这是对本包最常见的过度宣称路径。

---

## 2. 置信（Confidence）：我们有多敢放量

这里的置信不是模型 AUC，而是 **放量置信 = 证据是否够支撑当前 Cover**。

### 2.1 为什么要分层（核心可辩护论点）

CI 绿只能证明「管道不炸」，**不能**证明「审核同学觉得卡有用」。  
若用 C0 直接开大 Cover，等于用工程置信冒充人用置信——值班长无法对客诉负责。

| 层 | 符号 | 证据 | 允许的 Cover 天花板 | 缺它时的谎 |
|---|---|---|---|---|
| **C0 工程置信** | smoke/CI 绿 | 自动测试 | 预发可 r=1 **演练** | 「测过了」≠「人吃了」 |
| **C1 极性置信** | `dir_cover`、非假 flat | direction 落盘 | 才敢讲 S1/S2 分流 | 全 flat 当刷量 demo |
| **C2 人用置信** | \(n,\mathrm{UR}\) | 真人反馈 | 才敢按 UR 调 r | 小样本噪声驱动回滚抖动 |
| **C3 旁路置信** | G1–G7 | 解冻清单 | 才敢 ETA/规则 | 主链未稳就横向扩散 |

**硬规则：** 置信层级决定天花板 Cover，**不能用 C0 证明 C2**。

\[
r_{\max}
=
\begin{cases}
0 & \text{kill\_switch}\\
r_{\mathrm{trial}}\le 0.1 & \mathrm{C0}\land\neg\mathrm{C2}\quad\text{（仅预发试跑，不进现网全量）}\\
f(\mathrm{UR},n) & \mathrm{C2}\\
\end{cases}
\]

现状映射：C0≈有，C1≈无，C2≈无，C3≈停 U0/U1 → **现网正当 Cover 应保持低或演练态**；configs 里 `sample_rate=1.0` 是「闸门打开供导出」，不等于「已有 C2 放量许可」。

### 2.2 期望伤害：为什么 Cover 必须跟 UR 走

设一周暴露 \(N\) 张可打标工单，有用省时 \(b\)、无用浪费 \(c\)（人分钟）。粗期望：

\[
\begin{aligned}
\mathbb{E}[\text{benefit}] &\approx \mathrm{Cover}\cdot \mathrm{UR}\cdot N\cdot b \\
\mathbb{E}[\text{harm}] &\approx \mathrm{Cover}\cdot (1-\mathrm{UR})\cdot N\cdot c \\
\mathrm{Net} &\approx \mathrm{Cover}\cdot N\cdot \bigl(\mathrm{UR}\,b - (1-\mathrm{UR})\,c\bigr)
\end{aligned}
\]

取保守比 \(c=2b\)（误导比有用更贵——审核被带偏要纠错），\(N=100\)：

| Cover | UR | Net（相对 \(b\)） | 读法 |
|---|---|---|---|
| 0.1 | 0.2 | **−14** | 低有用还放一点也亏 |
| 0.1 | 0.7 | +1 | 小流量试探可正 |
| 0.5 | 0.2 | **−70** | 半量 × 差 UR = 大亏 |
| 0.5 | 0.7 | +5 | 健康才敢加 Cover |
| 1.0 | 0.2 | **−140** | 全量 × 差 UR 不可辩护 |
| 1.0 | 0.7 | +10 | 有 C2 才谈全量 |

**推论（可写进值班话术）：**  
1. Cover 是 **乘数**——UR 差时先砍 Cover，比改模型更快止损；  
2. 盈亏平衡约在 \(\mathrm{UR} > c/(b+c)=2/3\) 附近（本假设下）；故 **high_rate=0.7** 不是拍脑袋，是「误导更贵」下的谨慎线；  
3. 没有真人 UR 时，谈「把 r 拉满提效果」在期望上 **不可辩护**。

### 2.3 从有用率到放量置信（与 `suggest` 对齐）

`suggest_review_card_sample_rate`（**只建议、不自动写现网**）：

| 证据 | 动作建议 | 置信解释 | 脚本实测（例） |
|---|---|---|---|
| \(n<\min\_n(=5)\) | hold | **证据不足 → 不调 Cover** | n=0/3 → hold |
| \(\mathrm{UR}<\mathrm{kill\_rate}(0.2)\) | kill_switch 建议 on + r≤0.1 | 置信崩塌 → 硬回滚建议 | n=0+5, UR=0 → kill + r=0.1 |
| \(\mathrm{UR}<\mathrm{low}(0.4)\) | r ← max(0.1, 0.5r) | 置信低 → 收缩暴露 | n=1+4, UR=0.2 → r:1→0.5 |
| \(\mathrm{UR}\ge\mathrm{high}(0.7)\) | r 略升、关 kill | 置信高 → 谨慎放量 | n=7+3, UR=0.7, r0=0.4 → r=0.5 |
| 中间带 | hold | 先改话术再动 r | UR∈[0.4,0.7) → hold |

**非对称升量（实测路径）：** 在 UR≥0.7 持续成立时，从 r=0.1 爬到 1.0 约需 **8 步**（×1.25 与 +0.1 取大）：

```text
0.10 → 0.20 → 0.30 → 0.40 → 0.50 → 0.625 → 0.781 → 0.977 → 1.00
```

而降量：**一步减半**（或直接 kill）。  
→ **损人优先于抢覆盖**——这是对「审核被带偏成本 > 少看一张卡」的制度表达，不是同一步长的对称 A/B。

**可辩护点汇总：**  
1. 小样本不调（避免噪声驱动回滚抖动）；  
2. 降流量有下限 0.1，极端才 kill（保留观测窗，不全盲）；  
3. 升流量更慢；  
4. 输出是 `suggested.json` + disclaimer，**禁止自动 cp**——人在回路是置信的最后一道闸。

### 2.4 卡级「置信展示」（给审核看的）

| 字段 | 含义 | 回滚时可见性 |
|---|---|---|
| `gray_flags.allow` | 本卡是否建议写入 | false = 未注入 |
| `gray_flags.reason` | ok / not_sampled / kill_switch / queue_not_sampled / disabled | **可审计归因** |
| `suggested_action_level` | allow→L1_watch；否则 **L0_observe** | 动作级连带降级 |
| `queue_bucket=GRAY_DISABLED` | 回滚态显式可见 | 避免假装还在「刷量盯梢」 |

→ 回滚不是静默丢卡，而是 **降级为观察**；审核同学看到的是「系统主动降级」，不是「系统坏了」。这对客诉与值班复盘都可辩护。

---

## 3. 灰度回滚逻辑（Justifiable 设计）

### 3.1 回滚阶梯（由轻到重）

```text
L-roll-1  单队列 queue_sample_rates[q]=0     【外科手术】爆炸半径 = 一队
L-roll-2  全局 sample_rate 减半/降到 0.1   【收缩】爆炸半径 = Cover↓
L-roll-3  enabled=false                      【软关】非紧急全停
L-roll-4  kill_switch=true                   【硬关】优先应急；盖过一切
```

实测验收点：

| 动作 | allow | reason | 动作级 / 队列 |
|---|---|---|---|
| L-roll-4 | false | `kill_switch` | L0 |
| L-roll-1（队 r=0） | false | `queue_not_sampled` | `GRAY_DISABLED` |
| L-roll-3 | false | `disabled` | L0 |

**阶梯存在的理由：** 不是炫技，是让「坏一队」与「字段事故」用不同工具——避免一有客诉就只能全站 kill，也避免字段事故还在纠结调哪一队的 rate。

### 3.2 爆炸半径（为什么要队列级）

设队 \(q\) 日均工单 \(N_q\)，全局 Cover \(r\)，坏队占总量份额 \(s_q=N_q/\sum N\)。

| 策略 | 坏队曝光 | 好队误伤 |
|---|---|---|
| 只动全局 r→0 | 0 | **全部好队也被关** |
| L-roll-1 坏队 r=0 | 0 | **0**（好队仍 ≈ 原 Cover） |

本机对照：全局 r=0.5 + `bad_q=0` → bad Cover=0，good Cover≈0.50。  
**可辩护句：** 「局部差用局部闸」——风控值班的标准做法，迁移到卡注入层。

### 3.3 为什么这样可辩护（对照常见翻车）

| 原则 | 我们怎么做 | 反面（不可辩护） |
|---|---|---|
| **可预期覆盖** | md5 稳定采样，同 source 同结果 | 真随机导致复查不一致 |
| **可局部回滚** | 按队列覆盖 | 一队差就全局杀 |
| **硬开关独立** | kill_switch 盖过 enabled/rate | 和 sample_rate 缠在一起难关 |
| **人在回路** | suggested.json 不自动 cp | 反馈噪声自动关量 / 自动开量 |
| **状态可见** | reason + L0 + GRAY_DISABLED | 静默失败像系统坏了 |
| **不碰算法门** | 只闸注入，不改 Drill | 回滚误伤定位逻辑 / 难解释「模型变了」 |
| **非对称调参** | 降快升慢 + 小样本 hold | 有用率抖一下就全开 |
| **置信分层** | C0≠C2；r_max 跟证据走 | CI 绿就宣称业务成功 |

### 3.4 回滚决策表（值班用）

| 触发 | 置信判断 | 回滚动作 | 验收 |
|---|---|---|---|
| 单一队列 UR&lt;0.4，n≥5 | 局部差（C2 局部崩） | L-roll-1 该队 r=0 或 0.1 | 该队工单带卡↓，其它队不变 |
| 多队 UR 差 / 客诉「卡干扰」 | 系统差 | L-roll-2 或 L-roll-4 | 全局 allow≈0 或 Cover 腰斩 |
| 误 POST / 字段事故 | 接通事故（C0 之下） | **立即 L-roll-4** | 新导出全 L0；不先改话术 |
| Direction 大面积缺失 | C1 不可信 | 维持低 r；禁止解冻旁路 | dir_cover 修复前不升 r |
| UR&lt;0.2 且 n≥5 | C2 崩塌 | suggest→kill；人确认 L-roll-4 | 建议文件 + 现网 flags 一致 |
| 误触 kill | — | kill_switch=false 恢复 | 演练记录留痕 |

### 3.5 回滚演练（证明「能回得来」）

最低限度演练清单（解冻旁路前必做，对应解冻 G4/G7）：

1. r=1 导出 → `allow=true`，`reason=ok`  
2. `kill_switch=true` 再导出 → `allow=false` / L0 / `reason=kill_switch`  
3. kill=false，某队 r=0 → 仅该队 `GRAY_DISABLED` / `queue_not_sampled`  
4. 确认 **Drill/summary 未改**（只改 flags）  
5. （建议）跑一遍 `suggest` 在 n&lt;5 时输出 hold，证明不会「空数据自动杀」

演练留一份前后 `gray_flags` JSON，即「回滚可辩护」证据包。  
**没有演练记录就不谈解冻**——G4 不是形式主义，是证明 kill 路径在本环境真实可走。

### 3.6 一条因果链（从现象到动作）

```text
审核说「这卡在带偏」
  → 记反馈 useful/not_useful + queue_bucket
  → n(q)≥5？ 否 → 只改话术 / 继续收数（不调 r）
                 是 → UR(q)？
                      ≥0.7 → 可略升 r（人确认）
                      [0.4,0.7) → 先话术
                      <0.4  → L-roll-1 该队
                      <0.2 且多队 → L-roll-4
  → 验收：Cover(q)↓ 且 reason 可见；Drill 不变
```

这条链每一步都有**证据门槛**或**显式人确认**，所以回滚是「制度动作」而不是「感觉动作」。

---

## 4. 和「效果」怎么一起报（避免过度宣称）

推荐对外三句话：

1. **注入覆盖** 已由 `sample_rate` 线性可控（实测 Cover≈设定，|偏差|≲0.02 @n=1000）。  
2. **回滚** 支持队列级与 kill 硬关，状态对审核可见（L0 / GRAY_DISABLED）；升量约 8 步满量、降量一步减半。  
3. **业务效果** 仍待 Direction 落盘 + 真人 UR；在此之前 Cover↑ **不**等于人分钟↓。

禁止句式：「灰度开了所以审出加速成功了。」  
正当句式：「灰度闸门可调可回滚；人用与极性证据上来之前，Cover 保持与置信层级匹配。」

---

## 5. Takeaway

1. 落地强项在 **可控覆盖 + 可演练回滚**；弱项在 **极性与真人置信**。  
2. Coverage = Pr(allow)；Confidence = C0–C3 证据阶梯，决定 \(r_{\max}\)；**C0 不能冒充 C2**。  
3. 期望伤害里 Cover 是乘数：UR 差时先砍 Cover；high=0.7 / kill=0.2 与「误导更贵」一致。  
4. 回滚阶梯 L1 队列 → L4 kill；人确认；不改 Drill；状态可见。  
5. suggest 阈值（5 / 0.2 / 0.4 / 0.7）把「何时收缩/硬关」写死，升慢降快，可辩护。  
6. 下一步要「效果好看」：先补 Direction 与真人 UR，再让 Cover 跟着 C2 走——不是先把 r 拉满。
