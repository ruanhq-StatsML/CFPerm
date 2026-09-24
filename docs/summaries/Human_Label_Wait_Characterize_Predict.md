# 等人打标：怎么刻画 · 怎么预测

> 审出加速包当前停在 **human-gate NOT_READY**。  
> 本文只回答：等人打标这条逻辑 **怎么刻画、怎么预测、打标前后各干什么**。  
> 不改 Drill；不加新 feature；不解冻 ETA/规则。

相关：[`Review_Accel_Personas_Feedback_Loop.md`](Review_Accel_Personas_Feedback_Loop.md) · [`Direction_UR_Gray_Unfreeze.md`](Direction_UR_Gray_Unfreeze.md) · handoff [`artifacts/review_accel_human_handoff.md`](artifacts/review_accel_human_handoff.md) · gate [`artifacts/review_accel_human_gate.json`](artifacts/review_accel_human_gate.json)

---

## 0. 一句话

**等人打标 = 商业置信的门闩，不是工程缺陷。**  
管道（出卡/灰度/CLI）已通；缺的是真人「这卡有没有帮我少翻页」的二元票。没有票 → 有用率不可信 → 禁止用 smoke 冒充 ROI → feature freeze。

---

## 1. 状态机（等人打标逻辑）

```text
                    feature_freeze=True
                            │
     ┌──────────────────────▼──────────────────────┐
     │  WAIT_LABEL  (现在)                          │
     │  n_human=0 · smoke 可有可无 · 6h 巡检看板    │
     └──────────────┬───────────────────────────────┘
                    │  ≥1 条真人 useful|not_useful
                    │  (--reviewer 真名，非 smoke/bot/…)
                    ▼
              GATE READY
                    │
        ┌───────────┴───────────┐
        │ n=1..4                │ n≥5 by 队列
        │ 只允许：话术草稿/催收 │ 才允许：按 UR 调灰度建议
        └───────────┬───────────┘
                    │ G1–G7 …
                    ▼
              才议解冻 ETA/规则
```

| 状态 | 判据（刻画） | 允许 | 禁止 |
|---|---|---|---|
| **WAIT_LABEL** | `n_human < min_human(=1)` | 催打标、刷新看板、修接通 bug | 新 tip 模型 / Web UI / 自动 flags / ETA / 规则 |
| **READY_THIN** | `1≤n_human<5` | 读 note、改 1 句读法草稿 | 动 `sample_rate` / 解冻旁路 |
| **READY_UR** | `n≥5`（建议再按队切） | `suggest` 灰度建议 → **人确认** cp | 无确认自动覆盖现网 |
| **UNFREEZE_CAND** | G1–G7 全过 | 议 U2 ETA 软门 | 跳步解冻 |

现在：**WAIT_LABEL**。

---

## 2. 怎么刻画（指标，现在就能填）

### 2.1 Gate 层（是否还在等）

| 符号 | 定义 | 现值 | 好长相 |
|---|---|---|---|
| `n_human` | reviewer∉{smoke,bot,ci,test,ood_timer,agent,cursor} 且 label∈{useful,not_useful} | **0** | ≥1 → READY |
| `n_smoke` | 上述名单内反馈 | 1 | 仅连通性，**不计 UR** |
| `feature_freeze` | `n_human < min_human` | True | READY 后 False |
| `gate_age_h` | 自上次 NOT_READY 巡检起小时数 | 6h 一巡 | 看板 ts 刷新即可 |
| `handoff_sent` | 是否已把 §0 IM 发出 | 人工勾 | 1 |
| `owner_named` | OWNER / 审核同学是否落名 | 常空 | 1 |

脚本：`check_review_accel_human_gate.py` → `review_accel_human_gate.json`。

### 2.2 等待过程层（等人效率，不是模型 AUC）

| 符号 | 定义 | 用途 |
|---|---|---|
| `T_ask→label` | 发出 handoff → 首条真人标 的小时数 | 刻画「等人」瓶颈在不在 IM/责任人 |
| `n_ping` | 催打标次数 | 过高说明交接失败，不是卡内容问题 |
| `label_mix` | useful : not_useful | n=1 时只作极性提示，**不算 UR** |
| `note_rate` | 带 note 的真人标占比 | 预测「改话术有没有燃料」 |

### 2.3 打标后才启用的层（现在预测，不能当实测）

| 符号 | 定义 | 何时可信 |
|---|---|---|
| \(\mathrm{UR}\) | useful/(useful+not_useful) | n≥5；按队更好 |
| \(\mathrm{UR}(q)\) | 队列切片有用率 | 每队 n(q)≥5 才调灰度 |
| Score | 可点×(接通+有人用+ROI) | 有人用之前 Score 里「有人用」项=0 |

**刻画原则：** WAIT_LABEL 阶段只报 gate 层 + 等待过程层；**禁止**用 smoke 算 UR，禁止用「接通」冒充人分钟↓。

---

## 3. 怎么预测（等人期间 vs 打标之后）

### 3.1 等人期间（WAIT_LABEL）——能预测什么

| 预测问题 | 诚实答案 | 动作 |
|---|---|---|
| 再堆 feature 能否提高有用率？ | **预测无效**（无真人基线） | freeze；keep=false 的点子进 log 即可 |
| 把 r 拉满能否证明效果？ | Cover↑ ≠ 人用↑ | 维持演练态 Cover；不报业务成功 |
| 多久能 READY？ | \(\hat T \approx T_{\mathrm{ask}\to\mathrm{label}}\)；未知则「责任人未定」 | 填 OWNER；发 IM；6h 巡检刷新看板 |
| smoke useful 算不算？ | **不算** | gate 已排除 smoke reviewer |
| 现在解冻 ETA？ | **预测负 ROI**（G1 未过） | 禁 |

**预测句（等人中）：**  
唯一正期望动作 = **降低 \(T_{\mathrm{ask}\to\mathrm{label}}\)**（发 handoff、点名审核同学）；其它研发动作的人分钟 ROI 预测 = 0。

### 3.2 首条真人标落地之后（READY_THIN）——怎么预测下一步

设首票为 \(y_1\in\{\mathrm{useful},\mathrm{not\_useful}\}\)：

| \(y_1\) | 预测读法 | 下一步（仍 n&lt;5） |
|---|---|---|
| useful | 卡「至少有人吃」；**不能**外推 UR≥0.7 | 同队列再收 4 票；可改 1 句读法更贴 note |
| not_useful | 卡「至少有人嫌烦」；**不能**立刻 kill 全局 | 读 note：走错族→改词典；看不懂→改话术；干扰→考虑该队 r↓ 草稿（**不自动 cp**） |

**禁止预测：** 用 1 票推出「全站有用率」或「件均耗时已降」。

### 3.3 n≥5 之后（READY_UR）——灰度/话术预测

与 [`Direction_UR_Gray_Unfreeze.md`](Direction_UR_Gray_Unfreeze.md) §2 对齐：

\[
\begin{aligned}
\widehat{\mathrm{exposures}} &\propto r \\
\widehat{\mathrm{useful\ hits}} &\approx r\cdot\mathrm{UR}\cdot N \\
\widehat{\mathrm{harm}} &\propto r\cdot(1-\mathrm{UR})
\end{aligned}
\]

| UR 带 | 预测 | 动作建议（人确认） |
|---|---|---|
| ≥0.7 | 升 r 期望净收益为正（若误导成本≈2×有用） | 略升 sample_rate |
| [0.4,0.7) | 先改话术期望 ΔUR↑，再动 r | hold r |
| &lt;0.4 | 降 r 先减 harm | 队列 r↓ 或 kill 建议 |

### 3.4 和「等人」相关的失败预测

| 现象 | 预测根因 | 纠偏 |
|---|---|---|
| 看板刷了多轮仍 n_human=0 | handoff 未发出 / 无 OWNER | 停刷功能；只追 IM 与点名 |
| 连续 not_useful 且 note=空 | 打标无信息量 | 要求 note 必填一句 |
| 全是同一 reviewer | 有偏样本 | 再拉 1 人；再谈 UR |
| Direction 仍缺却收到 useful | 可能只是「有字总比没有好」 | 不解读为 S1 刷量已验证 |

---

## 4. 操作闭环（刻画 → 预测 → 动作）

```text
每 6h:
  check_review_accel_human_gate.py
  → 更新 gate.json + waiting_board
  → 若仍 WAIT_LABEL:
       刻画: n_human, handoff_sent, owner_named
       预测: 唯加速 READY 的杠杆 = 催打标
       动作: 发/重发 handoff §0；禁止新 feature
  → 若 READY:
       刻画: n, label_mix, note_rate
       预测: n<5 只改话术；n≥5 才 suggest 灰度
       动作: 恢复 10min 迭代仅限审出主链；解冻另走 G1–G7
```

打标命令（handoff §2）：

```bash
PYTHONPATH=. python3 scripts/tencent_gr/log_review_card_feedback.py \
  --card results/tencent_gr_review_agent_card/review_agent_card.json \
  --label useful \   # 或 not_useful
  --reviewer YOUR_REAL_NAME \
  --note "少翻了哪一页 / 哪里烦"
PYTHONPATH=. python3 scripts/tencent_gr/check_review_accel_human_gate.py
```

---

## 5. Takeaway

1. **刻画等人：** `n_human` / `feature_freeze` / `T_ask→label` / handoff 是否发出——不是 Acc，不是 Cover。  
2. **预测等人中：** 正期望动作只有催打标；堆 feature / 拉满 r / 解冻 ETA 的人用 ROI **预测为 0 或负**。  
3. **首票后：** useful/not_useful 只改「薄置信」读法；n&lt;5 不调灰度。  
4. **n≥5：** 才用 UR 预测升/降 r；人确认；不自动 cp。  
5. 现在：WAIT_LABEL → 看板已刷；要 READY，缺的是 **1 个真名 reviewer 的一票**。
