# 方法论落地 Roadmap：推理 · 对齐 · 客服增量账

> 目标：把 OnlineRFPerm / PO-risk / 轴拆分 **具体创新点** 接到推理监控与对齐门禁，
> 再用 **动作回本** 与 **周归因** 把「该不该动」变成可对账的 ¥。  
> 价值句式：**act-vs-ignore 量 × 单价**；AUROC / P@10 只决定路由，不直接等于贡献。

本原型已跑通数字（客服助手 / HaluEval hop）：

| 层 | 数字 |
|----|------|
| 跳变 | fire=1，ratio≈11.7；幻觉 quiet→fire ≈0.07→0.82 |
| 量 | 少工单 206 / 少退款 130 / 多承接 316 |
| 钱 | 毛 ¥16,656 · 净 ¥16,590 · 单价±20% 带 ¥13,259–¥19,921 |
| 周归因 | W1 ¥7,189（43.2%）+ W2 ¥9,467（56.8%），缺口 ¥0 |
| 动作回本 | 全部 **same_day_payback**；净ROI 29×–113× |

对照原型：`docs/biz/METHOD_BIZ_SCENARIO_PROTOTYPE.md`  
贡献账：`docs/biz/CS_ASSISTANT_CONTRIBUTION.md`  
**OnlineRFPerm → LLM 推理（low-hanging，可不接 ¥）**：`docs/biz/ONLINERFPERM_LLM_INFER_ROADMAP.md`

---

## 0. 一张总图（读完再往下）

```text
表征 / 特征库几何预检
        ↓
① OnlineRFPerm  制度火：P(Y|X) 或标签律是否 hop？     ← 方法论创新 A
        ↓
② 轴拆分        文风 P(X) vs 偏好/幻觉 P(Y|X) vs RAG缺口  ← 创新 B
        ↓
③ po_risk0      实例排序：审谁 / √PO 重加权谁             ← 创新 C
        ↓
④ 动作路由      audit_topk | retrieval_refresh | model_rollback | style-only
        ↓
⑤ 对照窗        fire 期内 acted vs ignored（不是 RCT）
        ↓
⑥ 周归因 ¥      哪一周贡献了多少（财务对账）               ← 运营创新 1
        ↓
⑦ 动作回本      日成本 / 日均毛¥ → 几天回本 + 净ROI         ← 运营创新 2
        ↓
⑧ 单位经济带    量固定扫单价 → 财务假设还站不站得住
```

**方法论回答「火了吗 / 审谁 / 动哪轴」；周归因+回本回答「这周值多少、这个动作几天回本」。两层缺一不可。**

---

## 1. 方法论层面的具体创新（要 elaborate 的部分）

下面每条都写成：**创新是什么 → 相对常规做法差在哪 → 落地接口 → 本原型证据**。

### 1.1 创新 A — OnlineRFPerm：在线制度跳变检测（不是 PSI，不是离线 AUROC）

| | 内容 |
|--|------|
| **创新** | 用上一批拟合的浅层 RF 作「制度探针」μ₀，对当前流批次算 **OOS 误差比** `e_now/e_prev`；越过 gate γ 则 **fire**，把「标签律 / 偏好映射 hop」变成可时间戳的事件。 |
| **vs 常规** | PSI/KS 只看边缘分布；离线 AUROC 是静态切片。OnlineRFPerm 盯的是 **预测性制度**是否变了（\(P(Y\mid X)\) 代理），可挂在线 FDR。 |
| **不是什么** | 不是事实核查器；fire ≠「这句话错了」。 |
| **落地接口** | `agod.run_rfperm_stream` / `fit_online_probe`；HF 原型 `scripts/agod/hf_landing_protos.py`。 |
| **本原型** | HaluEval cut：`fired=1`，ratio≈**11.7**；HH-RLHF：judge_err 0.147→0.539（ratio≈3.67）+ hop@cut。 |

**推理侧落地：** serving 日志按 batch（或按日）喂 OnlineRFPerm；Y = 幻觉标签 / 人工纠错 / 裁判分；X = hash 特征或现有 embedding。fire 打开「对照窗」——只有 fire 期内的 acted/ignored 才进贡献账。

**对齐侧落地：** DPO/RLHF 合并前，用黄金偏好批作 ref；候选合并批作 stream。fire → **拒合并 / 限量合并 / 强制 Top-k 人工复核**（见 `vw_alignment_merge_gate`）。

### 1.2 创新 B — 轴拆分：文风 \(P(X)\) ≠ 偏好/幻觉 \(P(Y\mid X)\) ≠ RAG 支持缺口

| 轴 | 统计对象 | 工具 | **正确动作** | **错误动作（禁）** |
|----|----------|------|--------------|-------------------|
| 文风 / register | 长度、正式度、嵌入几何 | RF-domain AUC / FSDS covariate | 刷风格塔、改 decoding、素材工单 | 当成对齐失败去动偏好头 |
| 偏好 / 幻觉映射 | chosen/rejected、忠实度标签 | OnlineRFPerm + PO | √PO 门控、回滚、偏好重训 | 用 PSI(长度)「解释」偏好掉点 |
| RAG 支持 | knowledge–answer overlap / rag_hit | Top-k rag_hit 均值 | `retrieval_refresh` | 一刀切回滚生成模型 |

| | 内容 |
|--|------|
| **创新** | 把「漂了」拆成 **可并行工单的轴**；每轴有独立 KPI 与钱路径（客服工单 vs CTR/品牌）。 |
| **本原型** | style_domain_auc≈0.999（文风可分）与幻觉 hop **分列**；客服主账只跟幻觉轴；文风走对照列。 |

**对齐创新点（具体）：** 合成偏好污染、文风域移、真实偏好漂移三者可同时发生——轴拆分让 CI 给出 **三张独立门禁票**，而不是一个模糊的「alignment score」。

### 1.3 创新 C — po_risk0：制度跳变下的实例优先级（排序器，不是因果分）

| | 内容 |
|--|------|
| **创新** | fire 后，对 T=1 样本用相对 μ₀ 的残差构造 `po_risk0_i`；**只排序**审计 / √PO 重加权对象。quiet 时权重退火回均匀。 |
| **vs 常规** | 随机抽检浪费人力；静态置信度在制度跳变后校准失效。po_risk0 绑定 **当前制度探针**。 |
| **不是什么** | 不是因果根因；P@10 ≠ 少退款数（少退款来自动作落地后的费率差）。 |
| **本原型** | po_risk0 P@10=1.0 → 驱动 `audit_topk`；审计件成本进入净贡献扣减。 |

**推理侧：** 在线队列按 po_risk0 截断 Top-k 进人工 / 二次生成 / 强制检索。  
**对齐侧：** 合并批内对高 po_risk0 对做人工偏好复核，再决定是否进 DPO 数据。

### 1.4 创新 D — 门控适应：fire → √PO / drop-old / 动作臂；quiet → 退火

| | 内容 |
|--|------|
| **创新** | 把检测与适应绑成 **有开关的闭环**：只有 fire 才开 √PO 或动作臂成本；quiet 强制退火，避免「永久重加权」漂成新偏置。 |
| **vs 常规** | 持续 importance weight 或永久 LoRA 补丁，缺少「何时关」。 |
| **本原型动作臂** | `audit_topk` / `retrieval_refresh` / `model_rollback`（按 rag_hit 路由）。 |

### 1.5 创新 E — 价值闭环：检测 → 对照窗 → 量 → ¥（方法必须接到账）

| | 内容 |
|--|------|
| **创新** | 规定：**方法论输出不得直接乘一个「玄学 AUROC¥」**。合法路径只有  
`fire 打开窗 → act vs ignore 费率差 → 量 → 单价 → 毛/净¥ → 周归因 / 回本`。 |
| **本原型** | 净 ¥16,590；周缺口 0；三动作全当天回本。 |

这是 **方法论文风上的创新约束**：没有对照臂就不谈贡献。

---

## 2. 推理侧落地（Inference Monitoring）

### Phase I — 探针上线（已具备原型）

1. 按日/按 batch 采 serving 特征 + 弱标签（纠错、rag_hit、退款/工单事后标）。
2. OnlineRFPerm gate（本原型 γ≈1.25 量级）→ `fct_shift_signal.fired`。
3. fire 日写入对照窗；默认一半流量可切 ignored 臂做对照（或历史同制度窗）。

### Phase II — 动作路由（与本账对齐）

```text
fired?
  no  → quiet：均匀权重，只记基线费率
  yes → rag_hit Top-k 低? → retrieval_refresh
        else              → model_rollback 或 audit_topk(po_risk0)
```

### Phase III — 推理预算门控（下一跳）

- fire 时：高 po_risk0 请求升配（更多检索 / 更大 draft / 人工接管）。
- quiet 时：降回标准路径（退火）。
- **KPI**：不是 AUROC，而是 **少工单/退款/承接** 与 **动作回本天数**。

### 推理侧方法论创新小结

| 点 | 落地物 |
|----|--------|
| 制度时钟 | OnlineRFPerm fire 时间戳 |
| 请求优先级 | po_risk0 Top-k |
| 路径切换 | RAG vs generation 轴 |
| 关断条件 | quiet 退火 + 回本不达标则停该臂 |

---

## 3. 对齐侧落地（Alignment / RLHF·DPO / Judge）

### 3.1 RFPerm-as-Judge（可包成调用）

```text
ref_batch  --fit-->  μ0（浅层 RF judge）
cand_batch         --OOS err-->  e_now/e_prev ≥ γ ? fire
fire: po_risk0 → Top-k 人工复核 / 拒合并
quiet: accept + 权重=1
```

**创新点：** 不必先训大裁判模型；用可解释浅层探针做 **合并门禁**，大模型裁判本身也可被 OnlineRFPerm 盯（裁判分作 Y vs 黄金集）。

### 3.2 三票门禁（对齐 CI）

| 票 | 失败含义 | 动作 |
|----|----------|------|
| 文风域票 | \(P(X)\) 污染 | 拒「风格漂移数据」进偏好集；开素材工单 |
| 映射票 | OnlineRFPerm fire | 拒合并或限量 + Top-k 复核 |
| RAG/任务票 | router_task_shift | 刷新对应专家/检索，不整模回滚 |

### 3.3 与客服账的边界（必须写死）

对齐失败的钱路径 ≠ 客服工单路径。  
HH 风格轴高 AUC **不要**写进客服贡献账；客服账只跟幻觉/纠错制度 hop。

---

## 4. 动作回本 — 为什么重要、怎么算、怎么用（重点 elaborate）

### 4.1 定义（可对账）

对每个动作臂 \(a\)：

\[
\text{payback\_days}(a)
= \frac{\text{日均动作成本}_a}{\text{日均毛贡献}_a}
= \frac{C_a / D_a}{G_a / D_a}
= \frac{C_a}{G_a}
\]

- \(G_a\)：该动作日内，相对 **ignored 基准费率** 的毛增量¥（少工单+少退款+多承接）。
- \(C_a\)：动作日成本合计（回滚工程、刷索引、审计人力等）。
- **\< 1 天** ⇒ `same_day_payback`（当天回本）。
- **净 ROI** \(= (G_a - C_a) / C_a\)。

SQL：`vw_cs_assist_action_payback`（`sql/biz_value/07_cs_week_attribution_payback.sql`）。

### 4.2 本原型数字（直观）

| 动作 | 天数 | 日均毛¥ | 日均成本¥ | 回本 | 分档 | 净ROI |
|------|------|---------|-----------|------|------|-------|
| audit_topk | 2 | ¥2,326 | ¥20 | **0.009 天** | same_day | **112.8×** |
| retrieval_refresh | 2 | ¥2,356 | ¥40 | **0.017 天** | same_day | **57.9×** |
| model_rollback | 3 | ¥2,431 | ¥80 | **0.033 天** | same_day | **29.4×** |

读法：

1. **三动作都当天回本** → 方法论 fire 后的动作在本单位经济下 **不是烧钱实验**。
2. **audit_topk 回本最快** → 人力件便宜、且 po_risk0 排序有效时，审计是「最先该开」的便宜开关。
3. **model_rollback 净¥最大（¥7,054）但回本略慢** → 贡献主力 ≠ 回本最快；排期要用 **双键**：先开快回本臂止血，再开高净额臂巩固。
4. 若某臂回本 > 3 天或净ROI \< 1 → **自动降级该臂**（方法论上等于「适应成本超过制度收益」）。

### 4.3 方法论创新含义

| 常规 | 本路线 |
|------|--------|
| 「模型更好了」用 AUROC/偏好准确率汇报 | 每个适应动作有 **回本时钟** |
| 适应成本沉没、事后讲故事 | 成本进分母，**同一周**就能否决贵臂 |
| 只比毛贡献 | 毛贡献排序 ≠ 回本排序；两张表一起看 |

**创新一句话：** 把「门控适应」从算法开关升级成 **带回收期的投资决策**；OnlineRFPerm 决定投资窗口，回本表决定仓位。

### 4.4 落地玩法（运营）

1. Oncall 看板：按 `payback_bucket` 着色；`slow_payback` 自动提工单「降本或停臂」。
2. 与 po_risk0 联动：P@10 掉到阈值以下 → audit 臂回本恶化 → 收紧 Top-k 或改路由。
3. 与单位经济联动：单价 −20% 后若某臂跌出 same_day → 该臂在保守财务口径下改 **试点** 而非默认。

---

## 5. 周归因 — 为什么重要、怎么算、怎么用（重点 elaborate）

### 5.1 定义（财务周账）

对每个 acted 日 \(t\)，用 **全样本 ignored 臂费率** 作反事实：

\[
\begin{aligned}
\Delta T_t &= (r^{\text{ign}}_{\text{ticket}} - r_t^{\text{ticket}}) \cdot n_t \\
\Delta R_t &= (r^{\text{ign}}_{\text{refund}} - r_t^{\text{refund}}) \cdot n_t \\
\Delta C_t &= (r_t^{\text{contain}} - r^{\text{ign}}_{\text{contain}}) \cdot n_t \\
G_t &= \Delta T_t \cdot p_T + \Delta R_t \cdot p_R + \Delta C_t \cdot p_C
\end{aligned}
\]

再按 ISO 周汇总 \(\sum_{t\in W} G_t\)，并要求：

\[
\sum_W G_W \;\approx\; G_{\text{ledger}}
\quad(\text{attribution\_gap} \approx 0)
\]

SQL：`vw_cs_assist_acted_day_value` → `vw_cs_assist_week_attribution` + `vw_cs_assist_attribution_check`。

### 5.2 本原型数字（直观）

| 周起始 | acted天 | 会话 | 少工单 | 少退款 | 多承接 | 毛¥ | 占总毛% |
|--------|---------|------|--------|--------|--------|-----|---------|
| 2026-09-14 | 3 | 600 | 86.1 | 57.0 | 135.7 | **¥7,189** | 43.2% |
| 2026-09-21 | 4 | 800 | 119.9 | 73.0 | 180.3 | **¥9,467** | 56.8% |

周归因合计 vs 总账缺口：**¥0（0%）**。

读法：

1. **缺口≈0** → 日账加总能对上总增量，财务可收。这是「可对账」的硬条件，不是展示。
2. **W2 占比 56.8%** → 跳变后期动作日更多/会话更多；周报应解释 **覆盖天数**，别夸成「模型又变强了」。
3. 周经营看板里 fire 周可能出现 **净 ops 为负**（ignored 臂质量成本暴涨）——那是 **对照叙事**：不动作更亏；归因表只计 acted 相对 ignored 的增量，两者要一起读。

### 5.3 方法论创新含义

| 常规 | 本路线 |
|------|--------|
| 月末拍一个「AI 节省估算」 | **周粒度**归因，可进经营周会 |
| 增量与总账对不上 | `attribution_check` 锁缺口 |
| 把流量增长算成算法贡献 | 日级用费率差×当日会话；周只是加总 |

**创新一句话：** 周归因把 OnlineRFPerm 的「火情周」翻译成 **财务可签字的增量周**；方法时钟与会计时钟对齐。

### 5.4 落地玩法（经营）

1. 周会模板：火情周 WoW + 本周归因¥ + 本周各动作回本（`CS_ASSISTANT_WEEKLY_OPS_BRIEF.md`）。
2. 若某周 gap \> 阈值 → 先修臂定义/费率，再谈方法论。
3. 与流量情景正交：周归因不改单价；扩量用 `vw_cs_assist_traffic_scenarios`。

---

## 6. 动作回本 × 周归因：怎么拼成一张经营决策表

| 决策问题 | 看哪张表 | 本原型答案 |
|----------|----------|------------|
| 这周算法侧贡献了多少可对账¥？ | 周归因 | W2 ¥9,467（56.8%） |
| 该先加哪类动作的排班？ | 回本（快）+ 净¥（大） | 先 audit 止血，rollback 巩固净额 |
| 财务单价更保守还站得住吗？ | 单位经济带 | −20% 仍净 ¥13,259，且回本仍应 ≪ 1 天 |
| 要不要扩到日 1 万会话？ | 流量情景 | prod_mid 月毛 ≈ ¥3.57M（线性外推） |
| 文风漂了要不要动客服预算？ | 轴拆分 | **不要**；开素材/CTR 路径 |

**拼法：** 周归因定「窗口价值」，回本定「窗口内仓位」，单位经济定「财务口径韧性」，方法论定「窗口从哪天开始算」。

---

## 7. 分阶段 Roadmap（方法 → 推理/对齐 → 账）

### Now（已具备 · 本分支）

- [x] HF 子集：HH 文风/Judge + HaluEval 幻觉/RAG/router
- [x] DuckDB 贡献账：毛/净、动作拆分、流量、周看板
- [x] **周归因 + 对账缺口锁 0**
- [x] **动作回本 + same_day 分档 + 净ROI**
- [x] 单位经济 ±20% 带 + 累计曲线
- [x] 方法×情景×单位经济对照原型

### Next（1–2 个小迭代）

1. **Ignored 机会成本**：fire 但未动作日的「桌上留白¥」→ 推动提高 acted 覆盖率。  
2. **回本×单价联合敏感度**：单价 −20% 后各臂 payback_bucket 是否降级。  
3. **对齐合并门禁样例**接入真实候选批（RFPerm-as-Judge CI）。  
4. **推理升配**：高 po_risk0 请求强制检索/人工的 A/B（仍用 act-vs-ignore 记账）。

### Later（产品化）

1. OnlineRFPerm 挂真实 serving 时钟 + online FDR。  
2. 多表面（客服 / 推荐文案 / RAG QA）共用轴拆分，**分账**不混账。  
3. 周归因自动进财务 CSV 接口；回本不达标自动关臂。  
4. 门控 √PO 与 LoRA/KV 预算调度打通（仍服从 quiet 退火）。

---

## 8. 方法论创新清单（对外可讲的短句）

1. **在线制度探针（OnlineRFPerm）**：把「模型/标签律 hop」变成可火警的时间戳事件。  
2. **轴拆分工单**：文风 / 偏好 / RAG 分路径，禁止一个分数打天下。  
3. **po_risk0 绑定制度**：跳变下的审计优先级，quiet 强制退火。  
4. **适应必带回本时钟**：每个动作臂有 payback_days 与净ROI。  
5. **火情周必对账**：周归因加总 ≡ 总增量（gap≈0），方法时钟对齐会计时钟。  
6. **价值路径立法**：只承认 act-vs-ignore × 单价；禁止 AUROC¥。

---

## 9. 怎么跑 / 读哪些文件

```bash
PYTHONPATH=. python3 scripts/agod/run_biz_value_sql_demo.py
python3 scripts/agod/cs_assist_hop_sensitivity.py
python3 scripts/agod/method_biz_scenario_prototype.py
pytest -q tests/test_cs_assist_contribution.py
```

| 文件 | 用途 |
|------|------|
| `docs/biz/METHOD_LANDING_ROADMAP.md` | 本文：创新 + 推理/对齐 + 回本/周归因 |
| `docs/biz/LANDING_CONTAIN_METHOD_ITER.md` | 落地场景 · 多承接 · 异同 · 迭代环 |
| `docs/biz/METHOD_BIZ_SCENARIO_PROTOTYPE.md` | 方法×情景×单位经济异同 |
| `docs/biz/CS_ASSISTANT_CONTRIBUTION.md` | 完整贡献数字 |
| `docs/biz/CS_ASSISTANT_WEEKLY_OPS_BRIEF.md` | 周会粘贴简报 |
| `docs/reports/HF_LANDING_ELABORATION.md` | HF 子集方法 elaborations |
| `sql/biz_value/07_cs_week_attribution_payback.sql` | 周归因 + 回本 SQL |
| `results/agod/biz_value_sql/cs_assist_action_payback.json` | 回本数字 |
| `results/agod/biz_value_sql/cs_assist_week_attribution.json` | 周归因数字 |

---

## 10. 边界（再写一次）

1. OnlineRFPerm / PO-risk **不是**事实核查，**不是**因果归因。  
2. 业务增量是 fire 期内 **对照代理**，不是 RCT。  
3. 回本用的是动作日成本假设；换你们财务口径时重跑 `07`/`08` 即可，不必改方法论。  
4. 周归因缺口必须先修到 ≈0，再对外讲周贡献。
