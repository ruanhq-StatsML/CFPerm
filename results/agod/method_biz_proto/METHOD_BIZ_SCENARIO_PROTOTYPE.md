# 方法论 × 业务情景 × 单位经济 — 全面对照原型

> 目标：把 OnlineRFPerm / PO-risk / 风格轴 与「客服助手全面情景 + 单位经济敏感度」放在同一张图上，
> 讲清**异同**与**relevance**。价值不来自 AUROC，来自 **act-vs-ignore 量 × 单价**。

## 对外一句

方法论确认幻觉制度跳变（ratio≈11.7）后，客服助手相对不动作少 206.0 工单 / 130.0 退款 / 多承接 316.0；净贡献 ¥16,590，单价±20% 净带 ¥13,259–¥19,921；优先动作 model_rollback。

## 1. 三层叠乘（relevance 主链）

```text
[方法论] OnlineRFPerm fire / po_risk0 / RAG轴
    → 打开「制度跳变对照窗」+ 动作路由
[业务情景] hop弱/基/强 × 流量 × 动作拆分
    → 少工单 / 少退款 / 多承接（量）
[单位经济] 工单¥ / 退款¥ / 承接¥ （±20%扫）
    → 毛/净贡献¥带（钱）
```

**一句话 relevance：**  
OnlineRFPerm 决定「哪段日子算跳变对照窗」；hop 情景把跳变强度变成量；单位经济把量变成可对账的 ¥ 带——三层叠乘，缺一不可。

## 2. 异同矩阵

### 相同（same）

- 都对「制度/画像是否变了」敏感
- 都需要对照臂（quiet / fire_ignored / fire_acted），不能只看绝对值
- 动作路由都依赖轴拆分（concept vs covariate / RAG vs generation）

### 不同（different）

- 方法论输出：fired / ratio / AUROC / P@10 / style_auc（诊断）
- 业务情景输出：少工单 / 少退款 / 多承接 / 动作净¥（运营）
- 单位经济输出：同一量下的 ¥ 带（财务假设）
- 方法论不保证因果归因；业务账是对照代理，不是 RCT

### 对照表（层 × 对象 × 信号 × 同/异 × relevance）

| 层 | 对象 | 信号（本原型） | 与业务的同 | 与业务的异 | relevance |
|----|------|----------------|------------|------------|-----------|
| 方法论 · OnlineRFPerm | P(Y|X) / P(X) 是否进入新制度 | `fired=1, ratio≈11.71` | 给「要不要动作」打时间戳；不直接给钱 | 不产出工单/退款/¥；AUROC≠贡献 | fire → 打开 act-vs-ignore 对照窗（制度跳变期） |
| 方法论 · po_risk0 | 跳变下谁该进审计队列 | `P@10=1.0, auroc=1.0` | 决定审计人力花在哪 | 排序器不是事实核查器；P@10≠少退款数 | 驱动 audit_topk 动作臂 + 审计件成本（净贡献扣减） |
| 方法论 · RAG / router | 检索缺口 vs 生成制度 vs 路由漂移 | `rag_top10=1.0, domain_auc≈0.81` | 决定「刷检索」还是「回滚模型」 | 不解释 GMV 因果 | 路由到 retrieval_refresh / model_rollback 动作拆分 |
| 方法论 · style_domain_auc（对照面） | P(X) 文风/画风画像 | `style_auc≈0.999, judge_err_ratio≈3.67` | 也是制度/画像监控 | 客服助手主账不走 CTR/品牌路径；勿与幻觉 hop 混账 | 证明要拆轴：风格工单 ≠ 客服工单；防误伤偏好头 |
| 业务情景 · hop 强度 | 同一 fire 下跳变强弱 → 量差 | `weak/base/strong net = ¥9,156 / ¥16,590 / ¥17,990` | 仍用 act-vs-ignore 对照 | 把 ratio/幻觉率翻译成工单·退款·承接量 | 方法论 hop 强度的业务可感知带（量） |
| 业务情景 · 流量缩放 | 每千会话单价 × 日会话 | `prod_mid 日10,000 → 月毛 ¥3,569,143` | 线性外推，不改费率 | 方法层没有「会话量」概念 | 把子集原型外推到生产量级 |
| 单位经济 · 单价敏感度 | 量固定，扫工单/退款/承接单价 | `净带 ¥13,259 → ¥16,590 → ¥19,921` | 仍基于同一套少工单/退款/承接量 | 与 fire/AUROC 无关；财务假设层 | 回答「单价变了贡献还站不站得住」——方法层答不了 |

## 3. 方法论侧（HF 子集落地）

### 3.1 幻觉制度（HaluEval）— 客服助手主账来源

| 项 | 值 |
|----|-----|
| dataset | `pminervini/HaluEval@qa_samples` |
| cut fire | `1`，ratio≈`11.71` |
| 幻觉率 quiet → fire | `0.07` → `0.82` |
| po_risk0 P@10 / AUROC | `1.0` / `1.0` |
| Top-10 rag_hit | `1.0` |

读法（方法层）：regime: OnlineRFPerm fire at cut => hallucination label law hopped；ranking: po_risk0 sorts audit candidates under shift (not a fact checker)；rag: low rag_hit among Top-k => retrieval-support gap, separate from generation hop；router: task_shift guides which expert/route to refresh next

**与业务的衔接：** `fired+ratio` → seed 的 `fire_halluc/hop_ratio`；`rag_hit` → `retrieval_refresh` vs `model_rollback`；`P@10` → 审计成本进入净贡献。

### 3.2 文风/偏好轴（HH-RLHF）— 对照面（不要混进客服账）

| 项 | 值 |
|----|-----|
| style_domain_auc | `0.9988925576901218` |
| judge_err_ratio | `3.670675300647549` |
| concept fire | `1` |

读法：style_axis: style_domain_auc high => register/文风 portrait moved (P(X))；preference_axis: judge_err_ratio / OnlineRFPerm fire => preference map P(Y|X) moved；action: quiet->uniform; fire->audit Top-k by po_risk0; optional sqrt(PO) on T=1；not_a_claim: not a moral judge; not causal feature attribution

**relevance：** 高 style_auc + 有/无 concept-fire 决定「只调素材」还是「动偏好」——**钱走 CTR/品牌路径，不是客服工单路径**。本原型把它放在对照列，防止方法论信号被当成同一本账。

## 4. 业务情景侧（全面场景）

### 4.1 Before → After（base hop）

| 费率 | quiet | fire+不动作 | fire+动作 |
|------|-------|-------------|-----------|
| 工单率 | 1.5% | 15.86% | 1.14% |
| 承接率提升 | — | — | **22.57 pp** |

### 4.2 Hop 强度情景（方法论 ratio 的业务翻译）

| 情景 | 净贡献¥ |
|------|---------|
| weak_hop | ¥9,156 |
| base_hop | **¥16,590** |
| strong_hop | ¥17,990 |

量（base）：少工单 **206.0** / 少退款 **130.0** / 多承接 **316.0** → 毛 **¥16,656**，优先动作 **`model_rollback`**。

生产中等流量月毛外推：**¥3,569,143**。

### 4.3 周归因（对账）

| 周 | 少工单 | 少退款 | 毛¥ | 占比 |
|----|--------|--------|-----|------|
| 2026-09-14 | 86.1 | 57.0 | **¥7,189** | 43.2% |
| 2026-09-21 | 119.9 | 73.0 | **¥9,467** | 56.8% |

### 4.4 动作回本

| 动作 | 回本 | 分档 | 净ROI |
|------|------|------|-------|
| audit_topk | 0.009 天 | same_day_payback | 112.8x |
| retrieval_refresh | 0.017 天 | same_day_payback | 57.9x |
| model_rollback | 0.033 天 | same_day_payback | 29.4x |

## 5. 单位经济敏感度（财务层）

量固定（少工单/退款/承接不变），只扫单价：

| 情景 | 工单¥ | 退款¥ | 承接¥ | 毛¥ | 净¥ |
|------|-------|-------|-------|-----|-----|
| all_minus20 | ¥20.0 | ¥64.0 | ¥2.8 | ¥13,325 | ¥13,259 |
| all_plus20 | ¥30.0 | ¥96.0 | ¥4.2 | ¥19,987 | ¥19,921 |
| base | ¥25.0 | ¥80.0 | ¥3.5 | ¥16,656 | ¥16,590 |
| refund_plus20 | ¥25.0 | ¥96.0 | ¥3.5 | ¥18,736 | ¥18,670 |
| ticket_plus20 | ¥30.0 | ¥80.0 | ¥3.5 | ¥17,686 | ¥17,620 |

**净贡献带（全单价 ±20%）：** ¥13,259 → **¥16,590** → ¥19,921  
（带宽 / base ≈ 0.4）

### 与方法论的异同（专指单位经济）

| | 内容 |
|--|------|
| **同** | 仍站在同一套 act-vs-ignore 量上；不另编故事 |
| **异** | 完全不依赖 fire/AUROC；回答的是财务假设，不是检测力 |
| **relevance** | 方法层证明「该动作」；情景层证明「有多少量」；单位经济证明「换成你们财务口径还值不值」 |

## 6. 一张总图：谁回答什么问题

| 问题 | 该问哪一层 | 本原型答案 |
|------|------------|------------|
| 模型/标签制度跳了吗？ | 方法论 OnlineRFPerm | fire=1，ratio≈11.7 |
| 该审哪些样本？ | 方法论 po_risk0 | P@10=1.0 |
| 刷检索还是回滚？ | 方法论 RAG + 业务动作拆分 | 优先 `model_rollback` |
| 动作相对不动作少多少单？ | 业务情景 | 206.0 工单 / 130.0 退款 |
| 值多少钱？ | 单位经济 | 净 ¥16,590；带 ¥13,259–¥19,921 |
| 文风漂了要不要动客服账？ | 方法论风格轴（对照） | **不要混账**；走素材/CTR 路径 |

## 7. 怎么跑

```bash
# 1) 业务账 + 单位经济（若尚未跑）
PYTHONPATH=. python3 scripts/agod/run_biz_value_sql_demo.py
python3 scripts/agod/cs_assist_hop_sensitivity.py

# 2) 本对照原型
python3 scripts/agod/method_biz_scenario_prototype.py
# → docs/biz/METHOD_BIZ_SCENARIO_PROTOTYPE.md
# → results/agod/method_biz_proto/summary.json
```

## 8. 边界（必须写死）

1. OnlineRFPerm / PO-risk **不是**事实核查器，也**不是**因果归因。
2. 业务增量是 **fire 期内 acted vs ignored 的对照代理**，不是 RCT。
3. 单位经济带只扫单价；流量情景只缩放会话量——两者正交，不要合成一个「玄学 AUROC¥」。
