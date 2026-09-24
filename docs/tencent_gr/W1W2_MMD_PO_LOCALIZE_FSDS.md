# Streamlined: Standardization → subset-MMD → FSDS

就三步：

```
1. StandardScaler(W1)
2. subset-level MMD²(item | W1 vs W2)
3. FSDS ranking   # StandardScaler → var → SelectKBest → HGB/LR
```

**定位：** 更快找到下钻归因的那部分客群 / subset 的方法论。

- **没有因果性** — 比较的是时段分布差异与排序，不声称 treatment effect
- **modality-agnostic** — 不绑图、不绑文本/图像；同一套 procedure 可在多层级、多颗粒度下钻（如 merchant → user → order，或任意实体粒度）

## 三步下钻（可选）：merchant → user → order

在漂移 subset 上再做定位：

```
Standardize → L1 merchant MMD → L2 user MMD → L3 order shift → FSDS
```

> **算法 elaborate（四联图逐步读法 + Drill 停钻规则）：**  
> [`../summaries/Graph_Localize_Algorithm_Elaboration.md`](../summaries/Graph_Localize_Algorithm_Elaboration.md)  
> / [`../summaries/Graph_Localize_Algorithm_Elaboration.tex`](../summaries/Graph_Localize_Algorithm_Elaboration.tex)
>
> **Tip / 正负向算法层（无行动项）：**  
> [`../summaries/Tip_Sign_Algorithm_Elaboration.md`](../summaries/Tip_Sign_Algorithm_Elaboration.md)
>
> **Tip / 正负向与反欺诈排查布控动作：**  
> [`../summaries/AntiFraud_Tip_Sign_Actions.md`](../summaries/AntiFraud_Tip_Sign_Actions.md)
>
> **Scope 收窄·先走审出加速卡：**  
> [`../summaries/Scope_Review_Agent_Card.md`](../summaries/Scope_Review_Agent_Card.md)
>
> **10分钟 OOD 迭代 + 机会串：**  
> [`../summaries/TenMin_OOD_Iteration_Opportunities.md`](../summaries/TenMin_OOD_Iteration_Opportunities.md)
>
> **如何 Motivate Agent 做商业价值：**  
> [`../summaries/Motivate_Agents_Commercial_Value.md`](../summaries/Motivate_Agents_Commercial_Value.md)
>
> **20-Agent 商业落地 procedure（均分发力，非 math-Eval）：**  
> [`../summaries/Twenty_Agent_Commercial_Landing.md`](../summaries/Twenty_Agent_Commercial_Landing.md)
>
> **反欺诈评估–刻画模型（CT特征+图谱变动+监测波次）：**  
> [`../summaries/AntiFraud_Eval_Char_Model.md`](../summaries/AntiFraud_Eval_Char_Model.md)
>
> **监测波次 → Drill 团伙支撑 → 连续图谱变动线索：**  
> [`../summaries/Wave_Drill_Clue_Pipeline.md`](../summaries/Wave_Drill_Clue_Pipeline.md)
>
> **基于 Graph Shift 的方法（不谈社区发现；cmean+MMD+PO）：**  
> [`../summaries/Graph_Shift_Method_Elaboration.md`](../summaries/Graph_Shift_Method_Elaboration.md)
>
> **连续监测可做事项 + 数据接口：**  
> [`../summaries/Graph_CT_AD_Interface.md`](../summaries/Graph_CT_AD_Interface.md)
>
> **连续时间团伙发现（Leiden 换钥匙 / Drill 不换锁）：**  
> [`../summaries/Graph_CT_Gang_Discovery.md`](../summaries/Graph_CT_Gang_Discovery.md)

```bash
PYTHONPATH=. python3 scripts/tencent_gr/run_three_step_subset_localize.py \
  --root data/tencent_subset --max-users 20000 --gap-days 30
```

商户 id 默认用 `item_feat.122`（加密 shop/advertiser 代理）；订单 = 终端成功边 `user_item_ts`（TencentGR 上即 click）。

## Run（扁平 item-MMD）
```bash
bash scripts/tencent_gr/run_behavior_shift_attribution.sh
```

## Outputs
- 扁平：`results/tencent_gr_standardize_mmd_fsds/`
- 三层：`results/tencent_gr_three_step_localize/`

## Data subset (Git LFS)

https://github.com/ruanhq-StatsML/CFPerm/tree/cursor/graph-localize-elab-abce/data/tencent_subset  
详见 [`data/tencent_subset/README.md`](../../data/tencent_subset/README.md)。
