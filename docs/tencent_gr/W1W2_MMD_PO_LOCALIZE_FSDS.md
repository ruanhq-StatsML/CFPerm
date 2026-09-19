# Streamlined: Standardization → subset-MMD → FSDS

就三步：

```
1. StandardScaler(W1)
2. subset-level MMD²(item | W1 vs W2)
3. FSDS ranking   # StandardScaler → var → SelectKBest → HGB/LR
```

**定位：** 更快找到下钻归因的那部分客群 / subset 的方法论。

- **没有因果性** — 比较的是时段分布差异与排序，不声称 treatment effect
- **modality-agnostic** — 不绑图算法；输入是 **图谱特征** \(X\)，同一套 cmean 可换颗粒度下钻
- **\(K^\star\)** — 只是图谱特征上 cmean 收完的最终支撑（跑一次 FSDS 的子集），不是图算法核
- **超纲不做** — community / ego / ULS·PPR·GraphScan / GNN；就这样

**停止 / 是否下钻：** 见
[`ATTRIBUTION_DRILL_STOP.md`](./ATTRIBUTION_DRILL_STOP.md)
（父集同一套 cmean：列 \(\delta\) = feature guidance / 是否本层 FSDS；行 \(r_u\) = drill；早停得 \(K^\star\) 再 FSDS）。

**特征×实体联合（LaTeX prototype）：**
[`Feature_Entity_Joint_CMean.tex`](./Feature_Entity_Joint_CMean.tex)
— conditional mean 原语、\(D\) 矩阵、F→E / E→F、非 subgroup justification。

**定稿 LaTeX prototype（整条闭环）：**
[`Graph_Feature_CMean_Localization.tex`](./Graph_Feature_CMean_Localization.tex)
— 图谱特征 → cmean → \(K^\star\) → FSDS；scope 锁死；可 compile 的 formulation。

**PO-risk for DS（讲武德）：**
[`PO_RISK_FOR_DS.md`](./PO_RISK_FOR_DS.md)
— period-PO `mean(τ̂²)` / VIMP → FSDS ranking prior；`po_help_select` one-liner；**不是** ATE。

**多层次下钻 × dissolve rollup（落地试探）：**
[`DRILL_DISSOLVE_LANDABILITY.md`](./DRILL_DISSOLVE_LANDABILITY.md)
— 叶子分向上聚合（非 PyPI dissolve）；商户层与 MMD 有重叠，不替换官方 drill。

**多层图调研 × 多模态归因逻辑：**
[`MULTILAYER_GRAPH_AND_MM_ATTR.md`](./MULTILAYER_GRAPH_AND_MM_ATTR.md)
— 多层图算法族取舍；模态=列块 / key；接到同一套 cmean 缩支撑闭环。

**Graph localization：**
[`GRAPH_LOCALIZATION_NO_COMMUNITY.md`](./GRAPH_LOCALIZATION_NO_COMMUNITY.md)
— **本阶段超纲，不做**（含 ULS / PPR-Nibble / GraphScan）；主协议只要图谱特征 + cmean。

## 三步下钻（可选）：merchant → user → order

在漂移 subset 上再做定位（**每层应过停止表，勿强制钻穿**）：

```
Standardize → L1 merchant MMD → L2 user MMD → L3 order shift → FSDS
```

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
