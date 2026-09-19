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

**停止 / 是否下钻：** 见
[`ATTRIBUTION_DRILL_STOP.md`](./ATTRIBUTION_DRILL_STOP.md)
（父集同一套 cmean：列 \(\delta\) = feature guidance / 是否本层 FSDS；行 \(r_u\) = drill；早停得 \(K^\star\) 再 FSDS）。

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
