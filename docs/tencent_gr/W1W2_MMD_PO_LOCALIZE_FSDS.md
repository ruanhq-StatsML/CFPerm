# Streamlined: Standardization → subset-MMD → FSDS

就三步：

```
1. StandardScaler(W1)
2. subset-level MMD²(item | W1 vs W2)
3. FSDS ranking   # StandardScaler → var → SelectKBest → HGB/LR
```

其他（PO-risk / conditional-mean / SMD 堆指标）不在主路径里。

## Run
```bash
bash scripts/tencent_gr/run_behavior_shift_attribution.sh

# optional GT
bash scripts/tencent_gr/run_behavior_shift_attribution.sh \
  --gt-items path/to/gt_items.csv --gt-orders path/to/gt_orders.csv
```

## Outputs (`results/tencent_gr_standardize_mmd_fsds/`)
| file | 含义 |
|---|---|
| `item_mmd_scores.csv` | item-level MMD² |
| `localized_subset_items.csv` | MMD top-k subset |
| `fsds_feature_ranking.csv` | FSDS 特征 ranking |
| `standardize_mmd_fsds.png` | subset-MMD + FSDS 图 |
| `STANDARDIZE_MMD_FSDS_REPORT.md` | 报告 |
