# 行为 / 购买欲变动归因（TencentGR）

一条 concise 落地链路：

```
Standardize(W1)
  → FE(W1) / FE(W2)          # gap ≥ 30d，无泄漏
  → subset = cmean ⊕ MMD ⊕ PO-risk
  → 可视化漂移商品
  → FSDS ranking             # StandardScaler → var → SelectKBest → HGB/LR
  → (可选) GT 订单/商品 evaluator
```

**读法：** 两窗之间用户行为与购买欲变了 → 哪些商品 subset 在动 → 哪些图谱特征在归因。

## Run
```bash
PYTHONPATH=. python3 scripts/tencent_gr/run_w1w2_mmd_po_localize_fsds.py \
  --root data/tencent_subset --max-users 20000 --gap-days 30 --localize-k 200

# 有 ground-truth 时直接挂上：
PYTHONPATH=. python3 scripts/tencent_gr/run_w1w2_mmd_po_localize_fsds.py \
  --root data/tencent_subset --gt-items path/to/gt_items.csv --gt-orders path/to/gt_orders.csv
```

GT CSV 列名兼容：`item_id` / `oid` / `sku_id` / `goods_id`；订单表可带 `order_id`。

## Outputs (`results/tencent_gr_w1w2_mmd_po_fsds/`)
| file | 含义 |
|---|---|
| `localized_subset_items.csv` | 漂移商品 subset |
| `w1w2_mmd_po_localize_fsds.png` | 归因可视化 |
| `fsds_feature_ranking.csv` | 特征 ranking |
| `gt_eval.json` | GT hit/P/R（有 `--gt-*` 时） |
| `W1W2_MMD_PO_LOCALIZE_FSDS_REPORT.md` | 报告 |

## 备注
- 旧的 share-linear localize→FSDS 原型见 `LOCALIZE_FSDS_TIME.md`（参考用）
- 主路径就是本文件这条：**MMD / PO / cmean → subset → viz → FSDS**
