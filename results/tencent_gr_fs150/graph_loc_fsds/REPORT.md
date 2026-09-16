# TencentGR：Graph localization → FSDS

数据：`user_feats_selected.parquet` n=6000 p_graph=91  
维大小：`{'user': 64, 'merchant': 11, 'product': 4, 'order': 12}`  Y=future_cnv  W=clock median

## 句式

```text
L1 图维 localization（user / merchant / product / order）
  → S* = Top-2
L2 仅在 S* 内 FSDS / LOGO / 特征排序
```

## Natural board（无注入）

- RF-domain AUC: **0.743**
- S*: **['product', 'user']**
- blended mass: `{'user': 0.336, 'merchant': 0.143, 'product': 0.435, 'order': 0.087}`
- LOGO share: `{'user': 0.0, 'merchant': 0.0, 'product': 1.0, 'order': 0.0}`

## Inject Hit@K（精确归因主考）

late 半窗平移目标维；`confound>0` 时其它维加干扰噪声。

| 注入维 | strength | confound | Graph S* | Flat Top-2 | H@2_g | H@2_f | H@1_g | H@1_f | rank_g | rank_f |
|--------|----------|----------|----------|------------|-------|-------|-------|-------|--------|--------|
| user | 0.8 | 0.0 | `['user', 'merchant']` | `['user', 'merchant']` | Y | Y | Y | Y | 1 | 1 |
| user | 0.8 | 0.35 | `['user', 'merchant']` | `['user', 'merchant']` | Y | Y | Y | Y | 1 | 1 |
| user | 1.5 | 0.0 | `['user', 'order']` | `['user', 'order']` | Y | Y | Y | Y | 1 | 1 |
| user | 1.5 | 0.35 | `['user', 'merchant']` | `['user', 'merchant']` | Y | Y | Y | Y | 1 | 1 |
| user | 2.5 | 0.0 | `['user', 'merchant']` | `['user', 'merchant']` | Y | Y | Y | Y | 1 | 1 |
| user | 2.5 | 0.35 | `['user', 'merchant']` | `['user', 'merchant']` | Y | Y | Y | Y | 1 | 1 |
| merchant | 0.8 | 0.0 | `['merchant', 'user']` | `['merchant', 'user']` | Y | Y | Y | Y | 1 | 1 |
| merchant | 0.8 | 0.35 | `['merchant', 'user']` | `['merchant', 'user']` | Y | Y | Y | Y | 1 | 1 |
| merchant | 1.5 | 0.0 | `['merchant', 'user']` | `['merchant', 'user']` | Y | Y | Y | Y | 1 | 1 |
| merchant | 1.5 | 0.35 | `['merchant', 'user']` | `['merchant', 'user']` | Y | Y | Y | Y | 1 | 1 |
| merchant | 2.5 | 0.0 | `['merchant', 'user']` | `['merchant', 'user']` | Y | Y | Y | Y | 1 | 1 |
| merchant | 2.5 | 0.35 | `['merchant', 'user']` | `['merchant', 'user']` | Y | Y | Y | Y | 1 | 1 |
| product | 0.8 | 0.0 | `['product', 'user']` | `['product', 'user']` | Y | Y | Y | Y | 1 | 1 |
| product | 0.8 | 0.35 | `['user', 'product']` | `['user', 'product']` | Y | Y | N | N | 2 | 2 |
| product | 1.5 | 0.0 | `['product', 'user']` | `['product', 'user']` | Y | Y | Y | Y | 1 | 1 |
| product | 1.5 | 0.35 | `['product', 'user']` | `['user', 'product']` | Y | Y | Y | N | 1 | 2 |
| product | 2.5 | 0.0 | `['product', 'user']` | `['product', 'user']` | Y | Y | Y | Y | 1 | 1 |
| product | 2.5 | 0.35 | `['product', 'user']` | `['product', 'user']` | Y | Y | Y | Y | 1 | 1 |
| order | 0.8 | 0.0 | `['order', 'user']` | `['order', 'user']` | Y | Y | Y | Y | 1 | 1 |
| order | 0.8 | 0.35 | `['order', 'user']` | `['order', 'user']` | Y | Y | Y | Y | 1 | 1 |
| order | 1.5 | 0.0 | `['order', 'user']` | `['order', 'user']` | Y | Y | Y | Y | 1 | 1 |
| order | 1.5 | 0.35 | `['order', 'user']` | `['order', 'user']` | Y | Y | Y | Y | 1 | 1 |
| order | 2.5 | 0.0 | `['order', 'user']` | `['order', 'user']` | Y | Y | Y | Y | 1 | 1 |
| order | 2.5 | 0.35 | `['order', 'user']` | `['order', 'user']` | Y | Y | Y | Y | 1 | 1 |

### Summary

| 指标 | 值 |
|------|-----|
| trials | 24 |
| **Hit@2 graph** | **1.000** |
| Hit@2 flat FSDS | 1.000 |
| **Hit@1 graph** | **0.958** |
| Hit@1 flat | 0.917 |
| mean rank graph / flat | 1.04 / 1.08 |
| hard (weak+confound) Hit@2 g/f | 1.0 / 1.0 |
| hard Hit@1 g/f | 0.75 / 0.75 |

对外一句：注入维 Hit@2（graph）= **100%**，Hit@1= **96%**
（flat @2/@1 = 100%/92%）；
自然窗 S*=['product', 'user']。

## 读法

- Hit@K 高 ⇒ 图维 L1 **能点到被污染的维**（维级精确归因）。
- hard 切片（弱注入+干扰）更接近线上噪窗。
- L2 只在 S* 内排特征；不宣称全空间唯一因果。

```bash
PYTHONPATH=. python3 scripts/tencent_gr/graph_loc_fsds.py
```
