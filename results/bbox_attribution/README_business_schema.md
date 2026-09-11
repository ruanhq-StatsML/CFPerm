# Business bbox feature schema (highly interpretable)

| 类别 | 特征 | 维度 |
|---|---|---:|
| 类别 | category histogram | 80 |
| 位置 | cx, cy mean/std | 4 |
| 尺寸 | w, h mean/std | 4 |
| 面积 | area mean/std | 2 |
| 宽高比 | aspect mean/std | 2 |
| 置信度 | score mean/std (GT placeholder) | 2 |
| 遮挡 | iscrowd 比例 | 1 |
| 空间分布 | densest quad + NN dist (+one-hot) | 6 |
| 空间关系 | pairwise IoU mean | 1 |
| ~~数量~~ | ~~count~~ | skipped |
| ~~RoI embedding~~ | ~~512~~ | omitted from business ranking |

Total named dims: **102**. RF AUC=0.600. mass image/text/bbox=0.525/0.449/0.026

## Top named ranks

| rank | feature | vimp |
|---:|---|---:|
| 1 | `cx_std` | 0.00501 |
| 2 | `aspect_std` | 0.00269 |
| 3 | `h_mean` | 0.00240 |
| 4 | `aspect_mean` | 0.00219 |
| 5 | `w_mean` | 0.00163 |
| 6 | `h_std` | 0.00152 |
| 7 | `area_mean` | 0.00136 |
| 8 | `cx_mean` | 0.00127 |
| 9 | `pairwise_iou_mean` | 0.00119 |
| 10 | `cat_hist_1` | 0.00108 |
| 11 | `w_std` | 0.00101 |
| 12 | `nn_center_dist_mean` | 0.00096 |
| 13 | `area_std` | 0.00093 |
| 14 | `cy_mean` | 0.00088 |
| 15 | `cy_std` | 0.00076 |
| 16 | `cat_hist_42` | 0.00019 |
| 17 | `iscrowd_ratio` | 0.00012 |
| 18 | `cat_hist_16` | 0.00010 |
| 19 | `cat_hist_15` | 0.00008 |
| 20 | `cat_hist_32` | 0.00007 |
