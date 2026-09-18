# 图谱特征 → 网格面板（无图算法）

主维：**user_id / item_id**。
图谱特征 = 漏斗 / 触点 credit / 共现强度 / 活跃度 rank 等聚合列。
网格 = `(user, item)` 行，左连 user 特征、右连 item 特征。

## 主产出
- `feature_grid.parquet` — 建模网格（label=`y_convert`）
- `feature_grid_sample.csv` / `feature_cols.json`
- `user_features.*` / `item_features.parquet` / `edge_ui_features.parquet`
- `localize_convert_path.csv` — convert 路径下钻
- `localize_item_covisit.parquet` — 共现 top-k 下钻（count only）

- grid rows: **286949**, feature cols: **48**

## Top items by linear credit share
```
 item_id  n_users  n_cnv  share_linear  share_first  share_last  n_covisit_neighbors  item_credit_rank
13340556       82      8      0.000693     0.000589    0.001768                 1196               1.0
 6685422       28      3      0.000681     0.000589    0.001179                  412               2.0
 5612917        1      1      0.000589     0.000589    0.000589                    8               7.0
13789107        1      1      0.000589     0.000589    0.000589                    8               7.0
14732041        1      1      0.000589     0.000589    0.000589                    8               7.0
 4612288        2      1      0.000589     0.000589    0.000589                   24               7.0
15428572        1      1      0.000589     0.000589    0.000589                    8               7.0
15099889        1      1      0.000589     0.000589    0.000589                    8               7.0
 9076491       11      1      0.000589     0.000589    0.000589                  165               7.0
 5339541        1      1      0.000589     0.000589    0.000589                    8               7.0
```
