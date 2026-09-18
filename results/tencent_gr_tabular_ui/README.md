# Tabular user / item features + localization drill-down

Main dims: **user_id**, **item_id**. No NetworkX / PageRank / Markov graph.

## Tables
- `user_features.*` — one row per user
- `item_features.parquet` — one row per item (funnel + first/last/linear credit shares)
- `edge_ui_features.parquet` — one row per (user, item)
- `localize_convert_path.csv` — drill-down: convert → path items + linear credit
- `localize_item_covisit.parquet` (+ `_sample.csv`) — item → top co-visited items (count)

meta: `{"n_users": 3000, "n_conv_users": 1697, "n_items": 194344, "n_edges_ui": 286949, "n_convert_path_rows": 90823, "co_window": 8, "top_covisit": 10, "main_dims": ["user_id", "item_id"], "note": "No graph/network algorithms \u2014 only groupby aggregations into tables."}`

## Top items by linear credit share
```
 item_id  n_users  n_cnv  share_linear  share_first  share_last  n_covisit_neighbors
13340556       82      8      0.000693     0.000589    0.001768                 1196
 6685422       28      3      0.000681     0.000589    0.001179                  412
 5612917        1      1      0.000589     0.000589    0.000589                    8
13789107        1      1      0.000589     0.000589    0.000589                    8
14732041        1      1      0.000589     0.000589    0.000589                    8
 4612288        2      1      0.000589     0.000589    0.000589                   24
15428572        1      1      0.000589     0.000589    0.000589                    8
15099889        1      1      0.000589     0.000589    0.000589                    8
 9076491       11      1      0.000589     0.000589    0.000589                  165
 5339541        1      1      0.000589     0.000589    0.000589                    8
```
