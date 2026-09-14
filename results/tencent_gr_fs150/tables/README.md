# 中间表

粒和 Spark 窗口见 `../EVAL_AND_SEQ.md` §7。

```
python3 scripts/tencent_gr/dump_block_tables.py --max-users 6000
```

| 文件 | 粒 | 行数（6000 user prefix） |
|---|---|---|
| `ev.parquet` | user × event | 423050 |
| `attr.parquet` | user × cnv（买前 asof） | 12866 |
| `post.parquet` | user × cnv（买后 + 满窗 y + `empty_any`/`lag_empty_any`） | 12866 |
| `user.parquet` | user × 1 | 6000 |
| `xy_left_right.parquet` | user × 1（左窗 X + 右窗 Y） | 6000 |

VIMP 管道用的矩阵：`xy_left_right.parquet`。Y = `y_ctr`（右窗 n_clk/n_exp，`n_exp<3` 为 NaN）、`y_n_clk`、`y_n_cnv`、`y_cvr`、`y_any_cnv`。X 含 `dec_hl7d_dec_clk`、`life_ctr`。
`python3 scripts/tencent_gr/feat_proto.py` 可重算。
