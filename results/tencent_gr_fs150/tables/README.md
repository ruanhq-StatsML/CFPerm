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
