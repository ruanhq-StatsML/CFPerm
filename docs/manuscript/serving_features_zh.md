# 盘上的 \(X\)：都可以直接 prototype

闸只读 `x_*`。问句、wiki 正文、gold 标记、HH chosen/rejected **都不是特征**。

```bash
PYTHONPATH=. python3 scripts/list_serving_features.py
PYTHONPATH=. python3 scripts/run_serving_gates.py
```

## 审核 / judge（13）

表：`results/manuscript/llm_audit/xy_*.csv`（HH / BeaverTails / WildGuard / ToxicChat 同一套 \(X\)）。\(Y=\) 过/不过。

| 列 | 含义 |
|---|---|
| `x_n_toks` | token 数 |
| `x_n_chars` | 字符数 |
| `x_avg_word` | 均词长 |
| `x_qmark` | 问号 |
| `x_bang` | 感叹号 |
| `x_hedge` | 含糊（maybe / I think / …） |
| `x_formal` | 正式连接词 |
| `x_i_count` | 第一人称 I |
| `x_newlines` | 换行 |
| `x_upper` | 大写比例 |
| `x_refuse` | 拒绝套话 |
| `x_please` | please |
| `x_thank` | thank |

多步表 `xy_hh_multistep_*.csv` 多一列 `x_step`（这一跳在线程里的位置）。一行一跳，Y 仍是过/不过。整段成功不是 Y。

## Graph-RAG 子图（7，窗上再聚合）

一问一行：`graph_rag_batches/xy_graph_query.csv`。闸打在这张 0/1 表上。

| 列 | 含义 |
|---|---|
| `x_n_nodes` | 标题节点 |
| `x_n_edges` | 共现边 |
| `x_mean_deg` | 平均度 |
| `x_n_cc` | 连通片 |
| `x_n_q_seeds` | 与问句 overlap 的 seed |
| `x_seed_frac` | seed 比例 |
| `x_lcc_frac` | 最大片占比 |

一窗一行：`xy_graph_batch.csv`，上述 7 个的均值 + 标准差（14 列）。这是 serving 表形，不是探针样本。

## 混合检索（14）

`hybrid_retrieval/xy_hotpot_hybrid.csv`。\(Y=\) 融合 top-\(k\) 盖住全部支撑标题。

| 列 | 含义 |
|---|---|
| `x_q_toks` / `x_q_chars` / `x_qmark` / `x_q_ents` | 问句大小与实体 |
| `x_n_cand` | 候选池大小（10） |
| `x_js_overlap` | 稀疏 vs 稠密 top-\(k\) Jaccard |
| `x_rank_corr` | 两路秩相关 |
| `x_sparse_margin` / `x_dense_margin` | 各路 top1−top2 |
| `x_rrf_top1_mass` | 融合 top-1 的 RRF 质量 |
| `x_fuse_uniq` | 融合 top-\(k\) 去重比例 |
| `x_mean_q_overlap` | 融合包与问句 overlap |
| `x_ks` / `x_kd` | 稀疏/稠密 \(k\) |

## 一问 × 10 标题快照（9）

`xy_hotpot_pairs.csv`。\(Y_j=\) 这个标题是不是支撑事实。`batch` 是题号。

| 列 | 含义 |
|---|---|
| `x_bm25` / `x_dense` / `x_rrf` | 三路分数（z-score） |
| `x_rank_sp` / `x_rank_de` | 逆秩 |
| `x_q_overlap` | 问句–段落 overlap |
| `x_title_seed` | 标题是不是 query seed |
| `x_title_deg` | 共现度 |
| `x_slot` | 10 槽位下标 |

稿件短表：`docs/manuscript/serving_features.tex`。
