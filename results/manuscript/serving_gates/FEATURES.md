# Serving-gate features (on disk)

These `x_*` columns are the prototype. No prompts, no wiki text, no gold flags, no HH chosen.
Gates: `PYTHONPATH=. python3 scripts/run_serving_gates.py`.

## 审核 / judge

`results/manuscript/llm_audit/xy_hh_helpful_consistent.csv` · n=1200 · 13 X · live gate

Y = pass / fail. Same 13 X on BeaverTails / WildGuard / ToxicChat. HH chosen is not Y.

| column | meaning |
|---|---|
| `x_n_toks` | token count (scaled) |
| `x_n_chars` | character count (scaled) |
| `x_avg_word` | mean word length |
| `x_qmark` | question marks |
| `x_bang` | exclamation marks |
| `x_hedge` | hedge phrases (maybe / I think / …) |
| `x_formal` | formal connectives (therefore / however / …) |
| `x_i_count` | first-person I |
| `x_newlines` | newline count |
| `x_upper` | uppercase fraction |
| `x_refuse` | refusal phrases (I can't / I'm sorry / …) |
| `x_please` | please |
| `x_thank` | thank |

## 审核 / judge（多步）

`results/manuscript/llm_audit/xy_hh_multistep_consistent.csv` · n=1200 · 14 X · live gate

Y = pass / fail on this hop. One row per assistant turn. X is this hop's geometry. HH chosen is not Y. Episode success is not Y.

| column | meaning |
|---|---|
| `x_n_toks` | token count (scaled) |
| `x_n_chars` | character count (scaled) |
| `x_avg_word` | mean word length |
| `x_qmark` | question marks |
| `x_bang` | exclamation marks |
| `x_hedge` | hedge phrases (maybe / I think / …) |
| `x_formal` | formal connectives (therefore / however / …) |
| `x_i_count` | first-person I |
| `x_newlines` | newline count |
| `x_upper` | uppercase fraction |
| `x_refuse` | refusal phrases (I can't / I'm sorry / …) |
| `x_please` | please |
| `x_thank` | thank |
| `x_step` | hop index in the episode (scaled) |

## Graph-RAG 子图（一问一行）

`results/manuscript/graph_rag_batches/xy_graph_query.csv` · n=1200 · 7 X · live gate

Y = supporting titles sit in the served pack. The probe uses these 7. Community hop rewires edges and serves the LCC.

| column | meaning |
|---|---|
| `x_n_nodes` | title nodes in the pool |
| `x_n_edges` | shared-token co-mention edges |
| `x_mean_deg` | mean degree |
| `x_n_cc` | connected components |
| `x_n_q_seeds` | titles overlapping the query |
| `x_seed_frac` | seed fraction |
| `x_lcc_frac` | largest-component fraction |

## Graph-RAG 子图（一窗一行）

`results/manuscript/graph_rag_batches/xy_graph_batch.csv` · n=15 · 14 X · table shape only

Y = pack-usable rate on the window. Serving-table shape: mean and std of the 7 graph coordinates. Not the probe sample.

| column | meaning |
|---|---|
| `x_n_nodes` | title nodes in the pool |
| `x_n_edges` | shared-token co-mention edges |
| `x_mean_deg` | mean degree |
| `x_n_cc` | connected components |
| `x_n_q_seeds` | titles overlapping the query |
| `x_seed_frac` | seed fraction |
| `x_lcc_frac` | largest-component fraction |
| `x_n_nodes_std` | window std of x_n_nodes |
| `x_n_edges_std` | window std of x_n_edges |
| `x_mean_deg_std` | window std of x_mean_deg |
| `x_n_cc_std` | window std of x_n_cc |
| `x_n_q_seeds_std` | window std of x_n_q_seeds |
| `x_seed_frac_std` | window std of x_seed_frac |
| `x_lcc_frac_std` | window std of x_lcc_frac |

## 混合检索

`results/manuscript/hybrid_retrieval/xy_hotpot_hybrid.csv` · n=1200 · 14 X · live gate

Y = all gold titles in fused top-k. Channel-agreement geometry. Gold hits are not X. Hop flips the dense channel.

| column | meaning |
|---|---|
| `x_q_toks` | query token count |
| `x_q_chars` | query character count |
| `x_qmark` | query has '?' |
| `x_q_ents` | capitalized query tokens |
| `x_n_cand` | candidate pool size (10) |
| `x_js_overlap` | Jaccard of sparse vs dense top-k titles |
| `x_rank_corr` | Spearman of sparse vs dense ranks |
| `x_sparse_margin` | BM25 top1 − top2 |
| `x_dense_margin` | dense top1 − top2 |
| `x_rrf_top1_mass` | RRF mass on fused top-1 |
| `x_fuse_uniq` | unique titles in fused top-k |
| `x_mean_q_overlap` | mean query–doc token overlap in fused pack |
| `x_ks` | sparse k |
| `x_kd` | dense k |

## Graph-RAG 快照（一问 × 10 标题）

`results/manuscript/hybrid_retrieval/xy_hotpot_pairs.csv` · n=12000 · 9 X · table shape only

Y = this title is a supporting fact. Native Hotpot pool shape. batch = query index, not wall-clock time.

| column | meaning |
|---|---|
| `x_bm25` | BM25 z-score of this title |
| `x_dense` | char-ngram cosine z-score of this title |
| `x_rrf` | RRF z-score of this title |
| `x_rank_sp` | inverse sparse rank |
| `x_rank_de` | inverse dense rank |
| `x_q_overlap` | query–paragraph token overlap |
| `x_title_seed` | title tokens overlap the query |
| `x_title_deg` | co-mention degree |
| `x_slot` | slot index in the 10-title pool |

Not features: raw prompt/question, wiki paragraph text, gold/supporting flags, HH chosen/rejected, single-channel Recall, hallucination rate.
