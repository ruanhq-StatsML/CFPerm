# Prediction table — hybrid retrieval / Graph-RAG

HotpotQA distractor validation. **Questions are not stored.** One row per query.

## How the row is made

Candidate pool is the example’s ~10 Wikipedia paras (already in Hotpot, no separate index).

1. **Sparse channel:** BM25(query, para text).
2. **Dense channel:** char-ngram hashing cosine (a second view, not a GPU embedder).
3. **Fuse:** RRF, \(k_s=k_d=k_{\mathrm{fuse}}=5\).
4. **\(Y_{\mathrm{fused}}=1\)** iff every gold supporting **title** is in the fused top-5. Not BM25 Recall, not EM.
5. **\(X\):** query geometry + how much the two channels agree. **Gold hit flags are not \(X\).**

Graph-RAG extra \(X\) on the same \(Y\): titles = nodes, shared-token edges, query-overlapping titles = seeds.

Hop overlay (`xy_hotpot_hybrid_hop.csv`): after `batch>=4`, flip the dense scores (embedding-pack swap). Pre-cut rows match native.

## Files

| File | n | Y | X |
|---|---:|---|---|
| `xy_hotpot_hybrid.csv` | 1200 | fused gold-title coverage | hybrid geometry |
| `xy_hotpot_hybrid_hop.csv` | 1200 | same, dense flipped after cut | recomputed after the flip |
| `xy_hotpot_graph.csv` | 1200 | same fused Y | hybrid + title-graph geometry |
| `xy_hotpot_two_stream.csv` | 2400 | `T=0` sparse-only usable, `T=1` dense-only usable | same hybrid X |

Schema of `xy_hotpot_hybrid.csv`:

`y,batch,x_q_toks,x_q_chars,x_qmark,x_q_ents,x_n_cand,x_js_overlap,x_rank_corr,x_sparse_margin,x_dense_margin,x_rrf_top1_mass,x_fuse_uniq,x_mean_q_overlap,x_ks,x_kd`

Rebuild (parquet under `data/hf_cache/retrieval/`, not committed):

```bash
python3 scripts/build_hybrid_retrieval_xy.py
```
