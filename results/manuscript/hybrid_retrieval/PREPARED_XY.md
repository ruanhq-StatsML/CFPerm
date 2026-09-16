# Prediction table — hybrid retrieval / Graph-RAG

HotpotQA distractor validation. **Questions and wiki text are not stored.**

## Native shape: one query × 10 wiki paras

Almost every distractor example is already a rectangular pool of 10 titles/paragraphs
and exactly 2 supporting titles. That tensor is the serving object::

    titles, docs          (10,)
    sparse BM25           (10,)
    dense char-ngram cos  (10,)
    RRF                   (10,)
    y_pair                (10,)   1 iff this title is supporting
    X_pair                (10, 9) channel scores, not gold flags
    batch                 query id, repeated 10 times

`xy_hotpot_pairs.csv` is this tensor flattened: **1200 queries × 10 rows**.
Print three examples (no paragraph text):

```bash
PYTHONPATH=. python3 scripts/prototype_hotpot_10para_shape.py
```

## Collapsed query-level table (optional readout)

Same pool, one row per query: \(Y=1\) iff **all** gold titles landed in fused top-5.
That is a collapse of the (10,) ranking, not the native shape.

Graph-RAG extra \(X\) on the collapsed \(Y\): titles = nodes, shared-token edges, query-overlapping titles = seeds.

Hop: after query `i // 80 >= 4`, flip dense scores (embedding-pack swap). Pair \(Y\) stays the gold mask; pair \(X\) moves. Collapsed fused \(Y\) can change.

## Files

| File | n | shape | Y |
|---|---:|---|---|
| `xy_hotpot_pairs.csv` | 12000 | query × 10 paras | this para is supporting |
| `xy_hotpot_pairs_hop.csv` | 12000 | same, dense flipped after cut | same gold mask |
| `xy_hotpot_hybrid.csv` | 1200 | 1 row / query | all gold titles in fused top-5 |
| `xy_hotpot_hybrid_hop.csv` | 1200 | collapsed + dense fracture | fused coverage |
| `xy_hotpot_graph.csv` | 1200 | collapsed + title-graph X | fused coverage |
| `xy_hotpot_two_stream.csv` | 2400 | `T=0` sparse-only, `T=1` dense-only | channel usable |

Pair schema:

`y,batch,x_bm25,x_dense,x_rrf,x_rank_sp,x_rank_de,x_q_overlap,x_title_seed,x_title_deg,x_slot`

Rebuild:

```bash
PYTHONPATH=. python3 scripts/prototype_hotpot_10para_shape.py
python3 scripts/build_hybrid_retrieval_xy.py
```
