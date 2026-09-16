# Graph-RAG pool snapshot (Hotpot)

HotpotQA distractor validation. **Questions and wiki text are not stored.**

This is a **graph snapshot**, not a continuous-time stream. Validation queries have no arrival order; the `batch` column is the query index, not wall-clock time. Do not feed this table to OnlineRFPerm as if it were serving traffic.

## Shape: one query × 10 wiki titles

Almost every distractor example is a rectangular pool of 10 titles and exactly 2 supporting titles::

    titles                 (10,)     nodes
    y_pair                 (10,)     1 iff this title is supporting
    title-graph X          seed / degree / slot
    channel scores         BM25 / cosine / RRF on the same 10 nodes

`xy_hotpot_pairs.csv` flattens that: **1200 queries × 10 rows**. \(X\) is graph geometry on those titles (plus the two ranking channels used to pack the subgraph). Nothing else.

Print three examples (no paragraph text):

```bash
PYTHONPATH=. python3 scripts/prototype_hotpot_10para_shape.py
```

Collapsed one-row-per-query files (`xy_hotpot_hybrid.csv`, `xy_hotpot_graph.csv`) are optional readouts of the same pool. A time-ordered Graph-RAG serving gate still needs a real request log.

Rebuild:

```bash
python3 scripts/build_hybrid_retrieval_xy.py
```
