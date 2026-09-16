# Graph-RAG serving windows (batch-aggregated graph features)

Hotpot distractor validation. Questions and wiki text are not stored.

This is a **table-shape prototype**, not a wall-clock stream. Validation queries have no arrival order; window index follows file order.

## Query table

`xy_graph_query.csv`: one row per query. X is title-graph geometry. Y=1 iff every supporting title sits in the served pack (local: seed ∪ 1-hop).

## Batch table

`xy_graph_batch.csv`: **one row per serving window**. X is mean and std of the seven graph coordinates on that window. Y is the pack-usable rate. A graph-pack / community refresh rewires edges and serves the largest connected component, then aggregates again.

Hop files apply that community pack after `batch >= 4`.

Rebuild:

```bash
PYTHONPATH=. python3 scripts/prototype_graph_pack_batch_agg.py
```
