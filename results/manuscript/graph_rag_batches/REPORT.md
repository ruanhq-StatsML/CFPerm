# Graph-RAG batch aggregate prototype

One serving batch = mean/std of title-graph features. Native Y: supporting titles sit in seed ∪ 1-hop. Hop Y: supporting titles sit in the largest CC after an edge rewire (community 换代). Hotpot file order is not wall-clock time.

| Regime | y pre | y post | mean Δ | hop@cut | n hop fires |
|---|---:|---:|---:|---|---:|
| native | 0.516 | 0.541 | 0.003 | no | 0 |
| community 换代 | 0.516 | 0.392 | -0.033 | no | 2 |

## Last-two around the cut (community 换代)

| abs batch | fire | e_prev | e_now | ratio |
|---:|---|---:|---:|---:|
| 2 | no | 0.350 | 0.387 | 1.107 |
| 3 | no | 0.387 | 0.438 | 1.129 |
| 4 | no | 0.438 | 0.500 | 1.143 |
| 5 | no | 0.500 | 0.113 | 0.225 |
| 6 | no | 0.113 | 0.100 | 0.889 |

Rebuild: `PYTHONPATH=. python3 scripts/prototype_graph_pack_batch_agg.py`.
