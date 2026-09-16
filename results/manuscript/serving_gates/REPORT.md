# Serving gates — same map, tables already on disk

One command. Two gates. Different `Y`. Fire is a hop of `P(Y|X)`, not a quality score.

```bash
PYTHONPATH=. python3 scripts/run_serving_gates.py
```

| Facet | Stream | Y | regime | pass | mean Δ | boot fires | hop@cut | n hop |
|---|---|---|---|---:|---:|---:|---|---:|
| 审核 / judge | HH helpful, consistent | 过 / 不过 | quiet | 0.491 | 0.037 | 1 | no | 0 |
| 审核 / judge | HH helpful, policy hop | 过 / 不过 | hop | 0.516 | 0.279 | 10 | yes | 2 |
| 审核 / judge | WildGuard native is_unharmful | 过 / 不过 | quiet | 0.917 | -0.027 | 0 | no | 0 |
| 审核 / judge | HH multi-step, consistent | 过 / 不过（这一跳） | quiet | 0.393 | 0.026 | 0 | no | 1 |
| 审核 / judge | HH multi-step, policy hop | 过 / 不过（这一跳） | hop | 0.527 | 0.605 | 10 | no | 3 |
| Graph-RAG 子图 | local pack (seed ∪ 1-hop) | 支撑节点在图包里 | quiet | 0.534 | 0.006 | 0 | no | 0 |
| Graph-RAG 子图 | rewire + largest-CC pack | 支撑节点在图包里 | hop | 0.425 | -0.029 | 0 | no | 2 |
| 混合检索 | RRF fused gold-in-topk | 融合后的答案成不成立 | quiet | 0.556 | 0.007 | 0 | no | 0 |
| 混合检索 | dense channel flipped after cut | 融合后的答案成不成立 | hop | 0.482 | 0.019 | 5 | no | 0 |

Quiet streams should stay near Δ = 0. Labeled hops should lift Δ and, when the adjacent-window ratio clears γ, fire last-two at the cut.
Multi-step audit is the same map: one row per assistant hop, Y = pass/fail for that hop, not the whole thread.
Graph-RAG and hybrid Hotpot tables are pool snapshots: `batch` is file/window index, not wall-clock time.
