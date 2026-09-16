# One map, two continuous-time jobs

The object is whether a frozen serving map \(P(Y\mid X)\) is still the same object.
What we actually run in continuous time is **moderation** and **agent reasoning**.
Graph-RAG and hybrid retrieval are ordinary serving maps under the same gates.
These are product use-cases. One hop:

```
1. observe the current serving pack
2. model takes the next step (reply / tool / cite subgraph / fuse channels)
3. write Y for that hop   (X is serving geometry, not raw text)
4. frozen f_ref predicts P(Y | X)
5. T = MSE - E_ref → rank p → online FDR
   quiet: this trajectory is gold
   fire:  current pack is not gold → cheapest swap → reset
```

Nothing else. Audit can be multi-step: one row per assistant hop, \(Y\) still pass/fail for that hop. Episode success is not \(Y\).

Drop-in LaTeX: `docs/manuscript/online_serving_gates.tex`.
Chinese: `online_serving_gates_zh.md`, `online_serving_gates_zh.html`.

| Facet | \(Y\) | Fire means | Does not mean |
|---|---|---|---|
| **Continuous-time audit** | pass / fail | policy-pack or judge refresh | refusal rate; HH chosen |
| **Agent reasoning** | this hop / path is legal (schema, tool, no loop) | the chain broke, or the agent is spinning | hallucination rate; eventual task success |
| Graph-RAG pack | citation on an edge; supporting nodes in the pack | graph-pack or community refresh | single-node relevance |
| Hybrid retrieval | fused answer is valid | one channel or the fusion broke | single-channel Recall |

Graph-RAG: one batch aggregates subgraph features (or a user-id list first). Nothing else.

The other object is group attribution: \(T\) is the queue / channel / hop bucket. Permute \(T\), not \(X\). See `docs/manuscript/group_attribution_zh.md`.

```bash
PYTHONPATH=. python3 scripts/run_serving_gates.py
PYTHONPATH=. python3 scripts/prototype_group_attribution.py
```
