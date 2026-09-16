# One map, two continuous-time jobs

The object is whether a frozen serving map \(P(Y\mid X)\) is still the same object.
What we actually run in continuous time is **moderation** and **agent reasoning**.
Graph-RAG and hybrid retrieval are ordinary serving maps under the same gates.
These are product use-cases. The LLM loop is one hop: observe the pack → take the next step → write \(Y\) → frozen \(f_{\mathrm{ref}}\) scores \(P(Y\mid X)\). Fire if that map has already moved. Nothing else.

Drop-in LaTeX: `docs/manuscript/online_serving_gates.tex`.
Chinese: `online_serving_gates_zh.md`, `online_serving_gates_zh.html`.

| Facet | \(Y\) | Fire means | Does not mean |
|---|---|---|---|
| **Continuous-time audit** | pass / fail | policy-pack or judge refresh | refusal rate; HH chosen |
| **Agent reasoning** | this hop / path is legal (schema, tool, no loop) | the chain broke, or the agent is spinning | hallucination rate; eventual task success |
| Graph-RAG pack | citation on an edge; supporting nodes in the pack | graph-pack or community refresh | single-node relevance |
| Hybrid retrieval | fused answer is valid | one channel or the fusion broke | single-channel Recall |

Graph-RAG: one batch aggregates subgraph features (or a user-id list first). Nothing else.

```bash
PYTHONPATH=. python3 scripts/run_serving_gates.py
```
