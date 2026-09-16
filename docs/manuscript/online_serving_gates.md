# One map, many labels

The object is whether a frozen serving map \(P(Y\mid X)\) is still the same object
on arrival windows `batch`. Write every stream as \((X,Y,\mathrm{batch})\)
(optional queue \(T\in\{0,1\}\)). Two gates, not interchangeable: frozen-ref
(lstsq MSE vs \(D_{\mathrm{ref}}\)) and last-two `hop_fires`. Facets share the
gates and differ only in \(Y\).

Drop-in LaTeX: `docs/manuscript/online_serving_gates.tex`.
Chinese writeup and slides: `online_serving_gates_zh.md`, `online_serving_gates_zh.html`.

| Facet | \(Y\) | Fire means | Does not mean |
|---|---|---|---|
| Reasoning smoothness | this hop / path is legal | the chain broke | hallucination rate |
| Auditor / judge | pass / fail | policy-pack or judge refresh | refusal rate; HH chosen |
| Graph-RAG pack | citation lands on an edge; supporting nodes sit in the pack | graph-pack or community refresh | single-node relevance |
| Hybrid retrieval | the fused answer is valid | one channel, or the fusion, broke | single-channel Recall |
| Agent next step | schema-ok, tool admitted, no loop | the protocol broke, or the agent is spinning | eventual task success |
| Synthetic gold | still usable as an observation | the synthetic distribution drifted | “looks human-written” |
| Serving refresh | current prompt + retrieval + generator still predicts success | time to swap the pack | whether to fine-tune the LLM |
| Experimental CUPED | covariate adjustment is still valid | re-fit the regression | which A/B arm won |

Graph-RAG is the only facet that needs an extra sentence: **one serving batch is the aggregate of subgraph features** (nodes, edges, mean degree, CCs, query seeds, seed fraction, LCC). \(Y=\) supporting nodes in the served pack. A community refresh rewires edges and serves the largest CC, then aggregates again.

Feature catalog (every `x_*` already on disk): `docs/manuscript/serving_features_zh.md`, `docs/manuscript/serving_features.tex`.

```bash
PYTHONPATH=. python3 scripts/list_serving_features.py
PYTHONPATH=. python3 scripts/run_serving_gates.py
```

Live tables already on disk: auditor, Graph-RAG pack, hybrid fusion. Same two gates. Rebuild the Graph-RAG windows with `scripts/prototype_graph_pack_batch_agg.py` if those csvs are missing.
