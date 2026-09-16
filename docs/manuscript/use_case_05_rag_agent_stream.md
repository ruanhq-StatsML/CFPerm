# Use-case 5 — Continuous serving streams (Graph-RAG / hybrid retrieval / agent next step)

Drop-in manuscript section. Same object as use-case 4: **one \(Y\), many \(X\), arrival windows**. Any stream whose data-generating map moves in continuous time can be written as \((X,Y,\mathrm{batch})\) and watched with the same two gates. Not attribution. Not a generation model.

The three scenes below are the same serving loop. Only the meaning of a hop changes.
The eight-facet table (smoothness, judge, Graph-RAG, hybrid, agent, synthetic gold, serving refresh, CUPED) is `docs/manuscript/online_serving_gates.tex`.

## Shared loop

Frozen serving policy \(f_{\mathrm{ref}}\) (prompt pack + retriever/graph/tools + generator config). Incoming traffic is batched. Probe error \(T_t=\mathrm{MSE}_t-E_{\mathrm{ref}}\) becomes a rank \(p\)-value; ADDIS/SAFFRON (async if windows overlap) decides fire. Last-two `hop_fires` is the adjacent-window gate, not the frozen-ref test.

```
load f_ref, E_ref, serving versions
burn-in (no decision)
for each batch:
    retrieve or act → generate next → align Y
    T = MSE - E_ref → p → online FDR
    quiet: keep gold; append T
    fire: stop treating current context/tools as gold
          → cheapest refresh that can recover
          → full reset of error pool, burn-in, FDR wealth
```

**Read-out.** Smooth stream: no fire. Fracture at a known cut: fire near onset; two consecutive fires are the operational trigger. After a refresh: reset, then quiet again. Hallucination rate, refusal rate, Recall@\(k\), and token cost are side gauges, not \(Y\).

---

## 5a. Graph-RAG serving refresh

**Business.** A Graph-RAG stack answers from a frozen knowledge graph: extract entities, pull a subgraph (local neighbours or global community summaries), pack them into a prompt, generate with citations. Corpus, community summaries, and prompt packs jump in versions. The question is *when that serving map is no longer the same object*, not whether a single chunk is relevant.

**Table.** One row per request, or one row per arrival window after aggregating graph features.

| | Graph-RAG |
|---|---|
| \(Y\) | citation lands on an edge / supporting nodes sit in the pack. Delayed labels are aligned by the time \(Y\) arrives |
| \(X\) | **window aggregate** of title-graph features: nodes, edges, mean degree, connected components, query seeds, seed fraction, LCC fraction |
| `batch` | arrival window. Overlapping community windows → async FDR with lag |
| \(T\) (optional) | \(0=\) local subgraph, \(1=\) global community. Do not average the two queues |

A hop of \(P(Y\mid X)\) is a graph pack / community recompute: drop or rewire edges, then aggregate again. Not single-node relevance. On fire: do not treat the current subgraph as gold; try refresh from cheap to expensive (prompt → retriever hops/community routing → re-extract graph / re-summarize → generator last). Promote a candidate only if a shadow copy of the *same* requests has lower window error **and** last-two is quiet; then reset.

**Prototype.** Quiet corpus + fixed graph: expect no fire. Injected fracture (rewired communities after a cut): fire at the cut. After a real refresh: new \(D_{\mathrm{ref}}\), subsequent windows quiet. Batch-aggregate tables: `scripts/prototype_graph_pack_batch_agg.py`.

**Hotpot on disk is a Graph-RAG snapshot, not a time stream.** Titles in the distractor pool are nodes; a shared-token edge is the cheap co-mention graph; seeds overlap the query; \(Y_j=1\) if title \(j\) is supporting. File: `xy_hotpot_pairs.csv` (query × 10). The `batch` column is query index, not arrival time. Use this table to learn the graph shape. Continuous-time refresh still needs a timestamped request log.

---

## 5b. Hybrid retrieval (vector + sparse + rerank)

**Business.** Production RAG is rarely one index. A request hits sparse (BM25 / keyword), dense (embedding \(k\)NN), optional graph or metadata filters, then a reranker, then the LLM. Indexes, embedding models, and rerankers refresh on different clocks. The object to monitor is whether *this mixture* still predicts a valid answer, not Recall@\(k\) of one retriever.

**Serving logic (one request).**

1. Query rewrite under prompt version \(v_{\mathrm{prompt}}\).
2. Parallel retrieve: sparse top-\(k_s\), dense top-\(k_d\), optional filters.
3. Fuse (RRF / weighted) → rerank → pack context.
4. Generate with citation offsets into the fused list.
5. Log \(X\) immediately; join \(Y\) when it lands.

**Table.**

| | Hybrid retrieval |
|---|---|
| \(Y\) | answer valid given the *fused* context (faithfulness, citation-in-list, task success). Not sparse-only or dense-only recall |
| \(X\) | query embedding; sparse/dense overlap; RRF ranks of cited chunks; rerank margin; \(k_s,k_d\); embedding-model id; reranker id; prompt id; filter-hit rate |
| `batch` | arrival window |
| \(T\) (optional) | \(0=\) sparse-heavy queries, \(1=\) dense-heavy, if the product already routes that way |

**What a hop is.** Embedding model swap, chunker change, reranker upgrade, or fusion weights drifting until the frozen pack no longer predicts \(Y\). Sparse still looking “fine” while dense/rerank jumped is exactly why \(Y\) must be on the fused answer, not on one channel.

**Next-step evaluation (the promotion test).** Fire is not auto-deploy.

- Live keeps the old hybrid pack. Shadow runs the candidate (new embedding, new \(k\), new rerank) on the same queries.
- Promote iff shadow window MSE is back near \(E_{\mathrm{ref}}\) *and* last-two does not fire on the first post-cut hops.
- Then full state reset. If only Recall@\(k\) moved and \(P(Y\mid X)\) did not recover, do not promote.

Cheapest refresh first: fusion weights / \(k\) → reranker → re-embed corpus → generator.

Hotpot’s 10-para pool is **not** this serving stream (no time order). It is the Graph-RAG snapshot in 5a.

---

## 5c. LLM / agent tool-trace next step

**Business.** An agent is a stream of hops, not one completion: observe → choose tool → call → observe → …. Each hop has a next-step map. The question is whether *this hop still predicts the next valid action*, not whether the final answer is factually true (that is a separate hallucination gauge).

**Serving logic (one hop).**

1. Pack observation \(o_t\) (user text, last tool payload, optional retrieved subgraph).
2. Policy emits tool name + arguments (or `stop`).
3. Tool returns \(o_{t+1}\) or a terminal answer.
4. When a hop label exists, write one row. Terminal success can be a second \(Y\) on a slower clock; do not collapse them.

**Table.**

| | Agent next step |
|---|---|
| \(Y_{\mathrm{hop}}\) | this hop is valid: schema-ok, tool admitted, arguments parse, non-loop. Available at hop end |
| \(Y_{\mathrm{task}}\) (optional, delayed) | episode success. Own stream; do not substitute for \(Y_{\mathrm{hop}}\) |
| \(X\) | observation embedding; last tool id; argument shape; hop index; retry count; retrieved-graph features if the agent is Graph-RAG-backed; prompt/tool-schema version |
| `batch` | hop index in the episode, or wall-clock window across episodes |
| \(T\) (optional) | tool family (search vs code vs DB) as two queues |

**What a hop is.** Tool schema change, broken auth, prompt-pack swap, or the policy starting to loop. On fire: do not auto-commit; do not write the trace into SFT/DPO gold; send high-error hops to review; optional freeze of the last memory write.

**Next-step evaluation.** The unit is the *next hop*, not the whole episode.

- Smooth traces (fixed tools, no schema cut): \(Y_{\mathrm{hop}}\) stream quiet.
- Fracture: flip tool names, permute argument keys, inject a loop: last-two should fire at that hop; frozen-ref fires within a short delay.
- After a tool-pack refresh: reset; following hops quiet on \(Y_{\mathrm{hop}}\).
- \(Y_{\mathrm{task}}\) may still fail (smooth wrong answers). That is hallucination/task eval, not this gate.

If the agent retrieves a graph or a hybrid index before the tool call, this scene *sits on top of* 5a/5b: retrieval fire means “do not trust context”; tool fire means “do not trust the next action.” Both can be on; they are not the same \(Y\).

---

## 中文

**凡是连续时间在变的数据，只要能写成 \((X,Y,\mathrm{batch})\)，都可以走这一套。** 冻住 serving 策略，看预测误差相对自己的历史池是否已经极端。Graph-RAG、混合检索、Agent 下一步是同一条流水线，只有 hop 的含义不同。

**Graph-RAG。** \(Y=\) 引用落在边上，或支撑节点在图包里。一个 serving batch 把窗内子图特征聚合起来（节点、边、度、连通片、query seed、最大片）。Fire = 图包 / community 换代（切点后 drop 或 rewire 边，再聚合），当前子图不当金标。不要盯单点相关性。Refresh 从模板 → 检索 → 重抽图，影子流量下一窗误差回来且 last-two quiet 才晋升，然后整池 reset。表形 prototype：`scripts/prototype_graph_pack_batch_agg.py`。八个面的总表：`docs/manuscript/online_serving_gates_zh.md`。

**混合检索。** 一次请求：改写 → 稀疏+稠密并行 → 融合/重排 → 生成。\(Y\) 必须打在 **融合后的答案** 上，不要打在单一通道 Recall。Hop = embedding / 切块 / 重排 / 融合权重把冻住的包用坏了。稀疏看起来还行、稠密已经跳，正是要盯融合 \(Y\) 的原因。下一步评估是影子对照 + last-two，不是 Recall@\(k\) 单独涨了就上线。

**Agent 下一步。** 每跳 observe → 选工具 → 调用。\(Y_{\mathrm{hop}}=\) 这一跳合不合法（schema、可调用、不打转），到跳结束就能写。任务成功是更慢的另一条 \(Y\)，不能替代。Fire = 工具协议/模板断了或开始 loop：停自动提交，轨迹不进金标。检索闸管「上下文还能不能信」，工具闸管「下一步还能不能做」；Graph-RAG agent 可以两闸叠上，不要合成一个分数。

幻觉率、拒绝率、单路 Recall 都不是这套 \(Y\)。读数仍是：光滑不火，断裂在 onset 火，refresh 后 reset 再 quiet。
