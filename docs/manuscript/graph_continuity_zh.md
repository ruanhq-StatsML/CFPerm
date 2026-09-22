# 图 continuity：user-id → list → 聚合，再上 OnlineRFPerm

这一面 **只走 OnlineRFPerm**（Sep 15 写稿 Algorithm 1）。其他面不必套这一套图编码器。

## 表怎么写

生产日志是：

```
user_id  →  这一窗里的图 list  [G1, G2, …]
         →  把这个 list aggregate 起来
```

每个到达窗 \(t\)、每个还在说话的 user：

\[
X_{u,t}
=
\bigl[\underbrace{\mathrm{mean}_{G\in L_{u,t}}\mathrm{GraphSAGE}(G)}_{\text{冻结的 mean-aggregator}}
\;\Vert\;
\underbrace{\mathrm{mean}_{G\in L_{u,t}} x_{\mathrm{hand}}(G)}_{\text{7 维人工图特征}}\bigr]
\]

\[
Y_{u,t}
=
\text{该 list 上图包可用率}
\quad
(\text{支撑节点落在 served pack 里})
\]

`batch` \(= t\)。探针打在 \((X_{u,t}, Y_{u,t}, \mathrm{batch})\) 上。再对 user 做一次均值，就是「每个 batch 一张图信息」的 serving 表形。

盘上 Hotpot **没有真 user-id**。Prototype 用 `query_index % 20` 只演示这条 log 形。文件顺序不是墙上时间。

## \(X\) 两截，拼起来就得了

| 截 | 维 | 是什么 | 不要做什么 |
|---|---:|---|---|
| GraphSAGE-mean | 16 | 节点特征 `(seed, deg, slot, in_lcc)`，两层 concat(self, mean-neighbors)，再 mean-pool 整张子图 | 不要在 trail 上继续训 SAGE。\(f_{\mathrm{ref}}\) 冻住，附录 B：在线更新会把 OOD 和模型漂移搅在一起 |
| 人工图特征 | 7 | `n_nodes, n_edges, mean_deg, n_cc, n_q_seeds, seed_frac, lcc_frac` | 不要再叠单点相关性、gold 标记、问句原文 |

Community / 图包换代：切点后 rewire 边，served pack 改成最大连通片。SAGE 邻居均值会动，人工连通片会动，\(Y\) 会动。这是 \(P(Y\mid X)\) 的 hop，不是「这个节点不相关」。

## OnlineRFPerm（只这一面）

对上面这张表：

1. \(D_{\mathrm{ref}}\) 上 fit 一次浅 RF，\(E_{\mathrm{ref}}\) 是参考窗 MSE。
2. 每窗 \(T_t=\mathrm{MSE}_t-E_{\mathrm{ref}}\)。
3. 秩 \(p_t=\#\{i<t: T_t\le T_i\}/t\)。小 \(p\) = 当前误差在历史池里已经极端。
4. 在线 FDR（ADDIS/SAFFRON；user 窗重叠则 async + lag）。Prototype 先画出秩 \(p\) 和 \(\alpha=0.05\)。
5. Last-two `hop_fires` 是相邻窗闸，和冻参考窗不是同一个问题。
6. Fire：当前图包不当金标；从最便宜的刷新做起（community 重切 / 再跑冻结 SAGE，而不是先微调大模型）；误差池、burn-in、FDR wealth **整池 reset**。

读数仍是：光滑不火；切点附近 \(T_t\) 抬起、\(p_t\) 变小；refresh 后 reset 再 quiet。

```bash
PYTHONPATH=. python3 scripts/prototype_graph_continuity.py
```

图：`results/manuscript/graph_continuity/pipeline.png`，`T_and_rank_p.png`。

## 另外几点建议（仅这一面）

1. **SAGE 冻在 \(D_{\mathrm{ref}}\)。** 切点后可以换图、换 pack，但不要把 SAGE 当在线 GNN 去 backprop。否则 \(T_t\) 分不清是图漂了还是编码器漂了。
2. **不要把 SAGE 和人工特征拆成两闸。** 拼成一个 \(X\)。两闸会把同一跳算两次。
3. **user 窗若滑动重叠，用 async ADDIS，lag = 重叠长度。** 这是稿子 §3.1，不是新方法。
4. **list 很长时用分位数聚合 MSE**（稿子 batch 大时 mean 不够）。现在每 user 每窗只有几跳，mean 够。
5. **\(Y\) 仍是图包可用，不是 embedding cosine。** cosine 是旁路。
6. **晋升：** 影子流量下一窗 \(T\) 回到 0 附近 **且** last-two quiet，再换包。
