# 在线 serving 闸：一张地图，多种标签

对象只有一个：**冻住的 serving 地图 \(P(Y\mid X)\) 还是不是同一个东西。**
不是归因，不是生成模型，也不是把旁路指标当成 \(Y\)。

凡是连续时间在变、能写成 \((X,Y,\mathrm{batch})\) 的流，都走这一套。
可选再加队列 \(T\in\{0,1\}\)（两路审核、local/global、sparse/dense），不要把两路平均成一个分数。

## 架构

```
              ┌──────────────────────────────────────┐
              │  P(Y | X) 还是同一张地图？            │
              │  表：(X, Y, batch)  可选 T            │
              │  闸：冻参考窗  +  last-two hop_fires  │
              └──────────────────────────────────────┘
                    │  只换 Y（和对应的 serving 几何 X）
     ┌────────┬─────┴──────┬─────────┬──────────┐
  smoothness judge  Graph-RAG  hybrid/agent  CUPED / 金标 / 刷新
```

**两个闸，不可互换。**

| | 冻参考窗 | Last-two `hop_fires` |
|---|---|---|
| 探针 | 在 \(D_{\mathrm{ref}}\) 上 fit 一次（lstsq MSE） | 每跳在 \(B_{t-1}\) 上重 fit（浅 RF） |
| 分数 | \(s_t=\mathrm{MSE}(B_t)\)，\(\Delta_t=s_t-\mu_{\mathrm{ref}}\) | 相邻窗 OOS \(e_{\mathrm{now}}/e_{\mathrm{prev}}\) |
| Fire | 相对参考窗已经显著变差（秩 \(p\) + 在线 FDR） | 相邻两窗的地图 hop 了 |
| 不问的问题 | 是不是刚刚那一跳 hop | 是不是已经偏离参考窗 |

Quiet：金标保留，误差入池。
Fire：当前包不当金标；从最便宜的刷新做起；误差池、burn-in、FDR wealth 整池 reset。
影子包晋升：下一窗误差回到 \(E_{\mathrm{ref}}\) 附近 **且** last-two quiet。

读数仍是三句：**光滑不火；断裂在 onset 火；refresh 后 reset 再 quiet。**

## 八个面：同一套闸，换 \(Y\)

Fire 是 \(P(Y\mid X)\) 的 hop，不是质量分。旁路指标可以挂在仪表盘上，不要写进这张表的 \(Y\)。

| 面 | \(Y\) | Fire 表示什么 | 不表示什么 |
|---|---|---|---|
| 推理 smoothness | 这一跳 / 路径是否合法 | 推理链断了 | 幻觉率 |
| 审核 / judge | 过 / 不过 | 政策包或 judge 换代 | 拒绝率；HH chosen |
| Graph-RAG 子图 | 引用落在边上；支撑节点在图包里 | 图包或 community 换代 | 单点相关性 |
| 混合检索 | 融合后的答案成不成立 | 某一路索引或融合坏了 | 单路 Recall |
| Agent 下一步 | schema 合法、工具可调用、不 loop | 协议断了，或开始空转 | 任务最终成功 |
| 合成数据金标 | 还能不能当观测样本 | 合成分布漂了 | 「像不像人写的」 |
| Serving 刷新时刻 | 当前模板+检索+生成还能否预测成功 | 该换包了 | 大模型该不该微调 |
| 实验 CUPED | 协变量调整还有效 | 该重做回归 | A/B 谁赢了 |

稿件可粘贴的表：

- 英文：`docs/manuscript/online_serving_gates.tex`（`tab:serving-facets`）
- 中文：`docs/manuscript/online_serving_gates_zh.tex`（`tab:serving-facets-zh`）

其余面的 \(X\) 用各面自己的 serving 几何即可（回复特征、融合几何、工具 schema、协变量）。
**只有 Graph-RAG 这一面需要多写一句：一个 batch 把图上的这些特征聚合起来，然后就够了。**
盘上每一列 `x_*` 的清单：`docs/manuscript/serving_features_zh.md`。

## Graph-RAG：一个 batch = 图特征聚合

子图几何（每条请求、每个 10 节点题包）：

| 特征 | 含义 |
|---|---|
| `n_nodes` | 标题节点数 |
| `n_edges` | 共现边（标题共享 token） |
| `mean_deg` | 平均度 |
| `n_cc` | 连通片数 |
| `n_q_seeds` | 与问句 overlap 的 seed 数 |
| `seed_frac` | seed 比例 |
| `lcc_frac` | 最大连通片占比 |

到达窗 \(t\) 上把这些特征做均值（可选再加标准差），得到这一窗的 \(X_t\)。
\(Y_t=\) 窗内「支撑节点落在 seed \(\cup\) 一跳邻居」的比例
（生产里也可以是引用落在边上 / 下游任务成功，同一列）。

**图包或 community 换代：** 切点后把边丢掉或按同样密度重接，再做同样的聚合。
不要另造一套奇怪的 \(X\)。

盘上 Hotpot distractor 是题库快照，没有到达顺序。`batch` 在 snapshot 表里是题号；
下面 prototype 用文件顺序只演示 **表形**，不声称墙上时间。连续时间 serving 仍要带时间戳的请求日志。

## 现在就能跑

盘上已经有三面的 `(X,Y,batch)` 表。同一条命令、同一套闸：

```bash
PYTHONPATH=. python3 scripts/list_serving_features.py
PYTHONPATH=. python3 scripts/run_serving_gates.py
```

| 面 | quiet 表 | hop 表 |
|---|---|---|
| 审核 / judge | `llm_audit/xy_hh_helpful_consistent.csv` | `xy_hh_helpful_hop.csv` |
| Graph-RAG 子图 | `graph_rag_batches/xy_graph_query.csv` | `xy_graph_query_hop.csv` |
| 混合检索 | `hybrid_retrieval/xy_hotpot_hybrid.csv` | `xy_hotpot_hybrid_hop.csv` |

读数写在 `results/manuscript/serving_gates/`。Graph-RAG 的闸打在一问一行的 0/1 表上（按窗分组）；一窗一行的聚合表是 serving 表形，不是探针样本。

图包窗若要重做：

```bash
PYTHONPATH=. python3 scripts/prototype_graph_pack_batch_agg.py
```

Hop：切点后 rewire 边，改用最大连通片当 community 图包。演示稿：`docs/manuscript/online_serving_gates_zh.html`。

还没有 csv 的面（smoothness / agent / 合成金标 / serving 刷新 / CUPED）同一套闸，换表即可。

