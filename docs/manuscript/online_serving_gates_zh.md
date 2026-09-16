# 在线 serving 闸：两件连续时间的事

对象只有一个：**冻住的 serving 地图 \(P(Y\mid X)\) 还是不是同一个东西。**

真正按连续时间做的是两件：**审核**，和 **agent 推理**。
Graph-RAG 图包、混合检索同一套闸，只换 \(Y\)，按普通 serving 写，不再往虚处扩。

```
              ┌──────────────────────────────────────┐
              │  P(Y | X) 还是同一张地图？            │
              │  表：(X, Y, batch)                    │
              │  闸：冻参考窗  +  last-two hop_fires  │
              └──────────────────────────────────────┘
                     │
          ┌──────────┴──────────┐
     连续时间审核            Agent 推理
          │
     普通 serving：Graph-RAG / 混合检索（同一套闸，换 Y）
```

**两个闸，不可互换。** 冻参考窗问「相对 \(D_{\mathrm{ref}}\) 是否已经变差」；last-two 问「相邻两窗有没有 hop」。
Quiet 留金标；fire 换包再整池 reset。读数：**光滑不火；切点火；refresh 后 reset 再 quiet。**

| 面 | \(Y\) | Fire 表示什么 | 不表示什么 |
|---|---|---|---|
| **连续时间审核** | 过 / 不过 | 政策包或 judge 换代 | 拒绝率；HH chosen |
| **Agent 推理** | 这一跳 / 路径是否合法（schema、工具、不 loop） | 推理断了，或开始空转 | 幻觉率；任务最终成功 |
| Graph-RAG 子图 | 引用落在边上；支撑节点在图包里 | 图包或 community 换代 | 单点相关性 |
| 混合检索 | 融合后的答案成不成立 | 某一路索引或融合坏了 | 单路 Recall |

稿件表：`docs/manuscript/online_serving_gates_zh.tex`。

Graph-RAG 仍只需一句：一个 batch 把图特征聚合起来（user-id → list 则先聚合 list）。见 `docs/manuscript/graph_continuity_zh.md`。

```bash
PYTHONPATH=. python3 scripts/run_serving_gates.py
```

连续时间审核的表已在 `results/manuscript/llm_audit/`。
