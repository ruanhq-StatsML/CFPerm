# 三塔票的在线 stacking：类目切换、泄漏、快照 GLS

上一张三塔把 `(cat_match, s_pos, s_neg)` **每步独立** GLS 一遍。y 全是 1 时，oracle 轴上的快照 `π` 几乎不动。缺的是 stacking 协议本身：

```
冻住画像 且 冻住 π_t
votes ← (snap, item, pool)     # 不准看 y
honest_step(stacker, votes, y) # 先计分再 update
sticky π → conf-damped LR      # FWD 开着
然后才 Welford / seen
```

这不是 MMoE。MMoE 是 `π(x)=softmax(W mean(U,I,N))`，门看输入。这里 `π_m(t)` 看的是 expert 的 one-step-ahead 损失。

---

## 票怎么造（类目频率 / 价格，不是原始塔向量）

| expert | 票（冻住的） | 类目一切换会怎样 |
|---|---|---|
| **user** | `⟨cat_hist, onehot(pos.cat)⟩` | 掉到 0，直到画像追上 |
| **item** | 价格亲和（和类目无关的 item 坐标） | 同价位时还活着 |
| **neg** | `1 − cos(U,N)` | 负样本池硬度 |

原始 `{U,I,N}` 不能当 GLS 库（会在 item 顶点塌掉）。在线 stacking 吃的是这些标量票，预测 y。

本流：前 18 步 Tools，之后 Sports，价格带相同。user 票 0.94 → 0.29。

---

## 谁在跟踪切换

| 方法 | π_item 前 → 后 | path TV | 延迟（π_item≥0.4） |
|---|---|---|---|
| Hedge + Fixed-Share（诚实） | 0.18 → 0.25 | 2.26 | 28 |
| 离散 OSL | 0.35 → 0.00 | 4.67 | None |
| 每步 GLS 快照 | TV=0.149 | 几乎不走 | — |

半程换人花的是单纯形 TV。快照 GLS 的 `π_user` 钉在约 2/3，**不是** `π(t)`。离散 OSL 是顶点（乱跳）。Hedge+share 在切换处把质量从 user 挪走。

---

## 两种泄漏

1. **画像泄漏**（他们 `CausalPosNeg` 的反面）：先 `update` 再计分。冷启动第一步 `cat_match` 从 0 变成 1。类目切换时泄漏只有 `1/n`，不是主戏。
2. **stacker 泄漏**：先 `π ← update` 再计分。in-sample 乐观。本流诚实 preq **0.086**，leaky stacker **0.066**（切换后 0.141 vs 0.108）。

终步诚实 Hedge `π`：user=0.09, item=0.11, neg=0.80。sticky actuator TV **1.46** < 快 π TV **2.26**（两套时间尺度，FWD 仍开）。

```bash
PYTHONPATH=. python3 -m tests.test_agod_tower_stack
PYTHONPATH=. python3 scripts/run_agod_tower_stack.py
```
