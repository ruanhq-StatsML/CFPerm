# 统计效率方法如何赋能 Tree-of-Thoughts（ToT）

> 把本栈已验证的尺子——**rank_eff / mse_eff / duty 算力账 / soft-burn / freeze Pareto**——
> 接到 ToT 的「展开 · 估值 · 剪枝 · 回溯」。  
> 不是另起一套 LLM 公式；是给 ToT **中间节点** 装上可检的统计门控。

相关：[`Soft_Weight_Burn_Logic.md`](./Soft_Weight_Burn_Logic.md) ·
[`RSI_PO_Efficiency_Pivot.md`](./RSI_PO_Efficiency_Pivot.md) ·
代码：`agod/tot_eff.py`

---

## 1. ToT 缺什么、我们有什么

经典 ToT（Yao et al.）：树搜索多条思维路径，用价值函数 \(V\) 评估中间 thought，再 BFS/DFS/beam 展开。

| ToT 组件 | 常见弱项 | 本栈统计赋能 |
|---|---|---|
| 分支生成 | 盲目扩宽 | **duty / n_control** 定预算 → 最多展开几叉 |
| 中间估值 \(V\) | LLM 自洽打分，不可复现 | **双目标**：\(V_{\mathrm{rank}}\) + \(V_{\mathrm{mse}}\)（排序≠下游） |
| 剪枝 | 阈值/启发式 | **burn 决策码** + Pareto（两轴都差就剪） |
| 算力 | 隐含 | **显式 FLOPs 账**：adapt vs weight；α 不占预算 |
| 终止 | 步数到了 | **SOFTEN_ONLY / KEEP_UNIFORM** 作安全叶 |

核心映射一句话：

> **Thought = 一条策略假设**（ref / probe / refit × √/∛/uniform × freeze…）；  
> **Value = 统计 scorecard**（不是模型自评）；  
> **Prune = 有条件 Pareto + soft-burn**。

---

## 2. 树长什么样（PO/AGOD 策略 ToT）

```text
root: OnlineRFPerm reject?
        │
        ├─ no ──► leaf KEEP_UNIFORM          [V_mse=baseline, flops≈0]
        │
        └─ yes ─┬─ adapt=ref     (flops≈ε)
                ├─ adapt=probe   (always-on)
                └─ adapt=refit   (duty·fit) ← 默认展开优先（rank_eff）
                        │
                        ├─ weight=uniform
                        ├─ weight=√     ── burn? ──► BURN_SQRT | prune
                        └─ weight=∛     ── burn? ──► BURN_CBRT | SOFTEN_ONLY
                                │
                                └─ (optional) freeze_early / freeze_low_share
                                         └── Pareto? keep : prune
```

每一节点携带统计状态：

```text
node = {
  adapt, alpha, freeze,
  rank_eff, mse_eff, rel_vs_uniform,
  expected_flops, burn_decision,
  children[]
}
```

---

## 3. 价值函数：为什么必须两维

ToT 若只用一个 \(V\)，会把「排得准」误当成「该加权」。

本栈强制：

\[
V_{\mathrm{rank}} = \mathrm{rank\_eff},\quad
V_{\mathrm{mse}} = \mathrm{mse\_eff},\quad
V_{\mathrm{cost}} = - \log(1+\mathrm{E[flops]})
\]

| 读法 | 动作 |
|---|---|
| \(V_{\mathrm{rank}}\) 高、\(V_{\mathrm{mse}}\) 低 | **展开分流枝**，**不**展开 IPTW 枝（排序≠下游） |
| 两者都高 | 可 `BURN_*` |
| 两者都低 / FLOPs 爆 | 剪枝或回退 `KEEP_UNIFORM` |
| 仅 FLOPs↓、MSE↑ | **假省钱**（synthetic freeze）→ 剪或标 `conditional_pareto` |

这就是赋能：ToT 的「哪个 thought 好」不再靠语感，而靠 **可复现的 scorecard**。

---

## 4. 剪枝规则（统计门控）

与 `soft_burn` / freeze 对齐：

1. **预算剪枝**：若 `E[flops] > budget` → 不展开（duty·n_control 超限）。  
2. **Burn 剪枝**：`rel_vs_uniform ≥ 1` → 不标 `BURN_*`；最多留 `SOFTEN_ONLY` 叶。  
3. **Pareto 剪枝**：相对 sibling，mse_rel≥1 且 flops_rel 未明显下降 → 剪。  
4. **Null/机会剪枝**（可接 transfer null）：excess≈0 的探针枝不展开。  
5. **校准剪枝**（可接 ECE）：高 excess + 高 ECE → 允许排序枝，禁止概率当权重。

---

## 5. 搜索算法（建议默认）

**Beam + 双目标**：

```text
beam_k = f(budget)           # 算力账定宽度
for depth in 1..D:
  for node in beam:
    expand legal children (adapt × α × freeze)
    score with (V_rank, V_mse, V_cost)
    prune by burn/Pareto/budget
  keep top-k on lexicographic:
    1) mse_eff > 0          # 下游不伤优先
    2) rank_eff             # 再比排序效率
    3) -expected_flops
```

这比「LLM 给 thought 打 1–10 分」更贴本栈证据：metro 可烧 √，其它 pack 多数 SOFTEN_ONLY。

---

## 6. 和 LLM-ToT 怎么接（赋能方式）

| 接法 | 做法 |
|---|---|
| **A. 策略 ToT（本仓默认）** | 节点 = 统计策略；LLM 只写解释文案，不参与 \(V\) |
| **B. LLM thought + 统计裁判** | LLM 提议「下一步试 refit+∛」；`soft_burn`/`po_eff` 判决 keep/prune |
| **C. 审出/广告卡 ToT** | S1/S2/S3 为枝；board metrics 为 \(V\)；软权重/限投永不由高 AUC 单独展开 |

推荐 **B**：发散用 LLM，收敛用统计——正好对应「排序可贵、下游要门控」。

---

## 7. 最小例子（本库 6 pack 读成 ToT）

```text
root ─reject─► adapt=refit          # V_rank 高
                 ├─ uniform         # 常为最优叶（5/6）
                 ├─ √ ─ metro only ► BURN_SQRT
                 └─ ∛ ─ else      ► SOFTEN_ONLY（不宣称胜利）
freeze 支路：electricity 进绿区留下；synthetic 标 conditional_pareto 后剪或降权
```

---

## 8. 不做什么

- 不用 ToT 给 image-OOD 的 PO 续命（已负结果）。  
- 不用单维 AUC/Spearman 当唯一 \(V\)。  
- 不用「更深的树」代替 duty 预算——深度也是算力账。

---

## 9. 一句话

> **ToT 负责想分叉；rank_eff / mse_eff / burn / duty 负责谁留下。**  
> 排序好只够把枝留在「分流」层；下游变好才允许烧软权重；算力由 duty×refit 记账，与 α 无关。
