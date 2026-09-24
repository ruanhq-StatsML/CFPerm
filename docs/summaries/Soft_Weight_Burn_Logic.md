# 软权重该不该烧 · 算力账（√ / ∛）

> 回答两件事：**(1) 软 IPTW 该不该烧？(2) 和算力账怎么分清？**

相关：[`RSI_PO_Efficiency_Pivot.md`](./RSI_PO_Efficiency_Pivot.md) ·
产物：`results/agod_po_power_eff/` · 代码：`agod/soft_burn.py`

---

## 1. 先把账拆开（算力账清楚）

| 账单 | 量级 | 谁决定 |
|---|---|---|
| **适应 FLOPs** | 大：RF / μ0 重拟合 | duty × n_control × fit（`po_eff`） |
| **加权 FLOPs** | 极小：\(w_i=\mathrm{PO}_i^\alpha\)，O(n) | α = ½ 或 ⅓ **几乎不动总账** |

所以：

> **选 √ 还是 ∛ ≠ 省算力问题。**  
> 真账单在 gate/refit；α 只改 IPTW 的风险形状。

---

## 2. 软权重该不该烧（决策树）

```text
gated_α 的 sig-MSE 比 uniform 更好？
        │
        ├─ 是 → BURN_α（√/∛ 取 rel 更低者）
        │
        └─ 否 → 还要走 reject 路径吗？
                  ├─ 否 → KEEP_UNIFORM（默认，5/6 pack）
                  └─ 是 → SOFTEN_ONLY：用 ∛ 减伤害（同 FLOPs）
```

| 决策码 | 含义 | 何时 |
|---|---|---|
| `KEEP_UNIFORM` | 不烧软权重 | 无证据 beat uniform |
| `BURN_SQRT` / `BURN_CBRT` | 烧，且 α 已选定 | gated rel&lt;1 |
| `SOFTEN_ONLY` | 不宣传 IPTW 胜利；若必须 gate → ∛ | gate 已开但伤 MSE |

**本库 6 pack 实证：** 多数 `SOFTEN_ONLY` / `KEEP_UNIFORM`；metro 偶发可 `BURN_SQRT`。

---

## 3. 和「排序 ≠ 下游」怎么拼

| 层 | 结论 |
|---|---|
| 排序 | refit 赢 rank_eff → **可以**为分流付适应 FLOPs |
| 下游 IPTW | 多数 mse_eff&lt;0 → **默认不烧** 软权重 |
| 软 α | 同 FLOPs；被迫 gate 时 ∛ ≤ √ 更常见 |

一句话：

> **算力可以花在 refit（为排序）；软权重默认不烧；非烧不可时选 ∛。**

---

## 4. 跑法

```bash
PYTHONPATH=. python3 scripts/run_soft_burn_scorecard.py \
  --summary results/agod_po_cbrt/summary.json \
  --out results/agod_soft_burn
```
