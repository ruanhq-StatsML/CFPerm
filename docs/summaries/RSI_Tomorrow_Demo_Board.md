# RSI 明日效果看板（高频探索）

> 节奏已从 **10 min → ~3 min**。方向 = PO/AGOD 适应效率（非广告案由、非 transfer-null）。  
> 分支：`cursor/rsi-model-stats-eff-abce`

---

## 一页结论（先看这个）

| # | 主张 | 证据 |
|---|---|---|
| 1 | **排序技能 ≠ IPTW 收益** | refit 赢 rank_eff 6/6，但多数 pack mse_eff&lt;0 |
| 2 | **低 duty 时 refit 更便宜** | mean duty≈0.24 ⇒ E[refit]/E[probe]≈0.24 |
| 3 | **√ vs ∛ 同 FLOPs** | 软权重是 bias–variance 旋钮；∛≤√ 多数 pack，但 **uniform 仍常最优** |

产物：
- [`../agod_po_eff/PO_EFF_SCORECARD.md`](../../results/agod_po_eff/PO_EFF_SCORECARD.md)
- [`../agod_po_power_eff/PO_POWER_EFF.md`](../../results/agod_po_power_eff/PO_POWER_EFF.md)（本轮）
- 账本：[`RSI_Iteration_Log.md`](./RSI_Iteration_Log.md)

---

## 高频 Protocol（3 min）

| 秒 | 动作 |
|---:|---|
| 0–20 | 抽 P2/P3/P4 或自修 |
| 20–120 | 落地 ≤1 模块 + 测试 |
| 120–160 | 刷新 scorecard / 明日看板一行 |
| 160–180 | commit · push · 等下一 tick |

下一优先：image-OOD 同尺 FLOPs · Grad-RFPerm freeze 并表 · 把三张卡合成一张 demo PNG。
