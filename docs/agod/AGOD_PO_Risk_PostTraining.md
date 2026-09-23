# PO-risk Metrics → 加速后训练（长短期融合 · 概念模态侧重）

> 主文（LaTeX）：[`AGOD_PO_Risk_PostTraining.tex`](AGOD_PO_Risk_PostTraining.tex) §\ref{sec:fuse}  
> 代码：`fuse_long_short` · `po_fuse` · `next_step_actuators_fused`（`agod/po_risk_train.py`）

**目标：** 加速后训练（FLOPs / step / wall-clock / T(Acc★)）；Acc 为约束。

---

## 关键发现：概念模态侧重要拆长/短

| 轨道 | 定义 | 管什么执行器 |
|---|---|---|
| **Long** \(L_m\) | \(\mathrm{EMA}(\mathrm{PO})\cdot(1+\mathrm{proto})+\eta\alpha_{t-1}\) | **freeze**（慢集合，防抖） |
| **Short** \(S_m\) | \(\Delta\mathrm{PO}\)（或 PO−EMA） | **step dump**（尖峰立刻给预算） |
| **Fused** \(\alpha\) | \(\omega_L z(L)+\omega_S z(S)-\lambda\mathrm{MMD}\) | **LR / stack prior** |

纯短期 → freeze 乱抖；纯长期 → 尖峰模态来不及倾倒 step → \(T(\mathrm{Acc}^\star)\) 下不来。

自适应：\(\mathrm{spike\_ratio}=\|S\|_\infty/\mathrm{MAD}(L)\)；平静偏 long，尖峰抬 \(\omega_S\)。

---

## Elaboration strategy（怎么迭代）

1. **先打点** — 每窗记 \((L,S,\omega,\mathrm{top\_concept},\mathrm{top\_spike},\mathrm{freeze},s,\lambda,\mathrm{FLOPs},\mathrm{Acc})\)  
2. **四路消融** — long-only / short-only / 单 Softmax 驱动全部旋钮 / **分层**（推荐）  
3. **调参顺序** — 先 freeze 分位（定 FLOPs）→ \(\omega,g\) → EMA \(\rho\) → \(\tau\) → 最后 \(\varepsilon\)  
4. **工况剧本** — 平静保 long；尖峰倾倒 short；概念交接跟 \(L\) 不追 \(S\)；\(M{=}2\) 饱和不冻  
5. **过线** — 分层在同 Acc 约束下 FLOPs 或 \(T(\mathrm{Acc}^\star)\) 优于单 α；freeze Jaccard 别因 \(S\) 乱抖  

默认：\(M\ge3\) 用不对称包上 `po_fuse`；只砍 FLOPs 用 `po_gated`；不能冻用 `po_budget`。

---

## 落地 P0（加速）

```text
Sense PO(/proto/MMD) → fuse_long_short → freeze←L, steps←S, LR←α
→ 报 FLOPs↓ 且 ΔAcc≥−ε
```
