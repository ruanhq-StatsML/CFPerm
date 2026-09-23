# PO-risk → 加速后训练：模态侧重 + 难样本 up-weight

> 主文：[`AGOD_PO_Risk_PostTraining.tex`](AGOD_PO_Risk_PostTraining.tex) §two-logics · §fuse  
> 代码：`fuse_long_short` / `po_fuse` · `po_iptw_weights`（√PO）

**目标：** 加速后训练；Acc 为约束。PO 不替换 loss，只重分预算。

---

## 符合预期的两条逻辑

| | **A. 模态侧重**（high-residual 概念） | **B. 难样本 up-weight**（reject 后） |
|---|---|---|
| 问什么 | 哪座塔在扛 \(P(Y\mid X)\) 残差概念漂？ | reject 批里哪些行仍然难？ |
| 分辨率 | 列 / 模态 \(m\) | 行 / 观测 \(i\) |
| 信号 | \(\mathrm{PO}_m\) → 长短期融合 \(\alpha\) | \(\mathrm{PO}_i=\|Y-\mu\|\)（可混 batch PO） |
| 执行器 | freeze←\(L\)；step dump←\(S\)；LR←\(\alpha\) | 仅 gate reject 后 \(w\propto\sqrt{\mathrm{PO}}\) |
| 加速方式 | 预算砸高残差塔，冻低残差 BWD → FLOPs↓ / \(T(\mathrm{Acc}^\star)\)↓ | 难行多吃梯度；平静窗 \(w=1\) 不白抬 |

```text
每个窗口 t:   PO_m  --fuse--> α --> {freeze, steps, LR}     ← 模态侧重
若 reject_t:  PO_i  --√----> w --> Fit_{t+1}                ← 难样本 up-weight
```

---

## A. 模态侧重（具体）

- **High-residual** = 该模态上控制拟合解释不了当前 \(Y\)（概念/残差，不是纯 MMD）。  
- 长期 \(L\)：慢性高残差 → 谁进概念集合（freeze 防抖）。  
- 短期 \(S\)：\(\Delta\mathrm{PO}\) 尖峰 → 本窗 step 倾倒给谁。  
- 同一 wall-clock 不再 \(1/M\) 均分，瞄准瓶颈塔。

## B. Reject 后难样本 up-weight（具体）

- **先 gate，后抬权**：OnlineRFPerm/CFPerm reject 才 \(w_i\propto\sqrt{\mathrm{PO}_i}\)；calm 保持 uniform。  
- 用 **√PO（soft）** 而非 raw \(\propto\mathrm{PO}\)（后者易过拟合当前 reject 批）。  
- 加速点：有限 step 少浪费在极易行上，难区域更快被下一 fit 纠正。

## 两者关系

- 同源（PO），正交旋钮：A 改 **哪颗头** 反传；B 改 **哪些行** 加权。  
- 不要混：\(\mathrm{PO}_m\) 高 ≠ 自动抬行权（没 reject 不抬）；reject ≠ 自动冻模态。

长短期融合细节见原文 §fuse；落地 KPI 仍是 FLOPs / \(T(\mathrm{Acc}^\star)\)，Acc 过线才算。

---

## ROI mapping（只维护 SQL）

与业务 ROI 的对接：**一张 mapping SQL 即可**，其它看板可有可无。

文件：[`po_posttrain_roi_map.sql`](po_posttrain_roi_map.sql)

| logic_id | cost_proxy | benefit_proxy | view |
|---|---|---|---|
| `modality_emphasis` | `flops_rel`, wall_clock | \(T(\mathrm{Acc}^\star)\), ΔAcc 约束 | `v_roi_modality_emphasis` |
| `hard_upweight` | reject 后 fit 耗时 | next-MSE drop vs uniform, P@20% | `v_roi_hard_upweight` |

过线看 `v_roi_ship_gate`：A 要求 FLOPs&lt;1 且 ΔAcc≥−0.5%；B 要求 reject 窗 next-MSE 相对 uniform 下降。
