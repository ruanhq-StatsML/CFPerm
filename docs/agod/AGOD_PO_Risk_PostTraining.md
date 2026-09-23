# PO-risk → ROI：下一步侧重哪一模态 / 哪一部分

> SQL（唯一要维护的 mapping）：[`po_posttrain_roi_map.sql`](po_posttrain_roi_map.sql)  
> 代码：`fuse_long_short` / `po_fuse` · `po_iptw_weights`

---

## 可以这样理解吗？—— **可以，但是「subset」有两层**

| 理解 | 对不对 | 具体是什么 subset |
|---|---|---|
| 限制 algorithm **update focus** 到某一块 | **对** | A：模态/塔子集；B：reject 批内难行子集 |
| 把整个 dataset **永久丢掉**只留一个 partition 再训 | **不对** | 其它模态 FWD 仍在；calm 行仍进 batch，只是 \(w=1\) |

一句话：**是「更新时侧重某个 subset」，不是「数据集切片后扔掉剩下的」。**

---

## 下一步应侧重什么（直接 map 到 ROI）

### A. 下一步侧重哪一个**模态**（`focus_id=next_modality`）

```text
L_m = EMA(PO_m)·(1+proto)     → 慢性高残差概念塔
S_m = ΔPO_m                   → 本窗尖峰
active = {m: L_m ≥ q30(L)}    → 模态 subset（谁还配拿 BWD）
next_chronic = argmax L       → 概念上该强调谁
next_spike   = argmax S       → 本窗 step 该砸给谁
ROI_A = flops_saved / acc_risk   (ΔAcc≥−0.5% 才算过线)
```

| 决策输出 | ROI 维 |
|---|---|
| `top_concept_mod` / `top_spike_mod` | 侧重对象（可审计） |
| `modality_subset_frac` = n_active/M | **成本**变小 → FLOPs↓ |
| `delta_acc` | **约束** |
| `t_to_acc_star` | **收益**（更快够到 Acc★） |

### B. 下一步强调数据的哪一**部分**（`focus_id=next_hard_rows`）

```text
仅当 rejected=1:
  PO_i = |Y−μ|
  hard subset ≈ top-20% by PO_i
  w_i ∝ √PO_i                 → 难行多吃梯度
ROI_B = next_mse_drop / fit_wall_clock
```

| 决策输出 | ROI 维 |
|---|---|
| `hard_row_subset_frac` | 强调了多大一块行 |
| `next_mse_drop` vs uniform | **收益** |
| `fit_wall_clock` | **成本** |

### 合在一起（`joint_focus`）

```text
每个 t:     先定模态 subset（A）→ 写进 ROI_A
若 reject:  再定难行 subset（B）→ 写进 ROI_B
正交：A 改哪颗头；B 改哪些行。过线看 v_roi_ship_gate。
```

---

## SQL 怎么读（efficient）

```sql
-- 本窗：下一步侧重谁？
SELECT window_t, next_modality_chronic, next_modality_spike,
       modality_subset, hard_row_subset_frac,
       roi_pass_modality, roi_pass_hard_rows
FROM v_roi_next_focus;

-- 跑级 ROI + 主导模态
SELECT * FROM v_roi_ship_gate;

-- subset 语义（防误解）
SELECT * FROM v_roi_subset_meaning;
```

其它看板可有可无；**focus → ROI 只维护这份 SQL。**
