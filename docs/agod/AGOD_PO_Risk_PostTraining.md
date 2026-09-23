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

---

## Sample-level weights：还有哪些？高残差 ≠「难」

可以这样理解：**未见得代表是难。**  
\(\mathrm{PO}_i = |Y_i-\hat\mu(X_i)|\) 量的是 **在当前控制拟合 \(\mu\) 下的残差大小**（reject 后才拿来加权），不是样本的内在难度、也不是标注质量标签。

| \(\mathrm{PO}_i\) 高可能是… | 若当「难」去猛抬权会怎样 |
|---|---|
| 概念/残差漂移（我们想纠正的） | 合理：下一 fit 更快对准这块 |
| 标注噪声 / 错标 | 危险：把噪声学进去 |
| 极端 outlier \(Y\) | 危险：梯度被单点绑架 |
| \(\mu\) 本身没拟合好（probe 差） | 误报：假「难」 |
| 稀有但其实好学 | 不一定需要高权 |

所以更准确的说法是 **high-residual rows（高残差行）**，口语说「难样本」只是 reject 门控下的操作简称。

### 其它 sample-level weight（不止 √PO）

| mode | \(w_i\) | 含义 | 何时想用 |
|---|---|---|---|
| `uniform` | \(1\) | 不侧重行 | calm 窗；或关掉 B |
| `sqrt`（默认） | \(\propto\sqrt{\mathrm{PO}}\) | **soft** 抬高残差 | reject 后主菜；少过拟合当前批 |
| `prop` | \(\propto\mathrm{PO}\) | 硬抬 | 易过拟合 reject 批，一般不如 sqrt |
| `cbrt` | \(\propto\mathrm{PO}^{1/3}\) | 更软 | 残差噪声大、想更保守 |
| `inv` | \(\propto 1/\mathrm{PO}\) | **压**高残差 | 假设高残差=噪声/污染时 |
| `dre` | \(\propto p_{\mathrm{cur}}/p_{\mathrm{ref}}\) | 纯协变量比 | 只信 \(X\)-shift、不信残差时 |
| top-\(k\) hard | 只保留残差 top 20% 其余 \(w{=}0\) 或很小 | 截断式 subset | 要比 soft 更狠的行子集 |
| gated vs always | reject 才非 uniform / 每窗都抬 | 门控 | **默认 gated**；always 易伤平静包 |

### 和 ROI 怎么读（仍挂同一维）

```text
PO_i 高  ≠  标签难
PO_i 高  =  在 μ 下残差大 → 候选「多给一点 sample weight」
是否给权     =  先看 reject（概念/分布门）再看 mode（sqrt/cbrt/…）
ROI_B        =  next_mse_drop / fit_time   （约束：别把噪声抬爆 Acc）
```

**操作句：** sample-level weights = 在 reject 批内按残差 **重分行预算**；默认 soft √PO；高残差只是「值得多看一眼」的排序键，**不自动等于难、更不自动等于该猛学。**
