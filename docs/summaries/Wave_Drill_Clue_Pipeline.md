# 监测波次 → Drill 团伙支撑 → 连续图谱分布变动线索

> 三条链合成一条反欺诈可用的 **连续图谱变动线索机**（算法层；行动项另表）。  
> 尺子锁定：MMD + PO + cmean；输出必带正负向（`direction` JSON）。

相关：[`Graph_Shift_Method_Elaboration.md`](Graph_Shift_Method_Elaboration.md) ·  
[`Tip_Sign_Algorithm_Elaboration.md`](Tip_Sign_Algorithm_Elaboration.md) ·  
[`AntiFraud_Tip_Sign_Actions.md`](AntiFraud_Tip_Sign_Actions.md) ·  
[`AntiFraud_Eval_Char_Model.md`](AntiFraud_Eval_Char_Model.md)

---

## 0. 一句话

```text
边流 (u,i,t)
  → 每个 Δt 一张图 G_t，算三路 T_t          【监测波次】
  → R_t=1 时 rank-avg(cmean,MMD,PO) 收核     【Drill 出团伙支撑 K*】
  → K* 上 FSDS + sign(δ), sign(Δȳ)          【连续图谱分布变动线索卡】
```

不是静态社区列表，也不是单点黑名单：是 **随时间冒出来的分布漂移波次 + 共变支撑块 + 可读 tip**。

---

## 1. 监测波次（WHEN）

### 1.1 波次是什么
时间轴按 \(\Delta t\)（日/时）切片，每片一张交互图 \(G_t\)（窗内 FE → \(X_t\)）。  
相对参考 \(G_{t-1}\) 或滚动 \(G_{t-h:t-1}\) 算标量：

\[
T_t^{\mathrm{MMD}}=\widehat{\mathrm{MMD}}^2(X_{t-1},X_t)-e_{\mathrm{MMD}},\quad
T_t^{\mathrm{PO}}=\mathrm{PO\text{-}gap}(t\!-\!1\to t)-e_{\mathrm{PO}},\quad
T_t^{\mathrm{cmean}}=\|\delta_t\|_2-e_{\mathrm{cm}}.
\]

OnlineRFPerm（EWMA-\(p\) + alpha-investing）→ \(R_t\in\{0,1\}\)。

### 1.2 「波次」的操作定义
- **一次波次** = 一段 \(R_t=1\) 的连续（或带短间隙的）报警区间 \([t_a,t_b]\)  
- 记录：起止、峰值 \(T\)、Lead vs 其它流、AR  
- \(R_t=0\)：**不** Drill、不产线索卡（省算力、防故事碎片化）

### 1.3 算法要点
- 三路并行、不互替；Fire 可用并集或工程规则（如 MMD∧(PO∨cmean)）  
- Burn-in 定 \(e_{\mathrm{ref}}\)；AR ≠ Type-I FAR  
- 两窗原型 = 单次人工波次（W1↔W2）；连续版 = 自动波次流

---

## 2. Drill 出团伙支撑（WHERE）

### 2.1 团伙支撑 = \(K^\star\)
在波次对应边集 \(S_t\)（通常 \(E_{t-h:t}\) 或 cur 片）上，对业务实体键  
（merchant → user → order；或 item）算：

\[
\mathrm{score}(v)=\mathrm{rank\text{-}average}(\mathrm{cmean},\mathrm{MMD},\mathrm{PO}).
\]

谱尖锐（Tail/CV）→ 提名 \(K\)；残差 \(R(K)=\|\delta(S\setminus K)\|^2/\|\delta(S)\|^2\) 小才接受。  
停钻得 **\(K^\star\)** = 这一波 **共同承担 graph shift 质量的实体块**。

### 2.2 为什么叫「团伙支撑」而不是「团伙认定」
- 算法目标：解释 \(\delta(S)\) 的检索支撑  
- 成员：在 score 尖端上共现/共变  
- **不是** ATE 子群、不是法律认定的欺诈组织  

读法：这一波异常 **落在这批实体上**；是否欺诈由人审/规则接力。

### 2.3 与监测的咬合
| | |
|---|---|
| 仅监测 | 知道「有波次」，不知道哪一块 |
| 仅 Drill | 每次全量下钻太贵、无时间结构 |
| **波次触发 Drill** | 只在 \(R_t=1\) 时定位 → \(K^\star(t)\) |

跨波次：对齐 \(K^\star\)（Jaccard）看加剧/消退；ID 漂了就当新支撑候选。

---

## 3. 连续图谱分布变动线索（WHAT 交付物）

### 3.1 线索卡 = 算法最小输出
一次波次 + 一次 Drill 后，冻结：

```json
{
  "wave": {"t_start": ..., "t_end": ..., "R": 1, "T_mmd": ..., "T_po": ..., "T_cmean": ...},
  "support": {"K_star": [...], "alpha": ..., "R_residual": ..., "key": "merchant|user|item"},
  "direction": {
    "Dy": ..., "sign_Dy": "pos|neg|flat",
    "tip_signs": {"u_span_sec": "+", "i_share_last": "+"},
    "tip_delta": {...},
    "report": "[Direction] Dy=... (pos); tip_signs={...}"
  },
  "tips": {"J_star": [...], "f_scores": {...}}
}
```

`direction` 由 `direction_report.build_direction_dict` 写入 `summary.json`。

### 3.2 「连续图谱分布变动」指什么
- **连续**：滑窗 \(G_t\)，波次沿时间排布  
- **图谱**：边特征来自交互图 FE（度、共现、share/credit…）  
- **分布变动**：比的是 \(P_{\mathrm{ref}}(X)\) vs \(P_{\mathrm{cur}}(X)\)，不是边集合对称差  
- **线索**：\(K^\star\) + tip + 正负向，给人审/规则用，不是自动定罪

### 3.3 Tip 与正负向（算法，极短）
1. **仅在 \(K^\star\)** 上 FSDS → \(J^\star\)（对 \(y\) 的判别 tip）  
2. \(\mathrm{sign}(\delta_j)\)：该 tip 在块上均值升/降  
3. \(\mathrm{sign}(\Delta\bar y)\)：块上 click 升/降  

列不改行门；无符号 = 输出不完整。

---

## 4. 三条链怎么合成一条故事

| 阶段 | 问题 | 输出 |
|---|---|---|
| 监测波次 | 何时图谱分布异常？ | \(R_t\)、波次区间、\(T_t\) 曲线 |
| Drill 团伙支撑 | 异常落在哪块共变实体？ | \(K^\star\)、\(\alpha\)、\(R(K)\) |
| 变动线索 | 往哪漂、哪些特征 tip？ | `direction` + \(J^\star\) |

时间线例子：
1. \(t=12\)：\(R=1\)，开波次  
2. Drill → merchants \(\{m_1,m_2\}\)、users 尖端  
3. 线索：`sign_Dy=pos`，`i_share_last=+` → 「末跳份额升且成功率升」的分布变动卡  
4. \(t=15\)：同支撑 score 仍尖 → 波次延续 / 加剧  
5. \(t=18\)：\(R=0\) 且 \(r\) 变平 → 波次结束  

---

## 5. 与两窗原型的关系

| 连续线索机 | 两窗原型 |
|---|---|
| 每个 \(\Delta t\) 监测 | 一次 W1↔W2 |
| 波次触发 Drill | 直接 localize |
| JSON 线索卡流 | 单次 `summary.json` |

公式相同：`rank-avg(cmean,MMD,PO)` + \(K^\star\) 上 FSDS + direction。  
实现入口：`run_w1w2_mmd_po_localize_fsds.py`（单波次）；连续 = 外层 `Δt` 循环 + OnlineRFPerm。

---

## 6. Takeaway

**监测波次**管何时；**Drill**管哪块团伙支撑；**线索卡**管分布往哪变、tip 正负向。  
三者顺序固定：\(R_t\to K^\star\to J^\star+\mathrm{direction}\)。这就是连续图谱分布变动线索的算法本体。
