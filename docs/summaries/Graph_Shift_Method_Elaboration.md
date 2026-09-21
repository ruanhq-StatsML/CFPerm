# 基于 Graph Shift 的方法 elaborate（不谈社区发现算法）

> 核心对象只有一个：**图上的分布漂移**（graph shift）。  
> 门仍是 **MMD + PO + cmean**；图算法（Leiden 等）本文不讨论。  
> 数据：`(user, item, timestamp)` 边 → 每个时间片一张图 \(G_t\)。

相关实现：`run_w1w2_mmd_po_localize_fsds.py`、`time_window_feats.feature_engineer`。

---

## 1. Graph shift 是什么

边流 \(\mathcal{E}=\{(u,i,a,t)\}\)。对时间片 / 窗 \(W\)：

\[
G_W=(V_W,E_W),\qquad
E_W=\{e\in\mathcal{E}:t_e\in W\}.
\]

每条边（或 `(u,i)` 聚合）带图派生特征 \(x_e\in\mathbb{R}^d\)  
（度、活跃跨度、共现、share/credit、pop mismatch… —— 全由窗内边算，**无网外泄漏**）。

**Graph shift** = 两个窗上边特征律的差异：

\[
P^{W_{\mathrm{ref}}}(X)
\;\;\text{vs}\;\;
P^{W_{\mathrm{cur}}}(X)
\qquad\text{（以及可选的 }Y=\text{click}\text{）}.
\]

不是「边集合对称差有多大」，而是：**同样一张「用户–物品 交互图」上，行为特征的分布漂了没有、漂在哪块支撑上**。

连续时间：\(W_{\mathrm{ref}}=[t-h,t)\)，\(W_{\mathrm{cur}}=[t,t+\Delta t)\)，或相邻片 \(G_{t-1},G_t\)。  
两窗原型：\(W_1\) early / \(W_2\) late（gap≥30d）——同一套定义。

---

## 2. 三把尺：怎么量 graph shift

全部在 **W1-fit 的 StandardScaler** 空间里算（尺度不抢戏）。

### 2.1 cmean（一阶 / tip 方向）

\[
\delta=\mu_{\mathrm{cur}}-\mu_{\mathrm{ref}},\qquad
\|\delta\|_2,\quad
\mathrm{sign}(\delta_j).
\]

回答：整体均值往哪移、哪维 tip 升/降（正向/负向）。

### 2.2 MMD（全分布）

\[
\widehat{\mathrm{MMD}}^2(X_{\mathrm{ref}},X_{\mathrm{cur}})
\quad\text{(RBF, median bandwidth)}.
\]

回答：不止均值，高阶/形状是否变。全局 \(T_t^{\mathrm{MMD}}\) 或实体条件 MMD。

### 2.3 PO-risk（风险窗 / 可预测差）

把窗指标记成 \(W\in\{0,1\}\)（ref vs cur），在 \((X,Y)\) 上拟合 period-PO：

\[
\tau(x)\ \text{相关量},\quad
\text{实体分}=\mathrm{mean}(\hat\tau^2\mid v),\quad
\text{列 VIMP}.
\]

回答：漂移是否落在「对结局/代理结局仍可预测」的方向（侧镜，进 rank-average）。

### 2.4 合成一把排序尺（Drill 用）

实体 \(v\)（merchant / user / item / order）：

\[
\mathrm{score}(v)=\mathrm{rank\text{-}average}\big(
  \|\mu_v^{\mathrm{cur}}-\mu_v^{\mathrm{ref}}\|_2,\ 
  \mathrm{MMD}^2_v,\ 
  \mathrm{mean}(\hat\tau^2\mid v)
\big).
\]

**这就是「基于 graph shift」的唯一实体排序**——不引入社区目标函数。

---

## 3. 方法骨架（when → where → what）

```text
边流 → 窗内 FE → X_t on G_t
         │
         ├─ 全局 T_t = {MMD, PO-gap, ‖δ‖} → OnlineRFPerm → R_t     # WHEN
         │
         └─ R_t=1（或离线两窗）:
              按实体键算 score(v) = rank-avg(cmean,MMD,PO)         # WHERE
              Tail/CV/R 决定下钻或停 → K*
              K* 上一次 FSDS + sign(δ), sign(Δȳ)                   # WHAT
```

### WHEN（连续）
\[
T_t^{\mathrm{MMD}}=\mathrm{MMD}^2(X_{t-1},X_t)-e_{\mathrm{ref}},\ \ldots
\]
拒识 \(R_t\)：graph shift **发生了**。

### WHERE（局部化）
父支撑 \(S=\) 当前边集（或 reject 窗并集）。  
对候选实体算 `score(v)`，谱尖锐则收成 \(K\subset S\)；残差

\[
R(K)=\frac{\|\delta(S\setminus K)\|^2}{\|\delta(S)\|^2}
\]

小才接受。嵌套业务键：`merchant → user → order`（表内已有，无需社区算法）。

### WHAT（特征 tip）
仅在 \(K^\star\)：`Scaler→Var→SelectKBest→HGB`，W1-train only。  
Graph shift 的「故事」落在 tip 特征 + 正向/负向。

---

## 4. 为什么这叫「基于 graph shift」而不是点异常

| | 点/边孤立异常 | **我们的 graph shift** |
|---|---|---|
| 对象 | 单点分数 | **分布** \(P_{\mathrm{ref}}\) vs \(P_{\mathrm{cur}}\) |
| 支撑 | 单边/单用户 | 实体块 \(K^\star\)（多边共撑 \(\delta\)） |
| 证据 | 阈值/规则 | MMD（形状）+ cmean（方向）+ PO（风险） |
| 时间 | 静态阈值 | 窗对比 / 连续 \(G_t\) |
| 输出 | 黑名单 ID | \(K^\star\) + tip + 符号（描述性） |

团伙/客群在本框架里 = **共同承担 graph shift 质量的支撑块** \(K^\star\)，  
由 shift 分数尖端自然形成，**不是**另跑一个社区目标。

---

## 5. 与两窗原型的一一对应

| 概念 | 代码 |
|---|---|
| \(G_{W}\) | `feature_engineer(root, t0, t1)` → `grid` / `edge_df` |
| Scaler | `fit_standardizer(W1)` |
| 实体 graph shift | `score_item_mmd` / 三列 cmean·MMD·PO → `rank_score` |
| \(K^\star\) | top-\(k\) items（或三层 merchant→user→order） |
| tip | `run_fsds` on localized edges |

```bash
PYTHONPATH=. python3 scripts/tencent_gr/run_w1w2_mmd_po_localize_fsds.py \
  --root data/tencent_subset --max-users 20000 --gap-days 30
```

连续版：把 `(W1,W2)` 换成滑动 `(t-h,t)` vs `(t,t+Δt)`，全局三路 \(T_t\) 接 OnlineRFPerm；  
**实体打分公式不变**。

---

## 6. 报告条（每次 shift 交付）

```text
[Shift]   ||δ||=...; MMD²=...; PO-gap=...; R_t=0|1
[Support] K*= {entities}; α=|K|/|S|; R(K)=...; eff=...
[Direction] Δȳ=... (pos|neg|flat); tip_signs={j:+/-}
```

禁止：只有 Top ID 列表、没有 shift 量与方向。

---

## 7. Takeaway

1. **Graph shift** = 图窗上边特征分布之差。  
2. 度量 = **cmean + MMD + PO**（标准化后 rank-average）。  
3. 流程 = 监测 shift → 在实体上分解 shift → 停在 \(K^\star\) → FSDS。  
4. 连续时间只改窗的滑法，**不改 shift 定义与三尺**。  
5. 不依赖社区发现；业务键足够做 WHERE。
