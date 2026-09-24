# 图上局部定位算法 elaborate

> 对应四联图：`docs/summaries/three_step_subset_localize.png`  
> 实现：`scripts/tencent_gr/run_three_step_subset_localize.py`  
> 决策理论：cmean Drill 产品式 → \(K^\star\) → 一次 FSDS（非因果）

```
Standardize(W1) → L1 merchant MMD → L2 user MMD → L3 order shift → FSDS
```

---

## 1. 在图上到底在定位什么

TencentGR 边 = `(user, item, ts)`；终端成功 = click（`y_convert=1`）。  
商户 = `item_feat.122`（加密 shop/advertiser）；未映射 item 退化为 singleton shop。

定位目标不是 ATE，而是：**两时段 \(P^{W_1}\) vs \(P^{W_2}\) 的分布漂移，落在哪一层实体子集上**。  
图结构只提供实体键（merchant / user / order），同一套 procedure 可换颗粒度。

| 层 | 实体 | 分数 | 图面板 |
|---|---|---|---|
| L1 | merchant | RBF-MMD²(W1,W2) | 左上蓝 |
| L2 | user（限制在 L1 商户内） | RBF-MMD² | 右上橙 |
| L3 | order = 成功边 `u\|i\|ts`（限制在 L2 用户内） | \(\|x-\mu_{\mathrm{user}}\|_2\) | 左下红 |
| stop | 特征 tip | FSDS F-score | 右下绿 |

---

## 2. 共用对象（cmean / Drill 语言）

父支撑 \(S\)（边集），下一实体键 \(v\)：

\[
\delta(S)=\mu_S^{(W_2)}-\mu_S^{(W_1)},\qquad
D_{v:}=\mu_v^{(W_2)}-\mu_v^{(W_1)},\qquad
r_v=\|D_{v:}\|_2.
\]

对齐度 \(a_v=\langle D_{v:},\delta\rangle/(\|D\|\|\delta\|+\varepsilon)\)。  
谱尖锐度：\(\mathrm{Tail}(r)=\) Top-\(K\) 质量占比；\(\mathrm{CV}(r)=\mathrm{sd}/\mathrm{mean}\)。  
尖端残差：\(R(K)=\|\delta(S\setminus K)\|^2/\|\delta(S)\|^2\)。

实践里 L1/L2 用 **MMD²** 代替 \(\|D\|_2\)（捕捉全分布差，不只均值）；L3 边太稀，退化为相对用户 W1 中心的 \(\ell_2\)。

---

## 3. 四步算法（对照四联图）

### Step 0 — 标准化（只 fit W1）

窗内 FE → 候选列（去 rank / 泄漏列）→ `StandardScaler` **仅在 W1** 上 fit，两窗共用。  
避免尺度主导 MMD；无泄漏。

### Step 1 — L1 merchant（图：左上）

对每个候选商户 \(m\)（两窗边数 \(\ge n_{\min}\)，按体积截断 candidates）：

\[
\mathrm{score}(m)=\widehat{\mathrm{MMD}}^2\big(X_m^{W_1},\,X_m^{W_2}\big)
\quad\text{(RBF, median bandwidth, unbiased)}.
\]

取 Top-\(k_m\) → 商户核 \(K_m\)。边集收窄到这些商户。

**读图：** MMD² ≈ 1.35–1.59，谱相对平缓——“若干商户一起漂”，不是单点爆炸。  
若 Tail/CV 低且 \(\|\delta\|\) 大 → 理论可 **D-stop** 在商户层；实践默认继续下钻到用户（业务可读）。

### Step 2 — L2 user（图：右上，**只在 \(K_m\) 内**）

同样 MMD²，但候选用户来自 L1 边集。取 Top-\(k_u\) → \(K_u\)。

**读图：** MMD² 跳到 30–70，谱尖锐——漂移集中在少数用户。这是典型 **D-drill**：行谱尖锐 → 继续收核，**暂不跑 FSDS**。

稀疏 fallback：MMD 候选不足时用 \(\|\mu_2-\mu_1\|_2\)；再不行按体积取 Top 用户。

### Step 3 — L3 order（图：左下，**只在 \(K_u\) 内**）

订单 ≈ 终端成功边；`order_id = user_item_ts`。  
对每条边（优先 W2）：

\[
\mathrm{shift}=\|x-\mu_{\mathrm{user}}^{(W_1)}\|_2
\quad\text{（标准化空间；无该用户 W1 均值则用 \(\|x\|\)）}.
\]

取 Top-\(k_o\) 作可视化尖端；FSDS 训练支撑取 **merchant×user 边集**（订单太稀、标签不够）。

**读图：** Top 订单与 Top 用户高度重合（如 `322226|448811`、`602140|448811`）——L3 是把 L2 尖端落到具体成交边上，不是另开一套故事。

### Step 4 — 一次 FSDS（图：右下）

仅在最终支撑 \(K^\star\)（实践 = L2 用户边；理论 = Drill 停层）上：

\[
\mathrm{Scaler}\to\mathrm{Var}\to\mathrm{SelectKBest}(F)\to\mathrm{HGB}/\mathrm{LogReg}.
\]

SelectKBest **只看 W1-train**；W2 仅 temporal holdout。

**读图：** `u_span_sec`、`i_credit_last` / `i_share_last` 主导 F-score —— 用户活跃跨度 + 物品末次归因份额是 tip；其余 covisit / mismatch 为次级。

---

## 4. Drill 产品式（何时停、何时继续）

\[
\mathrm{Drill}
=
\mathbf{1}\{\|\delta\|\ge\varepsilon_\delta\}
\cdot
\mathbf{1}\{\mathrm{Tail}(r)\ge\tau\lor\mathrm{CV}(r)\ge c^\star\}
\cdot
\mathbf{1}\{R(K)\le R^\star\}
\cdot
\mathbf{1}\{\mathrm{mass}/\pi\ \mathrm{ok}\}.
\]

| \(\|\delta\|\) | \(r\) 谱 | 动作 |
|---|---|---|
| 小 | 平 | Stop：无 FSDS、无下钻 |
| 大 | 平 | **停钻**；本层跑一次 FSDS |
| 大 | 尖 | **下钻**到 \(K\)；FSDS **推迟到 \(K^\star\)** |
| 尖但 \(R\approx 1\) | — | 假尖端，退回本层 FSDS |

**列（特征）从不授权下钻。** 大 \(\|\delta\|\) + 尖列份额只授权「若停在此则 FSDS」。  
嵌套默认：`merchant → user → order`；每层算 Drill；`Drill=0` 或触业务上限 → 定 \(K^\star\) → **唯一一次** FSDS。

---

## 5. 必须报的方向（正向 / 负向）

\(K^\star\) 定后：

- 结局：\(\Delta\bar y=\bar y_{W_2}-\bar y_{W_1}\)（convert 升/降/平）
- tip：\(\mathrm{sign}(\delta_j)=\mathrm{sign}(\mu_{W_2,j}-\mu_{W_1,j})\)
- 实体：\(s_v=\langle D_{v:},\delta\rangle\)（与父漂移同向 / 对冲）

禁止只交「Top merchants/users」而不给正向/负向。

---

## 6. 与图运行数字的对照（gap=30d 原型）

| | |
|---|---|
| k | merchant=80, user=47, order=52 |
| 局部边 | W1=324, W2=6743 |
| mapped_rate | ≈0.10（多数 item 用 singleton shop proxy） |
| L1 head MMD² | ≈1.59 … 1.35 |
| L2 head MMD² | ≈72 … 30（尖锐） |
| FSDS tips | `u_span_sec`, `i_credit_last`, `i_share_last`, … |

解释：L1 铺开候选商户 → L2 把质量收到少数用户 → L3 点名异常成交边 → FSDS 只在收窄后的支撑上讲特征故事。

---

## 7. 一句话 mnemonic

**列尖 + 行平 = 特征故事（停钻，FSDS）。**  
**行尖 = 实体故事（下钻）。**  
**FSDS 挂在最终支撑，不挂在每一探测层。**

```bash
PYTHONPATH=. python3 scripts/tencent_gr/run_three_step_subset_localize.py \
  --root data/tencent_subset --max-users 20000 --gap-days 30
```
