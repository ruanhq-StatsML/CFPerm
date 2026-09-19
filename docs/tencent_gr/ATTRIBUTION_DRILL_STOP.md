# 多模态归因：停止逻辑与是否下钻

> 非因果、modality-agnostic 的多颗粒度分布漂移定位协议。  
> 本文只钉 **停止 / 下钻** 口径；上游标准化与实体提名见
> [`W1W2_MMD_PO_LOCALIZE_FSDS.md`](./W1W2_MMD_PO_LOCALIZE_FSDS.md)。

**闭环一句话：**

```
标准化 → 实体核（mass + 稳定 π）
      → 父集同一套 cmean：列看 δ（feature / 是否本层 FSDS），行看 r_u（是否下钻）
      → 早停得 K* → 核上一次 FSDS
```

PO 边级分数作同口径旁证；全量 MMD leave-one-entity 重跑禁止；不谈因果 justification。

---

## 1. 这不是 subgroup analysis

| | Subgroup analysis | 本协议 |
|---|---|---|
| 问题 | 效应在子群是否不同 | 哪块 mass 的特征分布在窗间漂了 |
| 量 | ATE / 交互 | cmean 差、\(r_u\)、enrichment Lift、可选 MMD |
| 名单 | 推断对象 | 检索 / 交付核 |
| 下钻 | 再切子群估效应 | 换实体颗粒，看行谱是否还尖 |

Justification 换轨：要的是 **名单稳定性 + 业务交付口径 + enrichment/残余可操作性**，不是子群效应的多重校正。  
话术一旦滑成「这些商户效应更大」，就掉回 subgroup——协议禁止。

---

## 2. 业务顶 vs 统计下钻

两套东西，不要合成一个开关。

| | 业务顶 | 统计下钻（本文） |
|---|---|---|
| 由谁定 | 接结果的人（商户运营 / 人群 / case） | 父核内 cmean **行谱** |
| 含义 | **对外交付粒度 ≤ 顶** | 要不要建议再细一层 |
| 计算 | 顶以下仍可算 probe | `Drill=0` 则该支路 localization 终点 |

统计可以建议下钻；交付物仍截在业务顶。

---

## 3. 公共设定（列 / 行必须共用）

当前父集边 \(S\)（当前核或业务顶下的边），下一层实体键 `entity`：

```text
X      (n, d)   仅用 W1 fit 的 StandardScaler 之后
w      (n,)     0 = W1, 1 = W2
entity (n,)     下一层颗粒（merchant / user / order …）
```

**颗粒度（subset 原子）：**

| 层 | entity | 边集 |
|---|---|---|
| L1 | `merchant_id` | 商户内边 |
| L2 | `user_id` \| 父商户核 | 用户内边 |
| L3 | `order_id` \| 父用户核 | 单边 |

Mass 闸门（进谱才算）：

\[
U_{ok}=\{u:n_1(u)\ge n_{\min},\ n_2(u)\ge n_{\min}\}
\]

实体级与父级 conditional mean：

\[
\mu_u^{(w)}=\mathbb{E}[X\mid e=u,W=w],\quad
D_{u:}=\mu_u^{(1)}-\mu_u^{(0)}\in\mathbb{R}^d
\]

\[
\delta=\mu_S^{(1)}-\mu_S^{(0)}\in\mathbb{R}^d
\]

**矩阵口径：** \(D\in\mathbb{R}^{|U_{ok}|\times d}\) —— 一行一个实体，一列一个 feature。  
后面「看列 / 看行」都指这个 \(D\) 和同一个 \(\delta\)。Scaler、\(S\)、feature 列不允许两路各算各的。

---

## 4. Feature guidance = 看列（+ 父向量 \(\delta\)）

### 4.1 量

| 量 | 定义 | 含义 |
|---|---|---|
| \(\delta\) | \((d,)\) | 本层整体均值漂没漂、哪维在漂 |
| \(\|\delta\|_2\) | 标量 | mean-shift 总强度 |
| \(\mathrm{share}_j=\delta_j^2/\|\delta\|_2^2\) | \((d,)\) | 漂移能量落在哪些 feature |

### 4.2 判定（只指导特征侧，**不下钻**）

| 代号 | 条件 | 动作 |
|---|---|---|
| **G0** | \(\|\delta\|_2 < \varepsilon_\delta\) | 本层均值几乎无漂移 → **不跑 FSDS**；也不因 feature 去下钻 |
| **G1** | \(\|\delta\|_2 \ge \varepsilon_\delta\) 且 Top-\(k\) share 集中 | 本层值得跑 **一次** FSDS；盯 \(J^\star=\mathrm{TopShare}(\delta)\) |
| **G2** | \(\|\delta\|_2 \ge \varepsilon_\delta\) 但 share 很平 | 有漂移但维上不尖；FSDS 可选 / 预期弱 |

**允许的结论句式：**「本层漂/不漂」「盯哪些 feature」「本层要/不要 FSDS」。  
**禁止：** 从 G0–G2 直接写「因此下钻到 user」。

> Feature-level cmean difference 是好的 **guidance**（漂没漂、看哪维、要不要在本层 FSDS），  
> **不是**完整的下钻 signal。

---

## 5. Drill / 停止 = 看行谱（同一 \(D\)）

### 5.1 量

对每个 \(u\in U_{ok}\)：

\[
r_u=\|D_{u:}\|_2
\]

对齐（是否「大家都像父 \(\delta\)」；\(\delta\approx0\) 时跳过）：

\[
a_u=\frac{\langle D_{u:},\delta\rangle}{\|D_{u:}\|_2\,\|\delta\|_2+\varepsilon}
\]

谱汇总：

\[
\mathrm{Tail}(r)=\frac{\sum_{u\in\mathrm{Top}K}r_u}{\sum_u r_u},\quad
\mathrm{CV}(r)=\frac{\mathrm{sd}(r)}{\mathrm{mean}(r)+\varepsilon}
\]

### 5.2 判定（**唯一**下钻口径）

| 代号 | 条件 | 动作 |
|---|---|---|
| **D_stop** | Tail / CV 低（且多数 \(a_u\approx1\)） | 实体齐次 → **不下钻**；若 G1 则本层 FSDS 收工 |
| **D_drill** | Tail 或 CV 高 | 少数实体扛 shift → **下钻**；\(K=\mathrm{Top}K(r)\cap U_{ok}\)（再交稳定 \(\pi\)） |
| **D_check** | 钻后残差比 \(R=\|\delta_{S\setminus K}\|_2^2/(\|\delta_S\|_2^2+\varepsilon)\) | \(R\) 明显降 → 核成立；\(R\approx1\) → 假尖，**撤回** |

**允许的结论句式：**「\(r_u\) 谱尖 → 下钻到 \(K\)」「谱平 → 停在 \(S\)」。  
**禁止：** 用「\(\delta\) 某几维很大」代替「\(r\) 谱尖」。

### 5.3 为什么必须看行、不能只看 \(\delta\)

两个父集可以有相同的 \(\|\delta\|_2\)：

- 全体实体 \(D_{u:}\approx\delta\) → **停**（整层一起漂）  
- \(\delta\) 由少数实体撑起 → **钻**

只看 \(\delta\in\mathbb{R}^d\) 分不开这两种。行谱 \(r_u\) 这一下省不得；计算上只是按 entity 再一次 group mean，与只算 \(\delta\) 同级。

---

## 6. 决策表（口径合一）

| \(\|\delta\|\)（列 / 父） | \(r\) 谱（行） | 动作 |
|---|---|---|
| 小 | 平 | **停**；不 FSDS、不下钻 |
| 小 | 尖 | 少见；先查 mass / scaler；默认谨慎 |
| 大 | 平 | **停钻**；本层 **FSDS / 盯 \(J^\star\)** |
| 大 | 尖 | **下钻**到 \(K=\mathrm{Top}(r)\)；FSDS **等到 \(K^\star\) 固定再跑一次** |

口诀：

- **列尖、行平** = 特征故事，实体不用细挖  
- **行尖** = 实体故事，该下钻  
- **FSDS** 挂在最终交付支撑上，不挂在每一层探索上  

合成规则（描述性闸门，不是假设检验）：

\[
\mathrm{Drill}
=\mathbf{1}\{\mathrm{mass}/\pi\text{ ok}\}
\cdot\mathbf{1}\{\mathrm{Tail}(r)\ge\tau\ \lor\ \mathrm{CV}(r)\ge c^\star\}
\]

\(\mathrm{Drill}=0\) ⇒ 该支路 localization 终点 \(S\)；  
\(\mathrm{Drill}=1\) ⇒ 进入 \(K\)，再在 \(K\) 上重复 §3–§6（直到停或触业务顶）。

---

## 7. PO-risk `(n,)` 怎么进同一口径

边级分数 `s.shape == (n,)`（全父集 **fit 一次**）必须配同一 `entity`：

\[
g_u=\mathrm{mean}(s\mid e=u)
\quad\text{（或中心化后的实体份额）}
\]

\[
\eta^2=\frac{\mathrm{Var}_u(g_u)}{\mathrm{Var}_i(s)}
\]

- \(\eta^2\) 小且 Tail\((g)\) 低 → 与 **D_stop** 同向：分数变异不在下一层实体上  
- 与 \(r_u\) **同向尖** → 加信下钻  
- **反向** → 只报 divergence，默认 **不**用 \(g\) 推翻「行平则停」（除非显式改口径）

禁止：

- 用 \(\mathrm{Var}_x(\hat\tau(x))\) 当子集异质（轴不对）  
- 按商户 / 实体 **重 fit** PO 一千次  
- 无 `entity` 时宣称「从 `(n,)` 看出不用下钻」——没有划分就没有 drill 问题

---

## 8. 算力红线（和停止逻辑配套）

| 禁止 | 正法 |
|---|---|
| 对 \(\|U\|\) 个实体各跑一次全量 MMD / PO | 实体贡献用 **cmean / RFF 嵌入加减**；真 MMD 最多 Top‑几十抽检 |
| Observation-level LOO on MMD | 原子是实体；用行谱 \(r_u\) 或嵌入 LOEO |
| 每层探索都跑 FSDS | FSDS **一次**，在最终 \(K^\star\)（或停层的 \(S\)）上 |
| 固定强制 merchant→user→order | 每层过决策表；齐次则早停 |

稳定核（与行谱交）：

\[
\pi_e=\frac{1}{B}\sum_{b=1}^{B}\mathbf{1}\{e\in\mathrm{Top}K^{(b)}\},\quad
\mathcal{K}=\{e:\pi_e\ge\pi^\star\}
\]

交付认 \(\pi\) 名单，不认单次 TopK；单次 TopK 与 \(\mathcal{K}\) 的差集 = 不稳定噪声。

---

## 9. 最小落地清单

1. Scaler 只 fit W1；\(D,\delta,r\) 全用同一 \(X\)  
2. 先定 `entity` = 下一层颗粒，再算 \(D\)  
3. 报 \(\delta\) / share → **G0–G2**（feature + 是否本层 FSDS）  
4. 报 \(r\) / Tail / CV / \(a_u\) → **D_stop / D_drill**  
5. 钻则出 \(K\) + \(R\)；可选对一下 \(g_u\)  
6. FSDS 只在最终支撑上跑一次  
7. 报告分栏：**Feature guidance** | **Drill decision** | **Kernel \(K^\star\)** —— 禁止混写成一句

---

## 10. 和现有脚本的关系

| 组件 | 现状 | 本文角色 |
|---|---|---|
| `run_standardize_mmd_fsds.py` | 扁平 item-MMD → FSDS | 上游提名 / 终局 FSDS |
| `run_three_step_subset_localize.py` | 固定三层往下走 | 应包上 §6 早停，避免空转 |
| `run_cf_parsimonious_smoke.py` | period-W 上 CF \(\hat\tau^2\) | **并行透镜**；不进 \(\mathrm{Drill}\) 乘积 |
| 本文 | 口径 writeup | 停止 / 下钻的决策层 |

模态换特征矩阵 \(X\)、颗粒换 `entity` 键，§3–§6 不变。

---

## 参考阈值（工程默认，可校准）

非正式 α；按数据尺度调，改了要写进报告。

| 符号 | 建议起点 | 用途 |
|---|---|---|
| \(n_{\min}\) | 5–20（视层） | mass 闸门 |
| \(\varepsilon_\delta\) | 相对 \(\|\delta\|\) 历史分位或固定小阈值 | G0 |
| \(\rho_f\) / Top-\(k\) share | Top‑10% 维累计 ≥ 0.5 | G1 集中 |
| Tail / CV 门 | Tail@10% 或 CV 分位 | D_stop / D_drill |
| \(\pi^\star\) | 0.5–0.7 | 稳定核 |
| \(R^\star\) | 0.4–0.5 | D_check 核成立 |

---

**收束：**  
列 \(\delta\) = 本层漂不漂、看哪维、要不要 FSDS；  
行 \(r_u\) = 要不要下钻、钻到谁；  
停在行平或业务顶；FSDS 只打最终核。口径按决策表走，多模态归因的停止逻辑即闭环。
