# 多模态归因：停止逻辑与是否下钻

> 非因果、modality-agnostic 的多颗粒度分布漂移定位协议。  
> 本文只钉 **停止 / 下钻** 口径；上游标准化与实体提名见
> [`W1W2_MMD_PO_LOCALIZE_FSDS.md`](./W1W2_MMD_PO_LOCALIZE_FSDS.md)。

**闭环一句话：**

```
图谱特征 → X → 标准化
      → cmean：列 δ / 行 r_u（早停）
      → K* → 一次 FSDS
```

### \(K^\star\) 是什么（别神化）

**就是：在图谱特征 \(X\) 上，用 cmean 行谱收完之后剩下的那份边/实体支撑**——最终拿去跑一次 FSDS 的子集。

| \(K^\star\) 是 | \(K^\star\) 不是 |
|---|---|
| 图谱特征上的检索子集 | 图算法产物 |
| cmean（+\(\pi\)/mass）收支撑的结果 | community / ego / ULS / PPR / GraphScan |
| 交付用的边集或实体 id 集 | 因果核、效应子群 |

**其他都超纲：就这样。** 不做 community、ego、图上收核、GNN；输入侧就是图谱特征。

PO 边级分数作同口径旁证；全量 MMD leave-one-entity 重跑禁止。  
**不需要因果 justification，但必须写清 localization / 检索口径的 justification**（见 §1）——外形像 subgroup，最容易写歪。

---

## 1. 为什么像 subgroup、为什么不是、写的时候怎么 justify

这一节是整份协议里最 tricky、也最需要写进代码注释 / 报告前言的部分。

### 1.1 外形为什么像

程序长这样：

```text
先在父集上打分 → 取出高分实体核 → 只在核内再往下切一层再打分
```

和 post-hoc subgroup 的「先看主分析，再切子群看效应」**流程同构**。  
若报告再出现 Lift、Δ、CV、「更显著」之类词，读者会默认你在做子群效应——**长得像不是巧合，写的时候必须主动划界**。

### 1.2 目标函数不同（划界的硬核）

| | Subgroup analysis | 本协议（distribution-shift localization） |
|---|---|---|
| 科学问题 | 处理效应 \(\tau\) 是否随子群变 | 两窗特征分布 \(P^{W_1},P^{W_2}\) 的差异质量落在哪块 mass |
| 主对象 | \(\hat\tau(S)\)、交互项、CI | \(\delta\)、\(D\)、\(r_u\)、可选 \(\widehat{\mathrm{MMD}}^2\)、enrichment |
| 名单角色 | **推断对象**（常需多重比较 / 预登记） | **检索 / 交付候选**（要稳定 \(\pi\) + 业务顶） |
| 「更好」的含义 | 子群效应更大 / 更显著 | 子支撑上均值漂移更尖（\(r\) 谱）或残余 \(R\) 下降 |
| 下钻 | 再切子群继续估效应 → fishing | 换实体颗粒，问行谱是否还异质 → 可早停 |
| 需要的 justification | 因果识别 + 选择性推断 | **非因果声明** + 检索稳定性 + 交付口径 + 描述性闸门 |

所以：不是「免责所以不用 justify」，而是 **justify 换轨**。  
因果轨上的 pre-specification / FWER 不是本协议的主担保；本协议要担保的是——你交出去的核是可重复的漂移载体名单，不是 ATE 子群。

### 1.3 我们明确不声称什么（写进每份报告）

必须可逐条否定：

1. **不是 ATE / CATE**：\(W\) 是 period（W1/W2），不是 treatment。  
2. **不是「这些商户转化效应更大」**：\(r_u\) / Lift 只描述窗间特征（或分数）尖度。  
3. **Drill 门不是假设检验**：Tail/CV/\(R\) 阈值是工程闸门，不是 \(\alpha\)。  
4. **稳定 \(\pi_e\) 不是多重校正**：管的是 TopK 可重复，不管 Type I under no-effect null。  
5. **FSDS 特征重要性 ≠ 因果机制**。

少写一句「非因果、非 subgroup effect」，读者就会按 subgroup 读。

### 1.4 我们仍必须写清的 justification（具体条目）

写代码、写 PR、写结果报告时，按条勾：

| # | 要写明的 | 不写会怎样 |
|---|---|---|
| J1 | **任务**：W1/W2 特征分布漂移的多颗粒定位（retrieval） | 被读成效应估计 |
| J2 | **原子**：下一层 `entity` 的边集；\(D,\delta,r\) 同 scaler | 颗粒不清，和 subgroup 切片混谈 |
| J3 | **列 vs 行**：\(\delta\)=feature guidance；\(r_u\)=唯一下钻信号 | 「特征漂了所以下钻」——假 subgroup 话术 |
| J4 | **停止规则**：行平或业务顶或 \(R\)-reject；非「效应不显著故停」 | 停钻被读成「无异质效应」 |
| J5 | **名单**：交付认 \(\pi\) 核，不认单次 TopK；标 fragile | 事后钓鱼名单 |
| J6 | **Lift / enrichment**：若出现，定义 \(=\bar s_S/\bar s_C\)，禁止效应比 | 直接变 subgroup |
| J7 | **选择评估（可选加分）**：probe 定核、holdout 窗只评估 overlap/\(R\) | 完全 post-hoc 无诚实评估 |
| J8 | **业务顶**：统计建议与交付粒度分开写 | 「钻到 order」被当成必须对外口径 |

**一句话可放进模块 docstring：**

> This layer decides whether to refine a *localization support* by entity-level cmean row spectrum. It does **not** estimate subgroup treatment effects; \(W\) is period, gates are descriptive, and delivered kernels are stable retrieval sets under a business granularity cap.

### 1.5 最容易写歪的三句（禁止 → 替换）

| 禁止 | 替换 |
|---|---|
| 「下钻后效应更显著」 | 「下钻后 \(r\) 谱更尖 / \(R\) 下降，父集漂移由子核解释」 |
| 「这些用户是异质处理效应人群」 | 「这些用户是父核内 mean-shift 的主要载体（高 \(r_u\)）」 |
| 「CF/PO 证明应停止下钻」 | 「行谱齐次（+可选 \(g_u\) 同向）→ 描述性早停；CF 不进 Drill 乘积」 |

### 1.6 和 subgroup 的边界测试（自检）

写完一段结论，做三个替换测试：

1. 把文中所有「效应 / 显著 / uplift」删掉——若段落塌了，说明在靠 subgroup 语言撑着。  
2. 把「核」改成「检索名单 / 交付支撑」——若读不通，对象没立住。  
3. 问：若永远不下钻、只交本层 \(S\)+FSDS，协议是否仍自洽？——应自洽（行平即停是一等公民，不是失败）。

过不了这三条，就还在写假 subgroup，而不是 localization。

### 1.7 小结

外形像切子群，是因为 **都在数据依赖地缩小支撑**；  
不是 subgroup，是因为 **优化与声称的是漂移质量的局部化，不是 \(\tau\) 的切片推断**。  
因此：**不用 causality justification，但要用 §1.4 的 J1–J8 写满**——代码好写，话术和报告前言不能省。

---

## 2. 业务顶 vs 统计下钻

两套东西，不要合成一个开关。

| | 业务顶 | 统计下钻（本文） |
|---|---|---|
| 由谁定 | 接结果的人（商户运营 / 人群 / case） | 父核内 cmean **行谱** |
| 含义 | **对外交付粒度 ≤ 顶** | 要不要建议再细一层 |
| 计算 | 顶以下仍可算 probe | `Drill=0` 则该支路 localization 终点 |
| Justify 时 | 写清谁接结果、为何截在此粒 | 写清行谱门 + 非效应语言（§1.5） |

统计可以建议下钻；交付物仍截在业务顶。二者分栏写，避免「统计钻到哪交付到哪」被读成无约束的 subgroup fishing。

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
8. 报告 / 模块 docstring 勾完 **§1.4 J1–J8**；过一遍 **§1.6** 三句自检（效应词删除测试）

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

## 11. 逐步 elaborate：一整层怎么看

下面按「人在读报告」的顺序走一遍。假设业务顶在 user，当前停在 **某个商户核 \(S\)**，下一层 `entity = user_id`。

### Step A — 只算一次矩阵

```text
1. 边属于 S；X 用 W1-scaler
2. 按 (user, window) group mean → μ_u^(0), μ_u^(1)
3. D[u, :] = μ_u^(1) - μ_u^(0)          # (U_ok, d)
4. δ[:]    = μ_S^(1) - μ_S^(0)          # (d,)
5. r[u]    = ||D[u, :]||_2
6. share[j] = δ[j]^2 / ||δ||_2^2
```

到这里 **还没有** FSDS、没有千次 MMD。输出三块数：`δ`、`share`、`r`。

### Step B — 先读列（Feature guidance）

问三句，只许答特征侧：

1. \(\|\delta\|_2\) 算大吗？→ 否则 **G0**，本层故事很弱。  
2. 大的话，能量是否堆在少数维？→ `cumsum(sorted share)` 看 Top‑10% 维。  
3. 若集中 → 记下 \(J^\star\)（维名 / 图谱特征名），标记「本层 *若最终停在这里* 值得 FSDS」——**此刻先不跑**。

报告栏标题建议：`[Feature guidance] ||δ||=…; J*={…}; FSDS_if_stop=yes/no`

### Step C — 再读行（Drill）

问两句，只许答实体侧：

1. `r` 的 Tail@10% / CV 高不高？  
2. `a_u = cos(D_u, δ)` 是不是大多数靠近 1？

- 行平 + 多数对齐 → **D_stop**：商户内核里用户「一起漂」，再切 user 没信息。  
- 行尖 → **D_drill**：取出 `K = Top(r) ∩ {π≥π*}`（若还没有 π，先用 Top(r)+mass，下一轮 probe 补 π）。

报告栏：`[Drill] Tail(r)=…; CV=…; decision=STOP|DRILL; K=…`

### Step D — 若 DRILL，做 D_check（仍是 cmean）

```text
S' = S \ K 的边
δ_remain = mean(X|S',W2) - mean(X|S',W1)
R = ||δ_remain||² / ||δ||²
```

- \(R\) 降到位 → 承认 \(K\)，下一层父集换成 \(K\) 的边，**entity 换成 order**（或业务下一粒），回到 Step A。  
- \(R\approx1\) → 行谱假尖，**撤回**，按 D_stop 处理。

### Step E — 停层才 FSDS

当本层 `Drill=STOP`，或已触业务顶：

```text
K* = 当前交付支撑（停层的 S，或最后一次通过 D_check 的 K）
在 K* 的边上跑一次 FSDS
```

若 Step B 是 G0，可跳过 FSDS。  
若一路钻到业务顶仍行尖，交付截在顶，报告里可写「统计仍建议更细，但超业务顶，未交付」。

---

## 12. 四个格子的玩具数例

约定：\(d=2\)，\(U=4\) 个 user，数字是示意。

### 格 1：\(\|\delta\|\) 小、行平 → 全停

\[
D=\begin{bmatrix}0.05&0.02\\0.04&0.03\\0.06&0.01\\0.05&0.02\end{bmatrix},
\quad \delta\approx(0.05,0.02),\ \|\delta\|_2\text{ 小},\ r\approx\text{全相近}
\]

→ G0 + D_stop：不 FSDS、不下钻。  
话术：「本层均值几乎无漂移，用户间也无尖核。」

### 格 2：\(\|\delta\|\) 大、行平 → 停钻、本层特征

四个 user 都约等于 \(\delta=(2.0,\ 0.2)\)：

\[
r_u\approx\text{相同},\ a_u\approx1,\ \mathrm{share}_1\approx0.99
\]

→ G1 + D_stop。  
话术：「整层一起漂，且几乎全在 feature‑1；**不要下钻**；在当前 \(S\) 上跑 FSDS，盯 feature‑1。」  
这是「feature-level cmean 当 guidance」的 canonical 胜利场景。

### 格 3：\(\|\delta\|\) 大、行尖 → 下钻

\[
D=\begin{bmatrix}3.0&0.2\\ 2.8&0.1\\ 0.1&0.0\\ 0.1&0.0\end{bmatrix}
\quad\Rightarrow\quad
r\approx(3.0,\ 2.8,\ 0.1,\ 0.1),\ \mathrm{Tail@50\%}\text{ 很高}
\]

父 \(\delta\) 也会大（被前两行拉起来）。  
→ G1 信号有，但 **Drill 优先**：\(K=\{u_1,u_2\}\)，算 \(R\)；FSDS **先别跑**，等 \(K^\star\)。  
若只看 \(\delta\) / share，会误判成「本层做 FSDS 就收工」——漏掉实体局部化。

### 格 4：行尖但 \(R\not\downarrow\) → 撤回

Top(r) 取出的 \(K\) 去掉后 \(\delta_{S\setminus K}\approx\delta\) → \(R\approx1\)。  
常见原因：mass 极小的实体范数虚高、或 TopK 与真正载体不一致。  
→ 撤回下钻；收紧 mass / 改用份额加权 \(r_u\cdot n_u\) 再看一眼。

---

## 13. 和「只看 feature-level δ」的口径对照

| 你想回答的问题 | 该看 | 不该看 |
|---|---|---|
| 本层有没有均值漂移？ | \(\|\delta\|_2\) | 单次 MMD TopK |
| 盯哪些维 / 要不要本层 FSDS？ | \(\mathrm{share}(\delta)\) | \(r_u\) |
| 要不要再细一层实体？ | \(\mathrm{Tail}/\mathrm{CV}(r)\)、\(a_u\)、\(R\) | 单独的 \(\delta\) |
| 最终特征排序 | \(K^\star\) 上的 FSDS | 每一层探索都 FSDS |
| PO 旁证 | 同 `entity` 的 \(g_u,\eta^2\) | \(\mathrm{Var}(\hat\tau(x))\) |

**一句话对齐（重复强调）：**  
feature-level conditional-mean difference = 好 guidance；  
同一套 cmean 的 **行谱** = 下钻 / 停止的完整 signal；少算行谱会把格 2 和格 3 混掉。

---

## 14. 报告模板（强制分栏）

```text
## Layer L  (parent=S, drill_key=user_id, business_cap=user)

### Justification (non-causal; not subgroup effect)
- Task: period W1/W2 distribution-shift localization (retrieval)
- W = period, not treatment; r/Lift/R are descriptive gates
- Delivered kernel = stable π set under business cap (J1–J8)
- Joint order: full-dim | F→E (J* then r) | E→F  (declare one)
- Alt keys compared (if any): [category: eff=…; geo: eff=…]

### Feature guidance
- ||δ||_2 = …
- Top share features J* = […]
- Gate: G0 | G1 | G2
- FSDS_if_stop: yes/no

### Drill decision
- Tail(r)=…  CV(r)=…  mean(a_u)=…
- Gate: D_stop | D_drill
- K (if drill) = […]
- R (if drill) = …  → accept | reject

### Kernel / stop
- K* = …
- Reason: row-flat | business_cap | R-reject | …
- FSDS: run_once_on(K*) | skipped(G0)
```

三栏不能合成一句「因为特征漂了所以下钻」。  
Justification 栏不能省略——代码好写，这段是防读成 subgroup 的主担保。

---

## 15. 不同 entity key / value，与特征×实体联合

缩支撑的「刀」由 **划分键** 决定。同一父集 \(S\)、同一 \(X\)，换一把键就换一条 localization 路径。

### 15.1 Key–value 口径

| 词 | 含义 | 例 |
|---|---|---|
| **entity key** | 边表上的一列划分字段 | `merchant_id` / `user_id` / `category_id` / `item_id` / `geo` … |
| **entity value** | 该键下的一个水平 | 某个商户 id、某个用户 id |
| **划分** | \(S\to\{S_v:v\in\mathrm{values(key)}\}\) | 边按 value 分桶 |
| **行谱** | 对该 key 算的 \(D^{(key)},\ r^{(key)}\) | 只有指定 key 后才有 drill 问题 |

钉死：

- 没有 key，就没有「下一层」——`(n,)` 分数不能单独谈下钻。  
- 一次 Drill 决策 **只对一个 key**；多 key 是 **多条并列路径**，不是一个糊成的 \(r\)。  
- value 进核 \(K\) 的是 id 集合；交付写清 `(key=…, values=…)`。

### 15.2 两类 key 关系

**（A）嵌套键（树）** — 默认主路径  

```text
merchant → user → order
```

每层父集是上一层 \(K\) 的边；key 预先有包含关系。  
业务顶 = 树上允许的最深 key。  
早停 = 某一层 \(r\) 平，不再往子 key 走。

**（B）并列键（兄弟）** — 同一 \(S\) 上换刀  

```text
S 上分别用 key=category | key=geo | key=merchant
→ 三条 r 谱、三套 (α, γ)
```

用来回答：漂移更容易沿哪条业务语义收支撑。  
**不是**因果谁真；是哪把 key 的 \(\gamma/\alpha\)（捕获/收缩）更优、\(\pi\) 更稳、更好交。

### 15.3 多 key 怎么比（仍非因果）

对每个 key \(k\)，在同一 \(S\) 上：

\[
\alpha^{(k)}=\frac{|S_{K}|}{|S|},\quad
\gamma^{(k)}=1-R^{(k)},\quad
\mathrm{eff}^{(k)}=\frac{\gamma^{(k)}}{\alpha^{(k)}+\varepsilon}
\]

再加 \(\mathrm{Tail}(r^{(k)})\)、\(\bar\pi^{(k)}\)、业务可解释性。

| 结果 | 动作 |
|---|---|
| 只有一把 key 行尖且 \(\gamma\) 过关 | 沿该 key 下钻 |
| 多把都尖 | 取 \(\mathrm{eff}\) 高且 \(\pi\) 稳的主路径；其余写「并列路径」 |
| 嵌套路径与并列路径都出核 | 可交 **交集**（更稳、更小）或分栏交，禁止 silently 并成一个糊名单 |
| 全部行平 | 本层停；走 feature guidance / FSDS |

报告必须写 `drill_key=...`，禁止只写 Top ids 不写键名。

### 15.4 特征×实体联合（在同一张 \(D\) 上）

\(D\in\mathbb{R}^{U\times d}\) 本身就是联合对象：行=实体 value，列=feature。

三种用法，**顺序要写进口径**（换序会换核）：

**路径 F→E（先列后行，常用）**

```text
1) 父集上用 δ 得 J* = TopShare(δ)     # feature guidance
2) 限制列：D̃ = D[:, J*]
3) r_u = ||D̃_u||_2                   # 只在漂移维上比实体
4) 按 r 做 Drill / R-check
```

作用：噪声维稀释 \(r_u\) 时，先收特征再排实体，行谱更干净。  
注意：\(J^\star\) 来自父 \(\delta\)，不是 FSDS 重训；FSDS 仍只在最终 \(K^\star\)。

**路径 E→F（先行后列）**

```text
1) 全维 r_u 定 K
2) 只在 K 的边上重算 δ_K / share → J*_K
3) 解释与（停在 K 时的）FSDS
```

作用：先定位载体，再在核内看维——适合「实体很尖、父 \(\delta\) 被稀释」时。  
与 F→E 核不同要并列报告，不取silent并集冒充唯一真相。

**路径 Joint cell（读相互作用，一般不单决定 Drill）**

\[
\mathrm{cell}_{uj}=|D_{uj}|
\quad\text{或}\quad
\frac{|D_{uj}|}{\|D\|_*+\varepsilon}
\]

看热力：是否「少数 (entity, feature) 格子」撑起漂移。  

- 可出解释：`value=u* 在 feature=j* 上独漂`  
- **默认不单独当 Drill 门**（格子数 \(U\times d\)，更像事后读图）  
- 若要用：先 mass 闸门 + 只在 \(J^\star\times\mathrm{Top}(r)\) 子块上读，并写明 exploratory

### 15.5 联合逻辑的决策口诀

```text
δ / share     → 列：漂在哪维、本层 FSDS？
r (全维或 J*) → 行：沿当前 key 下不下钻？
多 key        → 哪把钥匙的 (γ, α, π) 更好？
E↔F 顺序      → 必须声明 F→E 或 E→F；两序不一致则分栏
cell 热力     → 解释用；进门需降范围并标注
```

和 subgroup 的关系：多 key / 联合只是 **多把检索刀**；每把刀仍读 \(\delta,r,\gamma\)，不读 \(\tau\)。Justify 时写清 key、顺序、是否联合，避免「又切又筛」被读成多重 subgroup fishing。

---

## 16. 还有其他的吗（同级短表）

| 块 | 要不要进默认协议 | 一句 |
|---|---|---|
| 嵌套 key 主路径 + 早停 | **要** | 主故事 |
| 并列 key 比 \(\mathrm{eff}\) | 建议有 | 同 \(S\) 换刀 |
| F→E 联合（\(J^\star\) 上再 \(r\)） | 建议有 | 降噪行谱 |
| E→F 联合 | 可选 | 父 \(\delta\) 被稀释时 |
| Cell 热力 | 解释可选 | 不默认当门 |
| 多 key 核交集 | 可选 | 更稳更小 |
| Holdout 上评 \(\gamma\)/overlap | 建议有 | 防定核窗虚 rate |
| 连续剥离曲线 | 可选增强 | TopK 的细版 |
| CF \(\hat\tau\) 进 Drill | **不要** | 轴不对 |

没有第三套并列大方法论；就是 **缩支撑 + 多钥匙 + 特征维约束** 三件事织在同一张 \(D\) 上。

**LaTeX prototype（公式 + justification）：**
[`Feature_Entity_Joint_CMean.tex`](./Feature_Entity_Joint_CMean.tex)。

**定稿整链 LaTeX：**
[`Graph_Feature_CMean_Localization.tex`](./Graph_Feature_CMean_Localization.tex)。

---

**收束：**  
列 \(\delta\) = 本层漂不漂、看哪维、要不要 FSDS；  
行 \(r_u\) = 沿**当前 entity key** 要不要下钻、钻到哪些 value；  
多 key = 多条路径比 \(\gamma/\alpha/\pi\)；特征×实体 = 声明 F→E 或 E→F 后再出核；  
停在行平或业务顶；FSDS 只打最终核。口径按决策表走，多模态归因的停止逻辑即闭环。
