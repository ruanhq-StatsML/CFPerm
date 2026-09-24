# Tip / 正负向：算法层 elaborate（不谈行动项）

> 行动项（L0–L3）已够用；本文只谈 **算法对象、估计顺序、与 Drill 的隔离**。  
> 代码锚点：`run_w1w2_mmd_po_localize_fsds.py`（cmean / MMD / PO / FSDS）。

---

## 1. 算法要解的三件事（彼此隔离）

| 层 | 对象 | 估计 | **不**做什么 |
|---|---|---|---|
| **Shift 监测** | 全局 \(T_t^{\mathrm{MMD}},T_t^{\mathrm{PO}},\|\delta\|\) | OnlineRFPerm → \(R_t\) | 不选 tip、不定 \(K^\star\) |
| **实体定位（行）** | \(r_v=\) rank-avg(cmean,MMD,PO) | Tail/CV/\(R(K)\) → \(K^\star\) | 不用 F-score 下钻 |
| **特征 tip（列）** | \(J^\star\) + \(\mathrm{sign}(\delta_j)\) | **仅在 \(K^\star\)** 上 FSDS + cmean 符号 | 不用 tip 改 \(K^\star\) |

**硬顺序：** \(R_t\)（可选）→ \(K^\star\) → \(J^\star\) / 符号。  
列尖 **从不**授权改行支撑；行尖 **推迟** FSDS 到停层。

---

## 2. 支撑上的 shift 向量（符号的母体）

窗 ref / cur（或 \(t\!-\!1,/\,t\)），W1-only scaler 后：

\[
\delta(S)=\mu_S^{\mathrm{cur}}-\mu_S^{\mathrm{ref}}\in\mathbb{R}^d.
\]

- **能量** \(\|\delta(S)\|_2^2=\sum_j\delta_j^2\)：有没有一阶 shift  
- **列份额** \(\mathrm{share}_j=\delta_j^2/\|\delta\|^2\)：哪些坐标贡献能量（**无符号**）  
- **符号** \(\mathrm{sign}(\delta_j)\)：该坐标升高还是降低  

实体层同构：\(D_{v:}=\mu_v^{\mathrm{cur}}-\mu_v^{\mathrm{ref}}\)，\(r_v=\|D_{v:}\|_2\)（实践可换 MMD² / PO）。

---

## 3. 实体分（WHERE）——算法定义

候选实体 \(v\)（item / merchant / user…），在标准化空间：

\begin{align}
\mathrm{cmean}(v)&=\|\mu_v^{\mathrm{cur}}-\mu_v^{\mathrm{ref}}\|_2,\\
\mathrm{MMD}(v)&=\widehat{\mathrm{MMD}}^2(X_v^{\mathrm{ref}},X_v^{\mathrm{cur}}),\\
\mathrm{PO}(v)&=\mathrm{mean}(\hat\tau(x)^2\mid v),
\end{align}

\[
\mathrm{score}(v)=\mathrm{rank\text{-}average}\big(\mathrm{cmean},\mathrm{MMD},\mathrm{PO}\big).
\]

PO 拟合（period \(W\in\{0,1\}\)）：交叉拟合 \(\hat m,\hat e\)，伪结果  
\((Y-\hat m)(W-\hat e)\)，再回归得 \(\hat\tau\)；VIMP 可作 **列 prior**，默认 **不进** `score(v)` 产品式以外的门。

停钻：

\[
\mathrm{Drill}=1\{\|\delta\|\ge\varepsilon\}
\cdot 1\{\mathrm{Tail}(r)\lor\mathrm{CV}(r)\}
\cdot 1\{R(K)\le R^\star\}
\cdot 1\{\mathrm{mass}/\pi\}.
\]

\[
R(K)=\|\delta(S\setminus K)\|^2/\|\delta(S)\|^2.
\]

\(K^\star=\) 最后接受的实体核。**到此为止还没有 tip。**

---

## 4. Tip 集 \(J^\star\)（WHAT）——算法定义

### 4.1 监督 tip（官方 FSDS）

数据：仅 \(K^\star\) 内边；标签 \(y=\) click（稀时可用 PO 代理 / PCA rank，原型里 PO 用连续代理）。

\[
\texttt{StandardScaler}\to\texttt{VarianceThreshold}
\to\texttt{SelectKBest}(F_{\mathrm{classif}},k)
\to\mathrm{HGB}/\mathrm{LogReg}.
\]

- SelectKBest / 模型 **只 fit W1-train（或 ref 内 holdout）**；cur / W2 仅评估  
- 入选特征集合 = \(J^\star_{\mathrm{FSDS}}\)  
- 排序键 = \(F\)-score（或模型后的重要度）；**无符号**

语义：在 **已定位的 shift 支撑上**，哪些图特征对 \(y\) 仍有判别力。  
≠「哪些特征造成了 shift」（那是 \(\delta\) / MMD-LOCO）。

### 4.2 Shift tip（描述性，与 FSDS 并列）

在同一 \(K^\star\)、同一 scaler 空间：

\begin{align}
j &\in J^\star_{\delta}:\ \text{Top by }\delta_j^2,\\
j &\in J^\star_{\mathrm{MMD\text{-}LOCO}}:\ \text{Top by }\Delta\mathrm{MMD}_{-j},\\
j &\in J^\star_{\mathrm{PO\text{-}VIMP}}:\ \text{Top by VIMP}.
\end{align}

可 **rank-average 三列** 得混合 tip 池，再交给官方 FSDS 做最终 \(k\)（`po_help_select` 路径）。  
默认交付：**FSDS \(J^\star\)** + 对每个 \(j\) 附上 \(\mathrm{sign}(\delta_j)\)。

**工程面板（交付）：** `feature_methods_panel.py` 把 cmean / MMD-LOCO / PO-VIMP / FSDS-F
四列 Top 与共识 / FSDS-only / shift-only 挂到 `summary.feature_methods` 与审出卡
（只读，不改 Drill）。FSDS-only = \(y\) 判别强但未必是 shift 坐标；shift-only = 结构漂了但对 \(y\) 未必尖。

### 4.3 为什么 FSDS 必须在 \(K^\star\) 之后

全图 FSDS：tip 被大众量子淹没，指向「全局可预测性」，不是「这一波 graph shift」。  
先收支撑再 FSDS：tip 条件在异常块上，与 Drill 目标一致。

---

## 5. 正负向 —— 纯算法，两套符号

### 5.1 结局符号
\[
\Delta\bar y(K^\star)=\mathbb{E}[y\mid K^\star,\mathrm{cur}]-\mathbb{E}[y\mid K^\star,\mathrm{ref}].
\]
\(\mathrm{sign}(\Delta\bar y)\in\{+, -, 0\}\)。  
与 tip **独立估计**；不进入 SelectKBest。

### 5.2 Tip 符号（挂在 \(J^\star\) 上）
对每个 \(j\in J^\star\)：

\[
\mathrm{sign}(\delta_j(K^\star))=\mathrm{sign}\big(
\mu_{K^\star,j}^{\mathrm{cur}}-\mu_{K^\star,j}^{\mathrm{ref}}
\big).
\]

注意：
- \(J^\star\) 由 **\(y\)-判别 / VIMP** 选出  
- 符号由 **\(X\)-边际 cmean** 给出  
- 二者可不一致（特征对 \(y\) 重要，但两窗均值几乎没动 → flat tip sign）——卡片仍报，研判时标 `sign=0`

### 5.3 算法卡片（最小充分统计）
```text
K*, α, R(K), score_method=rankavg(cmean,MMD,PO)
J* = [(j, F_j, sign(δ_j), δ_j² share)]
Δȳ, sign(Δȳ)
optional: PO-VIMP_j, MMD-LOCO_j
```

没有 \(\mathrm{sign}(\delta_j)\) / \(\Delta\bar y\) 的 tip 列表 = **算法输出不完整**（不是产品问题）。

---

## 6. 算法层「问题」清单（实现时要对的）

| 问题 | 处理 |
|---|---|
| Click 极稀，FSDS 不稳定 | PO 连续代理 \(Y\)；或 `po_help` 扩池再 FSDS；报 \(\sigma\) over seeds |
| Tip 无符号 / 符号与 F 不一致 | 强制附 cmean 符号；flat 显式标 0 |
| 先 FSDS 后 localize | **禁止**；打乱顺序则 tip 无条件含义 |
| 用 \(\|\delta_j\|\) 当下钻门 | **禁止**；列门 ≠ 行门 |
| MMD 与 cmean 打架 | rank-average 保留张力；可单列诊断，不单踢一门 |
| \(K^\star\) 太小 | mass 门；放宽 \(k\) 或停在上层再 FSDS |
| W2 泄漏进 SelectKBest | 管道只 fit ref；cur 仅 transform / 评估 |
| 把 \(\mathrm{sign}(\delta)\) 当 \(\mathrm{sign}(\tau)\) | 文档禁止；只描述 period 差 |

---

## 7. 与监测标量的接口（仍算法）

\[
T_t^{\mathrm{cmean}}=\|\delta(S_t)\|_2-e,\quad
T_t^{\mathrm{MMD}}=\mathrm{MMD}^2-e,\quad
T_t^{\mathrm{PO}}=\mathrm{PO\text{-}gap}-e.
\]

Reject 后 **冻结** 当步的 \(S_t\)，跑 §3–§5。  
Tip 不回写 OnlineRFPerm 的 \(T_t\)（避免循环：用 tip 再定义监测）。

可选下游（仍算法、非行动项）：tip-cmean 幅度  
\(\Delta_{\mathrm{tip}}(t)=\sum_{j\in J^\star}|\delta_j(t)|\)  
与 LOCO-confirm 做 **尖峰 timing**；不改变 \(K^\star\) 定义。

---

## 8. Takeaway（算法一句）

**行定位用 rank-avg(cmean,MMD,PO) 定 \(K^\star\)；列 tip 用 \(K^\star\) 上 FSDS 定 \(J^\star\)；正负向用同一支撑上的 \(\delta\) 与 \(\Delta\bar y\) 挂符号。三者估计目标不同，顺序固定，门不串。**
