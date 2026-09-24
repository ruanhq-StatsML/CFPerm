# 是否冻这一层：统计上能最优刻画什么

看板是产品启发式（PDF 也写了：not a methodology contribution）。下面只谈 **决策对象、能最优的部分、不能最优的部分**。T 仍是 batch 标签，不是处理。冻层不是「第 j 层导致了 MSE 下降」——那一层级的唯一分解做不到。

---

## 表 1 · 决策对象（先钉死，否则没有最优）

k 个 layer group → k+1 个嵌套策略，不是 k 个独立开关。

| 策略 | 可训练 | 冻住 |
|---|---|---|
| \(i=0\) | 无（或只 head，视实现） | 全部 |
| \(i\) | 只训 top \(i\) 组 | 底 \(k-i\) 组 |
| \(i=k\) | 全训 | 不冻 |

嵌套：trainable\((i)\subset\) trainable\((i+1)\)。要选的是 **切一刀的深度 \(i^\star\)**，不是对每一层做 k 次独立检验。

现在代码在 **两都崩** 之后才开 clones，\(i^\star=\arg\min_i \widehat{\mathrm{PO}}_i\)，另记 \(i^\star_{\mathrm{MSE}}=\arg\min_i \widehat{\mathrm{MSE}}_i\)（`PO_Dict` / `MSE_Dict`）。

---

## 表 2 · 两阶段（这是对的；最优在第二段，不在 2× 门槛）

| 段 | 问题 | 统计对象 | 不是什么 |
|---|---|---|---|
| 门 | 这一窗要不要进入「选冻深」 | 全局 PO（\(P(Y\mid X)\) hop 的距离）、serving MSE（经济损失）、MMD vs \(D_{\mathrm{ref}}\)（\(P(X)\)） | 不是 CATE(T) |
| 选 | 若进门，切在哪一刀 | 比较 k+1 个 **策略的风险** | 不是层的 Shapley / 唯一归因 |

门用 AND（PO 崩且 MSE 崩）是 **保守的交检验**：假冻（不该冻却冻）更贵。这不是 Neyman–Pearson 最优，是产品损失不对称。PO 崩 MSE 没崩 → watch，因为 RF PO 不该先塌。MSE 崩 PO 安静 + MMD 过线 → X-shift，**不要按 concept 去冻底训顶**（PDF 表里那一格写了 freeze/back-prop，统计上应读成系统级 \(P(X)\)，不是 `train_top`）。

2× 因果 MA 只是门槛，不是最优临界值。渐进里它可以永不破；那时 WHEN 是 \(T(\tau)\)/Brier 相对参考水平。

---

## 表 3 · 「最优」只能是策略选择，不能是层归因

不可能（和 talk 同一句）：观测到 \(\Delta\mathrm{MSE}\) 无法唯一拆成「第 j 层该冻」。无穷多种 drift 复制同一误差。

能谈最优的是 **有限策略类 \(\{0,\ldots,k\}\) 上的风险**：

\[
i^\star \in \arg\min_{i=0,\ldots,k}
\Big\{
\underbrace{R_i}_{\text{serving risk of clone }i}
+ \lambda\, C_i
\Big\}
\]

\(R_i=\mathbb{E}[\ell(Y,\hat f_i(X))]\)（MSE/Brier），\(C_i=\) 可训练参数量。\(\lambda\) 是算力/业务，不是从数据里估出来的因果权重。

三种可辩护的最优准则（由紧到松）：

| 准则 | 公式 | 何时用 | 和现在实现 |
|---|---|---|---|
| 最省算力、服务不差 | 最小的 \(i\) 使得 \(R_i \le R_k+\varepsilon\) | 冻得越多越好，只要房租不涨 | 还没做；应对 \(\mathrm{MSE}_i\) 做嵌套检验 |
| 服务最优 | \(\arg\min_i R_i\) | 只在乎下一批误差 | 代码里的 \(i^\star_{\mathrm{MSE}}\) |
| 机制距离最小 | \(\arg\min_i \mathrm{PO}_i\) | 把 clone 当 \(P(Y\mid X)\) 距离的探针 | **现在的 \(i^\star\)** |

PO 最小 ≠ MSE 最小。这就是层尺度上的 **服务误差 vs 份额误差**。最优产品规则不是选其中一个，而是：

**在 serving 可接受的 \(i\) 里，再按 PO 选（或反过来：PO 可接受里按 MSE 选）。**

不要用 PO 单独决定冻深（指针可以先脏）；不要用 MSE 单独当 concept 冻深（X-shift 时全训/冻底都会误伤）。

---

## 表 4 · 嵌套检验（比逐层独立检验更对口）

因为策略是嵌套的，对每一层做 k 个独立 2× 检验会 overlap、不可比。

从全训往下冻（推荐的最优形状）：

1. \(H_k\): 全训。始终可接受。
2. \(H_i\): \(R_i \le R_k+\varepsilon\)（再冻一层，服务不差）。
3. 找 **最小** 的 \(i\) 使得 \(H_i\) 不被拒 = 在服务约束下冻得最多。

闭式检验 / sequential from \(k\) down to \(0\)，FWER 沿这条链控。不必 online-bootstrap；batch 够大时用 batch 内对照或 \(D_{\mathrm{ref}}\) 上的配对差 \(\widehat{R}_i-\widehat{R}_k\)。

从底往上解冻：第一次 \(\Delta_i = R_i-R_{i+1}\) 显著为负的地方，是「再训一层开始划得来」的切点。\(\Delta_i\) **不必单调**（X-shift、遗忘），所以不能假设 isotonic 就停。实践：看整条 `MSE_Dict[0..k]`，取约束最优，而不是第一个 p<0.05。

---

## 表 5 · 层上的 PO 和 MSE 各回答什么

| 序列 | 问的问题 | 崩了意味着 | 单独拿来冻 |
|---|---|---|---|
| 全局 serving MSE | 现模型还付不付房租 | 要动策略 | 不够：不知切哪一刀，也不知是 X 还是 Y\|X |
| 全局 PO-risk | \(P(Y\mid X)\) 相对 \(D_{\mathrm{ref}}\) 的距离 | 机制 hop 的证据 | 不够：RF PO 不该先塌；且 T 不是处理 |
| 全局 MMD vs ref | \(P(X)\) | 协变量 | 崩了就不要按 concept 冻底 |
| \(\mathrm{MSE}_i\) | **策略 i** 的服务风险 | 这一刀的房租 | 选 \(i^\star_{\mathrm{MSE}}\) |
| \(\mathrm{PO}_i\) | 策略 i 训练后，机制距离还剩多少 | 这一刀有没有把 \(P(Y\mid X)\) 对上 | 选 \(i^\star_{\mathrm{PO}}\)；与 MSE 不一致时不要独断 |

Clone 的 \(\mathrm{PO}_i\) 用的是该 clone 的 \(\mu=\hat f_i\)，问的是「用这个冻深去拟合新 batch 之后，batch 指标 T 下的伪结局风险」。它仍是距离，不是「冻层 i 的因果效应」。

---

## 表 6 · 和渐进 / latency 的接法

冻深选择每窗一次，评估仍是 **服务误差面积**，不是 hop。

| 做法 | 统计后果 |
|---|---|
| 等 2× 两都崩再开 clones | 假冻少；渐进里可能永远不开 → 面积来自「该选 i 却一直全训」 |
| 每窗都开 k+1 clones 用 MSE 约束选 i | 更接近表 3 的最优；算力 = k+1 倍（PDF 的 10 层 + 全局就是这个） |
| 只在 serving 走起来之后比较 \(\mathrm{MSE}_i\) | 折中：门用服务误差，选用层风险。和「动手看服务、动塔看份额」同一句 |

算力够：滑动窗特征库（recon、T、Brier）当 clever covariate；冻深仍用 clones 的 \(\mathrm{MSE}_i,\mathrm{PO}_i\) 选，不要用 last-batch 层表征 MMD 当冻层开关。

---

## 一句

**最优刻画的是：在嵌套冻深 \(\{0,\ldots,k\}\) 上，于服务风险约束下选 \(i^\star\)。**  
2× 门、AND、\(\arg\min\mathrm{PO}\) 都是这个决策的启发式实现。不是层因果、不是唯一分解。MVP 保持启发式即可；若要把「冻这一层」写成统计程序，就做成表 4 的嵌套约束选择，并用 \(\mathrm{MSE}_i\) 与 \(\mathrm{PO}_i\) 不一致作为层上的服务/份额四格。
