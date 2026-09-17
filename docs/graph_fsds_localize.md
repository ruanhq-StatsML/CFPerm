# 图定位 → FSDS 收成一粒 → 两层 subset 归因

口径和冻层看板同一句：**定位，不是唯一分解**。T 是 batch 标签。份额是 ReLU 归一化，不是 Shapley，也不是 CATE。Y 只当 outcome，**不当特征、不当子集键**。

选了几个 feature 之后，在同一粒（默认 **订单**）上做 subset-localization，就能用 **MMD + PO-risk + Conditional Mean** 三支读数给出具体刻画：哪一个子群在往这三支上送、相对 other（补集）的 gap 有多大。

```
新 batch
  → 0. 冻结 D_ref：σ、分位数、实体画像、子集词表
  → 1. 各原生维度各自打三支（还不合表）
        订单行 / 商户节点(订单挂在点上) / 用户节点 / 商户表一阶矩
  → 2. FSDS：每个维度内给列打分，选出 loud / top-k（里面不能有 Y）
  → 3. 收到同一粒：订单（或商户）
        localize = 当前窗的 X；serve = 订单 X 用当前、实体列用 D_ref 画像
  → 4. 两层
        层 1  subset vs other：默认 **own-ref subset scan**（水平集 `{φ≥τ}` 或覆盖前缀）；
              lift 用 `merchant_id` / `user_id`（外键）。默认**不建图**。
        层 2  选中块的 LOGO（order / merchant / user）
  → 5. 对 loud subset 出三支肖像 + fingerprint
```

---

## 表 1 · 为什么必须先分维度、再 FSDS、再收粒

| 步 | 问什么 | 不先做会怎样 |
|---|---|---|
| 各维度三支 | 订单的 \(P(X)\)、商户节点的 \(P(X_{\mathrm{order}}\mid m)\)、商户表自己的 \(\mathbb{E}[X_m]\) 是不是在走 | 一上来 join 会把「商户目录变了」和「这张订单的 amount 变了」糊成一列，份额不可读 |
| FSDS | 在每个维度里，**哪些列**在驱动那三支 | 全表进 subset，噪声列把 MMD 稀释；也避免把 Y 或 Y 的分箱当特征 |
| 收成一粒 | 评估单位钉死（一张订单，或一个商户） | 订单指标和商户指标不能比份额；「谁贡献大」没有同一个分母 |
| 两层 | 哪一个 **subset**；这个 subset 里是哪一块 **feature block** | 只报全局 MMD 会漏掉「只有南区商户在走」 |

`merchant_id` / `user_id` 是外键，不是一张要跑社区发现的网络。变动 subset 是 Fast Subset Scan 的前缀（水平集 / 覆盖）。默认路径**不调用 networkx**。Louvain / kNN / SAGE 只在显式对照时才建。

---

## 表 0 · 多实体粒上的 subset scan（正经方法，不是图方法）

这是流行病学 / 异常定位里的 **spatial scan / Fast Subset Scan**（Kulldorff；Neill LTSS）：在一组实体上找最响的子集。势能可加时，最优解是按 \(\psi\) 排序的前缀，复杂度是一次排序，不是 \(2^N\)，也不是一张图的划分。

### 计算量

| 步 | 复杂度 | 贵不贵 |
|---|---|---|
| 每个实体打钟（MMD / \(\Delta\mu\) / 可选 PO） | \(O(N_{\mathrm{ent}}\cdot n_{\mathrm{bag}}^2)\) 量级（RBF-MMD 对袋） | **这是账单。** 16 个商户、每袋十几单可以忽略；上万商户、每袋上百单才需要抽袋 / 线性 MMD |
| 扫描本身（排序 + 取前缀） | \(O(N_{\mathrm{ent}}\log N_{\mathrm{ent}})\) | 可以忽略。\(N\) 是商户数或用户数，不是边数 |
| lift 到订单 | \(O(n_{\mathrm{order}})\) | `merchant_id ∈ loud_ids`，一次布尔 mask |
| 朴素穷举所有子集再对并袋打 MMD | \(2^{N}\) 次 MMD | **不做。** LTSS 用可加的 \(\sum\phi_i\)，所以前缀就是最优 |
| 图约束扫描（连通子图） / Louvain | 指数或迭代模块度 | **不做默认。** 那才是网络数据；这里外键不是边 |

所以 subset scan **作为切法几乎不花钱**。不要把它和「在图上搜连通异常块」混为一谈——后者才超纲、才贵。

Community detection 问「谁和谁连得密」。这里问的是：**哪一个实体粒上的哪一组 id，自己的 \(P(X)\) / \(P(Y\mid X)\) 在走。** 共单边的模块度切的是谁共享用户；南区用户跨商户随机挂，Louvain 跟种下的子集可以对齐也可以对不齐——对不齐不是调 resolution 能修的。图 / 网络数据不是这套 localization 的输入。

\[
\psi_i=\phi_i\quad\text{或}\quad\psi_i=\phi_i\cdot n_i,\qquad
\phi_i=\lVert\tilde v_i\rVert_1,\quad
v_i=(\mathrm{MMD}_i,\ \lVert\Delta\mu_X\rVert_i,\ \Delta\mu_Y_i,\ \mathrm{PO}_i)\ \text{vs own-ref}.
\]

两条前缀规则（同一条排序）：

\[
\begin{aligned}
\text{水平集}\quad
\hat S&=\{i:\psi_i\ge\tau\},\quad
\tau=\max(\mathrm{floor},\ \alpha\cdot\max_i\psi_i),\ \alpha=0.3;\\
\text{覆盖}\quad
\hat S&=\text{最短前缀，使 }\sum_{i\in\hat S}\psi_i\ge\beta\sum_i\psi_i,\ \beta=0.8.
\end{aligned}
\]

\(\tau\) 自适应：安静窗 \(\max\psi\) 小，floor 把噪声挡掉；onset 后 \(\max\psi\) 抬起来，前缀跟着走。覆盖前缀回答「偏移质量的 80% 落在哪些点」，水平集回答「谁过了响度门槛」。默认 subset 键用水平集；覆盖 / 质量加权（\(\phi\cdot n\)）写在同一粒上对照。

### 具体弄法（按实体粒扫，再 lift，再在订单粒上评）

粒 = 一种实体。这里两粒：商户、用户。`merchant_id` / `user_id` 只是订单挂在哪个实体上。**不要建成网络，也不要把两粒的 \(v_i\) 加进同一个 simplex。**

```
冻结 D_ref（σ、分位数、每个实体自己的 ref 袋）。新商户没有 own-ref → 跳过，不进扫描。

对每一粒 L ∈ {merchant, user} 独立：
  1. 每个实体 i 用自己的 D_ref 袋打捆
       v_i = (MMD_i, ‖Δμ_X‖_i, Δμ_Y_i, PO_i)     # own-ref
       φ_i = ‖ṽ_i‖_1
       ψ_i = φ_i            # 强度；要质量时改成 φ_i · n_i
  2. 排序 + 取前缀（Fast Subset Scan）
       水平集  Ŝ_L = { i : ψ_i ≥ τ }
       覆盖    Ŝ_L^cov = 累计 ψ 到 80% 的最短前缀
  3. lift 到订单（外键，不是图算法）
       订单 o ∈ Ŝ̂_L  ⇔  o 的该粒 id ∈ Ŝ_L

评估（唯一可加的单位 = 订单）：
       J(Ŝ̂_m, S*), J(Ŝ̂_u, S*), J(Ŝ̂_m, Ŝ̂_u)
       粒对了才回收：这个 DGP 变动在商户，用户粒 Jaccard 应接近 0

两层归因（默认用商户粒水平集 lift 出的标签）：
       订单标签 = loud / other（不是 C0/C1/C2）
       再打 π_MMD / π_PO / π_CMean 和 fingerprint
方向（捆 φ 混了就读切片，仍是扫描）：
       {i : MMD_i ≥ τ_MMD}、{i : |Δμ_Y|_i ≥ τ_Y}、…
```

| 步 | 做什么 | 不要做成 |
|---|---|---|
| 节点钟 | own-ref：这个商户/用户自己变了没有 | 用 full-ref 切（那是异质性） |
| 切 | 按 \(\psi\) 排序取前缀（水平集 / 覆盖 / \(\phi\cdot n\)） | Louvain / 模块度 / GraphRAG 社区 / 把一层切成 C0–C3 |
| 多层 | 每层自己扫，再 lift 到订单比 Jaccard | 把商户 MMD 和用户 MMD 丢进同一个份额；合成一张超图再切 |
| 图 | incidence：\(o\mapsto m(o),\ o\mapsto u(o)\) | 用共单连通当约束（用户跨南北随机挂，连通会把南北糊在一起） |
| 归因 | loud vs other，三支 + vs other gap | 把扫描结果再送进 Louvain 当 subset 键 |
| 方向 | metric slice：covariate 点亮 MMD / CMean_X；concept 点亮 CMean_Y / PO | 用社区标签代替「谁在变」 |

### 为什么这叫 localization，不叫 community detection

| | Subset scan（默认） | 共单 / bundled Louvain（对照，不是方法） |
|---|---|---|
| 问什么 | 哪些节点自己的钟响了，响的质量在哪一段前缀 | 谁和谁连成一块（模块度） |
| 可行集 | 排序的 \(N\) 个前缀（LTSS） | 节点的一个划分 \(c:V\to\{C_0,C_1,\ldots\}\) |
| 目标 | \(\sum_{i\in S}\psi_i\)（可加） | \(Q=\frac{1}{2m}\sum_{ij}(A_{ij}-k_i k_j/2m)\delta(c_i,c_j)\) |
| 种下 south、X 走 | 南区 \(\phi\) 高 → 进前缀 | 碰巧能收成一块，resolution 会碎 |
| 种下 south、X 不动（concept） | \(\Delta\mu_Y\) 把南区 \(\phi\) 抬起来 | 钟对不齐时 loud 社区 south_frac 掉到噪声 |
| 多层 | 每层一个扫描 | 共单图本来就是跨层的，更不能当切 |

求解器用 networkx 只是因为 incidence 图现成、无 torch。**默认切不用 `louvain_communities`。** 不要用 GraphRAG 合成图来代替这个目标。

### 同一套扫描上还能怎么取前缀（仍不是社区发现）

| 件 | 什么时候用 | 不是什么 |
|---|---|---|
| **水平集** `{φ≥τ}` | 默认 subset 键。门槛自适应 | 社区标签 |
| **覆盖前缀** 80% | 要看「偏移质量落在哪」；水平集太宽或太窄时对照 | 固定 top-k |
| **质量加权** \(\phi\cdot n\) | 小袋 \(\phi\) 虚高、大袋才是构成 | 只按 \(\phi\) 排序当贡献 |
| **metric slice** | 捆 \(\phi\) 把 MMD 和 \(\Delta\mu_Y\) 混在一起时，分开扫 | 第二套 subset 键；covariate 应点亮 MMD / CMean_X，concept 应点亮 CMean_Y / PO |
| 同方向拆开 | 只有 \(\cos(\tilde v_i,\tilde v_j)<0\) 的两拨 loud（X 走 vs \(Y\mid X\) hop）才拆 | 默认仍是一块 loud vs other |
| 同层邻接约束 | **仅当**你事先相信变动在该层的 geo / 类目图上连通。从 \(\arg\max\psi\) 往外长，目标仍是 \(\sum\psi\)，不是 \(Q\) | 共单图上做连通约束（会混层） |
| 结构 Louvain / SAGE-mean | 对照：模块度切不出 concept；SAGE 是读出 | 当变动 subset |

Y 进 \(\Delta\mu_Y\) / PO 是监控 outcome，**不进**节点特征、不进 \(Z\)。

冷启动（新商户没有 own-ref）：不把 full-ref 塞进扫描，否则又变成异质性。订单仍可挂 `other`。

---

## 表 0b · own-ref、同一维度、多层图怎么评

### 节点分数必须对「商户自己的 \(D_{\mathrm{ref}}\)」

| 钟 | 公式（节点 \(i\) = 一个商户的订单袋） | 问什么 | 会误伤什么 |
|---|---|---|---|
| **own-ref** | \(\mathrm{MMD}(X_{\mathrm{new}}(i), X_{\mathrm{ref}}(i))\)，\(\Delta\mu\) 同样用 \(i\) 自己的 ref | **这个商户自己变了没有** | 需要 \(n_{\mathrm{ref}}(i)\ge\min n\)；新商户没有自己的钟 → **跳过**，不进水平集 |
| **full-ref** | \(\mathrm{MMD}(X_{\mathrm{new}}(i), X_{\mathrm{ref}})\) | 这个商户和**全局混合**不像 | 从来就小众的商户永远 loud（异质性，不是 drift） |

变动 subset 的切只用 own-ref 水平集。full-ref 留着报「对整窗的贡献」（构成变化：南区单变多）。两口钟不要合成一个分数再切。

冷启动（新商户没有 own-ref）：不把 full-ref 塞进切，否则又变成异质性。订单仍可挂 `other`，两层归因时单独看。

σ 钉在**该节点自己的 ref 袋**上（own）或全局 \(D_{\mathrm{ref}}\)（full），都冻结，不在新窗重估。

### 同一维度怎么 localize

只在一个粒上切、只在这个粒上比份额：订单就订单，商户就商户。\(\pi_{\mathrm{MMD}}\) 的分母是这个粒上的 loud / other。不要把商户节点 MMD 和订单行 MMD 丢进同一个 simplex。

同一维度的评估：水平集 lift 成订单集合 \(\hat S\)，和种下的南区订单 \(S^\star\) 算 Jaccard；再在 \(\hat S\) 上打三支 + vs other。

### 多个维度 / 多层图怎么评估

图是分层的：商户层、用户层、（可选）订单 kNN 层。**每一层自己做 subset scan**，扫完 **全部 lift 到订单粒** 再比。这是唯一可加的评估单位。

```
商户层 scan {φ≥τ} / 覆盖前缀  →  订单集合 Ŝ_m
用户层 scan {φ≥τ} / 覆盖前缀  →  订单集合 Ŝ_u
订单粒 oracle（region）→  订单集合 S*
评估：J(Ŝ_m, S*)、J(Ŝ_u, S*)、J(Ŝ_m, Ŝ_u)
三支肖像只在 Ŝ_m / Ŝ_u 上打（已经是同一粒；标签是 loud / other）
```

| 看到 | 读法 |
|---|---|
| \(J(\hat S_m,S^\star)\) 高、\(J(\hat S_u,S^\star)\) 低 | 变动落在商户层；用户只是被共单带着走 |
| 两头 Jaccard 都高、层间 Jaccard 也高 | 两层指到同一批订单（仍不是唯一分解） |
| 层间 Jaccard 低、两头对 \(S^\star\) 都不高 | 切错层或 own-ref 样本不够 |
| 把两层的 \(v_i\) 直接加起来再取阈值 | **不要**。量纲和 n 都不同 |
| 把两层并成一张超图再 Louvain | **不要**。那又回到社区发现 |

多层不是 GraphRAG 合成一张超图。就是：每层一次 subset scan、lift 到订单、Jaccard。用户层在这个 DGP 里**不该**回收南区——用户是跨商户随机挂的。这正是评估：层对了才回收。

Louvain 对照仍按层切、再 lift，用来说明「切错目标会切碎」。不要把社区标签当默认 subset 键。不要在共单图上做连通约束——那会把南北用户糊在一起。

### 还有什么（同一套钟上的附件，不是新的面）

| 件 | 做什么 | 不要做成 |
|---|---|---|
| min_n | own-ref / 扫描 PO 的地板 | 小袋上打 RF PO |
| 强度 vs 质量 | \(\phi_i\) 高但 \(n_i\) 小 → 报 n 和 share | 只按 \(\phi\) 排序当贡献 |
| vs other | loud − other，两头对同一口 ref | 当成 Shapley |
| 结构 Louvain | 对照：模块度切不出 concept | 当变动 subset |
| 冻住 SAGE-mean | 结构读出 | 当切 |
| 新商户 | other | 用全局均值假造 own-ref 再切 |

---

## 表 2 · 三支读数（同一口 \(D_{\mathrm{ref}}\) 钟）

| 量 | 符号 | 问什么 | 口径 | 不是什么 |
|---|---|---|---|---|
| MMD | \(\mathrm{MMD}^2(Z_{\mathrm{new}}(S), Z_{\mathrm{ref}})\) 或 vs **own-ref** \(Z_{\mathrm{ref}}(S)\) | 这个子集的 \(P(X)\) 走了没有 | RBF，σ 钉在该空间的 \(D_{\mathrm{ref}}\) | 不是对历史 pairwise |
| PO-risk | \(\varphi=(Y-\mu)(T-e)\) 的风险 | 这个子集的 \(P(Y\mid X)\) 有没有 hop | 独立 RF nuisance；T=batch 标签 | 不是 CATE(T)，不是服务 MLP |
| Conditional Mean | \(\Delta\mu_Y(S)=\mathbb{E}[Y_{\mathrm{new}}\mid S]-\mathbb{E}[Y_{\mathrm{ref}}\mid S]\)，\(\Delta\mu_X(S)=\|\mathbb{E}[Z_{\mathrm{new}}\mid S]-\mathbb{E}[Z_{\mathrm{ref}}\mid S]\|\) | **一阶矩**怎么刻画这个子集 | 子集上用 \(\mathbb{E}[Y\mid S]\)；FSDS 列打分用冻结分箱的 \(\Delta(\mathbb{E}[Y\mid\mathrm{high}]-\mathbb{E}[Y\mid\mathrm{low}])\)，避免全局 \(\mathbb{E}[Y]\) 走的时候每一列都响 | \(\Delta\mu_Y(S)\) **不是** concept 检测器 |

三支一起读才有刻画：

| 看到 | 刻画（仍是定位，不是分解） |
|---|---|
| MMD↑、\(\|\Delta\mu_X\|\)↑、PO 安静；\(\Delta\mu_Y\) 可↑可 0 | 这个子集的 \(P(X)\) 在走。Y 均值若动，是因为 \(f\) 不是常数，**不是**单独的 concept |
| MMD 安静、PO↑、\(\Delta\mu_Y\)↑、\(\|\Delta\mu_X\|\) 小 | 这个子集的 \(P(Y\mid X)\) hop；Y 均值动是因为 \(f\) 变了 |
| 三支都响 | both：只报「这个子集贡献大」，不要拆成唯一的 concept/covariate 份额 |

**不要**用 \(\Delta\mu_Y\) 单独当 concept 开关。Covariate 把 amount 推走时，\(Y\mid X\) 没变，\(\mathbb{E}[Y\mid S]\) 也会走。

---

## 表 3 · 收到订单粒之后，评估单位是什么

订单粒的一行 = 一张订单。FSDS 选出的列拼成 \(Z\)：

\[
Z = \big[\;X^{\mathrm{order}}_{\mathrm{sel}}\;\big|\;X^{\mathrm{merchant}}_{\mathrm{sel}}\;\big|\;X^{\mathrm{user}}_{\mathrm{sel}}\;\big]
\]

两种 unify（泄漏规则不同）：

| mode | 订单列 | 商户/用户列 | 用来干什么 |
|---|---|---|---|
| `localize` | **当前窗** X | **当前窗** 实体 X（不含 Y） | 定位「现在谁在走」 |
| `serve` | 当前窗 X | **冻结的 \(D_{\mathrm{ref}}\) 画像** | 可进服务 / clever covariate；新商户用 ref 全局均值 |

评估只在 **同一个 \(Z\)** 上比。默认子集键是 **merchant-layer subset scan** lift 成的 `loud` / `other`。`region` 是种下的对照，不是切。`community` / `structural` 是 Louvain 对照，用来说明「切错目标会切不出 / 切碎 concept」。

两个对照，回答两个不同的问题：

| 对照 | 公式 | 评估什么 |
|---|---|---|
| vs own-ref | \(\mathrm{MMD}(Z_{\mathrm{new}}(S), Z_{\mathrm{ref}}(S))\)，\(\Delta\mu_Y(S)\) 用 \(S\) 自己的 ref | **这个子群自己变了没有**（南区商户的客是不是换了） |
| vs full \(D_{\mathrm{ref}}\) | \(\mathrm{MMD}(Z_{\mathrm{new}}(S), Z_{\mathrm{ref}})\) | **对整窗偏移的贡献**（含构成变化：南区订单变多） |
| vs other | \(\mathrm{metric}(S)-\mathrm{metric}(S^c)\)，两头都对同一口 \(D_{\mathrm{ref}}\) | **subset 和其他**：贡献是不是集中在 \(S\) |

层 1 份额（对 own-ref 的强度做 ReLU 归一化）：

\[
\pi^{\mathrm{MMD}}_S=\frac{\mathrm{ReLU}(\mathrm{MMD}_S)}{\sum_{S'}\mathrm{ReLU}(\mathrm{MMD}_{S'})}
\quad
\pi^{\mathrm{PO}}_S,\ 
\pi^{\mathrm{CMean}}_S
\text{ 同构。}
\]

CMean 份额用 \(\tfrac12|\Delta\mu_Y|+\tfrac12\|\Delta\mu_X\|\)（先各自非负）。贡献

\[
\mathrm{contrib}_S=\pi^{\mathrm{MMD}}_S+\pi^{\mathrm{PO}}_S+\pi^{\mathrm{CMean}}_S.
\]

\(\arg\max_S \mathrm{contrib}_S\) 就是 loud subset。某一支如果只是噪声（两边都接近 0，或相对差太小）**不进份额**——否则 concept 时 MMD 的 0.010 vs 0.011 会被归一化成 0/1，把南区的 \(\Delta\mu_Y\) 盖掉。再看 \(\mathrm{mix}^{\mathrm{MMD}}_S=\pi^{\mathrm{MMD}}_S/\mathrm{contrib}_S\) 等，给 fingerprint。

层 2：在选中的 \(Z\) 上对 `{order, merchant, user}` 做 LOGO（和模态那套同一套 \(\pi_{\mathrm{loss}},\pi_{\mathrm{PO}},\pi_{\mathrm{MMD}}\)）。全 batch 一次，loud subset 再一次。读法：层 1 指出 **南区订单**；层 2 指出是 **订单块的 amount** 还是 **商户块的 gmv** 在送。

---

## 表 4 · 选了 feature 之后，肖像长什么样

对 loud subset \(S^\star\) 必须同时报这些，才叫「具体刻画」：

| 字段 | 含义 |
|---|---|
| `n`, `share_of_batch` | 质量。25% 的单吃了 70% 的 \(\pi\) → 强度集中 |
| `mmd`, `po`, `cmean_x`, `cmean_y` | 三支 + 一阶矩（X 和 Y） |
| \(\pi\) 和 mix | 相对其他子集，哪一支在送 |
| `gap_vs_other` | subset − other；负的说明贡献其实在补集 |
| `fingerprint` | `x_shift` / `concept` / `both` / `quiet`（启发式，不是识别） |
| `read` | 一句话：X 走了 / \(Y\mid X\) hop / 两头都响 |
| 层 2 LOGO | 选中块里谁 loud |

评估（玩具里可回收）按 **订单粒** 写死：

1. **特征回收**：FSDS 集合含种下的列，且 **不含 Y**。
2. **扫描回收**：merchant-layer `{φ≥τ}`（以及覆盖前缀）lift 后 `south_frac` 高（covariate 应接近 1）。结构 / bundled Louvain **不要求**回收。
3. **构成 vs 强度**：扫描可以不含全部南区商户（冷启动、n 不够）；看 loud 的 \(\pi\) 和 gap vs other。覆盖前缀看质量落在哪；质量加权 \(\phi\cdot n\) 压小袋。
4. **指纹**：loud 上的 mix 偏 MMD / PO / CMean；切片对不上时先读 metric slice，不要换回 Louvain。

---

## 表 5 · 泄漏清单（硬约束）

| 冻结在 \(D_{\mathrm{ref}}\) | 新窗可以用 | 永远不进 \(Z\) / 子集键 |
|---|---|---|
| σ、订单列的 median/std | 当前订单 X、当前实体 X（localize） | **Y**、label、残差分箱、未来 batch |
| 实体画像（serve 模式） | PO / CMean_Y 的 **outcome** | 用新窗 Y 重切的 HH / 小时桶去当 subset |
| 子集词表（south/north 来自商户） | 服务误差（冻住的 \(\mu_{\mathrm{ref}}\)） | 把 PO 的 \(\mu\) 当下一轮特征再选一遍 |

HH（划分）选的是图上的 merchant 扫描 / region，**不是 Y**。FSDS 的 CMean_Y 用 Y 当左边的均值，列仍是 X。

---

## 表 6 · 玩具 DGP（`make_order_graph_stream`）

8 个商户，北 0–3 / 南 4–7。Y = 转化（amount + merchant_cat + user_tenure 的 logit）。onset 后 **只污染南区订单**：

| kind | 种下的事 | 订单粒该看到 |
|---|---|---|
| `covariate_south` | amount、channel、merchant_gmv 走；\(f\) 不变 | FSDS 含 amount / channel / merchant_gmv；loud=south；π_MMD、\(\|\Delta\mu_X\|\)、gap_MMD>0；\(\Delta\mu_Y\) 可能跟着走；fingerprint `x_shift` |
| `concept_south` | 南区 amount 系数翻号并下移截距；\(P(X)\) 固定 | FSDS 含 amount（靠 CMean 关联差，不是靠 MMD）；loud=south；π_PO 和 \(|\Delta\mu_Y|\) 南>北；MMD 份额被门控掉；fingerprint `concept` |
| `both` | 两件同时 | 后期 MMD 往往盖过 PO（X 走得更响）；肖像里仍要读 PO 和 \(\Delta\mu_Y\)，不要只看 mix |

实测表：`results/graph_fsds_localize/TABLES.md`。特征库看板：`results/graph_fsds_localize/library.html`（`library_heatmap.png` / `library_time.png` / `merchant_scan.png`）。

---

## 表 8 · 特征库（这个 dataset 上的可视化）

九列、三粒、不含 Y。看板把 FSDS 分数贴回目录，而不是画一张网络。

| grain | 列 | 角色 | covariate 种下 | concept 种下 |
|---|---|---|---|---|
| order | amount | 进 logit | \(P(X)\) | \(P(Y\mid X)\) |
| order | channel | 渠道 | \(P(X)\) | — |
| order | hour, n_items | 噪声 | — | — |
| merchant | merchant_gmv | GMV | \(P(X)\) | — |
| merchant | merchant_cat | 进 logit，不种 shift | — | — |
| merchant | n_skus | 噪声 | — | — |
| user | user_tenure | 进 logit，不种 shift | — | — |
| user | user_hist_freq | 噪声 | — | — |

读法：种下的列应变红、进 FSDS；噪声列保持淡。商户扫描条形图红 = 南区，● = 进 \(\{φ\geτ\}\)。**没有边、没有社区。**

---

## 表 7 · 接到冻层看板

| 这里看到 | 看板格子 | 不要做什么 |
|---|---|---|
| 南区 MMD / CMean_X loud，全局 PO×MSE 还 keep | 全局可 keep，**子集已经响** | 不要等全窗 2× 才报南区 |
| 南区 PO loud，服务 MSE 还没崩 | watch | 不要按 concept 去冻底 |
| 南区 MMD loud + 层 2 mix_MMD 在 order 块 | X-shift 落在订单特征 | 不要冻顶；可 stem / 重加权南区 |
| serve 模式画像 | clever covariate | 不要把 localize 的当前商户均值写进生产特征再训 Y |

代码：

- `Python/src/stream_dgps.py` — `make_order_graph_stream` + `FEATURE_LIBRARY`
- `Python/src/graph_fsds_localize.py` — 特征库、FSDS、unify、subset vs other、肖像
- `Python/src/order_graph_nx.py` — subset scan（水平集 / 覆盖 / \(\phi\cdot n\)）；Louvain 只做对照
- `scripts/run_graph_fsds_localize.py` — `library.html` 看板
- `tests/test_graph_fsds_localize.py`

```bash
PYTHONPATH=Python/src:. python3 -m unittest tests.test_graph_fsds_localize
PYTHONPATH=Python/src:. python3 scripts/run_graph_fsds_localize.py
```
