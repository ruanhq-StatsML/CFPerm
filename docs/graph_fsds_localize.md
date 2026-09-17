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
        层 1  subset vs other：默认 **bundled Louvain 社区**（变动子集）；region / structural 只做对照
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

图在这里就是实体图：`order → merchant`、`order → user`。包是 **networkx**（不是 GraphRAG / PyG / DGL）。社区发现有两套目标，不能混。

---

## 表 0 · 用哪个包、社区发现到底在优化什么

| | 结构 Louvain | **Bundled-shift Louvain（默认）** |
|---|---|---|
| 包 | `networkx==3.6.x` | 同一个包 |
| 求解器 | `nx.community.louvain_communities` | **同一个求解器** |
| 图 | 共单（merchant–user）+ 当前窗 kNN(X) | 商户节点；边 = 三支捆的**同方向亲和** |
| 目标 \(Q\) | 模块度：谁和谁连得密 | **不是**模块度语义。边是「一样在变」，切出来的社区是变动 subset 的代理 |
| 种下 south、X 走（covariate） | kNN(X) 可能碰巧把南区聚在一起 | 应用目标：loud 社区的 south_frac 高 |
| 种下 south、X 不动（concept） | 结构图几乎不变，**切不出变动子集** | 捆里的 \(\Delta\mu_Y\) / PO 把南区商户连在一起 |
| GraphSAGE | 冻住的 mean-aggregator 读出，不是切 | 不用它来切 |
| GraphRAG | 不用 | 不用 |

Vanilla Louvain 最大化

\[
Q=\frac{1}{2m}\sum_{ij}\Big(A_{ij}-\frac{k_i k_j}{2m}\Big)\delta(c_i,c_j).
\]

这回答「稠密子图」。共单图的稠密块 = 谁共享用户，**不是**谁的 \(P(X)\) / \(P(Y\mid X)\) 在走。

Bundled 图把多指标捆进边：

\[
v_i=\big(\mathrm{MMD}_i,\;\|\Delta\mu_X\|_i,\;|\Delta\mu_Y|_i,\;\mathrm{PO}_i\big)
\quad\text{（相对 }D_{\mathrm{ref}}\text{，节点 = 商户的订单袋）}
\]

\[
w_{ij}=\mathrm{ReLU}\big(\cos(\tilde v_i,\tilde v_j)\big)\cdot\mathbf{1}[\phi_i\ge\tau\text{ 或 }\phi_j\ge\tau],
\quad \phi_i=\mathbf{1}^\top \tilde v_i.
\]

安静的商户孤立。Louvain 在这张图上切 = 把**同方向 loud** 的点收成一块。然后两层归因仍在订单粒上对社区打三支。

Y 进 \(\Delta\mu_Y\) / PO 是监控 outcome，**不进**节点特征、不进 \(Z\)。这不是泄漏进服务模型；用 Y 去找「哪一块转化在 hop」就是 localization 本身。不要用 GraphRAG 合成图来代替这个目标。

**Multi-bundled 三层：**

| 层 | 对象 | 捆 |
|---|---|---|
| 节点 | 每个商户的订单袋 vs \(D_{\mathrm{ref}}\) | \(v_i\) |
| 边 | 同方向亲和 \(w_{ij}\) | \(\cos(\tilde v_i,\tilde v_j)\) |
| 社区 | 社区内订单并起来，再打三支 + \(\pi\) vs other | 层 1 归因 |

求解器用 networkx Louvain，是因为它现成、无 torch、和冻层看板一样不引入可训练 GNN（可训练 GNN 会把 encoder drift 和 graph drift 糊在一起）。目标换了，求解器没换。

---

## 表 0b · 怎么定位变动 subset：own-ref、同一维度、多层图怎么评

### 节点分数必须对「商户自己的 \(D_{\mathrm{ref}}\)」

| 钟 | 公式（节点 \(i\) = 一个商户的订单袋） | 问什么 | 会误伤什么 |
|---|---|---|---|
| **own-ref** | \(\mathrm{MMD}(X_{\mathrm{new}}(i), X_{\mathrm{ref}}(i))\)，\(\Delta\mu\) 同样用 \(i\) 自己的 ref | **这个商户自己变了没有** | 需要 \(n_{\mathrm{ref}}(i)\ge\min n\)；新商户没有自己的钟 → **跳过**，不进 bundled 图 |
| **full-ref** | \(\mathrm{MMD}(X_{\mathrm{new}}(i), X_{\mathrm{ref}})\) | 这个商户和**全局混合**不像 | 从来就小众的商户永远 loud（异质性，不是 drift） |

变动 subset 的切图只用 own-ref。full-ref 留着报「对整窗的贡献」（构成变化：南区单变多）。两口钟不要合成一个分数再切。

冷启动（新商户没有 own-ref）：不把 full-ref 塞进切图，否则又变成异质性。订单仍可挂 `C_quiet`，两层归因时单独看。

σ 钉在**该节点自己的 ref 袋**上（own）或全局 \(D_{\mathrm{ref}}\)（full），都冻结，不在新窗重估。

### 同一维度怎么 localize

只在一个粒上切、只在这个粒上比份额：订单就订单，商户就商户。\(\pi_{\mathrm{MMD}}\) 的分母是这个粒上的社区。不要把商户节点 MMD 和订单行 MMD 丢进同一个 simplex。

同一维度的评估：loud 社区 lift 成订单集合 \(\hat S\)，和种下的南区订单 \(S^\star\) 算 Jaccard；再在 \(\hat S\) 上打三支 + vs other。

### 多个维度 / 多层图怎么评估

图是分层的：商户层、用户层、（可选）订单 kNN 层。**每一层自己切**，切完 **全部 lift 到订单粒** 再比。这是唯一可加的评估单位。

```
商户层 bundled Louvain  →  订单集合 Ŝ_m
用户层 bundled Louvain  →  订单集合 Ŝ_u
订单粒 oracle（region）→  订单集合 S*
评估：J(Ŝ_m, S*)、J(Ŝ_u, S*)、J(Ŝ_m, Ŝ_u)
三支肖像只在 Ŝ_m / Ŝ_u 上打（已经是同一粒）
```

| 看到 | 读法 |
|---|---|
| \(J(\hat S_m,S^\star)\) 高、\(J(\hat S_u,S^\star)\) 低 | 变动落在商户层；用户只是被共单带着走 |
| 两头 Jaccard 都高、层间 Jaccard 也高 | 两层指到同一批订单（仍不是唯一分解） |
| 层间 Jaccard 低、两头对 \(S^\star\) 都不高 | 切错层或 own-ref 样本不够 |
| 把两层的 \(v_i\) 直接加起来再 Louvain | **不要**。量纲和 n 都不同 |

多层不是 GraphRAG 合成一张超图。就是：每层一个 bundled 图、一个 Louvain、lift 到订单、Jaccard。用户层在这个 DGP 里**不该**回收南区——用户是跨商户随机挂的。这正是评估：层对了才回收。

### 还有什么（同一套钟上的附件，不是新的面）

| 件 | 做什么 | 不要做成 |
|---|---|---|
| min_n | own-ref / 社区 PO 的地板 | 小袋上打 RF PO |
| 强度 vs 质量 | \(\phi_i\) 高但 \(n_i\) 小 → 报 n 和 share | 只按 \(\phi\) 排序当贡献 |
| vs other | 社区 − 补集，两头对同一口 ref | 当成 Shapley |
| 结构 Louvain | 对照：模块度切不出 concept | 当变动 subset |
| 冻住 SAGE-mean | 结构读出 | 当切 |
| 新商户 | C_quiet | 用全局均值假造 own-ref 再切 |

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

评估只在 **同一个 \(Z\)** 上比。默认子集键是 **bundled Louvain 社区**（变动 subset 的代理）。`region` 是种下的对照，不是切图目标。`structural` 是共单模块度，用来说明「切错目标会切不出 concept」。

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
2. **社区回收**：bundled Louvain 的 loud 社区 `south_frac` 高（covariate 应接近 1）。结构 Louvain 对 concept **不要求**回收。
3. **构成 vs 强度**：社区可以碎成几块；看 loud 那一块的 \(\pi\) 和 gap vs other，不要求全部南区商户在同一社区。
4. **指纹**：社区上的 mix 偏 MMD / PO / CMean。

---

## 表 5 · 泄漏清单（硬约束）

| 冻结在 \(D_{\mathrm{ref}}\) | 新窗可以用 | 永远不进 \(Z\) / 子集键 |
|---|---|---|
| σ、订单列的 median/std | 当前订单 X、当前实体 X（localize） | **Y**、label、残差分箱、未来 batch |
| 实体画像（serve 模式） | PO / CMean_Y 的 **outcome** | 用新窗 Y 重切的 HH / 小时桶去当 subset |
| 子集词表（south/north 来自商户） | 服务误差（冻住的 \(\mu_{\mathrm{ref}}\)） | 把 PO 的 \(\mu\) 当下一轮特征再选一遍 |

HH（划分）选的是图上的 merchant/region，**不是 Y**。FSDS 的 CMean_Y 用 Y 当左边的均值，列仍是 X。

---

## 表 6 · 玩具 DGP（`make_order_graph_stream`）

8 个商户，北 0–3 / 南 4–7。Y = 转化（amount + merchant_cat + user_tenure 的 logit）。onset 后 **只污染南区订单**：

| kind | 种下的事 | 订单粒该看到 |
|---|---|---|
| `covariate_south` | amount、channel、merchant_gmv 走；\(f\) 不变 | FSDS 含 amount / channel / merchant_gmv；loud=south；π_MMD、\(\|\Delta\mu_X\|\)、gap_MMD>0；\(\Delta\mu_Y\) 可能跟着走；fingerprint `x_shift` |
| `concept_south` | 南区 amount 系数翻号并下移截距；\(P(X)\) 固定 | FSDS 含 amount（靠 CMean 关联差，不是靠 MMD）；loud=south；π_PO 和 \(|\Delta\mu_Y|\) 南>北；MMD 份额被门控掉；fingerprint `concept` |
| `both` | 两件同时 | 后期 MMD 往往盖过 PO（X 走得更响）；肖像里仍要读 PO 和 \(\Delta\mu_Y\)，不要只看 mix |

实测表：`results/graph_fsds_localize/TABLES.md`。图：`subset_shares.png`。

---

## 表 7 · 接到冻层看板

| 这里看到 | 看板格子 | 不要做什么 |
|---|---|---|
| 南区 MMD / CMean_X loud，全局 PO×MSE 还 keep | 全局可 keep，**子集已经响** | 不要等全窗 2× 才报南区 |
| 南区 PO loud，服务 MSE 还没崩 | watch | 不要按 concept 去冻底 |
| 南区 MMD loud + 层 2 mix_MMD 在 order 块 | X-shift 落在订单特征 | 不要冻顶；可 stem / 重加权南区 |
| serve 模式画像 | clever covariate | 不要把 localize 的当前商户均值写进生产特征再训 Y |

代码：

- `Python/src/order_graph_nx.py` — networkx 图、结构 Louvain、bundled-shift Louvain、冻住 SAGE-mean
- `Python/src/graph_fsds_localize.py` — 维度定位、FSDS、unify、subset vs other、肖像
- `Python/src/stream_dgps.py` — `make_order_graph_stream`（接着用这批合成订单）
- `scripts/run_graph_fsds_localize.py`
- `tests/test_graph_fsds_localize.py`

```bash
PYTHONPATH=Python/src:. python3 -m unittest tests.test_graph_fsds_localize
PYTHONPATH=Python/src:. python3 scripts/run_graph_fsds_localize.py
```
