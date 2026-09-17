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
        层 1  subset vs other：south/north（图上的割）
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

图在这里就是实体图：`order → merchant`、`order → user`。**节点定位** = 挂在这个点上的订单 vs \(D_{\mathrm{ref}}\)（以及 vs 该点自己的 ref）。不是 GraphSAGE，不把 GNN 当新的面。

---

## 表 2 · 三支读数（同一口 \(D_{\mathrm{ref}}\) 钟）

| 量 | 符号 | 问什么 | 口径 | 不是什么 |
|---|---|---|---|---|
| MMD | \(\mathrm{MMD}^2(Z_{\mathrm{new}}(S), Z_{\mathrm{ref}})\) 或 vs **own-ref** \(Z_{\mathrm{ref}}(S)\) | 这个子集的 \(P(X)\) 走了没有 | RBF，σ 钉在该空间的 \(D_{\mathrm{ref}}\) | 不是对历史 pairwise |
| PO-risk | \(\varphi=(Y-\mu)(T-e)\) 的风险 | 这个子集的 \(P(Y\mid X)\) 有没有 hop | 独立 RF nuisance；T=batch 标签 | 不是 CATE(T)，不是服务 MLP |
| Conditional Mean | \(\Delta\mu_Y(S)=\mathbb{E}[Y_{\mathrm{new}}\mid S]-\mathbb{E}[Y_{\mathrm{ref}}\mid S]\)，\(\Delta\mu_X(S)=\|\mathbb{E}[Z_{\mathrm{new}}\mid S]-\mathbb{E}[Z_{\mathrm{ref}}\mid S]\|\) | **一阶矩**怎么刻画这个子集 | 分箱用 ref 的 median，不在新窗重切 | \(\Delta\mu_Y\) **不是** concept 检测器 |

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

评估只在 **同一个 \(Z\)** 上比。子集键是图上的割（`region` = south/north，来自 merchant），**不是** Y，也不是 \(Y\) 的残差分箱。

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

1. **特征回收**：FSDS 集合含种下的列（covariate → `amount`,`channel`,`merchant_gmv`；concept → `amount` 走 CMean_Y），且 **不含 Y**。
2. **子集回收**：`loud_subset = south`；\(\pi^{\mathrm{MMD}}_{\mathrm{south}}>\pi^{\mathrm{MMD}}_{\mathrm{north}}\)（covariate）或 \(|\Delta\mu_Y|\) 南区更大（concept）。
3. **构成 vs 强度**：`share_of_batch` 接近 ½ 而 gap>0 → 不是因为南区单更多，是南区内部在走。
4. **指纹**：covariate 的 mix 偏 MMD；concept 的 mix 偏 PO（PO 够分、n 够时）。对不上也不强行拆。

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
| `covariate_south` | amount、channel、merchant_gmv 走；\(f\) 不变 | FSDS 含 amount；loud=south；π_MMD、\(\|\Delta\mu_X\|\)、gap_MMD>0；\(\Delta\mu_Y\) 可能跟着走 |
| `concept_south` | 南区 amount 系数翻号；\(P(X)\) 固定 | loud=south；\(|\Delta\mu_Y|\) 南>北；MMD 不大 |
| `both` | 两件同时 | MMD 仍在南区；PO / CMean_Y 也在南区 |

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

- `Python/src/graph_fsds_localize.py` — 维度定位、FSDS、unify、subset vs other、肖像
- `Python/src/stream_dgps.py` — `make_order_graph_stream`
- `scripts/run_graph_fsds_localize.py`
- `tests/test_graph_fsds_localize.py`

```bash
PYTHONPATH=Python/src:. python3 -m unittest tests.test_graph_fsds_localize
PYTHONPATH=Python/src:. python3 scripts/run_graph_fsds_localize.py
```
