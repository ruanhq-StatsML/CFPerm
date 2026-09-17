# 特征库监控：哪几列在走、怎么走、落在哪些商户

给运营 / 风控的一页说明。**定位，不是唯一分解。** Y 只当转化结果，不当特征。T 是 batch 标签。份额是相对强弱，不是 Shapley，也不是因果效应。

全局看板（PO × MSE × MMD）回答「整窗要不要冻层」。特征库回答下一句：**具体是哪几列、落在哪些商户。** 南区已经在走的时候，全窗仍可能 keep training——库就是为了把这件事写出来。

数据：合成订单流。冻结参考窗 480 单，之后每批 160 单 × 3 批，onset = 第 1 批。16 个商户，后一半种成「南区」。三套场景：只动 \(P(X)\)（covariate）、只动 \(P(Y\mid X)\)（concept）、两头一起（both）。目录九列，分订单 / 商户 / 用户三粒。**不是图方法**：`merchant_id` 是外键，用来把「哪些商户」对回订单。

看板：`library.html`。

---

## 三张图各回答什么

**图 1 · `library_heatmap.png` —— 各 feature 的 shift 程度（就是这个）**

行 = 目录里的列，列 = 三支读数 + 综合分。颜色越红，相对参考窗越响。

| 格子 | 问什么 |
|---|---|
| MMD | 这一列的分布 \(P(X)\) 走了没有 |
| \(\|\Delta\mu_X\|\) | 这一列的均值走了多少 |
| \(\|\Delta\mu_Y\|\) | 按参考窗切开后，\(Y\) 的均值差（concept 会亮；covariate 也可能跟着亮，因为 \(f\) 不是常数） |
| score | FSDS 综合分；蓝框 / `selected` = 这一窗真正拿去刻画 subset 的列 |

`*` = 种下的列。covariate 最后一批：amount 6.74、channel 2.34、merchant_gmv 3.77，噪声列 < 0.8。concept：amount 靠 \(|\Delta\mu_Y|=0.17\) 亮，MMD 接近 0。

**图 2 · `library_time.png` —— 程度随 batch 怎么抬**

横轴 batch，纵轴 FSDS score。粗线 = 种下的列。onset 前应接近安静；onset 后 covariate / both 上金额、渠道、GMV 一起抬，concept 只有 amount 缓升。用来看**何时开始、是否还在走**。

**图 3 · `merchant_scan.png` —— 不是 feature 的 shift 程度。是哪些商户在走。**

每根条是一个商户相对**自己**参考袋的响度 \(\phi\)。红 = 南区，蓝 = 北区。虚线 = 门槛 \(\tau\)，● = 过线（loud）。covariate 最后一批：过线的几乎全是南区（订单 south_frac = 0.85）。concept 更散，和「\(X\) 没动、只有 \(Y\mid X\) 在 hop」一致。

前两张管 **哪几列**；第三张管 **哪几个商户**。不要把第三张读成 feature importance。

---

## 特征库怎么刻画

目录先钉死，再打分。Y 不进目录。

| 粒 | 列 | 角色 |
|---|---|---|
| 订单 | amount, channel | 金额（进转化）、渠道 |
| 订单 | hour, n_items | 噪声 |
| 商户 | merchant_gmv | GMV |
| 商户 | merchant_cat, n_skus | 类目进转化但不种 shift；SKU 噪声 |
| 用户 | user_tenure, user_hist_freq | tenure 进转化但不种 shift；频次噪声 |

每一窗、每一列报四件事：三支读数、是否 loud、是否被 FSDS 选中、（玩具里）是否种下。选中的列收到**同一粒（订单）**再打 subset 肖像：loud vs other 的 \(\pi_{\mathrm{MMD}},\pi_{\mathrm{PO}},\pi_{\mathrm{CMean}}\)。不要把商户 MMD 和订单 MMD 丢进同一个份额。

读法（仍是定位）：

| 看到 | 刻画 |
|---|---|
| MMD / \(\|\Delta\mu_X\|\) 亮，PO 安静 | 这列的客群 / 金额结构在走 |
| MMD 安静，\(\|\Delta\mu_Y\|\) / PO 亮 | 这列上 \(Y\mid X\) 在 hop |
| 三支都响 | both：报「这列贡献大」，不要拆成唯一的 concept/covariate 份额 |

---

## 怎么评估

玩具里种下了答案，所以可以回收。线上没有种下，就看「选中列对不对得上业务直觉」和「loud 商户的 gap vs other」。

| 检查 | 通过长什么样 | 这批玩具 |
|---|---|---|
| 泄漏 | 目录和选中列都没有 Y | 通过 |
| 特征回收 | 种下的列进 FSDS；噪声列不进 | covariate：amount, channel, merchant_gmv。hour / n_skus / hist_freq 未选 |
| 方向 | covariate 亮 MMD；concept 亮 \(\Delta\mu_Y\)，MMD 接近 0 | 见图 1 |
| 时机 | onset 前安静，之后种下列抬升 | 见图 2 |
| 商户回收 | loud 订单的南区比例高；用户侧不应回收（用户跨商户随机挂） | covariate \(J=0.77\)，south_frac \(=0.85\)；用户 Jaccard \(=0\) |
| 补集 | loud − other 的 MMD / \(\Delta\mu\) 应为正 | covariate 最后一批 loud MMD \(=0.38\)，other \(\approx 0\) |

线上用法：参考窗冻结 σ 和分箱；新商户没有自己的历史 → 不进商户排名（冷启动），订单仍可挂 other。不要用「和全局平均不像」当「这个商户变了」——那是异质性。

接到冻层看板：库已经响、全局还 keep → **子集先报**，不要等全窗 2×。南区 MMD loud 且选中的是订单金额 / 渠道 → 更像 X-shift，可重加权或 stem，不要按 concept 去冻底。

---

## 为什么有用

1. **全窗会漏子集。** 一半商户在走，混合 MMD 可能还没到冻层阈值。库按列、按商户切开。
2. **列和商户是两件不同的事。** 金额走了，和南区走了，要分开报。图 1/2 vs 图 3。
3. **口径和现有看板同一口钟。** 还是 MMD / PO / Conditional Mean vs \(D_{\mathrm{ref}}\)。没有新的因果故事，没有图社区。
4. **算得起。** 九列各打一次一维 MMD；商户排名是 16 个数排个序。

局限（写进交付）：这是定位，不是「金额贡献了 60% 的漂移」那种唯一分解；concept 比 covariate 难收（\(X\) 不动，只能靠 \(Y\) 的一阶矩）；用户粒在这个 DGP 里本来就不该回收南区。

数字细表：`TABLES.md`。
