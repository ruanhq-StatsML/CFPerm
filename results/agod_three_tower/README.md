# 三塔原型：类目频率、评论次数、负样本池 → AGOD

不做特征系统。在线状态只有因果画像：`cat_hist`、`review_cnt`、价格 μ/σ。三塔是

| 塔 | 输入（冻住的） | 向量 |
|---|---|---|
| **U user** | `cat_hist` + `review_cnt` | 类目单纯形向均匀先验收缩，收缩量由 n 决定 |
| **I item** | 当前商品类目 + 价格 | one-hot + 价格，目录不是评分 |
| **N neg** | 同一 item 塔，seen-excluded 池 | 他们 `CausalPosNeg` 的负样本 |

分数 `s = cos(U,I) − mean_k cos(U,N_k)`。AGOD **不**把原始 {U,I,N} 当 GLS 库（holdout 在它们张成的空间里会在 item 顶点塌掉）。票是正交打包的标量，holdout 只看 y，再 `π→LR`。顺序仍是 **先冻住画像再计分，最后才 update Welford**。

---

## 类目频率怎么刻画

`cat_hist` 是用户塔在类目单纯形上的**方向**。和当前商品的对齐就是

```
match = ⟨cat_hist, onehot(item.cat)⟩
```

也就是「过去有多大比例的评论落在这个类目」。正样本 match 应高于负样本；same-cat hard 负样本会把这个差压小——这就是负样本池硬度。

终盘画像：`{'Tools': 0.888889, 'Sports': 0.111111}`，`review_cnt=9`。

---

## 评论次数怎么刻画

`review_cnt` 不是又一个稀疏 ID，是用户塔的**置信度**：

```
conf = 1 − exp(−n / n0)
U = (1−conf)·uniform + conf·cat_hist
LR_damped = 1 + conf · (LR − 1)
```

n→0：U 塌成均匀，AGOD 不许猛调 LR（冷启动）。n 大：相信类目方向，才允许 π 离开 1/3、LR 离开 1。这和 stacking 里 `N_eff` / Kish 是同一类量。

---

## 负样本池

| pool | pos−neg gap | in-seen | cat_match neg |
|---|---:|---:|---:|
| random global | 0.448 | 0.00 | 0.19 |
| seen-excluded | 0.386 | 0.00 | 0.31 |
| same-cat hard | 0.021 | 0.00 | 0.89 |

他们流水线是 seen-excluded + 全局 10k 池（上架时间不过滤）。hard 负样本 gap 更小，才是在考 `cat_hist` 而不是考「类目不同」。

---

## 回到 AGOD 对齐

三塔向量负责 **打分**。AGOD 的 expert 票是打到正交轴上的标量（否则 GLS 会在 item 顶点塌掉）。holdout 只看 y：

```
π = GLS / direction-match(vecs, g_hold)
LR = π_to_lr(π)   # 再乘 conf
```

holdout 是 **y 决定的 oracle 轴**（正样本要 user/item 轴，负样本要 neg 轴），不是把 {U,I,N} 再投影回去——那会在 item 上塌成顶点。expert 票是正交打包的标量：user=`⟨cat_hist, pos.cat⟩`，item=`cos(U,I)`，neg=`cos(U,N)`。所以 **s_user 跟着类目频率走**，π 是 GLS 组合，LR 再乘 `conf(review_cnt)`。终步 `π`：user=0.71, item=0.29, neg=0.00。这不是 MMoE：门不看 x，看 holdout。MMoE 的 `softmax(W mean(experts))` 只对照，不进这条环。

```bash
PYTHONPATH=. python3 -m tests.test_agod_three_tower
PYTHONPATH=. python3 scripts/run_agod_three_tower.py
```
