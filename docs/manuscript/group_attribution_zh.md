# 多组别归因（CFPerm）：permute \(T\)，不是 permute \(X\)

闸问的是 **时间上** \(P(Y\mid X)\) 还在不在。这一件是另一张表：

\[
T \in \{0,1,\ldots,K-1\}
\quad\text{是组别（队列 / 检索路 / 第几跳）}
\]

问的是哪些 \(X_j\) 带着 **组间** 的 \(P(Y\mid X)\) 差，不是哪一维造成了时间 hop。13 个审核 \(x_*\) 是组别对比的归因坐标，不是时间闸的归因维。HH chosen 不是 \(Y\)。

## 逻辑

```
1. 观察当前包，写下这一跳的 Y（过/不过，或这一路是否可用）
2. 组别 T 是队列 / 通道 / hop 桶，不是时间 batch
3. 冻住的组间图：importance_j = 组 × X_j 交互有多大
4. permute T（打乱组别，不打乱 X），得到 null importance
5. 特征 p = #{null ≥ 观察} / B
   跨特征阈值 = 各组 null 上尾的再上尾
   至少 top_k 维超过阈值 → reject
```

和 `R/CFPerm_vimp.R` 同一套决定。Prototype 用线性组×X 交互，只为把 permute-\(T\) 的循环看清楚。换成 permuCATE / LOCO / GRF 是同一个循环、换 vimp。

## 盘上的组别

| 流 | \(T\) | \(Y\) |
|---|---|---|
| HH 两路队列 | helpful / harmless | 审核过/不过 |
| 真标签两路 | BeaverTails / ToxicChat | 真审核标签 |
| 混合检索两路 | sparse 可用 / dense 可用 | 这一路 gold 是否覆盖 |
| 多步审核三组 | hop 0 / hop 1 / hop 2+ | 这一跳过/不过 |

```bash
PYTHONPATH=. python3 scripts/prototype_group_attribution.py
```

Post-hoc 子群定位：把 subset indices 拽出来，pairwise 打 **MMD** 和 **PO-risk**（不只是条件均值）。见 `scripts/posthoc_localization.py`。

Planted：三组里只有组 2 依赖 `x1`，应 reject 且 top 是 `x1`。Null：\(Y\) 只靠 \(X\)、不靠 \(T\)，应 quiet。

## 和其他对象怎么分开

| 对象 | 问什么 | 动哪一列 |
|---|---|---|
| 在线 serving 闸 | \(P(Y\mid X)\) 是否随时间 hop | 时间 `batch` |
| 组别归因 CFPerm | 哪些 \(X_j\) 带着组间差 | 组别 \(T\) |

不要把 CFPerm 的命中维读成「所以时间闸是这维引起的」。
