# Attribution / one-pass / FSDS 选特征

## 1. 为什么该用 `merge_asof`

现在的 dict 扫一遍 = last-touch as-of join。事件表上就该写成：

```python
# 同 item last-click
pd.merge_asof(
    cnv.sort_values("ts"), clk.sort_values("ts"),
    by=["user_id", "item_id"], left_on="ts", right_on="ts",
    direction="backward",
)
# 任意 last-click：by=["user_id"] 即可
# clk_ts 空 → cnv_wo_prior_clk
```

嵌套 `seq` 里 per-user for-loop 等价；全表 / Spark 必须 asof，不要 Python 状态机。

## 2. One-pass（不要 9 次 funnel）

`t_end` 已知后 **一趟**：

```
age = t_end - ts
if 0 <= age < W:  count[W][act] += 1          # 所有窗
dec[hl][act] += exp(-log(2)*age/hl)
session / last_act→trans / iid 频次 / last_clk[iid]
```

交叉 `x_a__b` 在 pass 后再乘。别再对每个窗口扫全序列。

## 3. Attribution 刻画什么

不是因果效应。是两个 last-touch 时延：

| 量 | 定义 | 业务 |
|---|---|---|
| `clk2cnv` | 同 item 上次 clk → cnv（分钟） | 这件商品的决策时长 |
| `anyclk2cnv` | 全局上次任意 clk → cnv | 跨品浏览后成交 |
| `cnv_wo_prior_clk` | 该 item 从未 clk | 漏点 / 直达 / 曝光即买 |

没有 ignorability，不能当 CATE。只是路径描述。

## 4. 怎么评估（FSDS，不是 F-score）

F-score = 同样本 Cor(X,Y)，转化用户的计数特征必然头名。

该用 CFPerm/PO-risk 选 **漂移驱动**：

- 时钟切两批 `W`（`seq_t_end` 中位数）
- `Y = future_cnv`
- **RF-domain**：P(W\|X) VIMP → P(X) / covariate
- **PO-risk**：φ=(Y-μ)(W-e)，R=E[τ̂(X)²]，LOGO Δ_g=R-R^{(-g)} → 概念相关
- Kendall τ(F-score rank, PO-VIMP rank)：低 ⇒ 150 板不是 OOD 板

Attribution 家族：LOGO Δ 大 = 决策时长/漏点结构在批间变了；Δ≈0 = 只是水平关联，F-score 抬它没用。

## 5. 这 150 上的结果（时钟 = seq_t_end 中位切半）

- RF-domain AUC **0.753**：P(X) 有漂（后来的用户行为画像不同）
- PO-risk R **0.0018**：P(Y|X) 几乎没跳（和 rfperm fire=0 一致）
- Kendall τ(PO-VIMP, F-score) **0.05**：F 板 ≠ OOD 板
- attr LOGO Δ **负**：决策时延不是漂移源；F-score 把它排前面只是「转化过的人还会转化」

下一步：family LOGO → 组内再 CFPerm/LOCO 到单特征；别再堆 1000 维 F 选。
