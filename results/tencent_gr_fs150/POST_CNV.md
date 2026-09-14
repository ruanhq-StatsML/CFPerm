# 买完会不会带动后续点击

计算（状态机 ≡ asof）：`results/tencent_gr_fs150/POST_CALC.md`，扫 `scripts/tencent_gr/onepass_post.py`。

last-touch（backward asof）= 买之前怎么点过来。
这块（forward asof）= **买完还会不会点**。`trans_cnv_to_exp` 只是邻接 Markov，不够。

## 中间表

```
ev
 ├ attr  = merge_asof(cnv, clk, backward)   # 同品/任意 last-click
 ├ post  = merge_asof(cnv, clk, forward)    # 下一次点击 + 同长窗 before/after 计数
 │          cum asof: n_after = C(t+W)-C(t), n_before = C(t)-C(t-W)
 └ user  = 未删失转化上的 rate / lift p50，再和 funnel/decay/sess merge
```

转化粒。删失：窗内已点到 → 1；跟满 W 没点 → 0；跟不满 → 不当 0。

## 刻画什么（不是因果）

| 量 | 公式 | 在问什么 |
|---|---|---|
| `y_post_clk_T` | 1{任意下一次点击 ≤ T} | 买完人还在不在点（平台活跃） |
| `y_post_same_T` | 1{同品下一次点击 ≤ T} | 这件商品买完还看不看（复访/晒单/后悔） |
| `lift_T` | n_after / (n_before+1) | 相对买前同长窗，点击量抬没抬 |
| `next_clk_same_sess` | dt_next ≤ 30min | 当场续点，还是隔场回访 |

任意点击的 lift 会被「本来就在逛」污染：连着买两单时，后一单的预热点会算进前一单的 after。
**同品 lift / y_post_same** 干净得多。当场续点用 5m/30min，回访用 1d/7d。
这不是 CATE：没有对照、没有 ignorability。只是路径描述。要预测的是「这一单之后会不会点」，不是「买导致多点」。

## 预测：Y 与 X 必须切开

```
Y = y_post_clk_1d   # 或 y_post_same_1d；只用跟满 1d 的转化
X 只能是 cnv_ts 已知：
  路径  wo_prior_clk, dt_item/any, item_within_5m/1h     # 已有 attr 表
  买前量 n_clk_before_{1h,1d,7d}, n_clk_same_before      # 基线活跃，不是 Y
  当场  sess_pos, sess_clk_before                         # 热场续点的主混杂
  钱    log1p_price
  滞后  n_prior_cnv, lag_post_clk_1d_rate                 # 此前各单的买后点击率，shift(1)
不准进 X：n_clk_after_* / lift_* / dt_next / y_post_* / 原始 price（和 log1p 共线）
```

为什么要这些：

- **买前量**：人本来就爱点，买后也会点。这是必须先控的基线。
- **路径**：冲动（within 5m / wo_prior_clk）vs 长犹豫。前者更容易当场续点，后者更像买完走人。
- **当场深度**：30min 场还没关，下一击几乎是续逛，不是「购买带动」。
- **滞后买后率**：这个人以前买完爱不爱点——用户倾向，给下一单用。
- **价格**：贵的可能回去反复看；便宜的买完即走。要数据说话，不先验锁方向。

用户表只并 **历史倾向**（`post_clk_1d_rate` 等），给 `future_cnv` 那类用户粒任务。
预测「这一单之后」必须停在转化粒。

## 这批 prefix（n_cnv=12866, n_user=3273）

| 窗 | 未删失 | 任意点击率 | 同品点击率 |
|---|---:|---:|---:|
| 5m | 12767 | 0.005 | 0.000 |
| 1h | 12746 | 0.023 | 0.000 |
| 1d | 12564 | 0.113 | 0.001 |
| 7d | 11822 | 0.412 | 0.001 |

1d lift p50=0.000 mean=0.114；P(lift>1)=0.014；P(n_after>n_before)=0.089。
任意：before 0.136 → after 0.135；
同品：before 0.001 → after 0.001。
有下一次点击时，P(仍在当场 30min)=0.023；dt_next p50=5.8 d。

读这批：买完 **不抬** 后续点击。1d after≈before，lift 中位数 0，P(after>before)=0.09。
5m/1h 几乎不点；同品 1d 只有约千分之一（这件商品买完基本不再点它）。
7d 任意 0.41、dt_next 中位约 6 天 → 那是人还在平台上逛，不是购买带动的余热。
同品 last-click 在这批对得上的极少，所以 attr 同品桶对「买后任意点击」几乎没信息；
有用的是任意点击间隔、7d 买前量和滞后买后率。

## 预测 Y=y_post_clk_1d（按用户 70/30）

| X | HGB AUC | LogReg AUC | HGB AP |
|---|---:|---:|---:|
| volume | 0.678 | 0.678 | 0.274 |
| +path | 0.689 | 0.670 | 0.278 |
| +lag | 0.728 | 0.708 | 0.304 |
| +sess_money_lag | 0.745 | 0.713 | 0.371 |
| Y=clk_1h allX | 0.708 | 0.670 | 0.113 |
| Y=same_1d allX | nan | nan | nan |

LogReg 系数（标准化后，Y=任意 1d 点击；>0 更像买完还点）：

| feat | coef |
|---|---:|
| `lag_post_clk_1d_rate` | +0.370 |
| `n_clk_before_7d` | +0.334 |
| `dt_any_min` | -0.111 |
| `item_within_1h` | -0.053 |
| `item_within_5m` | -0.052 |
| `n_clk_before_1d` | +0.048 |
| `n_prior_cnv` | +0.039 |
| `sess_pos` | -0.033 |
| `sess_clk_before` | +0.030 |
| `dt_item_min` | -0.021 |
| `log1p_price` | +0.012 |
| `n_clk_before_1h` | -0.005 |
| `wo_prior_clk` | +0.000 |
| `n_clk_same_before` | -0.000 |

volume-only 已经能到 ~0.68：1d 任意点击主要是「本来就爱点的人还在点」。
+lag 再涨，说明「这个人以前买完爱不爱点」是下一单的主信号。
path/sess/price 几乎不再涨：同品买后点击近乎 0，1d Y 也不是当场续点。
同品 Y 这批只有十几正例，不够建模——要刻画「买完还看这件」先承认事件极稀。
