# 买完会不会带动后续点击

last-touch（backward asof）= 买之前怎么点过来。
这块（forward asof）= **买完还会不会点**。`trans_cnv_to_exp` 只是邻接 Markov，不够。

跑：`python3 scripts/tencent_gr/run_post_cnv.py`（写本文件的数字段）。

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
  钱    price, log1p_price
  滞后  n_prior_cnv, lag_post_clk_1d_rate                 # 此前各单的买后点击率，shift(1)
不准进 X：n_clk_after_* / lift_* / dt_next / y_post_*
```

为什么要这些：

- **买前量**：人本来就爱点，买后也会点。这是必须先控的基线。
- **路径**：冲动（within 5m / wo_prior_clk）vs 长犹豫。前者更容易当场续点，后者更像买完走人。
- **当场深度**：30min 场还没关，下一击几乎是续逛，不是「购买带动」。
- **滞后买后率**：这个人以前买完爱不爱点——用户倾向，给下一单用。
- **价格**：贵的可能回去反复看；便宜的买完即走。要数据说话，不先验锁方向。

用户表只并 **历史倾向**（`post_clk_1d_rate` 等），给 `future_cnv` 那类用户粒任务。
预测「这一单之后」必须停在转化粒。

## 这批 prefix

先跑 `run_post_cnv.py` 填表。
