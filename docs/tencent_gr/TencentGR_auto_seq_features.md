# TencentGR-10M automated sequence features

- users: **12000**
- feature dim: **138**
- funnel: expose=0 / click=1 / conversion=2
- price proxy: item `115` (global p50=370.0)

## ARPU / 付费次数 / 优惠敏感度

| 业务概念 | 特征 | 公式 | 说明 |
|---|---|---|---|
| 付费次数 | `pay_cnt` | `#conversion` | action=2 作 pay intensity |
| ARPU | `arpu_sum/mean/p50_proxy` | 转化上 `item.115` | 客单排序强度，非字面人民币 |
| 日均 ARPU | `arpu_per_active_day_proxy` | `arpu_sum / active_days` | 去掉「活得久」混杂 |
| 优惠敏感 | `deal_cheap_cnv_share` | 转化价≤全局p50 占比 | 爱买相对便宜货 |
| 优惠敏感 | `deal_vs_self_price_share` | 低于自身成交中位价占比 | 等更便宜再买 |
| 优惠敏感 | `deal_hesitation_exp_before_cnv` | 转化前同item曝光次数 | 比价/等券 |
| 点击归因 | `attr_clk2cnv_min_*` | 同item上次点击→转化分钟 | 你要的归因时长 |

## 为什么 user 侧要进模

1. `103–110` 是官方用户先验，丢掉白扔信号。
2. 标量桶做人口/等级分桶，利于 CTR/CVR 校准。
3. list → 长度+熵 = 兴趣广度，无需明文。
4. `*_known` 缺失本身可预测（新/低活）。
5. 侧信息=静态长期，seq=动态行为，正交互补。

## Key stats

```json
{
  "pay_cnt": {
    "mean": 2.4153333333333333,
    "p50": 1.0,
    "p90": 7.0,
    "nonzero_rate": 0.5691666666666667
  },
  "arpu_sum_proxy": {
    "mean": 262.44983333333334,
    "p50": 0.0,
    "p90": 740.0,
    "nonzero_rate": 0.17258333333333334
  },
  "arpu_mean_proxy": {
    "mean": 103.58202087569337,
    "p50": 0.0,
    "p90": 534.0,
    "nonzero_rate": 0.17258333333333334
  },
  "arpu_per_active_day_proxy": {
    "mean": 6.4598759800881576,
    "p50": 0.0,
    "p90": 17.235918114143956,
    "nonzero_rate": 0.17258333333333334
  },
  "deal_cheap_cnv_share": {
    "mean": 0.0684857634951385,
    "p50": 0.0,
    "p90": 0.0,
    "nonzero_rate": 0.08366666666666667
  },
  "deal_hesitation_exp_before_cnv": {
    "mean": 0.002072271065242234,
    "p50": 0.0,
    "p90": 0.0,
    "nonzero_rate": 0.0045
  },
  "deal_vs_self_price_share": {
    "mean": 0.01240859637984638,
    "p50": 0.0,
    "p90": 0.0,
    "nonzero_rate": 0.03175
  },
  "attr_clk2cnv_min_p50": {
    "mean": 648.5228070175439,
    "p50": 114.9,
    "p90": 1117.57,
    "nonzero_rate": 0.0015833333333333333
  },
  "attr_clk2cnv_min_mean": {
    "mean": 652.8216374269006,
    "p50": 114.9,
    "p90": 1117.57,
    "nonzero_rate": 0.0015833333333333333
  },
  "attr_exp2clk_min_p50": {
    "mean": 6602.626612903226,
    "p50": 2859.266666666667,
    "p90": 19152.266666666666,
    "nonzero_rate": 0.0025833333333333333
  },
  "1d_ctr": {
    "mean": 0.007669146825396826,
    "p50": 0.0,
    "p90": 0.0,
    "nonzero_rate": 0.015916666666666666
  },
  "7d_cvr": {
    "mean": 0.0556665663040663,
    "p50": 0.0,
    "p90": 0.0,
    "nonzero_rate": 0.05425
  },
  "life_ctcvr": {
    "mean": 0.030247676980717883,
    "p50": 0.01020408163265306,
    "p90": 0.07777777777777778,
    "nonzero_rate": 0.5690833333333334
  },
  "sess_n": {
    "mean": 77.12625,
    "p50": 80.0,
    "p90": 89.0,
    "nonzero_rate": 1.0
  },
  "sess_bounce_rate": {
    "mean": 0.7782744395703453,
    "p50": 0.7901234567901234,
    "p90": 0.8888888888888888,
    "nonzero_rate": 0.9995
  }
}
```

Catalog: `results/tencent_gr/feature_catalog.json` (138 dims).
