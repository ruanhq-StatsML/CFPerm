# 转化粒模型 v1

Y = 满窗 `y_post_clk_1d`。X 停在 `cnv_ts`。按用户切。SKU 支路看每天构成率。

```
post 转化粒                    daily mix: sku_clk=0.0011  empty=0.484  sess0=0.984
        |                      SKU_GATE=OFF (构成率<1%，只监控)
   满窗滤 y_post_clk_1d
   按 user_id 70/30
        |
  ┌─────────┬──────────┬───────────┬──────────┬────────────┬─────────────────┐
  │ heat    │ empty    │ any_path  │ sess     │ lag        │ sku × GATE      │
  │ 7d/1d/1h│ empty_any│ dt_any    │ pos      │ post_clk_1d│ wo_prior_clk    │
  │ n_clk   │ lag_empty│ (+miss)   │ clk_before│ n_prior    │ dt_item/within  │
  └─────────┴──────────┴───────────┴──────────┴────────────┴─────────────────┘
        | concat + User×Ctx 交叉（SKU 交叉 × GATE）
        +-- LogReg  (scale → 线性)           系数可读
        +-- HGB     (depth=3, 80 iter)       吃 NaN
        +-- MLP     (16 → 8 → 1, ReLU)       浅层非线性
        |
  Ŷ = P(满窗后 1d 内任意点击 | cnv_ts 已知)
  不是 P(哪次点击导致购买)
```
每天：dump 中间表 → 打构成 → SKU 占比跨过 1% 才打开 sku 支路 → 按同一协议重训。

## daily mix（这批 prefix）

- n_cnv=12866  sku_clk_share=**0.0011**  empty_any=0.484  sess0=0.984
- SKU_GATE=关（阈值 1%）

## 满窗 1d 预报（user 70/30）

| 消融 | LogReg AUC | HGB AUC | MLP AUC | HGB AP |
|---|---:|---:|---:|---:|
| heat | 0.677 | 0.678 | 0.677 | 0.274 |
| +empty | 0.683 | 0.687 | 0.682 | 0.283 |
| +any_path | 0.688 | 0.686 | 0.653 | 0.276 |
| +sess | 0.687 | 0.685 | 0.686 | 0.279 |
| +lag | 0.728 | 0.746 | 0.730 | 0.374 |
| +cross | 0.740 | 0.748 | 0.729 | 0.386 |

LogReg 系数（+cross，标准化后）：

| feat | coef |
|---|---:|
| `uc_lag_post_clk_1d_rate__x__empty_any` | +1.977 |
| `n_clk_before_7d` | +0.411 |
| `dt_any_min` | -0.311 |
| `uc_lag_post_clk_1d_rate__x__dt_any_log` | +0.236 |
| `sess_clk_before` | +0.178 |
| `empty_any` | -0.158 |
| `dt_any_min_miss` | -0.158 |
| `uc_n_clk_7d_log__x__sess_clk_before` | -0.142 |
| `uc_n_clk_7d_log__x__dt_any_log` | -0.141 |
| `lag_empty_any` | -0.066 |
| `sess_pos` | -0.049 |
| `uc_lag_empty_any__x__empty_any` | -0.044 |
| `n_clk_before_1d` | +0.036 |
| `n_prior_cnv` | +0.034 |
| `n_clk_before_1h` | -0.029 |
| `lag_post_clk_1d_rate` | +0.024 |

AP/AUC 只说明满窗后还会不会点预报得怎样。SKU 漏斗闭没闭合看 daily mix，不看这张表。
