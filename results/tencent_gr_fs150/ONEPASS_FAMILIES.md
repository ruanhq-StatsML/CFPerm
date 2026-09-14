# Session / Markov / Cross — one-pass

150 里真正用到的原子很少。其余窗口漏斗可以后补。

## Session（gap=1800s）

切 session = 相邻 ts 差 >1800。**一趟**关上一场、开下一场：

```
if ts - prev_ts > 1800: close(); open()
cur_len++; cur_clk += (act==clk); cur_cnv += (act==cnv)
close last
```

关场时记：len、dur_min=(end-start)/60、clk条数、cnv条数、bounce=(len<=1)。

150 用到的：

| 名 | 计算 |
|---|---|
| sess_n | 场次数 |
| sess_bounce_rate | #len≤1 / sess_n |
| sess_depth_clk_mean | 每场 clk 条数均值 |
| sess_depth_cnv_mean | 每场 cnv 条数均值 |
| sess_clk_sess_rate | 含 clk 的场次比例 |
| sess_cnv_sess_rate | 含 cnv 的场次比例 |

其余 `len_p50/last_*` 生成了但没进 150。别算。

## Transitions（1 阶 3×3 + 8 个三元组）

状态机：记住 prev、prev2。

```
t1[prev, act] += 1
t2[prev2, prev, act] += 1
trans_a_to_b = t1[a,b] / (#pairs + 1e-6)
trans2_x_y_z = t2[x,y,z] / (#triples + 1e-6)
```

不是 P(b\|a)，是联合占比。150 用到：`exp→clk`, `exp→cnv`, `cnv→exp`, `cnv→cnv`, `exp→exp`, `clk→cnv`, `trans2_exp_cnv_exp`。

业务：`cnv→exp` = 买完又逛；`exp→clk` = 漏斗在走；`clk→cnv` = 点了就买。

## Cross（过完原子再乘，O(k²)，k=22）

固定键：

```
life_ctr, life_cvr, life_ctcvr, 7d_ctr, 7d_cvr, 1d_ctr, 1d_cvr,
pay_cnt, arpu_sum_proxy, arpu_mean_proxy, deal_cheap_cnv_share,
sess_bounce_rate, sess_n, active_days, hist_len, item_entropy_clk,
hours_since_last_clk, hours_since_last_cnv,
attr_clk2cnv_min_p50, attr_anyclk2cnv_min_p50,
dec_hl1d_dec_ctr, dec_hl1d_dec_cvr
```

```
sq_a = a²
log1p_abs_a = log1p(|a|)
x_a__b = a*b     # i<j，上三角
```

这是强度 × 犹豫 / 强度 × 广度。PO 板爱用它，因为单特征 F 看不到「转化率高但 bounce 也高」。

## 同一趟还能带上的（别再扫）

| 家族 | 状态 | 150 里要的 |
|---|---|---|
| funnel | `if 0≤age<W: n[W][act]++` | n/share/ctr/cvr/ctcvr、trend=短-长 |
| decay | `Σ 2^{-age/hl}` | hl1d/3d/7d 的 clk、cnv、arpu |
| attr | last_clk[iid] / last_any_clk | asof，见上条 |
| money | cnv 且 price≠None | sum/mean/p50 |
| diversity | iid 计数 | H_clk、n_uniq_*、hist_len、active_days |
| ui | 按 iid 聚合 | only_exp_share、exp_per_item_mean |

`ui` 等价 `groupby(item_id)`，不是第二趟时间扫描。
