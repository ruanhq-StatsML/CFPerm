# 150 特征：可落地计算

跑的是 `prefix_frac=0.75`。特征只看 prefix，label 是 suffix 里有没有 cnv（你这边可以不管 label）。

## 0. 事件

```
seq 一条: {item_id, action_type, timestamp}
action_type: 0=exp  1=clk  2=cnv
price = item_feat["115"][item_id]   # 缺/非有限/负 → None
按 timestamp 升序
e = (iid, act, ts, price)
```

## 1. 切 prefix（必须先切，再算）

```
if n < 5:
    prefix = all
else:
    t_cut = t0 + int((t1 - t0) * 0.75)      # t0=首ts, t1=末ts
    prefix = [e | ts <= t_cut]
    if prefix 空: prefix = evs[:max(1, min(n-1, int(n*0.75)))]
```

下面全部在 prefix 上。`t_end = prefix[-1].ts`。

## 2. 三个公用函数

```
rate(n, d) = n/d if d>0 else 0

entropy(counts):                         # EPS=1e-6
    tot = sum(counts.values())
    if tot<=0: return 0
    return -sum( p * log(p + 1e-6)  for c in counts.values() if c>0
                 for p in [c/tot] )

stats(xs, prefix):                       # 空列表
    n=0, min=p50=mean=p90=std=-1
    非空: n, min, median, mean, np.percentile(xs,90), np.std(xs)  # ddof=0
```

## 3. 漏斗窗口 funnel(W)

`sub = { e | t_end-W < ts <= t_end }`  
life 用 `W=10**12`（等于全史）。

```
n_exp, n_clk, n_cnv, n_tot
n_uniq_item = distinct iid
ctr  = rate(n_clk, n_exp)
cvr  = rate(n_cnv, n_clk)
ctcvr= rate(n_cnv, n_exp)
clk_share = rate(n_clk, n_tot)
cnv_share = rate(n_cnv, n_tot)
cnv_price_mean = mean(price | act=cnv, price≠None)   # 空→0
cnv_price_sum  = sum (...)                            # 空→0
log1p_n_clk = log1p(n_clk)
log1p_n_cnv = log1p(n_cnv)
```

窗口秒：`1d=86400, 3d=3*86400, 7d=7*86400, 14d=14*86400, 30d=30*86400`。

字段名：`life_*` / `7d_*` / … 就是上面这些键。

trend（150 里用到的）：

```
trend_n_cnv_7d_minus_30d = 7d_n_cnv - 30d_n_cnv     # 短减长，常≤0
trend_n_cnv_3d_minus_14d = 3d_n_cnv - 14d_n_cnv
trend_n_cnv_1d_minus_7d  = 1d_n_cnv - 7d_n_cnv
trend_n_clk_1d_minus_7d  = 1d_n_clk - 7d_n_clk
trend_cvr_3d_minus_14d   = 3d_cvr - 14d_cvr
trend_cvr_7d_minus_30d   = 7d_cvr - 30d_cvr
trend_n_cnv_3d_over_14d  = rate(3d_n_cnv, 14d_n_cnv if 14d_n_cnv!=0 else 1e-6)
trend_n_cnv_7d_over_30d  = rate(7d_n_cnv, 30d_n_cnv if 30d_n_cnv!=0 else 1e-6)
```

## 4. 指数衰减 decay(hl)

```
λ = log(2) / hl
w(e) = exp(-λ * max(t_end - ts, 0))
dec_clk = Σ w  over act=clk
dec_cnv = Σ w  over act=cnv
dec_arpu= Σ w*price  over act=cnv and price≠None
```

150 用到：`hl1d=86400, hl3d=3*86400, hl7d=7*86400`。  
字段：`dec_hl7d_dec_cnv` = 半衰期 7d 的 `dec_cnv`（中间那个 `dec_` 是 decay 函数返回键）。

## 5. session（SESS_GAP=1800s）

扫一遍：`ts - 上一条.ts > 1800` → 新 session。

```
sess_n
sess_bounce_rate     = (#session with len<=1) / sess_n
sess_depth_cnv_mean  = mean(每个 session 的 cnv 条数)
sess_depth_clk_mean  = mean(每个 session 的 clk 条数)
sess_cnv_sess_rate   = (#session 含≥1 cnv) / sess_n
sess_clk_sess_rate   = (#session 含≥1 clk) / sess_n
```

## 6. 多样性

```
hist_len     = len(prefix)
active_days  = #{ ts // 86400 }
n_uniq_clk_item = #{iid | 有 clk}
n_uniq_cnv_item = #{iid | 有 cnv}
item_entropy_clk = entropy(clk 的 iid 频次)
item_entropy_cnv = entropy(cnv 的 iid 频次)
```

## 7. 归因（扫一遍，按时间）

状态：`last_exp[iid], first_exp, last_clk[iid], first_clk, last_any_clk, exp_cnt[iid]`

```
exp: last_exp[iid]=ts; first_exp.setdefault; exp_cnt[iid]++
clk: last_clk[iid]=ts; first_clk.setdefault; last_any_clk=ts
     （同 item last_exp 存在时记 exp2clk，150 没用）
cnv:
    if last_clk[iid] 存在:  clk2cnv_min.append( (ts-last_clk[iid])/60 )
    else:                   cnv_wo_clk += 1
    if last_any_clk 存在:   anyclk2cnv_min.append( (ts-last_any_clk)/60 )
```

```
attr_cnv_wo_prior_clk_cnt  = cnv_wo_clk
attr_cnv_wo_prior_clk_rate = rate(cnv_wo_clk, n_cnv)
attr_anyclk2cnv_min_n      = len(anyclk2cnv_min)
attr_anyclk2cnv_min_p50    = median(anyclk2cnv_min)   # 空→-1
attr_anyclk2cnv_min_p90    = p90(...)                 # 空→-1
attr_anyclk2cnv_min_std    = std(...)                 # 空→-1
attr_clk2cnv_min_p50       = median(clk2cnv_min)      # 同 item；空→-1
```

`anyclk` = 全局上次任意点击；`clk2cnv` = **同 item** 上次点击。单位分钟。

## 8. 变现 / UI / 转移

```
pay_cnt            = n_cnv(life)
pay_user           = 1 if pay_cnt>0 else 0
prices             = [price | act=cnv, price≠None]
arpu_sum_proxy     = sum(prices)          # 空→0   ≡ life_cnv_price_sum
arpu_mean_proxy    = mean(prices)         # 空→0   ≡ life_cnv_price_mean
arpu_p50_proxy     = median(prices)       # 空→0
log1p_arpu_sum     = log1p(arpu_sum_proxy)
arpu_per_active_day_proxy = arpu_sum / active_days if active_days>0 else 0
cvr_pay_per_click  = rate(pay_cnt, n_clk) ≡ life_cvr
ctvr_pay_per_expose= rate(pay_cnt, n_exp) ≡ life_ctcvr
30d_cnv_price_mean = funnel(30d).cnv_price_mean
```

按 iid 分组（组内已按时间）：

```
ui_exp_per_item_mean = mean(每个 item 的 n_exp)
ui_only_exp_share    = #{item: n_exp>0, n_clk=0, n_cnv=0} / #item
```

相邻转移（分母是对数 + 1e-6，**不是**条件概率）：

```
tot1 = (#相邻对) + 1e-6
trans_a_to_b = count(act_i=a, act_{i+1}=b) / tot1

tot2 = (#三元组) + 1e-6
trans2_exp_cnv_exp = count(exp,cnv,exp) / tot2
```

## 9. 交叉（生成后再选进 150 的）

在已经算完的原子上：

```
sq_a        = a * a
log1p_abs_a = log1p(|a|)
x_a__b      = a * b
```

## 10. 150 个 = 上面原子的名字

| # | 名称 | 就是 |
|---|---|---|
| 1 | life_log1p_n_cnv | log1p(n_cnv life) |
| 2 | log1p_pay_cnt | =1 |
| 3 | log1p_abs_pay_cnt | =1（pay_cnt≥0） |
| 4 | trans_cnv_to_exp | §8 |
| 5 | trans_exp_to_cnv | §8 |
| 6 | 30d_log1p_n_cnv | log1p(n_cnv 30d) |
| 7 | item_entropy_cnv | §6 |
| 8 | sess_cnv_sess_rate | §5 |
| 9 | dec_hl7d_dec_cnv | §4 hl=7d |
| 10 | sess_depth_cnv_mean | §5 |
| 11 | life_cnv_share | n_cnv/n_tot life |
| 12 | 30d_n_cnv | n_cnv 30d |
| 13 | n_uniq_cnv_item | §6 |
| 14 | life_n_cnv | n_cnv life |
| 15 | pay_cnt | =14 |
| 16 | attr_cnv_wo_prior_clk_cnt | §7 |
| 17 | ui_only_exp_share | §8 |
| 18 | ui_exp_per_item_mean | §8 |
| 19 | log1p_abs_life_ctcvr | log1p(\|life_ctcvr\|) |
| 20 | x_pay_cnt__sess_bounce_rate | pay_cnt × bounce |
| 21 | 30d_cnv_share | n_cnv/n_tot 30d |
| 22 | 14d_log1p_n_cnv | log1p(n_cnv 14d) |
| 23 | trans_exp_to_exp | §8 |
| 24 | 14d_n_cnv | n_cnv 14d |
| 25 | trans2_exp_cnv_exp | §8 |
| 26 | trend_n_cnv_7d_minus_30d | §3 |
| 27 | dec_hl3d_dec_cnv | §4 hl=3d |
| 28 | x_pay_cnt__hist_len | pay_cnt × hist_len |
| 29 | x_pay_cnt__sess_n | pay_cnt × sess_n |
| 30 | x_life_ctcvr__hist_len | life_ctcvr × hist_len |
| 31 | 7d_log1p_n_cnv | log1p(n_cnv 7d) |
| 32 | x_life_ctcvr__sess_n | life_ctcvr × sess_n |
| 33 | trend_n_cnv_3d_minus_14d | §3 |
| 34 | life_ctcvr | n_cnv/n_exp life |
| 35 | ctvr_pay_per_expose | =34 |
| 36 | 7d_n_cnv | n_cnv 7d |
| 37 | x_life_ctr__life_cvr | life_ctr × life_cvr |
| 38 | x_pay_cnt__active_days | pay_cnt × active_days |
| 39 | 14d_cnv_share | n_cnv/n_tot 14d |
| 40 | attr_anyclk2cnv_min_n | §7 |
| 41 | x_life_cvr__item_entropy_clk | life_cvr × H_clk |
| 42 | attr_cnv_wo_prior_clk_rate | §7 |
| 43 | pay_user | pay_cnt>0 |
| 44 | log1p_abs_life_cvr | log1p(\|life_cvr\|) |
| 45 | x_life_ctcvr__active_days | life_ctcvr × active_days |
| 46 | x_life_ctcvr__sess_bounce_rate | life_ctcvr × bounce |
| 47 | 30d_ctcvr | n_cnv/n_exp 30d |
| 48 | x_pay_cnt__item_entropy_clk | pay_cnt × H_clk |
| 49 | trend_n_cnv_1d_minus_7d | §3 |
| 50 | log1p_abs_attr_anyclk2cnv_min_p50 | log1p(\|p50\|)；空 p50=-1 |
| 51 | 3d_log1p_n_cnv | log1p(n_cnv 3d) |
| 52 | x_life_ctcvr__item_entropy_clk | life_ctcvr × H_clk |
| 53 | 3d_n_cnv | n_cnv 3d |
| 54 | 14d_ctcvr | n_cnv/n_exp 14d |
| 55 | dec_hl1d_dec_cnv | §4 hl=1d |
| 56 | 7d_cnv_share | n_cnv/n_tot 7d |
| 57 | trans_cnv_to_cnv | §8 |
| 58 | life_cvr | n_cnv/n_clk life |
| 59 | cvr_pay_per_click | =58 |
| 60 | x_life_ctr__pay_cnt | life_ctr × pay_cnt |
| 61 | 30d_cvr | n_cnv/n_clk 30d |
| 62 | x_life_cvr__sess_bounce_rate | life_cvr × bounce |
| 63 | 7d_log1p_n_clk | log1p(n_clk 7d) |
| 64 | 14d_cvr | n_cnv/n_clk 14d |
| 65 | x_life_ctcvr__arpu_mean_proxy | life_ctcvr × arpu_mean |
| 66 | 14d_log1p_n_clk | log1p(n_clk 14d) |
| 67 | x_pay_cnt__arpu_mean_proxy | pay_cnt × arpu_mean |
| 68 | sq_pay_cnt | pay_cnt² |
| 69 | x_life_cvr__hist_len | life_cvr × hist_len |
| 70 | life_log1p_n_clk | log1p(n_clk life) |
| 71 | 7d_ctcvr | n_cnv/n_exp 7d |
| 72 | x_life_cvr__sess_n | life_cvr × sess_n |
| 73 | 30d_log1p_n_clk | log1p(n_clk 30d) |
| 74 | x_life_ctr__7d_cvr | life_ctr × 7d_cvr |
| 75 | log1p_abs_7d_cvr | log1p(\|7d_cvr\|) |
| 76 | x_7d_cvr__item_entropy_clk | 7d_cvr × H_clk |
| 77 | dec_hl3d_dec_clk | §4 |
| 78 | dec_hl7d_dec_clk | §4 |
| 79 | x_7d_ctr__pay_cnt | 7d_ctr × pay_cnt |
| 80 | 7d_n_clk | n_clk 7d |
| 81 | log1p_abs_item_entropy_clk | log1p(\|H_clk\|) |
| 82 | 3d_cnv_share | n_cnv/n_tot 3d |
| 83 | x_life_cvr__active_days | life_cvr × active_days |
| 84 | 7d_clk_share | n_clk/n_tot 7d |
| 85 | 3d_log1p_n_clk | log1p(n_clk 3d) |
| 86 | item_entropy_clk | §6 |
| 87 | x_life_ctcvr__attr_anyclk2cnv_min_p50 | life_ctcvr × anyclk p50 |
| 88 | attr_anyclk2cnv_min_p90 | §7 |
| 89 | log1p_abs_7d_ctr | log1p(\|7d_ctr\|) |
| 90 | x_life_ctr__attr_anyclk2cnv_min_p50 | life_ctr × anyclk p50 |
| 91 | life_cnv_price_mean | mean(cnv price) life |
| 92 | arpu_mean_proxy | =91 |
| 93 | log1p_arpu_sum | log1p(Σ cnv price) |
| 94 | log1p_abs_arpu_sum_proxy | ≈93 |
| 95 | x_sess_bounce_rate__item_entropy_clk | bounce × H_clk |
| 96 | arpu_p50_proxy | median(cnv price) |
| 97 | x_arpu_mean_proxy__sess_bounce_rate | arpu_mean × bounce |
| 98 | x_life_cvr__7d_ctr | life_cvr × 7d_ctr |
| 99 | trend_n_cnv_3d_over_14d | §3 |
| 100 | 14d_clk_share | n_clk/n_tot 14d |
| 101 | 14d_n_clk | n_clk 14d |
| 102 | 30d_cnv_price_mean | mean(cnv price) 30d |
| 103 | log1p_abs_arpu_mean_proxy | log1p(\|arpu_mean\|) |
| 104 | log1p_abs_life_ctr | log1p(\|life_ctr\|) |
| 105 | x_life_ctr__hist_len | life_ctr × hist_len |
| 106 | sq_arpu_mean_proxy | arpu_mean² |
| 107 | attr_anyclk2cnv_min_std | §7 |
| 108 | trend_n_cnv_7d_over_30d | §3 |
| 109 | x_life_ctr__sess_n | life_ctr × sess_n |
| 110 | x_pay_cnt__attr_anyclk2cnv_min_p50 | pay_cnt × anyclk p50 |
| 111 | 3d_ctcvr | n_cnv/n_exp 3d |
| 112 | x_life_ctr__active_days | life_ctr × active_days |
| 113 | x_7d_ctr__7d_cvr | 7d_ctr × 7d_cvr |
| 114 | arpu_per_active_day_proxy | Σprice / active_days |
| 115 | sess_clk_sess_rate | §5 |
| 116 | 7d_cvr | n_cnv/n_clk 7d |
| 117 | dec_hl1d_dec_clk | §4 hl=1d |
| 118 | trend_n_clk_1d_minus_7d | §3 |
| 119 | life_clk_share | n_clk/n_tot life |
| 120 | x_7d_cvr__sess_bounce_rate | 7d_cvr × bounce |
| 121 | 30d_n_clk | n_clk 30d |
| 122 | x_7d_cvr__attr_clk2cnv_min_p50 | 7d_cvr × **同item** clk2cnv p50 |
| 123 | n_uniq_clk_item | §6 |
| 124 | 1d_log1p_n_cnv | log1p(n_cnv 1d) |
| 125 | life_n_clk | n_clk life |
| 126 | x_life_cvr__arpu_mean_proxy | life_cvr × arpu_mean |
| 127 | x_arpu_mean_proxy__hist_len | arpu_mean × hist_len |
| 128 | life_cnv_price_sum | Σ cnv price life |
| 129 | arpu_sum_proxy | =128 |
| 130 | x_life_ctcvr__7d_ctr | life_ctcvr × 7d_ctr |
| 131 | x_arpu_sum_proxy__sess_bounce_rate | arpu_sum × bounce |
| 132 | sess_depth_clk_mean | §5 |
| 133 | 3d_n_clk | n_clk 3d |
| 134 | 30d_clk_share | n_clk/n_tot 30d |
| 135 | trend_cvr_3d_minus_14d | §3 |
| 136 | x_arpu_mean_proxy__sess_n | arpu_mean × sess_n |
| 137 | x_hist_len__item_entropy_clk | hist_len × H_clk |
| 138 | trans_exp_to_clk | §8 |
| 139 | trend_cvr_7d_minus_30d | §3 |
| 140 | 1d_n_cnv | n_cnv 1d |
| 141 | dec_hl7d_dec_arpu | §4 |
| 142 | x_7d_cvr__hist_len | 7d_cvr × hist_len |
| 143 | sq_item_entropy_clk | H_clk² |
| 144 | x_7d_cvr__sess_n | 7d_cvr × sess_n |
| 145 | x_life_ctcvr__pay_cnt | life_ctcvr × pay_cnt |
| 146 | x_sess_n__item_entropy_clk | sess_n × H_clk |
| 147 | trans_clk_to_cnv | §8 |
| 148 | x_arpu_mean_proxy__active_days | arpu_mean × active_days |
| 149 | 3d_clk_share | n_clk/n_tot 3d |
| 150 | x_arpu_sum_proxy__arpu_mean_proxy | arpu_sum × arpu_mean |

## 对拍时容易错的点

1. 窗口左开右闭：`(t_end-W, t_end]`。落在 `t_end-W` 上的点不算。
2. 转移分母 `+#pairs + 1e-6`，不是 `count(from=a)`。
3. 归因 p50 空值是 **-1**，不是 0；所以 `log1p_abs_*` 空时 = `log(2)`。
4. `np.std` 默认 ddof=0。
5. price=None 的 cnv **不进** arpu 的 sum/mean/median，但 **进** pay_cnt。
6. 150 里没有 user_feat 侧列（`u_side_*`）。
7. 对拍脚本：`scripts/tencent_gr/auto_feats_select_train.py` 里对应函数。
