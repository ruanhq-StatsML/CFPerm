# TencentGR 150 特征：名称 → 计算

全在用户序列 **prefix** 上算。事件 `(item, act, ts, price)`，`act∈{exp,clk,cnv}`。
`t_end`=最后一条 ts。窗口 `Wd` = `(t_end-d, t_end]`。
`rate(n,d)=n/d`（d=0→0）。session 切分：相邻 ts 差 >30min。
price = item_feat `115`。全局 p50 只用于没进这 150 的 deal bin。

**记号**
- `n_a(W)`：窗口内 act=a 的条数；`life` = 全史（W=∞）
- `share_a(W)=n_a(W)/n_tot(W)`
- `ctr=n_clk/n_exp`，`cvr=n_cnv/n_clk`，`ctcvr=n_cnv/n_exp`
- `dec_hlT(a)=Σ_e 2^{-age/T}`，age=`t_end-ts`，T∈{1d,3d,7d}
- `H_a`：对该 act 的 item 频次做 `−Σ p log(p+ε)`
- `x_a__b = a·b`，`sq_a=a²`，`log1p_abs_a=log(1+|a|)`

同值别名（F 分相同）：`pay_cnt≡life_n_cnv`，`log1p_pay_cnt≡life_log1p_n_cnv≡log1p_abs_pay_cnt`，`ctvr_pay_per_expose≡life_ctcvr`，`cvr_pay_per_click≡life_cvr`，`arpu_sum_proxy≡life_cnv_price_sum`，`arpu_mean_proxy≡life_cnv_price_mean`。

| # | 名称 | 计算 |
|---|---|---|
| 1 | life_log1p_n_cnv | log1p(n_cnv(life)) |
| 2 | log1p_pay_cnt | log1p(pay_cnt)，=1 |
| 3 | log1p_abs_pay_cnt | log1p(\|pay_cnt\|)，=1 |
| 4 | trans_cnv_to_exp | 相邻对里 (cnv→exp) 占比 |
| 5 | trans_exp_to_cnv | 相邻对 (exp→cnv) 占比 |
| 6 | 30d_log1p_n_cnv | log1p(n_cnv(30d)) |
| 7 | item_entropy_cnv | H_cnv |
| 8 | sess_cnv_sess_rate | 含≥1次 cnv 的 session 比例 |
| 9 | dec_hl7d_dec_cnv | dec_hl7d(cnv) |
| 10 | sess_depth_cnv_mean | 每 session cnv 次数的均值 |
| 11 | life_cnv_share | share_cnv(life) |
| 12 | 30d_n_cnv | n_cnv(30d) |
| 13 | n_uniq_cnv_item | 转化过的 distinct item 数 |
| 14 | life_n_cnv | n_cnv(life) |
| 15 | pay_cnt | 转化条数，=14 |
| 16 | attr_cnv_wo_prior_clk_cnt | 该 item 转化前没有 clk 的 cnv 条数 |
| 17 | ui_only_exp_share | 只曝光、无 clk/cnv 的 item 占比 |
| 18 | ui_exp_per_item_mean | 每 item 曝光次数均值 |
| 19 | log1p_abs_life_ctcvr | log1p(\|life_ctcvr\|) |
| 20 | x_pay_cnt__sess_bounce_rate | pay_cnt × sess_bounce_rate |
| 21 | 30d_cnv_share | share_cnv(30d) |
| 22 | 14d_log1p_n_cnv | log1p(n_cnv(14d)) |
| 23 | trans_exp_to_exp | 相邻对 (exp→exp) 占比 |
| 24 | 14d_n_cnv | n_cnv(14d) |
| 25 | trans2_exp_cnv_exp | 相邻三元组 (exp,cnv,exp) 占比 |
| 26 | trend_n_cnv_7d_minus_30d | n_cnv(7d) − n_cnv(30d) |
| 27 | dec_hl3d_dec_cnv | dec_hl3d(cnv) |
| 28 | x_pay_cnt__hist_len | pay_cnt × hist_len |
| 29 | x_pay_cnt__sess_n | pay_cnt × sess_n |
| 30 | x_life_ctcvr__hist_len | life_ctcvr × hist_len |
| 31 | 7d_log1p_n_cnv | log1p(n_cnv(7d)) |
| 32 | x_life_ctcvr__sess_n | life_ctcvr × sess_n |
| 33 | trend_n_cnv_3d_minus_14d | n_cnv(3d) − n_cnv(14d) |
| 34 | life_ctcvr | n_cnv/n_exp（life） |
| 35 | ctvr_pay_per_expose | pay_cnt / n_exp，=34 |
| 36 | 7d_n_cnv | n_cnv(7d) |
| 37 | x_life_ctr__life_cvr | life_ctr × life_cvr |
| 38 | x_pay_cnt__active_days | pay_cnt × active_days |
| 39 | 14d_cnv_share | share_cnv(14d) |
| 40 | attr_anyclk2cnv_min_n | 有「任意 item 上次 clk」可归因的 cnv 条数 |
| 41 | x_life_cvr__item_entropy_clk | life_cvr × H_clk |
| 42 | attr_cnv_wo_prior_clk_rate | #16 / n_cnv |
| 43 | pay_user | 1 if pay_cnt>0 else 0 |
| 44 | log1p_abs_life_cvr | log1p(\|life_cvr\|) |
| 45 | x_life_ctcvr__active_days | life_ctcvr × active_days |
| 46 | x_life_ctcvr__sess_bounce_rate | life_ctcvr × sess_bounce_rate |
| 47 | 30d_ctcvr | ctcvr(30d) |
| 48 | x_pay_cnt__item_entropy_clk | pay_cnt × H_clk |
| 49 | trend_n_cnv_1d_minus_7d | n_cnv(1d) − n_cnv(7d) |
| 50 | log1p_abs_attr_anyclk2cnv_min_p50 | log1p(\|median(Δt_min: 上次任意clk→cnv)\|) |
| 51 | 3d_log1p_n_cnv | log1p(n_cnv(3d)) |
| 52 | x_life_ctcvr__item_entropy_clk | life_ctcvr × H_clk |
| 53 | 3d_n_cnv | n_cnv(3d) |
| 54 | 14d_ctcvr | ctcvr(14d) |
| 55 | dec_hl1d_dec_cnv | dec_hl1d(cnv) |
| 56 | 7d_cnv_share | share_cnv(7d) |
| 57 | trans_cnv_to_cnv | 相邻对 (cnv→cnv) 占比 |
| 58 | life_cvr | n_cnv/n_clk（life） |
| 59 | cvr_pay_per_click | pay_cnt / n_clk，=58 |
| 60 | x_life_ctr__pay_cnt | life_ctr × pay_cnt |
| 61 | 30d_cvr | cvr(30d) |
| 62 | x_life_cvr__sess_bounce_rate | life_cvr × sess_bounce_rate |
| 63 | 7d_log1p_n_clk | log1p(n_clk(7d)) |
| 64 | 14d_cvr | cvr(14d) |
| 65 | x_life_ctcvr__arpu_mean_proxy | life_ctcvr × arpu_mean_proxy |
| 66 | 14d_log1p_n_clk | log1p(n_clk(14d)) |
| 67 | x_pay_cnt__arpu_mean_proxy | pay_cnt × arpu_mean_proxy |
| 68 | sq_pay_cnt | pay_cnt² |
| 69 | x_life_cvr__hist_len | life_cvr × hist_len |
| 70 | life_log1p_n_clk | log1p(n_clk(life)) |
| 71 | 7d_ctcvr | ctcvr(7d) |
| 72 | x_life_cvr__sess_n | life_cvr × sess_n |
| 73 | 30d_log1p_n_clk | log1p(n_clk(30d)) |
| 74 | x_life_ctr__7d_cvr | life_ctr × cvr(7d) |
| 75 | log1p_abs_7d_cvr | log1p(\|cvr(7d)\|) |
| 76 | x_7d_cvr__item_entropy_clk | cvr(7d) × H_clk |
| 77 | dec_hl3d_dec_clk | dec_hl3d(clk) |
| 78 | dec_hl7d_dec_clk | dec_hl7d(clk) |
| 79 | x_7d_ctr__pay_cnt | ctr(7d) × pay_cnt |
| 80 | 7d_n_clk | n_clk(7d) |
| 81 | log1p_abs_item_entropy_clk | log1p(\|H_clk\|) |
| 82 | 3d_cnv_share | share_cnv(3d) |
| 83 | x_life_cvr__active_days | life_cvr × active_days |
| 84 | 7d_clk_share | share_clk(7d) |
| 85 | 3d_log1p_n_clk | log1p(n_clk(3d)) |
| 86 | item_entropy_clk | H_clk |
| 87 | x_life_ctcvr__attr_anyclk2cnv_min_p50 | life_ctcvr × median(任意clk→cnv 分钟) |
| 88 | attr_anyclk2cnv_min_p90 | p90(任意clk→cnv 分钟) |
| 89 | log1p_abs_7d_ctr | log1p(\|ctr(7d)\|) |
| 90 | x_life_ctr__attr_anyclk2cnv_min_p50 | life_ctr × median(任意clk→cnv 分钟) |
| 91 | life_cnv_price_mean | 全史 cnv 的 price 均值 |
| 92 | arpu_mean_proxy | 同上，=91 |
| 93 | log1p_arpu_sum | log1p(Σ cnv price) |
| 94 | log1p_abs_arpu_sum_proxy | log1p(\|Σ cnv price\|)，≈93 |
| 95 | x_sess_bounce_rate__item_entropy_clk | sess_bounce_rate × H_clk |
| 96 | arpu_p50_proxy | 全史 cnv price 中位数 |
| 97 | x_arpu_mean_proxy__sess_bounce_rate | arpu_mean × sess_bounce_rate |
| 98 | x_life_cvr__7d_ctr | life_cvr × ctr(7d) |
| 99 | trend_n_cnv_3d_over_14d | n_cnv(3d) / n_cnv(14d) |
| 100 | 14d_clk_share | share_clk(14d) |
| 101 | 14d_n_clk | n_clk(14d) |
| 102 | 30d_cnv_price_mean | 30d 内 cnv 的 price 均值 |
| 103 | log1p_abs_arpu_mean_proxy | log1p(\|arpu_mean\|) |
| 104 | log1p_abs_life_ctr | log1p(\|life_ctr\|) |
| 105 | x_life_ctr__hist_len | life_ctr × hist_len |
| 106 | sq_arpu_mean_proxy | arpu_mean² |
| 107 | attr_anyclk2cnv_min_std | std(任意clk→cnv 分钟) |
| 108 | trend_n_cnv_7d_over_30d | n_cnv(7d) / n_cnv(30d) |
| 109 | x_life_ctr__sess_n | life_ctr × sess_n |
| 110 | x_pay_cnt__attr_anyclk2cnv_min_p50 | pay_cnt × median(任意clk→cnv 分钟) |
| 111 | 3d_ctcvr | ctcvr(3d) |
| 112 | x_life_ctr__active_days | life_ctr × active_days |
| 113 | x_7d_ctr__7d_cvr | ctr(7d) × cvr(7d) |
| 114 | arpu_per_active_day_proxy | (Σ cnv price) / active_days |
| 115 | sess_clk_sess_rate | 含≥1次 clk 的 session 比例 |
| 116 | 7d_cvr | cvr(7d) |
| 117 | dec_hl1d_dec_clk | dec_hl1d(clk) |
| 118 | trend_n_clk_1d_minus_7d | n_clk(1d) − n_clk(7d) |
| 119 | life_clk_share | share_clk(life) |
| 120 | x_7d_cvr__sess_bounce_rate | cvr(7d) × sess_bounce_rate |
| 121 | 30d_n_clk | n_clk(30d) |
| 122 | x_7d_cvr__attr_clk2cnv_min_p50 | cvr(7d) × median(同 item 上次clk→cnv 分钟) |
| 123 | n_uniq_clk_item | 点击过的 distinct item 数 |
| 124 | 1d_log1p_n_cnv | log1p(n_cnv(1d)) |
| 125 | life_n_clk | n_clk(life) |
| 126 | x_life_cvr__arpu_mean_proxy | life_cvr × arpu_mean |
| 127 | x_arpu_mean_proxy__hist_len | arpu_mean × hist_len |
| 128 | life_cnv_price_sum | Σ cnv price（life） |
| 129 | arpu_sum_proxy | 同上，=128 |
| 130 | x_life_ctcvr__7d_ctr | life_ctcvr × ctr(7d) |
| 131 | x_arpu_sum_proxy__sess_bounce_rate | arpu_sum × sess_bounce_rate |
| 132 | sess_depth_clk_mean | 每 session clk 次数均值 |
| 133 | 3d_n_clk | n_clk(3d) |
| 134 | 30d_clk_share | share_clk(30d) |
| 135 | trend_cvr_3d_minus_14d | cvr(3d) − cvr(14d) |
| 136 | x_arpu_mean_proxy__sess_n | arpu_mean × sess_n |
| 137 | x_hist_len__item_entropy_clk | hist_len × H_clk |
| 138 | trans_exp_to_clk | 相邻对 (exp→clk) 占比 |
| 139 | trend_cvr_7d_minus_30d | cvr(7d) − cvr(30d) |
| 140 | 1d_n_cnv | n_cnv(1d) |
| 141 | dec_hl7d_dec_arpu | Σ_{cnv} 2^{-age/7d} · price |
| 142 | x_7d_cvr__hist_len | cvr(7d) × hist_len |
| 143 | sq_item_entropy_clk | H_clk² |
| 144 | x_7d_cvr__sess_n | cvr(7d) × sess_n |
| 145 | x_life_ctcvr__pay_cnt | life_ctcvr × pay_cnt |
| 146 | x_sess_n__item_entropy_clk | sess_n × H_clk |
| 147 | trans_clk_to_cnv | 相邻对 (clk→cnv) 占比 |
| 148 | x_arpu_mean_proxy__active_days | arpu_mean × active_days |
| 149 | 3d_clk_share | share_clk(3d) |
| 150 | x_arpu_sum_proxy__arpu_mean_proxy | arpu_sum × arpu_mean |

**依赖的未单列原子**
- `hist_len`：prefix 事件条数
- `active_days`：distinct `ts//86400`
- `sess_n`：session 数；`sess_bounce_rate`：长度≤1 的 session 占比
- `life_ctr`：n_clk/n_exp（life）
- `attr_clk2cnv_min_p50`：同 item 上次 clk→cnv 的分钟中位数（#122 用）
- `attr_anyclk2cnv_min_*`：全局上次任意 clk→该 cnv 的分钟
