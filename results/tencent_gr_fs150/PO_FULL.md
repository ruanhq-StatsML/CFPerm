# Direct PO on full one-pass X

Y=`future_cnv`. W=`1{t_end > median}`. 不裁列、不换 Y、不做 LOGO。

## `onepass_full`  n=6000 p=107  pos=0.157  RF-domain **0.788**  PO-risk **0.001699**

family mass | RF-domain (P(X)) | PO-VIMP |
|---|---:|---:|
| funnel | 0.423 | 0.435 |
| cross | 0.267 | 0.228 |
| decay | 0.058 | 0.126 |
| session | 0.098 | 0.085 |
| markov | 0.122 | 0.075 |
| post | 0.010 | 0.035 |
| attr | 0.013 | 0.015 |
| diversity | 0.007 | 0.002 |
| money | 0.001 | 0.000 |

| rank | feat | family | PO-VIMP | RF-VIMP |
|---|---|---|---:|---:|
| 1 | `14d_n_exp` | funnel | 0.0525 | 0.0481 |
| 2 | `x_life_ctcvr__sess_bounce_rate` | cross | 0.0448 | 0.0299 |
| 3 | `30d_n_exp` | funnel | 0.0418 | 0.0366 |
| 4 | `life_cnv_share` | funnel | 0.0401 | 0.0055 |
| 5 | `7d_n_exp` | funnel | 0.0360 | 0.0368 |
| 6 | `hist_len` | funnel | 0.0350 | 0.0147 |
| 7 | `active_days` | funnel | 0.0271 | 0.0248 |
| 8 | `dec_hl7d_dec_clk` | decay | 0.0269 | 0.0107 |
| 9 | `dec_hl7d_dec_cnv` | decay | 0.0262 | 0.0146 |
| 10 | `dec_hl1d_dec_clk` | decay | 0.0261 | 0.0064 |
| 11 | `sess_bounce_rate` | session | 0.0222 | 0.0056 |
| 12 | `sess_n` | session | 0.0221 | 0.0103 |
| 13 | `x_life_cvr__sess_bounce_rate` | cross | 0.0221 | 0.0026 |
| 14 | `x_life_ctcvr__active_days` | cross | 0.0218 | 0.0521 |
| 15 | `x_life_ctr__attr_anyclk2cnv_min_p50` | cross | 0.0214 | 0.0083 |
| 16 | `dec_hl3d_dec_cnv` | decay | 0.0210 | 0.0129 |
| 17 | `trans_cnv_to_exp` | markov | 0.0206 | 0.0114 |
| 18 | `dec_hl3d_dec_clk` | decay | 0.0196 | 0.0088 |
| 19 | `sess_cnv_sess_rate` | session | 0.0182 | 0.0250 |
| 20 | `3d_n_exp` | funnel | 0.0174 | 0.0251 |
| 21 | `trans_exp_to_cnv` | markov | 0.0162 | 0.0077 |
| 22 | `life_n_exp` | funnel | 0.0161 | 0.0168 |
| 23 | `trend_cvr_3d_minus_14d` | funnel | 0.0158 | 0.0011 |
| 24 | `x_life_ctcvr__sess_n` | cross | 0.0155 | 0.0322 |
| 25 | `trans2_exp_cnv_exp` | markov | 0.0143 | 0.0180 |
| 26 | `x_life_cvr__7d_ctr` | cross | 0.0135 | 0.0085 |
| 27 | `log1p_abs_life_ctcvr` | cross | 0.0134 | 0.0301 |
| 28 | `post_clk_dt_p50` | post | 0.0127 | 0.0020 |
| 29 | `x_hist_len__item_entropy_clk` | cross | 0.0123 | 0.0061 |
| 30 | `sess_depth_cnv_mean` | session | 0.0123 | 0.0269 |
| 31 | `x_life_ctr__active_days` | cross | 0.0120 | 0.0178 |
| 32 | `x_post_clk_1d_rate__life_ctcvr` | cross | 0.0119 | 0.0018 |
| 33 | `14d_cnv_share` | funnel | 0.0116 | 0.0093 |
| 34 | `trans_exp_to_exp` | markov | 0.0113 | 0.0682 |
| 35 | `30d_ctr` | funnel | 0.0108 | 0.0061 |
| 36 | `trans_exp_to_clk` | markov | 0.0105 | 0.0165 |
| 37 | `attr_anyclk2cnv_min_p50` | attr | 0.0104 | 0.0038 |
| 38 | `x_life_ctr__sess_n` | cross | 0.0104 | 0.0038 |
| 39 | `x_sess_n__item_entropy_clk` | cross | 0.0103 | 0.0038 |
| 40 | `7d_ctcvr` | funnel | 0.0102 | 0.0026 |
| 41 | `x_pay_cnt__active_days` | cross | 0.0102 | 0.0380 |
| 42 | `14d_clk_share` | funnel | 0.0099 | 0.0048 |
| 43 | `7d_clk_share` | funnel | 0.0098 | 0.0125 |
| 44 | `30d_cnv_share` | funnel | 0.0095 | 0.0108 |
| 45 | `14d_ctr` | funnel | 0.0087 | 0.0065 |
| 46 | `3d_cnv_share` | funnel | 0.0077 | 0.0039 |
| 47 | `post_clk_7d_rate` | post | 0.0075 | 0.0005 |
| 48 | `dec_hl1d_dec_cnv` | decay | 0.0064 | 0.0042 |
| 49 | `sess_depth_clk_mean` | session | 0.0060 | 0.0202 |
| 50 | `14d_ctcvr` | funnel | 0.0055 | 0.0129 |
| 51 | `30d_ctcvr` | funnel | 0.0055 | 0.0139 |
| 52 | `14d_cvr` | funnel | 0.0053 | 0.0009 |
| 53 | `7d_n_cnv` | funnel | 0.0052 | 0.0005 |
| 54 | `14d_n_clk` | funnel | 0.0048 | 0.0040 |
| 55 | `x_life_ctcvr__pay_cnt` | cross | 0.0047 | 0.0296 |
| 56 | `7d_cnv_share` | funnel | 0.0047 | 0.0045 |
| 57 | `3d_ctr` | funnel | 0.0046 | 0.0008 |
| 58 | `attr_anyclk2cnv_min_std` | attr | 0.0045 | 0.0016 |
| 59 | `post_n_clk_before_1d_mean` | post | 0.0045 | 0.0002 |
| 60 | `3d_clk_share` | funnel | 0.0042 | 0.0052 |
| 61 | `post_delta_1d_p50` | post | 0.0041 | 0.0010 |
| 62 | `sess_clk_sess_rate` | session | 0.0038 | 0.0100 |
| 63 | `7d_cvr` | funnel | 0.0034 | 0.0002 |
| 64 | `life_cvr` | funnel | 0.0033 | 0.0018 |
| 65 | `trend_n_cnv_3d_minus_14d` | funnel | 0.0031 | 0.0018 |
| 66 | `14d_n_cnv` | funnel | 0.0028 | 0.0008 |
| 67 | `3d_ctcvr` | funnel | 0.0026 | 0.0012 |
| 68 | `30d_cvr` | funnel | 0.0025 | 0.0026 |
| 69 | `7d_ctr` | funnel | 0.0024 | 0.0159 |
| 70 | `post_n_clk_after_1d_mean` | post | 0.0024 | 0.0019 |
| 71 | `life_ctcvr` | funnel | 0.0024 | 0.0201 |
| 72 | `7d_n_clk` | funnel | 0.0021 | 0.0016 |
| 73 | `3d_n_clk` | funnel | 0.0021 | 0.0013 |
| 74 | `x_post_lift_1d_p50__life_cvr` | cross | 0.0019 | 0.0006 |
| 75 | `x_post_clk_1d_rate__attr_wo_clk_rate` | cross | 0.0019 | 0.0021 |
| 76 | `30d_n_clk` | funnel | 0.0017 | 0.0022 |
| 77 | `trans_clk_to_cnv` | markov | 0.0016 | 0.0001 |
| 78 | `life_clk_share` | funnel | 0.0015 | 0.0229 |
| 79 | `3d_cvr` | funnel | 0.0014 | 0.0001 |
| 80 | `30d_clk_share` | funnel | 0.0013 | 0.0138 |
| 81 | `life_ctr` | funnel | 0.0013 | 0.0060 |
| 82 | `post_clk_same_sess_rate` | post | 0.0013 | 0.0001 |
| 83 | `post_clk_1d_rate` | post | 0.0012 | 0.0022 |
| 84 | `item_entropy_clk` | diversity | 0.0012 | 0.0014 |
| 85 | `30d_n_cnv` | funnel | 0.0011 | 0.0008 |
| 86 | `post_lift_1d_p50` | post | 0.0009 | 0.0004 |
| 87 | `post_clk_1h_rate` | post | 0.0005 | 0.0000 |
| 88 | `n_uniq_clk_item` | diversity | 0.0005 | 0.0056 |
| 89 | `post_n_obs_1d` | post | 0.0001 | 0.0014 |
| 90 | `3d_n_cnv` | funnel | 0.0000 | 0.0004 |
| 91 | `attr_wo_clk_cnt` | attr | 0.0000 | 0.0009 |
| 92 | `attr_wo_clk_rate` | attr | 0.0000 | 0.0002 |
| 93 | `attr_cnv_n` | attr | 0.0000 | 0.0008 |
| 94 | `pay_cnt` | money | 0.0000 | 0.0015 |
| 95 | `life_n_clk` | funnel | 0.0000 | 0.0024 |
| 96 | `life_n_cnv` | funnel | 0.0000 | 0.0187 |
| 97 | `attr_firstclk2cnv_min_p50` | attr | 0.0000 | 0.0061 |
| 98 | `attr_clk2cnv_min_p50` | attr | 0.0000 | 0.0000 |
| 99 | `attr_clk2cnv_within_1h_rate` | attr | 0.0000 | 0.0000 |
| 100 | `attr_clk2cnv_within_1d_rate` | attr | 0.0000 | 0.0000 |
| 101 | `attr_clk2cnv_within_5m_rate` | attr | 0.0000 | 0.0000 |
| 102 | `post_same_1h_rate` | post | 0.0000 | 0.0000 |
| 103 | `post_same_1d_rate` | post | 0.0000 | 0.0000 |
| 104 | `post_lift_same_1d_p50` | post | 0.0000 | 0.0000 |
| 105 | `post_clk_5m_rate` | post | 0.0000 | 0.0000 |
| 106 | `attr_clk2cnv_within_30m_rate` | attr | 0.0000 | 0.0000 |
| 107 | `x_post_same_1d_rate__life_cvr` | cross | 0.0000 | 0.0000 |

## `f150`  n=6000 p=150  pos=0.157  RF-domain **0.753**  PO-risk **0.001782**

family mass | RF-domain (P(X)) | PO-VIMP |
|---|---:|---:|
| cross | 0.345 | 0.420 |
| funnel | 0.372 | 0.284 |
| decay | 0.049 | 0.133 |
| session | 0.066 | 0.069 |
| markov | 0.092 | 0.065 |
| attr | 0.010 | 0.028 |
| money | 0.043 | 0.002 |
| diversity | 0.024 | 0.000 |

| rank | feat | family | PO-VIMP | RF-VIMP |
|---|---|---|---:|---:|
| 1 | `life_cnv_share` | funnel | 0.0448 | 0.0136 |
| 2 | `x_life_ctcvr__sess_bounce_rate` | cross | 0.0372 | 0.0186 |
| 3 | `sess_depth_cnv_mean` | session | 0.0334 | 0.0229 |
| 4 | `dec_hl3d_dec_clk` | decay | 0.0306 | 0.0117 |
| 5 | `3d_clk_share` | funnel | 0.0298 | 0.0065 |
| 6 | `x_life_ctr__active_days` | cross | 0.0281 | 0.0235 |
| 7 | `dec_hl1d_dec_clk` | decay | 0.0263 | 0.0089 |
| 8 | `x_life_cvr__7d_ctr` | cross | 0.0259 | 0.0044 |
| 9 | `sess_cnv_sess_rate` | session | 0.0257 | 0.0226 |
| 10 | `dec_hl7d_dec_clk` | decay | 0.0255 | 0.0105 |
| 11 | `dec_hl3d_dec_cnv` | decay | 0.0232 | 0.0067 |
| 12 | `attr_anyclk2cnv_min_std` | attr | 0.0226 | 0.0024 |
| 13 | `x_life_ctr__sess_n` | cross | 0.0204 | 0.0074 |
| 14 | `x_pay_cnt__active_days` | cross | 0.0195 | 0.0268 |
| 15 | `x_life_ctcvr__sess_n` | cross | 0.0178 | 0.0229 |
| 16 | `30d_cnv_share` | funnel | 0.0173 | 0.0217 |
| 17 | `x_life_ctr__attr_anyclk2cnv_min_p50` | cross | 0.0170 | 0.0155 |
| 18 | `trans_cnv_to_exp` | markov | 0.0153 | 0.0236 |
| 19 | `dec_hl7d_dec_cnv` | decay | 0.0148 | 0.0066 |
| 20 | `trend_n_cnv_3d_minus_14d` | funnel | 0.0145 | 0.0016 |
| 21 | `30d_clk_share` | funnel | 0.0144 | 0.0066 |
| 22 | `trans_exp_to_clk` | markov | 0.0142 | 0.0078 |
| 23 | `14d_clk_share` | funnel | 0.0140 | 0.0104 |
| 24 | `trend_cvr_3d_minus_14d` | funnel | 0.0139 | 0.0016 |
| 25 | `x_sess_n__item_entropy_clk` | cross | 0.0135 | 0.0040 |
| 26 | `trans2_exp_cnv_exp` | markov | 0.0134 | 0.0121 |
| 27 | `x_life_cvr__sess_bounce_rate` | cross | 0.0133 | 0.0022 |
| 28 | `trans_exp_to_exp` | markov | 0.0131 | 0.0390 |
| 29 | `x_life_ctcvr__active_days` | cross | 0.0130 | 0.0298 |
| 30 | `log1p_abs_life_ctcvr` | cross | 0.0123 | 0.0314 |
| 31 | `x_life_ctcvr__pay_cnt` | cross | 0.0120 | 0.0152 |
| 32 | `x_hist_len__item_entropy_clk` | cross | 0.0115 | 0.0059 |
| 33 | `3d_cnv_share` | funnel | 0.0111 | 0.0004 |
| 34 | `x_life_cvr__hist_len` | cross | 0.0111 | 0.0015 |
| 35 | `x_life_ctcvr__item_entropy_clk` | cross | 0.0105 | 0.0021 |
| 36 | `x_life_cvr__item_entropy_clk` | cross | 0.0099 | 0.0013 |
| 37 | `x_sess_bounce_rate__item_entropy_clk` | cross | 0.0094 | 0.0096 |
| 38 | `ctvr_pay_per_expose` | funnel | 0.0093 | 0.0257 |
| 39 | `x_pay_cnt__sess_bounce_rate` | cross | 0.0089 | 0.0146 |
| 40 | `life_ctcvr` | funnel | 0.0089 | 0.0126 |
| 41 | `x_pay_cnt__hist_len` | cross | 0.0084 | 0.0057 |
| 42 | `dec_hl1d_dec_cnv` | decay | 0.0084 | 0.0023 |
| 43 | `trend_n_clk_1d_minus_7d` | funnel | 0.0083 | 0.0010 |
| 44 | `14d_ctcvr` | funnel | 0.0082 | 0.0083 |
| 45 | `x_life_ctcvr__attr_anyclk2cnv_min_p50` | cross | 0.0079 | 0.0053 |
| 46 | `14d_cvr` | funnel | 0.0079 | 0.0002 |
| 47 | `x_life_ctr__pay_cnt` | cross | 0.0078 | 0.0027 |
| 48 | `log1p_abs_life_ctr` | cross | 0.0074 | 0.0026 |
| 49 | `x_life_ctr__7d_cvr` | cross | 0.0071 | 0.0011 |
| 50 | `x_life_ctcvr__7d_ctr` | cross | 0.0071 | 0.0036 |
| 51 | `x_life_ctr__life_cvr` | cross | 0.0071 | 0.0014 |
| 52 | `x_life_cvr__sess_n` | cross | 0.0065 | 0.0010 |
| 53 | `7d_cnv_share` | funnel | 0.0061 | 0.0022 |
| 54 | `sess_depth_clk_mean` | session | 0.0061 | 0.0134 |
| 55 | `x_pay_cnt__attr_anyclk2cnv_min_p50` | cross | 0.0060 | 0.0011 |
| 56 | `7d_ctcvr` | funnel | 0.0059 | 0.0028 |
| 57 | `log1p_abs_7d_ctr` | cross | 0.0058 | 0.0096 |
| 58 | `x_7d_cvr__sess_n` | cross | 0.0057 | 0.0014 |
| 59 | `trend_n_cnv_7d_over_30d` | funnel | 0.0056 | 0.0025 |
| 60 | `ui_only_exp_share` | funnel | 0.0055 | 0.0773 |
| 61 | `ui_exp_per_item_mean` | funnel | 0.0052 | 0.0823 |
| 62 | `log1p_abs_attr_anyclk2cnv_min_p50` | cross | 0.0052 | 0.0030 |
| 63 | `x_pay_cnt__item_entropy_clk` | cross | 0.0049 | 0.0003 |
| 64 | `dec_hl7d_dec_arpu` | decay | 0.0046 | 0.0019 |
| 65 | `trend_cvr_7d_minus_30d` | funnel | 0.0045 | 0.0015 |
| 66 | `7d_clk_share` | funnel | 0.0045 | 0.0189 |
| 67 | `x_7d_cvr__sess_bounce_rate` | cross | 0.0041 | 0.0016 |
| 68 | `x_7d_ctr__7d_cvr` | cross | 0.0041 | 0.0004 |
| 69 | `x_pay_cnt__sess_n` | cross | 0.0040 | 0.0047 |
| 70 | `x_life_cvr__active_days` | cross | 0.0040 | 0.0005 |
| 71 | `trans_exp_to_cnv` | markov | 0.0040 | 0.0066 |
| 72 | `x_7d_ctr__pay_cnt` | cross | 0.0038 | 0.0039 |
| 73 | `trend_n_cnv_7d_minus_30d` | funnel | 0.0035 | 0.0005 |
| 74 | `30d_ctcvr` | funnel | 0.0035 | 0.0107 |
| 75 | `sess_clk_sess_rate` | session | 0.0035 | 0.0066 |
| 76 | `trans_clk_to_cnv` | markov | 0.0034 | 0.0005 |
| 77 | `14d_n_clk` | funnel | 0.0034 | 0.0012 |
| 78 | `x_life_ctcvr__hist_len` | cross | 0.0034 | 0.0046 |
| 79 | `life_clk_share` | funnel | 0.0032 | 0.0120 |
| 80 | `14d_log1p_n_clk` | funnel | 0.0031 | 0.0017 |
| 81 | `14d_n_cnv` | funnel | 0.0031 | 0.0010 |
| 82 | `7d_n_clk` | funnel | 0.0030 | 0.0026 |
| 83 | `30d_cvr` | funnel | 0.0028 | 0.0011 |
| 84 | `3d_ctcvr` | funnel | 0.0028 | 0.0020 |
| 85 | `x_life_ctcvr__arpu_mean_proxy` | cross | 0.0027 | 0.0031 |
| 86 | `trend_n_cnv_1d_minus_7d` | funnel | 0.0026 | 0.0003 |
| 87 | `14d_cnv_share` | funnel | 0.0025 | 0.0081 |
| 88 | `attr_anyclk2cnv_min_p90` | attr | 0.0023 | 0.0052 |
| 89 | `x_7d_cvr__item_entropy_clk` | cross | 0.0023 | 0.0010 |
| 90 | `x_arpu_mean_proxy__hist_len` | cross | 0.0022 | 0.0030 |
| 91 | `x_life_ctr__hist_len` | cross | 0.0020 | 0.0188 |
| 92 | `x_7d_cvr__hist_len` | cross | 0.0019 | 0.0000 |
| 93 | `trans_cnv_to_cnv` | markov | 0.0019 | 0.0024 |
| 94 | `3d_log1p_n_clk` | funnel | 0.0018 | 0.0020 |
| 95 | `attr_anyclk2cnv_min_n` | attr | 0.0017 | 0.0004 |
| 96 | `7d_log1p_n_clk` | funnel | 0.0016 | 0.0004 |
| 97 | `x_arpu_mean_proxy__sess_n` | cross | 0.0015 | 0.0020 |
| 98 | `7d_log1p_n_cnv` | funnel | 0.0014 | 0.0013 |
| 99 | `x_life_cvr__arpu_mean_proxy` | cross | 0.0013 | 0.0002 |
| 100 | `x_pay_cnt__arpu_mean_proxy` | cross | 0.0013 | 0.0008 |
| 101 | `x_arpu_sum_proxy__sess_bounce_rate` | cross | 0.0012 | 0.0031 |
| 102 | `arpu_per_active_day_proxy` | money | 0.0011 | 0.0025 |
| 103 | `14d_log1p_n_cnv` | funnel | 0.0011 | 0.0006 |
| 104 | `life_cvr` | funnel | 0.0010 | 0.0009 |
| 105 | `attr_cnv_wo_prior_clk_cnt` | attr | 0.0010 | 0.0013 |
| 106 | `x_arpu_mean_proxy__active_days` | cross | 0.0010 | 0.0035 |
| 107 | `life_n_clk` | funnel | 0.0009 | 0.0009 |
| 108 | `30d_n_cnv` | funnel | 0.0008 | 0.0016 |
| 109 | `30d_n_clk` | funnel | 0.0007 | 0.0024 |
| 110 | `log1p_abs_life_cvr` | cross | 0.0006 | 0.0006 |
| 111 | `trend_n_cnv_3d_over_14d` | funnel | 0.0005 | 0.0004 |
| 112 | `sq_pay_cnt` | cross | 0.0005 | 0.0007 |
| 113 | `life_cnv_price_sum` | money | 0.0002 | 0.0015 |
| 114 | `3d_n_clk` | funnel | 0.0001 | 0.0013 |
| 115 | `log1p_arpu_sum` | money | 0.0001 | 0.0015 |
| 116 | `x_arpu_mean_proxy__sess_bounce_rate` | cross | 0.0001 | 0.0022 |
| 117 | `cvr_pay_per_click` | funnel | 0.0000 | 0.0003 |
| 118 | `3d_n_cnv` | funnel | 0.0000 | 0.0001 |
| 119 | `7d_n_cnv` | funnel | 0.0000 | 0.0026 |
| 120 | `3d_log1p_n_cnv` | funnel | 0.0000 | 0.0004 |
| 121 | `pay_cnt` | money | 0.0000 | 0.0123 |
| 122 | `life_n_cnv` | funnel | 0.0000 | 0.0092 |
| 123 | `n_uniq_cnv_item` | diversity | 0.0000 | 0.0089 |
| 124 | `life_log1p_n_cnv` | funnel | 0.0000 | 0.0002 |
| 125 | `log1p_abs_pay_cnt` | cross | 0.0000 | 0.0002 |
| 126 | `log1p_pay_cnt` | money | 0.0000 | 0.0203 |
| 127 | `pay_user` | money | 0.0000 | 0.0000 |
| 128 | `attr_cnv_wo_prior_clk_rate` | attr | 0.0000 | 0.0002 |
| 129 | `n_uniq_clk_item` | diversity | 0.0000 | 0.0020 |
| 130 | `x_7d_cvr__attr_clk2cnv_min_p50` | cross | 0.0000 | 0.0012 |
| 131 | `1d_log1p_n_cnv` | funnel | 0.0000 | 0.0000 |
| 132 | `7d_cvr` | funnel | 0.0000 | 0.0009 |
| 133 | `log1p_abs_arpu_mean_proxy` | cross | 0.0000 | 0.0034 |
| 134 | `30d_cnv_price_mean` | money | 0.0000 | 0.0000 |
| 135 | `sq_arpu_mean_proxy` | cross | 0.0000 | 0.0014 |
| 136 | `life_log1p_n_clk` | funnel | 0.0000 | 0.0044 |
| 137 | `30d_log1p_n_clk` | funnel | 0.0000 | 0.0017 |
| 138 | `log1p_abs_7d_cvr` | cross | 0.0000 | 0.0004 |
| 139 | `item_entropy_clk` | diversity | 0.0000 | 0.0029 |
| 140 | `log1p_abs_item_entropy_clk` | cross | 0.0000 | 0.0005 |
| 141 | `arpu_mean_proxy` | money | 0.0000 | 0.0028 |
| 142 | `life_cnv_price_mean` | money | 0.0000 | 0.0009 |
| 143 | `log1p_abs_arpu_sum_proxy` | cross | 0.0000 | 0.0013 |
| 144 | `arpu_p50_proxy` | money | 0.0000 | 0.0014 |
| 145 | `item_entropy_cnv` | diversity | 0.0000 | 0.0103 |
| 146 | `30d_log1p_n_cnv` | funnel | 0.0000 | 0.0005 |
| 147 | `arpu_sum_proxy` | money | 0.0000 | 0.0004 |
| 148 | `1d_n_cnv` | funnel | 0.0000 | 0.0007 |
| 149 | `sq_item_entropy_clk` | cross | 0.0000 | 0.0043 |
| 150 | `x_arpu_sum_proxy__arpu_mean_proxy` | cross | 0.0000 | 0.0018 |
