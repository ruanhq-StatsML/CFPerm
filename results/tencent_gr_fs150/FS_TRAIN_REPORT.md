# TencentGR auto-1000 → feature-select → train

- users: **6000**
- generated dims: **1000** (target 1000)
- used after leak-drop: **1000**
- selected: **150** via `f`
- label: `future_cnv` (pos rate=0.157)
- vs previous 128: overlap **128**, only-new **22**, dropped **0**

## Metrics

```json
{
  "hgb": {
    "auc": 0.8289681271549912,
    "ap": 0.4913645915722099,
    "acc": 0.8546666666666667,
    "sec": 0.19318842887878418,
    "n_selected": 150
  },
  "logreg": {
    "auc": 0.8140795559666975,
    "ap": 0.49058922761210694,
    "acc": 0.8553333333333333,
    "sec": 2.5088274478912354,
    "n_selected": 150
  }
}
```

## Top-20 selected

```json
[
  {
    "name": "life_log1p_n_cnv",
    "score": 1090.419620914414
  },
  {
    "name": "log1p_pay_cnt",
    "score": 1090.419620914414
  },
  {
    "name": "log1p_abs_pay_cnt",
    "score": 1090.419620914414
  },
  {
    "name": "trans_cnv_to_exp",
    "score": 1086.6931851891777
  },
  {
    "name": "trans_exp_to_cnv",
    "score": 1074.140827566234
  },
  {
    "name": "30d_log1p_n_cnv",
    "score": 1054.2762717589294
  },
  {
    "name": "item_entropy_cnv",
    "score": 1052.248090203517
  },
  {
    "name": "sess_cnv_sess_rate",
    "score": 1001.0841463781658
  },
  {
    "name": "dec_hl7d_dec_cnv",
    "score": 1001.0097728556357
  },
  {
    "name": "sess_depth_cnv_mean",
    "score": 973.4289935569008
  },
  {
    "name": "life_cnv_share",
    "score": 955.1276348195488
  },
  {
    "name": "30d_n_cnv",
    "score": 942.3901172986189
  },
  {
    "name": "n_uniq_cnv_item",
    "score": 903.9045001562869
  },
  {
    "name": "life_n_cnv",
    "score": 902.891119820815
  },
  {
    "name": "pay_cnt",
    "score": 902.891119820815
  },
  {
    "name": "attr_cnv_wo_prior_clk_cnt",
    "score": 902.4805085175046
  },
  {
    "name": "ui_only_exp_share",
    "score": 894.205593113959
  },
  {
    "name": "ui_exp_per_item_mean",
    "score": 891.2403609264442
  },
  {
    "name": "log1p_abs_life_ctcvr",
    "score": 873.3539457009244
  },
  {
    "name": "x_pay_cnt__sess_bounce_rate",
    "score": 864.1658129790109
  }
]
```

## 150-feature list

1. `life_log1p_n_cnv`
2. `log1p_pay_cnt`
3. `log1p_abs_pay_cnt`
4. `trans_cnv_to_exp`
5. `trans_exp_to_cnv`
6. `30d_log1p_n_cnv`
7. `item_entropy_cnv`
8. `sess_cnv_sess_rate`
9. `dec_hl7d_dec_cnv`
10. `sess_depth_cnv_mean`
11. `life_cnv_share`
12. `30d_n_cnv`
13. `n_uniq_cnv_item`
14. `life_n_cnv`
15. `pay_cnt`
16. `attr_cnv_wo_prior_clk_cnt`
17. `ui_only_exp_share`
18. `ui_exp_per_item_mean`
19. `log1p_abs_life_ctcvr`
20. `x_pay_cnt__sess_bounce_rate`
21. `30d_cnv_share`
22. `14d_log1p_n_cnv`
23. `trans_exp_to_exp`
24. `14d_n_cnv`
25. `trans2_exp_cnv_exp`
26. `trend_n_cnv_7d_minus_30d`
27. `dec_hl3d_dec_cnv`
28. `x_pay_cnt__hist_len`
29. `x_pay_cnt__sess_n`
30. `x_life_ctcvr__hist_len`
31. `7d_log1p_n_cnv`
32. `x_life_ctcvr__sess_n`
33. `trend_n_cnv_3d_minus_14d`
34. `life_ctcvr`
35. `ctvr_pay_per_expose`
36. `7d_n_cnv`
37. `x_life_ctr__life_cvr`
38. `x_pay_cnt__active_days`
39. `14d_cnv_share`
40. `attr_anyclk2cnv_min_n`
41. `x_life_cvr__item_entropy_clk`
42. `attr_cnv_wo_prior_clk_rate`
43. `pay_user`
44. `log1p_abs_life_cvr`
45. `x_life_ctcvr__active_days`
46. `x_life_ctcvr__sess_bounce_rate`
47. `30d_ctcvr`
48. `x_pay_cnt__item_entropy_clk`
49. `trend_n_cnv_1d_minus_7d`
50. `log1p_abs_attr_anyclk2cnv_min_p50`
51. `3d_log1p_n_cnv`
52. `x_life_ctcvr__item_entropy_clk`
53. `3d_n_cnv`
54. `14d_ctcvr`
55. `dec_hl1d_dec_cnv`
56. `7d_cnv_share`
57. `trans_cnv_to_cnv`
58. `life_cvr`
59. `cvr_pay_per_click`
60. `x_life_ctr__pay_cnt`
61. `30d_cvr`
62. `x_life_cvr__sess_bounce_rate`
63. `7d_log1p_n_clk`
64. `14d_cvr`
65. `x_life_ctcvr__arpu_mean_proxy`
66. `14d_log1p_n_clk`
67. `x_pay_cnt__arpu_mean_proxy`
68. `sq_pay_cnt`
69. `x_life_cvr__hist_len`
70. `life_log1p_n_clk`
71. `7d_ctcvr`
72. `x_life_cvr__sess_n`
73. `30d_log1p_n_clk`
74. `x_life_ctr__7d_cvr`
75. `log1p_abs_7d_cvr`
76. `x_7d_cvr__item_entropy_clk`
77. `dec_hl3d_dec_clk`
78. `dec_hl7d_dec_clk`
79. `x_7d_ctr__pay_cnt`
80. `7d_n_clk`
81. `log1p_abs_item_entropy_clk`
82. `3d_cnv_share`
83. `x_life_cvr__active_days`
84. `7d_clk_share`
85. `3d_log1p_n_clk`
86. `item_entropy_clk`
87. `x_life_ctcvr__attr_anyclk2cnv_min_p50`
88. `attr_anyclk2cnv_min_p90`
89. `log1p_abs_7d_ctr`
90. `x_life_ctr__attr_anyclk2cnv_min_p50`
91. `life_cnv_price_mean`
92. `arpu_mean_proxy`
93. `log1p_arpu_sum`
94. `log1p_abs_arpu_sum_proxy`
95. `x_sess_bounce_rate__item_entropy_clk`
96. `arpu_p50_proxy`
97. `x_arpu_mean_proxy__sess_bounce_rate`
98. `x_life_cvr__7d_ctr`
99. `trend_n_cnv_3d_over_14d`
100. `14d_clk_share`
101. `14d_n_clk`
102. `30d_cnv_price_mean`
103. `log1p_abs_arpu_mean_proxy`
104. `log1p_abs_life_ctr`
105. `x_life_ctr__hist_len`
106. `sq_arpu_mean_proxy`
107. `attr_anyclk2cnv_min_std`
108. `trend_n_cnv_7d_over_30d`
109. `x_life_ctr__sess_n`
110. `x_pay_cnt__attr_anyclk2cnv_min_p50`
111. `3d_ctcvr`
112. `x_life_ctr__active_days`
113. `x_7d_ctr__7d_cvr`
114. `arpu_per_active_day_proxy`
115. `sess_clk_sess_rate`
116. `7d_cvr`
117. `dec_hl1d_dec_clk`
118. `trend_n_clk_1d_minus_7d`
119. `life_clk_share`
120. `x_7d_cvr__sess_bounce_rate`
121. `30d_n_clk`
122. `x_7d_cvr__attr_clk2cnv_min_p50`
123. `n_uniq_clk_item`
124. `1d_log1p_n_cnv`
125. `life_n_clk`
126. `x_life_cvr__arpu_mean_proxy`
127. `x_arpu_mean_proxy__hist_len`
128. `life_cnv_price_sum`
129. `arpu_sum_proxy`
130. `x_life_ctcvr__7d_ctr`
131. `x_arpu_sum_proxy__sess_bounce_rate`
132. `sess_depth_clk_mean`
133. `3d_n_clk`
134. `30d_clk_share`
135. `trend_cvr_3d_minus_14d`
136. `x_arpu_mean_proxy__sess_n`
137. `x_hist_len__item_entropy_clk`
138. `trans_exp_to_clk`
139. `trend_cvr_7d_minus_30d`
140. `1d_n_cnv`
141. `dec_hl7d_dec_arpu`
142. `x_7d_cvr__hist_len`
143. `sq_item_entropy_clk`
144. `x_7d_cvr__sess_n`
145. `x_life_ctcvr__pay_cnt`
146. `x_sess_n__item_entropy_clk`
147. `trans_clk_to_cnv`
148. `x_arpu_mean_proxy__active_days`
149. `3d_clk_share`
150. `x_arpu_sum_proxy__arpu_mean_proxy`

## vs 128 board

- overlap: 128
- new (not in 128): `life_cnv_price_sum`, `x_life_ctcvr__7d_ctr`, `x_arpu_sum_proxy__sess_bounce_rate`, `sess_depth_clk_mean`, `3d_n_clk`, `30d_clk_share`, `trend_cvr_3d_minus_14d`, `x_arpu_mean_proxy__sess_n`, `x_hist_len__item_entropy_clk`, `trans_exp_to_clk`, `trend_cvr_7d_minus_30d`, `1d_n_cnv`, `dec_hl7d_dec_arpu`, `x_7d_cvr__hist_len`, `sq_item_entropy_clk`, `x_7d_cvr__sess_n`, `x_life_ctcvr__pay_cnt`, `x_sess_n__item_entropy_clk`, `trans_clk_to_cnv`, `x_arpu_mean_proxy__active_days`, `3d_clk_share`, `x_arpu_sum_proxy__arpu_mean_proxy`
- dropped from 128: (none)

## Recipe

1. **Generate**: combinatorial windows × funnel × decay × session × attribution × ARPU/deal × Markov × TOD × crosses ≈ 1000.
2. **Select**: VarianceThreshold → SelectKBest(F / MI) → top-k.
3. **Train**: HistGradientBoosting + LogisticRegression holdout AUC/AP.

Leakage note: for `pay_user`, raw `pay_cnt` / `life_n_cnv` / `arpu_sum` are dropped before selection.
Label `future_cnv` uses prefix features and suffix conversion, so conversion-count features are not a direct leak.
