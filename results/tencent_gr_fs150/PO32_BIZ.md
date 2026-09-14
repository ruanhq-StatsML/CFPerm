# PO-32 业务解读（对着 one-pass / asof 看）

选出来的不是「买过多少次」，是 **路径长什么样的人后面还会买**。
`pay_cnt` 三件套没进板。进板的是 share / 衰减点击 / 场碎 × 会买 / 转移 / 决策时长散度。

---

## 1. 漏斗：买在行为里占多大（不是次数）

循环里 `n[w][cnv] / n[w][tot]`、`n_cnv/n_exp`。

| 特征 | 在代码里 | 业务 |
|---|---|---|
| life_cnv_share | 全史转化 / 总事件 | 这个人的轨迹里「买」是常态还是偶尔一次 |
| 30d_cnv_share | 近 30 天同上 | 最近还在不在买 |
| life_ctcvr / log1p_abs_life_ctcvr | 看见→买 | 曝光有没有变成单 |
| 3d/14d/30d_clk_share | 窗内点击占比 | 最近手热不热（点，还没到买） |
| trend_n_cnv_3d_minus_14d | 3d 买次数 − 14d | 近 3 天买在掉（常 ≤0） |
| trend_cvr_3d_minus_14d | 3d cvr − 14d cvr | 点了更不爱买了 |

## 2. 衰减：还在不在点/买

循环里 `2^{-age/hl}`。比硬截 7 天滑。

| 特征 | 业务 |
|---|---|
| dec_hl1d/3d/7d_dec_clk | 昨天 / 三天 / 一周内点击热度还在不在 |
| dec_hl3d/7d_dec_cnv | 近期成交热度 |

板子更吃 **clk 衰减**（还在逛），成交衰减次之。人走了先是不点。

## 3. 场次：这场逛完没有、买了没有

30min 没动 `close_sess()`。

| 特征 | 业务 |
|---|---|
| sess_depth_cnv_mean | 一场里成交几下（深度买家 vs 点一下走） |
| sess_cnv_sess_rate | 多少场真正成交（不是来了就买） |
| sess_n / bounce | 来了多少场、看一条就走的比例（多在 cross 里） |

决策 >30min 的转化会被切成新场，还可能算 bounce——这是定义，不是 bug。

## 4. 归因：点完多久买、稳不稳

asof 中间表：每笔转化一行。`anyclk` = 全局上次点；同商品 last-click 在表里但 PO-32 只用了任意点击的 **标准差**。

| 特征 | 业务 |
|---|---|
| attr_anyclk2cnv_min_std | 有时秒下、有时纠结很久 → 决策不稳定 |
| x_life_ctr__attr_anyclk2cnv_min_p50 | 爱点 × 下单慢：手勤但犹豫 |

没选 `wo_prior_clk`：漏点/直达不是这块的主故事。

## 5. 转移：上一步到这一步

`prev_act → act`。

| 特征 | 业务 |
|---|---|
| trans_exp_to_clk | 看见就点，漏斗在走 |
| trans_exp_to_exp | 连着看、不点 |
| trans_cnv_to_exp | 买完又逛（复购/连带） |
| trans2_exp_cnv_exp | 看见→买→又看 |

## 6. Cross：两件事同时成立（循环后再乘）

不是新扫描。把「会买」和「逛的样子」焊在一起。

| 特征 | 读法 |
|---|---|
| x_life_ctcvr__sess_bounce_rate | 会买，但场很碎（高意向、低耐心） |
| x_life_cvr__sess_bounce_rate | 点了会买 × 看一眼走 |
| x_life_ctr__active_days | 爱点 × 来的天数多 |
| x_life_ctcvr__active_days | 会买 × 活跃天数 |
| x_life_ctr / ctcvr __ sess_n | 漏斗好 × 来的场次多 |
| x_pay_cnt__active_days | 买过的次数摊在活跃天上（强度密度） |
| x_life_cvr__7d_ctr | 会买 × 近 7 天还在点 |
| x_sess_n__item_entropy_clk | 场次多 × 点得杂（广撒网） |
| x_hist_len__item_entropy_clk | 历史长 × 点得杂 |
| x_life_ctcvr__pay_cnt | 看见就买 × 买过的单量 |

---

一句话：这 32 个在刻画 **「还会买的人，最近还在点、场里能成交、买完还逛、决策时长不稳定」**，不是「历史单量最大的人」。
