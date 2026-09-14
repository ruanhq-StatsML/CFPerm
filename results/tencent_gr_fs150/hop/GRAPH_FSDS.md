# Funnel FSDS + networkx graph features

PO-risk **没有因果性**。φ=(Y−μ)(W−e) 是早/晚对买后 Y 的距离，不是谁导致转化。
图谱：左窗 clk/cnv 的 user—item、user—merchant 二部图，店投影 PageRank/聚类，同场共点。
LOGO 把 graph 整包拿掉。不是 GNN，不是 DFS。

## 订单粒 漏斗 ∪ graph

n=12564 p=18 pos=0.113 W1=0.500  RF-domain **0.781**  PO-risk **0.000175**

| hop | n | RF-domain | PO-VIMP | LOGO Δ | LOGO share |
|---|---:|---:|---:|---:|---:|
| funnel_order | 5 | 0.162 | 0.161 | +0.00004 | 0.617 |
| graph | 8 | 0.176 | 0.297 | +0.00002 | 0.383 |
| funnel_user | 5 | 0.662 | 0.541 | -0.00025 | 0.000 |

| rank | feat | hop | PO-VIMP | RF-VIMP |
|---|---|---|---:|---:|
| 1 | `dec_hl7d_dec_clk` | funnel_user | 0.1706 | 0.1685 |
| 2 | `life_ctr` | funnel_user | 0.1231 | 0.0697 |
| 3 | `n_prior_cnv` | funnel_order | 0.1063 | 0.1093 |
| 4 | `g_m_clust` | graph | 0.0920 | 0.0048 |
| 5 | `life_cnv_share` | funnel_user | 0.0907 | 0.3350 |
| 6 | `sess_bounce_rate` | funnel_user | 0.0857 | 0.0782 |
| 7 | `post_clk_same_sess_rate` | funnel_user | 0.0713 | 0.0105 |
| 8 | `g_m_pr` | graph | 0.0583 | 0.0063 |
| 9 | `g_m_user_deg` | graph | 0.0398 | 0.0053 |
| 10 | `g_m_proj_deg` | graph | 0.0315 | 0.0032 |
| 11 | `n_clk_before_1d` | funnel_order | 0.0298 | 0.0379 |
| 12 | `g_i_user_deg` | graph | 0.0291 | 0.0040 |
| 13 | `g_u_item_deg` | graph | 0.0226 | 0.0773 |
| 14 | `g_u_merch_deg` | graph | 0.0215 | 0.0751 |
| 15 | `log1p_price` | funnel_order | 0.0132 | 0.0115 |
| 16 | `sess_clk_before` | funnel_order | 0.0070 | 0.0005 |

LOCO ΔR（>0 才是拿掉后 R 下降）：

| feat | hop | impurity | LOCO ΔR |
|---|---|---:|---:|
| `n_prior_cnv` | funnel_order | 0.1063 | +3.382469e-05 |
| `life_ctr` | funnel_user | 0.1231 | +3.038211e-05 |
| `log1p_price` | funnel_order | 0.0132 | +2.206327e-05 |
| `g_m_clust` | graph | 0.0920 | +2.192832e-05 |
| `g_m_pr` | graph | 0.0583 | +2.084229e-05 |
| `sess_clk_before` | funnel_order | 0.0070 | +2.063166e-05 |
| `g_u_item_deg` | graph | 0.0226 | +2.005091e-05 |
| `n_clk_before_1d` | funnel_order | 0.0298 | +1.959166e-05 |
| `post_clk_same_sess_rate` | funnel_user | 0.0713 | +1.785751e-05 |
| `g_i_user_deg` | graph | 0.0291 | +1.740579e-05 |
| `g_u_merch_deg` | graph | 0.0215 | +1.463403e-05 |
| `g_m_proj_deg` | graph | 0.0315 | +1.231647e-05 |
