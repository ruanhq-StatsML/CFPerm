# Hop FSDS prototype（user / order / item / merchant）

不是图谱 GNN。跳：`item → merchant`（上架），`cnv → order → item → merchant`。
店名是抽的英文（或 `--messy` / `--vllm`）。订单号 `randint` 在转化行上。
PO-risk φ=(Y−μ)(W−e)；LOGO=整跳拿掉；LOCO=一列拿掉。R 不当检验。

店粒：该店自己的时间 50% 切开，左 X 右 Y=`ctr`（右窗 n_clk/n_exp），W=`1{t_end>中位}`。
订单粒：Y=满窗 `y_post_clk_1d`，X=用户历史 + 店聚合 + 当单 ctx。

## 店粒

n=400 p=18 pos=0.034 W1=0.500  RF-domain **0.639**  PO-risk **0.000001**

| hop | n | RF-domain | PO-VIMP | LOGO Δ | LOGO share |
|---|---:|---:|---:|---:|---:|
| structure | 3 | 0.110 | 0.171 | +0.00000 | 0.474 |
| volume | 6 | 0.445 | 0.221 | +0.00000 | 0.392 |
| user_hop | 5 | 0.281 | 0.283 | +0.00000 | 0.133 |
| recency | 2 | 0.065 | 0.073 | -0.00000 | 0.000 |
| shop_draw | 2 | 0.099 | 0.253 | -0.00000 | 0.000 |

| rank | feat | hop | PO-VIMP | RF-VIMP |
|---|---|---|---:|---:|
| 1 | `m_star` | shop_draw | 0.1352 | 0.0680 |
| 2 | `m_open_days` | shop_draw | 0.1174 | 0.0308 |
| 3 | `m_n_item` | volume | 0.1073 | 0.1618 |
| 4 | `m_user_dec_clk` | user_hop | 0.0861 | 0.0885 |
| 5 | `m_ctr` | structure | 0.0677 | 0.0356 |
| 6 | `m_cnv_share` | structure | 0.0672 | 0.0581 |
| 7 | `m_n_user` | volume | 0.0647 | 0.1946 |
| 8 | `m_user_life_ctr` | user_hop | 0.0635 | 0.0256 |
| 9 | `m_user_bounce` | user_hop | 0.0539 | 0.0547 |
| 10 | `m_dec_cnv` | recency | 0.0478 | 0.0450 |
| 11 | `m_user_same_sess` | user_hop | 0.0432 | 0.0535 |
| 12 | `m_user_life_cnv_share` | user_hop | 0.0361 | 0.0588 |
| 13 | `m_cvr` | structure | 0.0359 | 0.0158 |
| 14 | `m_dec_clk` | recency | 0.0249 | 0.0204 |
| 15 | `m_n_exp` | volume | 0.0184 | 0.0621 |
| 16 | `m_n_cnv` | volume | 0.0136 | 0.0066 |

LOCO ΔR（>0 才是拿掉后 R 下降）：

| feat | hop | impurity | LOCO ΔR |
|---|---|---:|---:|
| `m_user_same_sess` | user_hop | 0.0432 | +1.978696e-07 |
| `m_open_days` | shop_draw | 0.1174 | +1.017525e-07 |
| `m_n_clk` | volume | 0.0054 | +9.528812e-08 |
| `m_n_exp` | volume | 0.0184 | +5.617686e-08 |
| `m_cvr` | structure | 0.0359 | +5.347744e-08 |
| `m_n_user` | volume | 0.0647 | +4.561230e-08 |
| `m_user_dec_clk` | user_hop | 0.0861 | +2.867648e-08 |
| `m_n_cnv` | volume | 0.0136 | +2.260965e-08 |
| `m_dec_clk` | recency | 0.0249 | +1.928431e-08 |
| `m_ctr` | structure | 0.0677 | +1.805297e-08 |
| `m_n_order` | volume | 0.0117 | +8.538457e-09 |
| `m_star` | shop_draw | 0.1352 | +8.263659e-09 |

## 订单粒（三跳）

n=12564 p=28 pos=0.113 W1=0.500  RF-domain **0.763**  PO-risk **0.000192**

| hop | n | RF-domain | PO-VIMP | LOGO Δ | LOGO share |
|---|---:|---:|---:|---:|---:|
| merchant | 18 | 0.079 | 0.443 | +0.00004 | 1.000 |
| order | 5 | 0.235 | 0.120 | -0.00003 | 0.000 |
| user | 5 | 0.686 | 0.436 | -0.00024 | 0.000 |

| rank | feat | hop | PO-VIMP | RF-VIMP |
|---|---|---|---:|---:|
| 1 | `dec_hl7d_dec_clk` | user | 0.1337 | 0.1284 |
| 2 | `life_ctr` | user | 0.1165 | 0.0624 |
| 3 | `n_prior_cnv` | order | 0.0727 | 0.1760 |
| 4 | `sess_bounce_rate` | user | 0.0671 | 0.1258 |
| 5 | `life_cnv_share` | user | 0.0607 | 0.3481 |
| 6 | `post_clk_same_sess_rate` | user | 0.0581 | 0.0217 |
| 7 | `m_open_days` | merchant | 0.0575 | 0.0029 |
| 8 | `m_user_life_ctr` | merchant | 0.0472 | 0.0070 |
| 9 | `m_user_dec_clk` | merchant | 0.0444 | 0.0052 |
| 10 | `m_star` | merchant | 0.0372 | 0.0041 |
| 11 | `m_user_life_cnv_share` | merchant | 0.0371 | 0.0047 |
| 12 | `m_user_bounce` | merchant | 0.0334 | 0.0087 |
| 13 | `m_dec_cnv` | merchant | 0.0330 | 0.0065 |
| 14 | `m_ctr` | merchant | 0.0321 | 0.0050 |
| 15 | `m_cnv_share` | merchant | 0.0301 | 0.0044 |
| 16 | `n_clk_before_1d` | order | 0.0267 | 0.0222 |

LOCO ΔR（>0 才是拿掉后 R 下降）：

| feat | hop | impurity | LOCO ΔR |
|---|---|---:|---:|
| `n_prior_cnv` | order | 0.0727 | +3.490445e-05 |
| `m_n_cnv` | merchant | 0.0058 | +2.561340e-05 |
| `m_user_life_ctr` | merchant | 0.0472 | +2.523569e-05 |
| `m_n_order` | merchant | 0.0048 | +1.962866e-05 |
| `m_n_clk` | merchant | 0.0027 | +1.856517e-05 |
| `m_n_user` | merchant | 0.0136 | +1.743558e-05 |
| `empty_any` | order | 0.0087 | +1.718817e-05 |
| `m_user_dec_clk` | merchant | 0.0444 | +1.597380e-05 |
| `life_ctr` | user | 0.1165 | +1.480483e-05 |
| `m_cnv_share` | merchant | 0.0301 | +1.410068e-05 |
| `m_user_life_cnv_share` | merchant | 0.0371 | +1.116564e-05 |
| `m_n_exp` | merchant | 0.0159 | +9.197833e-06 |

样例店名： Delta Bazaar Market, Silver Goods Studio, Maple Emporium Co, Jade Depot Ltd

`python3 scripts/tencent_gr/hop_fsds_proto.py`
