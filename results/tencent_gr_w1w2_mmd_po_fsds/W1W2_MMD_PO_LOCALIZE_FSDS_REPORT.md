# 行为 / 购买欲变动归因：MMD + PO-risk + conditional-mean → FSDS

直观 concise 链路（**先 Standardization**）：
1. W1/W2 独立 FE → standardized X
2. **Subset localization**：conditional-mean / MMD² / PO-risk
3. **可视化** 漂移商品 subset（用户行为 & 购买欲变动）
4. **FSDS** 标准流 → 归因特征 ranking
5. （可选）GT 订单/商品 list → hit / precision / recall@k

## Protocol
0. **StandardScaler** fit on W1（pipeline 最上面）
1. FE ×2，gap = **30**d
2. Item subset rank-average(cmean, MMD, PO) → top-**200**
3. Viz subset
4. FSDS: **StandardScaler** → var → SelectKBest(k=15) → HGB/LogReg
5. GT evaluator（`--gt-items` / `--gt-orders`）

## Localized subset (head)
| rank | item_id | cmean | MMD² | PO τ² | n_W1 | n_W2 |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 226288 | 40.606 | 1.7227 | 0.0000 | 65 | 170 |
| 2 | 3885368 | 16.317 | 1.5826 | 0.0000 | 23 | 49 |
| 3 | 749302 | 21.844 | 1.4055 | 0.0000 | 8 | 60 |
| 4 | 53521 | 21.736 | 1.4754 | 0.0000 | 9 | 68 |
| 5 | 560077 | 41.236 | 1.3424 | 0.0000 | 4 | 120 |
| 6 | 2189611 | 12.806 | 1.5225 | 0.0000 | 5 | 31 |
| 7 | 3757633 | 16.021 | 1.6812 | 0.0000 | 18 | 37 |
| 8 | 3576787 | 20.600 | 1.5551 | 0.0000 | 12 | 68 |
| 9 | 345842 | 12.252 | 1.5350 | 0.0000 | 22 | 11 |
| 10 | 1510013 | 66.433 | 1.3910 | 0.0000 | 6 | 167 |
| 11 | 1971100 | 21.451 | 1.2927 | 0.0000 | 9 | 65 |
| 12 | 2535188 | 29.231 | 1.3583 | 0.0000 | 14 | 88 |

global PO-risk = **0.000001**

## FSDS feature ranking
| rank | feature | F-score | selected |
|---:|---|---:|---:|
| 1 | `i_share_linear` | 1.061 | 1 |
| 2 | `i_credit_linear` | 1.061 | 1 |
| 3 | `i_share_first` | 1.057 | 1 |
| 4 | `i_credit_first` | 1.057 | 1 |
| 5 | `i_log1p_n_users` | 0.840 | 1 |
| 6 | `i_log1p_n_covisit` | 0.773 | 1 |
| 7 | `i_n_covisit_neighbors` | 0.697 | 1 |
| 8 | `i_log1p_n_exp` | 0.666 | 1 |
| 9 | `i_n_users` | 0.652 | 1 |
| 10 | `i_n_exp` | 0.520 | 1 |
| 11 | `ui_pop_mismatch` | 0.291 | 1 |
| 12 | `u_log1p_n_events` | 0.221 | 1 |
| 13 | `u_n_events` | 0.063 | 1 |
| 14 | `u_log1p_n_uniq` | 0.040 | 1 |
| 15 | `i_share_last` | 0.034 | 1 |

## Holdout
- W1 user-holdout: n=962/314 | ranking only | selected: `u_n_events`, `u_log1p_n_events`, `u_log1p_n_uniq`, `i_n_exp`, `i_n_users`, `i_credit_first`, `i_credit_linear`, `i_share_first`
- W2 temporal: n=962/8815 | hgb AUC=0.445 | logreg AUC=0.358 | selected: `u_n_events`, `u_log1p_n_events`, `u_log1p_n_uniq`, `i_n_exp`, `i_n_users`, `i_credit_first`, `i_credit_linear`, `i_share_first`

## GT evaluator
- items: n_gt=42 | P@100=0.400 | R@100=0.952
- orders: n=31 | coverage=0.968
