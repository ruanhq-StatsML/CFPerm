# Edge grain schema (auto)

- rows=286949 · users=3000 · items=194344
- grain: `one_row_equals_one_user_item_edge`
- sort: `e_last_ts` · Y: `y_convert`
- ops |X|=35 · content |X|=14

## Equal-count calendar spans (smoke)

| N | chunk | span_h | y_rate |
|---:|---:|---:|---:|
| 1000 | 0 | 1132.6 | 0.0520 |
| 1000 | 1 | 717.1 | 0.0420 |
| 1000 | 2 | 418.8 | 0.0510 |
| 1000 | 3 | 268.6 | 0.0700 |
| 1000 | 4 | 204.2 | 0.0630 |
| 2000 | 0 | 1850.9 | 0.0470 |
| 2000 | 1 | 687.5 | 0.0605 |
| 2000 | 2 | 378.5 | 0.0620 |
| 2000 | 3 | 278.2 | 0.0605 |
| 2000 | 4 | 225.9 | 0.0680 |

## ops_X

`e_n_exp`, `e_n_clk`, `e_ctr`, `u_n_events`, `u_n_exp`, `u_n_clk`, `u_n_uniq_items`, `u_ctr`, `u_span_sec`, `u_log1p_n_events`, `u_log1p_n_uniq_items`, `u_log1p_n_exp`, `u_log1p_n_clk`, `u_user_activity_rank`, `i_n_exp`, `i_n_clk`, `i_n_users`, `i_ctr`, `i_credit_first`, `i_credit_last`, `i_credit_linear`, `i_share_first`, `i_share_last`, `i_share_linear`, `i_n_covisit_neighbors`, `i_log1p_n_exp`, `i_log1p_n_clk`, `i_log1p_n_users`, `i_log1p_n_covisit_neighbors`, `i_log1p_credit_linear`, `i_item_pop_rank`, `i_item_credit_rank`, `e_log1p_exp`, `e_log1p_clk`, `ui_pop_mismatch`

## content_X

`u_n_uniq_items`, `u_span_sec`, `u_log1p_n_uniq_items`, `i_credit_first`, `i_credit_last`, `i_credit_linear`, `i_share_first`, `i_share_last`, `i_share_linear`, `i_n_covisit_neighbors`, `i_log1p_n_covisit_neighbors`, `i_log1p_credit_linear`, `i_item_credit_rank`, `ui_pop_mismatch`

