# CausalForest parsimonious empirical smoke

> **Not a causal ATE claim.** `W` = period (W1/W2). CF is a lean empirical tool
> for heterogeneity / feature rank next to **MMD localization**.

- gap=30d | rows=8000 | feats=21 | CF trees=60
- global mean(τ̂²) ≈ **452.815690** (PO-risk-like scalar)
- merchant Spearman(CF τ̂², MMD) = **-0.330**
- merchant top-30 overlap = **0.000**
- user Spearman = **0.077** | overlap@30 = **0.000**

## CF feature importance (head)
| rank | feature | importance |
|---:|---|---:|
| 1 | `i_n_covisit_neighbors` | 0.2557 |
| 2 | `i_log1p_n_covisit` | 0.2337 |
| 3 | `i_share_linear` | 0.1890 |
| 4 | `i_credit_linear` | 0.1768 |
| 5 | `ui_pop_mismatch` | 0.0471 |
| 6 | `u_n_exp` | 0.0306 |
| 7 | `u_n_uniq_items` | 0.0163 |
| 8 | `i_share_first` | 0.0084 |
| 9 | `u_log1p_n_events` | 0.0078 |
| 10 | `u_span_sec` | 0.0073 |
| 11 | `i_n_exp` | 0.0057 |
| 12 | `u_n_events` | 0.0056 |

## Role in the pipeline
- **MMD path**: subset localization (who drifted)
- **CF path**: parsimonious tau-hat / feature_importances_ (alternate drift score + ranking)
- **Stop-drill / business logic**: unchanged — CF does not replace drill triggers

## Takeaway from this smoke
CF merchant/user rankings **do not track MMD** here (Spearman merchant≈-0.33, top-30 overlap=0).
So CF is a useful **parallel empirical lens** (esp. feature_importances_), not a drop-in
replacement for MMD localization. Keep both; do not justify stop-drill from CF alone.
