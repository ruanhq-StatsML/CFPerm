# Streamlined: Standardization → subset-MMD → FSDS

```
StandardScaler(W1) → item MMD²(W1,W2) → FSDS ranking
```

## Protocol
1. **Standardization** fit on W1
2. **Subset-level MMD** → top-**200** items (gap=30d)
3. **FSDS**: StandardScaler → var → SelectKBest(k=15) → HGB/LogReg

## Localized subset (by MMD²)
| rank | item_id | MMD² | n_W1 | n_W2 |
|---:|---:|---:|---:|---:|
| 1 | 226288 | 1.7191 | 65 | 170 |
| 2 | 3757633 | 1.6812 | 18 | 37 |
| 3 | 86522 | 1.6342 | 4 | 13 |
| 4 | 2708068 | 1.6155 | 10 | 29 |
| 5 | 543572 | 1.6041 | 13 | 4 |
| 6 | 498321 | 1.5967 | 10 | 50 |
| 7 | 3885368 | 1.5826 | 23 | 49 |
| 8 | 4486643 | 1.5759 | 3 | 29 |
| 9 | 2160858 | 1.5730 | 6 | 11 |
| 10 | 259284 | 1.5708 | 6 | 29 |
| 11 | 1608867 | 1.5702 | 8 | 30 |
| 12 | 107021 | 1.5676 | 5 | 26 |

## FSDS feature ranking
| rank | feature | F-score | selected |
|---:|---|---:|---:|
| 1 | `i_credit_last` | 1.007 | 1 |
| 2 | `i_share_last` | 1.007 | 1 |
| 3 | `u_span_sec` | 0.799 | 1 |
| 4 | `u_n_exp` | 0.702 | 1 |
| 5 | `i_share_first` | 0.671 | 1 |
| 6 | `i_credit_first` | 0.671 | 1 |
| 7 | `i_credit_linear` | 0.480 | 1 |
| 8 | `i_share_linear` | 0.480 | 1 |
| 9 | `u_n_uniq_items` | 0.361 | 1 |
| 10 | `i_n_covisit_neighbors` | 0.279 | 1 |
| 11 | `i_n_users` | 0.255 | 1 |
| 12 | `i_n_exp` | 0.209 | 1 |
| 13 | `ui_pop_mismatch` | 0.164 | 1 |
| 14 | `u_log1p_n_uniq` | 0.134 | 1 |
| 15 | `i_log1p_n_covisit` | 0.118 | 1 |

## Holdout
- W1 user-holdout: n=937/298 | ranking only
- W2 temporal: n=937/7223 | hgb AUC=0.295 | logreg AUC=0.189
