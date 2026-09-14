# Hard-sample ranking: reference PO vs re-fit PO-learner

Primary claim = **PO ranks hard OOD rows**. Downstream MSE is secondary.

On each OnlineRFPerm reject batch:

1. `truth_i = |Y_i − μ_oracle(X_i)|` with μ_oracle fit on recent∪current (diagnostic).
2. Score `ref_po` / `probe_po` / `refit_po` on the same rows.
3. Measure ranking quality (Spearman, Precision@20%, Lift, NDCG, AUROC).

## Hard-row ranking (reject batches only) — primary

| dataset | spearman ref→probe→**refit** | P@20% ref→probe→**refit** | AUROC ref→probe→**refit** | Lift@20% **refit** | NDCG **refit** |
|---|---|---|---|---:|---:|
| `metro_interstate` | 0.30→0.32→**0.32** | 0.29→0.36→**0.35** | 0.64→0.67→**0.67** | 1.75 | 0.54 |
| `beijing_pm25` | 0.08→0.39→**0.40** | 0.35→0.52→**0.47** | 0.59→0.77→**0.77** | 2.35 | 0.74 |
| `stocks_AAPL` | 0.29→0.67→**0.67** | 0.52→0.71→**0.69** | 0.72→0.93→**0.92** | 3.46 | 0.90 |
| `waymo_proxy` | 0.11→0.55→**0.53** | 0.26→0.60→**0.57** | 0.57→0.84→**0.83** | 2.83 | 0.79 |
| `stocks_MSFT` | 0.38→0.75→**0.72** | 0.49→0.70→**0.63** | 0.75→0.92→**0.90** | 3.17 | 0.89 |
| `stocks_IWM` | 0.61→0.73→**0.69** | 0.65→0.71→**0.71** | 0.86→0.92→**0.91** | 3.57 | 0.90 |

**Best Spearman wins:** ref=0, probe=4, **refit=2**
**Best Precision@20% wins:** ref=0, probe=5, **refit=1**

### Metric definitions (top 20% = hard)

- **Spearman**: `corr(rank(PO), rank(truth))` — full order concordance.
- **Precision@k**: `|Top_k(PO) ∩ Top_k(truth)| / k` — hard-set recovery.
- **Lift@k**: Precision@k / (k/n) — vs random (1.0 = chance).
- **NDCG@k**: graded by truth hardness — rewards ordering the *hardest* first.
- **AUROC**: truth top-20% as positive class — threshold-free hard detection.

## Downstream next-MSE (sig-only) — secondary

| dataset | n_sig | unif | ref_po | probe_po | **refit_po** | dre | best |
|---|---:|---:|---:|---:|---:|---:|---|
| `metro_interstate` | 6 | 8.602e+05 | 7.925e+05 | 7.839e+05 | **7.57e+05** | 1.642e+06 | `refit_po` |
| `beijing_pm25` | 5 | 2364 | 2261 | 2639 | **2678** | 3577 | `ref_po` |
| `stocks_AAPL` | 6 | 0.0007784 | 0.0008489 | 0.0008779 | **0.0009023** | 0.0008401 | `uniform` |
| `waymo_proxy` | 22 | 0.009992 | 0.01039 | 0.01112 | **0.01078** | 0.01154 | `uniform` |
| `stocks_MSFT` | 8 | 0.0006164 | 0.0006494 | 0.0007019 | **0.0007034** | 0.0005892 | `dre` |
| `stocks_IWM` | 7 | 0.0005025 | 0.0005518 | 0.0005562 | **0.0005362** | 0.0005117 | `uniform` |

**MSE wins (sig-only):** `uniform`=3, `ref_po`=1, `probe_po`=0, `refit_po`=1, `dre`=1
**refit_po < ref_po (MSE):** `2/6`
**refit_po < probe_po (MSE):** `3/6`

See `docs/agod/AGOD_hard_rank_eval.md` for the full protocol.
