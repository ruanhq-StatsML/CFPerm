# FSDS time split (TencentGR 150)

Clock = `_seq_t_end`. No shuffle.
Pos rate by quartile: {'Q1': 0.31066666666666665, 'Q2': 0.27666666666666667, 'Q3': 0.04, 'Q4': 0.0}

| clock | n | pos | RF-domain AUC | PO-risk | τ(PO, F-score) |
|---|---:|---:|---:|---:|---:|
| median | 6000 | 0.157 | 0.753 | 0.00178 | 0.051 |
| q1q4 | 3000 | 0.155 | 0.881 | 0.00110 | 0.082 |

Time holdout HGB (by quartile, not a dumb last-25% — Q4 has no positives):

| split | set | AUC | AP | Acc | pos_tr | pos_te |
|---|---|---:|---:|---:|---:|---:|
| Q1→Q2 | f150 | 0.7252 | 0.5140 | 0.7487 | 0.311 | 0.277 |
| Q1→Q2 | po32 | 0.7265 | 0.5299 | 0.7447 | 0.311 | 0.277 |
| Q1→Q2 | rfdomain32 | 0.7148 | 0.5042 | 0.7447 | 0.311 | 0.277 |
| Q12→Q3 | f150 | 0.7274 | 0.1187 | 0.9153 | 0.294 | 0.040 |
| Q12→Q3 | po32 | 0.7164 | 0.1289 | 0.9140 | 0.294 | 0.040 |
| Q12→Q3 | rfdomain32 | 0.7127 | 0.1336 | 0.9007 | 0.294 | 0.040 |
| Q12→Q4 | f150 | nan | 0.0000 | 1.0000 | 0.294 | 0.000 |
| Q12→Q4 | po32 | nan | 0.0000 | 1.0000 | 0.294 | 0.000 |
| Q12→Q4 | rfdomain32 | nan | 0.0000 | 1.0000 | 0.294 | 0.000 |

RF-domain = covariate / P(X). PO-risk = φ=(Y-μ)(W-e).
Q4 pos=0 is right-censoring on this clock, not a concept hop.
Q1→Q2 is the comparable time OOS. Q12→Q3 is the label-rate drop.
τ(PO, F) stays ~0: F board is not FSDS.
