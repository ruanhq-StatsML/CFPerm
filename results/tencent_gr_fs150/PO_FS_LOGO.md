# PO-risk feature selection vs F-score (TencentGR 150)

Clock = `_seq_t_end` median. W=later half. Y=`future_cnv`.
RF-domain AUC **0.753**. PO-risk R=E[τ̂²] **0.001782**.

Kendall τ: PO-VIMP vs F-score **0.050**; RF-domain vs F-score **0.115**; PO vs RF-domain **0.381**.

## Family mass

| family | n | RF-domain (P(X)) | PO-VIMP | LOGO Δ | LOGO share |
|---|---:|---:|---:|---:|---:|
| funnel | 53 | 0.372 | 0.288 | +0.00001 | 0.517 |
| session | 4 | 0.066 | 0.070 | +0.00001 | 0.483 |
| attr | 5 | 0.010 | 0.027 | -0.00004 | 0.000 |
| cross | 59 | 0.345 | 0.413 | -0.00012 | 0.000 |
| decay | 7 | 0.049 | 0.134 | -0.00007 | 0.000 |
| diversity | 4 | 0.024 | 0.000 | -0.00003 | 0.000 |
| markov | 7 | 0.092 | 0.066 | -0.00006 | 0.000 |
| money | 11 | 0.043 | 0.002 | -0.00008 | 0.000 |

F-score ranks in-sample Y association. PO-LOGO ranks which family carries the *batch map* φ=(Y-μ)(W-e). Low Kendall ⇒ the 150 F-board is not an OOD board.
