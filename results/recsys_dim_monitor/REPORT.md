# Recsys grain monitor — conclusions in OnlineRFPerm / CFPerm form

Same two questions as the papers. **OnlineRFPerm** (Sep 15 Algorithm 1): when did predictive error leave the reference pool? First rejection, delay, FAR. **CFPerm / RFPerm VIMP** (MetaLearner Algorithm 1): which columns contribute, recovered by Kendall-τ against the planted magnitude amount ≻ merchant_gmv ≻ channel. **FSDS + south/north** is the post-hoc localization on this table — not a graph cut.

Data is the synthetic recommendation-shaped order stream: user, merchant, order, conversion Y. Slice X by grain (order / merchant / user / all). Y is never a feature. T is the batch label. Shares are localization, not Shapley or CATE.

n_ref=400, n_new=120, batches=6, labeled onset=2. Planted accounts = south merchants. Planted columns: covariate/both = amount, channel, merchant_gmv; concept = amount in Y|X.

## Dataset (paper §5 style)

Each incoming batch is 120 orders. Frozen RF is fit once on D_ref (400 orders). After onset, only **south** merchants are shifted. User features never walk. Concatenating every grain into `all` is the anti-pattern: more columns, optimistic in-sample E_ref, rank-p FAR at t=0.

## WHEN (OnlineRFPerm Algorithm 1)

Frozen RF on D_ref. Each batch T_t = MSE_t − E_ref. Hop = last-two 1.5× (jump detector). Rank-p vs a ref-null pool; ADDIS primary, SAFFRON contrast. Quiet T≤0 is fed as p=1. Delay = first rejection − onset_true.

| kind | grain | ADDIS | rank-p | hop | last T | last MMD |
|---|---|---|---|---|---|---|
| covariate_south | order | 2 (d=0) | 2 (d=0) | miss | 0.067 | 0.0952 |
| covariate_south | user | 2 (d=0) | FAR@0 | miss | 0.0569 | 0.000627 |
| covariate_south | all | FAR@0 | FAR@0 | miss | 0.0531 | 0.0773 |
| concept_south | order | 2 (d=0) | 2 (d=0) | miss | 0.113 | -0.000532 |
| both | order | 2 (d=0) | 2 (d=0) | miss | 0.157 | 0.0952 |

Findings (OnlineRFPerm §5.2 style):

1. Order grain + ADDIS hits at the labeled onset on covariate (delay=0), concept, and both. That is the slice that actually moved.
2. Last-two hop never fires on this six-batch walk. Gradual south-only amount/channel walk is a trend in T / MMD, not a 1.5× jump. Rank-p and ADDIS are the marks, matching the paper: hop is a jump detector.
3. User grain MMD stays quiet (0.000627). ADDIS can still mark t=2 because Y moved (amount is in the logit) while user X did not — performance-relevant shift with quiet X, the concept fingerprint on the wrong grain.
4. Concatenated `all` grain rank-p / ADDIS fire at t=0 (FAR). Nine-column in-sample E_ref is optimistic. Slice the log the way it is written.
5. Concept order-grain MMD is quiet while T rises — Y|X moved, P(X) did not.

Figure: `online_rfperm_T.png`. Sequential T_t / p_t: Table 4 in `TABLES.md`.

## WHICH columns (RFPerm / CFPerm / FSDS, MetaLearner Algorithm 1)

Nuisances on φ=(Y−μ)(T−e) fit once. CFPerm VIMP = extra MSE of predicting φ after permuting a column. RFPerm ΔMSE = extra MSE of the frozen f_ref after permuting a column of the new batch. FSDS is the univariate MMD / CMean catalog. Kendall-τ vs planted magnitude is the ranking metric in the MetaLearner paper. CFPerm global reject = max VIMP vs 95% of T-permuted nulls (B=12, not 500).

| kind | grain | FSDS recovered | τ_FSDS | RFPerm ΔMSE top | τ_RFPerm | CFPerm φ top | τ_CFPerm | reject |
|---|---|---|---|---|---|---|---|---|
| covariate | order | amount,channel | 0.913 | amount,n_items,channel | 0.548 | n_items,channel,amount | -0.183 |  |
| covariate | all | amount,channel,merchant_gmv | 0.691 | amount,merchant_cat,user_hist_freq | 0.182 | n_items,channel,user_hist_freq | -0.182 |  |
| concept | order | amount | 0.236 | channel,amount,hour | 0.236 | amount,hour,channel | 0.707 |  |
| both | all | amount,channel,merchant_gmv | 0.691 | merchant_cat,merchant_gmv,channel | 0.0364 | channel,user_hist_freq,amount | 0.255 |  |
| covariate | user | — |  | user_tenure,user_hist_freq |  | user_hist_freq,user_tenure |  |  |

Findings (MetaLearner ranking style):

1. FSDS on the native grains recovers the planted columns: covariate order → amount, channel; covariate all → amount, channel, merchant_gmv; concept order → amount (amount via CMean_Y, not MMD).
2. User grain recovers nothing of amount/channel/gmv — those columns are not in that slice. A quiet user catalog is the correct negative control.
3. CFPerm global reject at B=12 does not fire. Ranking, not the max-vs-null test, is the readout here (paper uses B=500).
4. φ-VIMP can put a noise column (n_items) first on covariate; FSDS and RFPerm ΔMSE are the methods that track the planted X-walk. Concept is the reverse: amount leads φ-VIMP because Y|X moved.

Figures: `rfperm_vimp.png` (red = planted), `fsds_rank.png`.

## WHICH accounts (post-hoc FSDS localization)

Subset key = region (south / north), not Y. Own-ref clock. Three readouts together: MMD (P(X)), CMean (||ΔE[X]|| and ΔE[Y]), PO-risk (P(Y|X)).

Covariate order-grain south: MMD=0.35, CMean_X=2.55, PO=7.59e-05, CMean_Y=0.349. North MMD=0.00109.
Covariate user-grain south MMD=-0.00133 — users mix across merchants.
Concept order-grain south: MMD=-0.00137, PO=0.00307, ΔE[Y]=-0.0412.

## What to tell a production recsys

1. Slice the serving table the way the log is written (order / merchant / user). Do not dump every id embedding into one simplex — `all` is the FAR grain.
2. OnlineRFPerm + ADDIS on the slice that actually moved. A quiet user MMD does not mean the order grain is quiet, and a user-grain T mark can just be Y walking through another grain.
3. After a mark, FSDS names the columns (Kendall-τ); RFPerm ΔMSE is the frozen-model companion; CFPerm φ-VIMP is for Y|X. South/north (or any frozen account key) names the accounts.
4. Reset the error pool after a refresh (paper Appendix A.2).

Fine tables: `TABLES.md`. One-pager: `JUSTIFY.md`.
