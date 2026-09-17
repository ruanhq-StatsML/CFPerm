# Justify: 特征库 — 哪几列在走、落在哪些商户

全局 PO × MSE × MMD 看板管冻层。特征库管定位：目录九列（订单 / 商户 / 用户），Y 不当特征。图 1/2 = 各列 shift 程度；图 3 = 商户排名，**不是** feature 程度。`merchant_id` 是外键，不是图方法。

n_ref=480, n_new=160, merchants=16, onset=1. 种下南区。定位，不是唯一分解。

| kind | 种下的列 | FSDS 选中（最后一批） | loud 订单 south_frac | J(商户, south) | J(用户, south) | 指纹 |
|---|---|---|---:|---:|---:|---|
| covariate_south | amount, channel, merchant_gmv | amount, channel, merchant_gmv | 0.85 | 0.77 | 0 | x_shift |
| concept_south | amount（\(Y\mid X\)） | amount, channel | 0.52 | 0.46 | 0 | 弱；MMD 安静 |
| both | amount, channel, merchant_gmv | amount, channel, merchant_gmv | 0.85 | 0.77 | 0 | x_shift 盖过 PO |

通过：无 Y 泄漏；covariate 噪声列未选；onset 后种下列抬升；用户粒 Jaccard=0（用户跨商户随机挂，不该回收）。

接到冻层：库响、全局 keep → 子集先报。南区金额/渠道 MMD loud → 不要按 concept 冻底。
