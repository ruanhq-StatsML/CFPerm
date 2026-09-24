# Board reason codes

Reason codes auto-generated from adjacent-board summary. Routing / claim-control only; allows_ship_model=false always. Includes TencentGR advertising funnel scenario when pack present.

catalog=board_rc_v2_ad · n=32

## Advertising funnel scenario (TencentGR)

- family: **灌入族** (`S2_inject`)
- sub: 转化下行/差流 (`ad_convert_dip`)
- sign_Dy(board)=neg · mean_Δȳ=-0.007833333333333333
- read: 相邻窗转化掉：先查差流灌入、落地劣化、账户被打压残留——勿直接判商户/创意「变差」
- ops−content gap: 0.1607264618473112

```
【广告漏斗场景卡 · TencentGR】
族: 灌入族 (`S2_inject`) / 子类: 转化下行/差流 (`ad_convert_dip`)
读法: 相邻窗转化掉：先查差流灌入、落地劣化、账户被打压残留——勿直接判商户/创意「变差」
sign_Dy(board)=neg · mean_Δȳ=-0.007833333333333333 · N=1000
ops−content gap=0.1607264618473112
ops tops: ['e_log1p_exp', 'e_n_exp', 'i_log1p_n_exp']
content tops: ['i_credit_last', 'i_share_last', 'i_item_credit_rank']
动作默认: L1_watch · auto_ban=false · 禁止 HGB 直接上线
强度 = 买量基线；credit/share = 末跳候选 shortlist，非真驱动证书
```

| code | sev | pack@N | title |
|---|---|---|---|
| `RC_BOARD_NOT_SHIP` | block_claim | — | 看板禁止直接上线 HGB |
| `RC_SHORTLIST_NOT_DRIVER` | block_claim | — | shortlist ≠ 真驱动 |
| `RC_AD_LAST_TOUCH_CANDIDATE` | investigate | tencent_gr_content@1000 | 末跳/路径归因候选（落地侧） |
| `RC_CONTENT_CREDIT_CANDIDATE` | investigate | tencent_gr_content@1000 | content 面板 credit/share 候选 |
| `RC_CONTENT_CREDIT_CANDIDATE` | investigate | tencent_gr_content@2000 | content 面板 credit/share 候选 |
| `RC_AD_CONVERT_DIP` | watch | tencent_gr@1000 | 相邻窗转化率下行（灌入/差流倾向） |
| `RC_DIFFDB_TOKEN_NO_TRAVEL` | watch | diffusiondb@1000 | DiffDB prompt token 对 NSFW 不稳传 |
| `RC_DIFFDB_TOKEN_NO_TRAVEL` | watch | diffusiondb@2000 | DiffDB prompt token 对 NSFW 不稳传 |
| `RC_OPS_CONTENT_GAP` | watch | tencent_gr@1000 | ops−content 缺口：强度解释大半 transfer |
| `RC_PROBE_DISAGREE` | watch | metro_interstate@1000 | HGB 与 LogReg transfer 分歧 |
| `RC_PROBE_DISAGREE` | watch | metro_interstate@2000 | HGB 与 LogReg transfer 分歧 |
| `RC_PROBE_DISAGREE` | watch | tencent_gr@1000 | HGB 与 LogReg transfer 分歧 |
| `RC_PROBE_DISAGREE` | watch | tencent_gr@2000 | HGB 与 LogReg transfer 分歧 |
| `RC_SHIFTING_DRIVERS` | watch | metro_interstate@1000 | 高 transfer 但驱动名单在换 |
| `RC_AD_BUY_INTENSITY` | info | tencent_gr@1000 | 投放/曝光强度基线（买量侧） |
| `RC_AD_FUNNEL_CONTEXT` | info | tencent_gr@1000 | 广告漏斗语境：曝→点→转化 |
| `RC_INTENSITY_BASELINE` | info | tencent_gr@1000 | 强度基线（曝光量在传） |
| `RC_INTENSITY_BASELINE` | info | tencent_gr@2000 | 强度基线（曝光量在传） |
| `RC_PERSISTENT_ASSOCIATION` | info | beijing_pm25@1000 | 关联稳传且 top 稳定 |
| `RC_PERSISTENT_ASSOCIATION` | info | beijing_pm25@2000 | 关联稳传且 top 稳定 |
| `RC_PERSISTENT_ASSOCIATION` | info | metro_interstate@2000 | 关联稳传且 top 稳定 |
| `RC_PERSISTENT_ASSOCIATION` | info | tencent_gr@1000 | 关联稳传且 top 稳定 |
| `RC_PERSISTENT_ASSOCIATION` | info | waymo_proxy@1000 | 关联稳传且 top 稳定 |
| `RC_PERSISTENT_ASSOCIATION` | info | waymo_proxy@2000 | 关联稳传且 top 稳定 |
| `RC_PLANTED_SANITY_OK` | info | waymo_proxy@1000 | 种植漂移包探针正常 |
| `RC_PLANTED_SANITY_OK` | info | waymo_proxy@2000 | 种植漂移包探针正常 |
| `RC_RANK_MEAN_SPLIT` | info | diffusiondb@1000 | FSDS 排名与 cmean 漂移不一致 |
| `RC_RANK_MEAN_SPLIT` | info | diffusiondb@2000 | FSDS 排名与 cmean 漂移不一致 |
| `RC_RANK_MEAN_SPLIT` | info | tencent_gr@1000 | FSDS 排名与 cmean 漂移不一致 |
| `RC_RANK_MEAN_SPLIT` | info | tencent_gr@2000 | FSDS 排名与 cmean 漂移不一致 |
| `RC_WEAK_TRANSFER` | info | diffusiondb@1000 | 关联下一段传不过去 |
| `RC_WEAK_TRANSFER` | info | diffusiondb@2000 | 关联下一段传不过去 |

## Paste for agent

```
【board reason-codes】
来源: sample-chunk adjacent transfer board（非定罪 / 非上线门）

—— 广告漏斗场景 ——
【广告漏斗场景卡 · TencentGR】
族: 灌入族 (`S2_inject`) / 子类: 转化下行/差流 (`ad_convert_dip`)
读法: 相邻窗转化掉：先查差流灌入、落地劣化、账户被打压残留——勿直接判商户/创意「变差」
sign_Dy(board)=neg · mean_Δȳ=-0.007833333333333333 · N=1000
ops−content gap=0.1607264618473112
ops tops: ['e_log1p_exp', 'e_n_exp', 'i_log1p_n_exp']
content tops: ['i_credit_last', 'i_share_last', 'i_item_credit_rank']
动作默认: L1_watch · auto_ban=false · 禁止 HGB 直接上线
强度 = 买量基线；credit/share = 末跳候选 shortlist，非真驱动证书

—— 全量 reason-codes ——

- `RC_BOARD_NOT_SHIP` [block_claim]: 本卡来自相邻切窗 transfer 探针，不是线上打分器。禁止因 AUC 高直接上线 HGB；缺 label delay / 校准 / Acc·时延 / 灰发回滚。
- `RC_SHORTLIST_NOT_DRIVER` [block_claim]: 看板只做特征族 shortlist，不做真驱动认定。真驱动需 ablation / PO-risk / 干预实验。
- `RC_AD_LAST_TOUCH_CANDIDATE` [investigate] `tencent_gr_content`@N=1000: 去量后 credit/share 进 top：广告路径末跳/份额 **候选**。对齐审出卡末跳桶；仍 L1，不做自动限投。
    · i_credit_last: 末次归因偏高（末跳嫌疑）
    · i_share_last: 末次份额偏高（末跳嫌疑）
    · i_item_credit_rank: 商户归因秩变动
- `RC_CONTENT_CREDIT_CANDIDATE` [investigate] `tencent_gr_content`@N=1000: Tencent content（已去量）：credit/share 进入 top——作为末跳/路径归因 **候选** 进审出卡，不做真驱动/L2。
    · i_credit_last: 末次归因偏高（末跳嫌疑）
    · i_share_last: 末次份额偏高（末跳嫌疑）
    · i_item_credit_rank: 商户归因秩变动
- `RC_CONTENT_CREDIT_CANDIDATE` [investigate] `tencent_gr_content`@N=2000: Tencent content（已去量）：credit/share 进入 top——作为末跳/路径归因 **候选** 进审出卡，不做真驱动/L2。
    · i_share_last: 末次份额偏高（末跳嫌疑）
    · i_credit_last: 末次归因偏高（末跳嫌疑）
    · i_item_credit_rank: 商户归因秩变动
- `RC_AD_CONVERT_DIP` [watch] `tencent_gr`@N=1000: 相邻窗 mean_Δȳ≈-0.0078（neg）：转化下行倾向。广告读法=灌入/差流/落地劣化候选，先 L1 盯梢。
    · i_credit_last: 末次归因偏高（末跳嫌疑）
    · i_share_last: 末次份额偏高（末跳嫌疑）
    · i_item_credit_rank: 商户归因秩变动
- `RC_DIFFDB_TOKEN_NO_TRAVEL` [watch] `diffusiondb`@N=1000: DiffusionDB：prompt TF-IDF 对 image_nsfw 弱传。风格 token 常进 SelectKBest，但不构成稳定 NSFW 过滤。
    · mucha: 特征 `mucha`
    · alphonse mucha: token `alphonse mucha`
    · alphonse: 特征 `alphonse`
- `RC_DIFFDB_TOKEN_NO_TRAVEL` [watch] `diffusiondb`@N=2000: DiffusionDB：prompt TF-IDF 对 image_nsfw 弱传。风格 token 常进 SelectKBest，但不构成稳定 NSFW 过滤。
    · alphonse mucha: token `alphonse mucha`
    · alphonse: 特征 `alphonse`
    · mucha: 特征 `mucha`
- `RC_OPS_CONTENT_GAP` [watch] `tencent_gr`@N=1000: ops−content gap≈0.16 @N=1000：大半 transfer 是强度；content 残差才是归因候选空间。
    · i_credit_last: 末次归因偏高（末跳嫌疑）
    · i_share_last: 末次份额偏高（末跳嫌疑）
    · i_item_credit_rank: 商户归因秩变动
- `RC_PROBE_DISAGREE` [watch] `metro_interstate`@N=1000: `metro_interstate` @N=1000: HGB=0.97 vs LogReg=0.61 分歧≥0.15。别把非线性探针优势误写成业务必然；对照 content/线性面板。
    · m5: 特征 `m5`
    · m0: 特征 `m0`
    · m14: 特征 `m14`
- `RC_PROBE_DISAGREE` [watch] `metro_interstate`@N=2000: `metro_interstate` @N=2000: HGB=0.98 vs LogReg=0.64 分歧≥0.15。别把非线性探针优势误写成业务必然；对照 content/线性面板。
    · m5: 特征 `m5`
    · m0: 特征 `m0`
    · m6: 特征 `m6`
- `RC_PROBE_DISAGREE` [watch] `tencent_gr`@N=1000: `tencent_gr` @N=1000: HGB=1.00 vs LogReg=0.83 分歧≥0.15。别把非线性探针优势误写成业务必然；对照 content/线性面板。
    · e_log1p_exp: 边曝光强度(log)
    · e_n_exp: 边曝光次数
    · i_log1p_n_exp: 商户曝光规模(log)
- `RC_PROBE_DISAGREE` [watch] `tencent_gr`@N=2000: `tencent_gr` @N=2000: HGB=1.00 vs LogReg=0.66 分歧≥0.15。别把非线性探针优势误写成业务必然；对照 content/线性面板。
    · e_log1p_exp: 边曝光强度(log)
    · e_n_exp: 边曝光次数
    · i_log1p_n_exp: 商户曝光规模(log)
- `RC_SHIFTING_DRIVERS` [watch] `metro_interstate`@N=1000: `metro_interstate` @N=1000: AUC 高但驱动名单 Jaccard≈0.23，预测还在、名单在换——按 regime/composition 读，勿锁死单一 tip。
    · m5: 特征 `m5`
    · m0: 特征 `m0`
    · m14: 特征 `m14`
- `RC_AD_BUY_INTENSITY` [info] `tencent_gr`@N=1000: 买量/曝光强度在相邻窗稳传：这是投放量基线。可解释「量在不在」，不能单独解释「创意/落地好不好」。
    · e_log1p_exp: 边曝光强度(log)
    · e_n_exp: 边曝光次数
    · i_log1p_n_exp: 商户曝光规模(log)
- `RC_AD_FUNNEL_CONTEXT` [info] `tencent_gr`@N=1000: 广告场景：边=曝光/点击流，Y=y_convert，商户键≈广告主。看板按每 N 条边切窗做 transfer 探针，服务审出分流而非出价模型。
    · e_log1p_exp: 边曝光强度(log)
    · e_n_exp: 边曝光次数
- `RC_INTENSITY_BASELINE` [info] `tencent_gr`@N=1000: Tencent ops：transfer 几乎由曝光/点击强度贡献。这是强度基线，不是 content tip；勿写成商户内容问题。
    · e_log1p_exp: 边曝光强度(log)
    · e_n_exp: 边曝光次数
    · i_log1p_n_exp: 商户曝光规模(log)
- `RC_INTENSITY_BASELINE` [info] `tencent_gr`@N=2000: Tencent ops：transfer 几乎由曝光/点击强度贡献。这是强度基线，不是 content tip；勿写成商户内容问题。
    · e_log1p_exp: 边曝光强度(log)
    · e_n_exp: 边曝光次数
    · i_log1p_n_exp: 商户曝光规模(log)
- `RC_PERSISTENT_ASSOCIATION` [info] `beijing_pm25`@N=1000: `beijing_pm25` @N=1000: 高 transfer 且 top Jaccard≈0.56，关联稳传——先当 persistent correlate，勿直接当因果驱动。
    · p8: 特征 `p8`
    · p0: 特征 `p0`
    · p2: 特征 `p2`
- `RC_PERSISTENT_ASSOCIATION` [info] `beijing_pm25`@N=2000: `beijing_pm25` @N=2000: 高 transfer 且 top Jaccard≈0.55，关联稳传——先当 persistent correlate，勿直接当因果驱动。
    · p8: 特征 `p8`
    · p0: 特征 `p0`
    · p2: 特征 `p2`
- `RC_PERSISTENT_ASSOCIATION` [info] `metro_interstate`@N=2000: `metro_interstate` @N=2000: 高 transfer 且 top Jaccard≈0.55，关联稳传——先当 persistent correlate，勿直接当因果驱动。
    · m5: 特征 `m5`
    · m0: 特征 `m0`
    · m6: 特征 `m6`
- `RC_PERSISTENT_ASSOCIATION` [info] `tencent_gr`@N=1000: `tencent_gr` @N=1000: 高 transfer 且 top Jaccard≈0.70，关联稳传——先当 persistent correlate，勿直接当因果驱动。
    · e_log1p_exp: 边曝光强度(log)
    · e_n_exp: 边曝光次数
    · i_log1p_n_exp: 商户曝光规模(log)
- `RC_PERSISTENT_ASSOCIATION` [info] `waymo_proxy`@N=1000: `waymo_proxy` @N=1000: 高 transfer 且 top Jaccard≈0.72，关联稳传——先当 persistent correlate，勿直接当因果驱动。
    · x0: 特征 `x0`
    · x8: 特征 `x8`
    · x3: 特征 `x3`
- `RC_PERSISTENT_ASSOCIATION` [info] `waymo_proxy`@N=2000: `waymo_proxy` @N=2000: 高 transfer 且 top Jaccard≈1.00，关联稳传——先当 persistent correlate，勿直接当因果驱动。
    · x0: 特征 `x0`
    · x8: 特征 `x8`
    · x5: 特征 `x5`
- `RC_PLANTED_SANITY_OK` [info] `waymo_proxy`@N=1000: Waymo proxy 种植漂移上探针高 transfer——看板 sanity 通过；若此处失败，先修探针再读业务包。
    · x0: 特征 `x0`
    · x8: 特征 `x8`
    · x3: 特征 `x3`
- `RC_PLANTED_SANITY_OK` [info] `waymo_proxy`@N=2000: Waymo proxy 种植漂移上探针高 transfer——看板 sanity 通过；若此处失败，先修探针再读业务包。
    · x0: 特征 `x0`
    · x8: 特征 `x8`
    · x5: 特征 `x5`
- `RC_RANK_MEAN_SPLIT` [info] `diffusiondb`@N=1000: `diffusiondb` @N=1000: FSDS∩cmean≈0.00，排名选中的特征和均值漂移特征不是一回事——查法分开写。
    · mucha: 特征 `mucha`
    · alphonse mucha: token `alphonse mucha`
    · alphonse: 特征 `alphonse`
- `RC_RANK_MEAN_SPLIT` [info] `diffusiondb`@N=2000: `diffusiondb` @N=2000: FSDS∩cmean≈0.00，排名选中的特征和均值漂移特征不是一回事——查法分开写。
    · alphonse mucha: token `alphonse mucha`
    · alphonse: 特征 `alphonse`
    · mucha: 特征 `mucha`
- `RC_RANK_MEAN_SPLIT` [info] `tencent_gr`@N=1000: `tencent_gr` @N=1000: FSDS∩cmean≈0.04，排名选中的特征和均值漂移特征不是一回事——查法分开写。
    · e_log1p_exp: 边曝光强度(log)
    · e_n_exp: 边曝光次数
    · i_log1p_n_exp: 商户曝光规模(log)
- `RC_RANK_MEAN_SPLIT` [info] `tencent_gr`@N=2000: `tencent_gr` @N=2000: FSDS∩cmean≈0.04，排名选中的特征和均值漂移特征不是一回事——查法分开写。
    · e_log1p_exp: 边曝光强度(log)
    · e_n_exp: 边曝光次数
    · i_log1p_n_exp: 商户曝光规模(log)
- `RC_WEAK_TRANSFER` [info] `diffusiondb`@N=1000: `diffusiondb` @N=1000: 下一段 AUC≈0.60，关联基本不传。不要把本窗 top 特征当成稳定过滤器。
    · mucha: 特征 `mucha`
    · alphonse mucha: token `alphonse mucha`
    · alphonse: 特征 `alphonse`
- `RC_WEAK_TRANSFER` [info] `diffusiondb`@N=2000: `diffusiondb` @N=2000: 下一段 AUC≈0.61，关联基本不传。不要把本窗 top 特征当成稳定过滤器。
    · alphonse mucha: token `alphonse mucha`
    · alphonse: 特征 `alphonse`
    · mucha: 特征 `mucha`
```


## Detail

### `RC_BOARD_NOT_SHIP` — 看板禁止直接上线 HGB
- severity: **block_claim** · family: ship_gate
- pack: `None` · chunk: `None`
- copy: 本卡来自相邻切窗 transfer 探针，不是线上打分器。禁止因 AUC 高直接上线 HGB；缺 label delay / 校准 / Acc·时延 / 灰发回滚。
- next: 保持 promote_HGB_to_production=false, 若要上线另走 holdout Acc + calibration + gray rollback
- evidence: `{"promote_HGB_to_production": false, "reason": "Adjacent-chunk transfer probe ≠ online scorer. Missing: label delay, calibration under serve skew, Acc/latency budget, abstain/rollback, and causal/PO checks for 'true drivers'."}`

### `RC_SHORTLIST_NOT_DRIVER` — shortlist ≠ 真驱动
- severity: **block_claim** · family: driver_gate
- pack: `None` · chunk: `None`
- copy: 看板只做特征族 shortlist，不做真驱动认定。真驱动需 ablation / PO-risk / 干预实验。
- next: shortlist 进 localize/PO 队列, 禁止把 top_fsds 写成定罪话术
- evidence: `{"note": "transfer shortlist only"}`

### `RC_AD_LAST_TOUCH_CANDIDATE` — 末跳/路径归因候选（落地侧）
- severity: **investigate** · family: ad_funnel
- pack: `tencent_gr_content` · chunk: `1000`
- copy: 去量后 credit/share 进 top：广告路径末跳/份额 **候选**。对齐审出卡末跳桶；仍 L1，不做自动限投。
- next: 映射 tip 桶 i_share_last / i_credit_last, 人审辨末跳操控 vs 正常回收
- evidence: `{"credit_hits": ["i_credit_last", "i_share_last", "i_item_credit_rank", "i_log1p_credit_linear", "i_share_linear"]}`

### `RC_CONTENT_CREDIT_CANDIDATE` — content 面板 credit/share 候选
- severity: **investigate** · family: tencent_content
- pack: `tencent_gr_content` · chunk: `1000`
- copy: Tencent content（已去量）：credit/share 进入 top——作为末跳/路径归因 **候选** 进审出卡，不做真驱动/L2。
- next: 映射 TIP_BUCKETS 末跳/份额话术, 进 localize / PO shortlist, 人审 useful 票后再谈升桶
- evidence: `{"mean_auc": 0.8350687106424594, "mean_logreg_auc": 0.8445150735163063, "mean_jaccard": 0.5714285714285714, "mean_fsds_cmean_jaccard": 0.29298941798941797, "mean_delta_Y": -0.007833333333333333, "reading": "strong transfer + stable top feat`

### `RC_CONTENT_CREDIT_CANDIDATE` — content 面板 credit/share 候选
- severity: **investigate** · family: tencent_content
- pack: `tencent_gr_content` · chunk: `2000`
- copy: Tencent content（已去量）：credit/share 进入 top——作为末跳/路径归因 **候选** 进审出卡，不做真驱动/L2。
- next: 映射 TIP_BUCKETS 末跳/份额话术, 进 localize / PO shortlist, 人审 useful 票后再谈升桶
- evidence: `{"mean_auc": 0.9003182536545639, "mean_logreg_auc": 0.8982588120690723, "mean_jaccard": 0.42857142857142855, "mean_fsds_cmean_jaccard": 0.1574074074074074, "mean_delta_Y": -0.017333333333333336, "reading": "transfer holds; feature set parti`

### `RC_AD_CONVERT_DIP` — 相邻窗转化率下行（灌入/差流倾向）
- severity: **watch** · family: ad_funnel
- pack: `tencent_gr` · chunk: `1000`
- copy: 相邻窗 mean_Δȳ≈-0.0078（neg）：转化下行倾向。广告读法=灌入/差流/落地劣化候选，先 L1 盯梢。
- next: 对 S2_inject 队列, 核对定向变更与打压残留
- evidence: `{"sign_Dy_board": "neg", "mean_delta_Y": -0.007833333333333333}`

### `RC_DIFFDB_TOKEN_NO_TRAVEL` — DiffDB prompt token 对 NSFW 不稳传
- severity: **watch** · family: diffusiondb
- pack: `diffusiondb` · chunk: `1000`
- copy: DiffusionDB：prompt TF-IDF 对 image_nsfw 弱传。风格 token 常进 SelectKBest，但不构成稳定 NSFW 过滤。
- next: NSFW 故事不要押在单一风格 token, 可选：concat CFG/step 或换 embedding 后再跑探针
- evidence: `{"mean_auc": 0.6048427572700964, "mean_logreg_auc": 0.6161109153514835, "mean_jaccard": 0.2400793650793651, "mean_fsds_cmean_jaccard": 0.0, "mean_delta_Y": 0.014843780122963442, "reading": "weak transfer (association does not travel)", "top`

### `RC_DIFFDB_TOKEN_NO_TRAVEL` — DiffDB prompt token 对 NSFW 不稳传
- severity: **watch** · family: diffusiondb
- pack: `diffusiondb` · chunk: `2000`
- copy: DiffusionDB：prompt TF-IDF 对 image_nsfw 弱传。风格 token 常进 SelectKBest，但不构成稳定 NSFW 过滤。
- next: NSFW 故事不要押在单一风格 token, 可选：concat CFG/step 或换 embedding 后再跑探针
- evidence: `{"mean_auc": 0.6091376939024363, "mean_logreg_auc": 0.6187648813337822, "mean_jaccard": 0.3392857142857143, "mean_fsds_cmean_jaccard": 0.0, "mean_delta_Y": 0.029074908243802683, "reading": "weak transfer (association does not travel)", "top`

### `RC_OPS_CONTENT_GAP` — ops−content 缺口：强度解释大半 transfer
- severity: **watch** · family: tencent_gap
- pack: `tencent_gr` · chunk: `1000`
- copy: ops−content gap≈0.16 @N=1000：大半 transfer 是强度；content 残差才是归因候选空间。
- next: 汇报时先报 gap 再报 content tops, 强度基线与 content 候选分两行写进工单
- evidence: `{"chunk": 1000, "ops_auc": 0.9957951724897706, "content_auc": 0.8350687106424594, "auc_gap_ops_minus_content": 0.1607264618473112, "ops_top0": ["e_log1p_exp", "e_n_exp", "i_log1p_n_exp"], "content_top0": ["i_credit_last", "i_share_last", "i`

### `RC_PROBE_DISAGREE` — HGB 与 LogReg transfer 分歧
- severity: **watch** · family: probe
- pack: `metro_interstate` · chunk: `1000`
- copy: `metro_interstate` @N=1000: HGB=0.97 vs LogReg=0.61 分歧≥0.15。别把非线性探针优势误写成业务必然；对照 content/线性面板。
- next: 并列报告双探针, 优先信任两探针都同意的特征族
- evidence: `{"mean_auc": 0.9718807034328297, "mean_logreg_auc": 0.6090001912223876, "mean_jaccard": 0.2334656084656085, "mean_fsds_cmean_jaccard": 0.2556689342403628, "mean_delta_Y": -20.366857142857175, "reading": "strong transfer but shifting drivers`

### `RC_PROBE_DISAGREE` — HGB 与 LogReg transfer 分歧
- severity: **watch** · family: probe
- pack: `metro_interstate` · chunk: `2000`
- copy: `metro_interstate` @N=2000: HGB=0.98 vs LogReg=0.64 分歧≥0.15。别把非线性探针优势误写成业务必然；对照 content/线性面板。
- next: 并列报告双探针, 优先信任两探针都同意的特征族
- evidence: `{"mean_auc": 0.9753765691624277, "mean_logreg_auc": 0.6396613442498914, "mean_jaccard": 0.5476190476190476, "mean_fsds_cmean_jaccard": 0.30952380952380953, "mean_delta_Y": 15.445333333333414, "reading": "strong transfer + stable top feats (`

### `RC_PROBE_DISAGREE` — HGB 与 LogReg transfer 分歧
- severity: **watch** · family: probe
- pack: `tencent_gr` · chunk: `1000`
- copy: `tencent_gr` @N=1000: HGB=1.00 vs LogReg=0.83 分歧≥0.15。别把非线性探针优势误写成业务必然；对照 content/线性面板。
- next: 并列报告双探针, 优先信任两探针都同意的特征族
- evidence: `{"mean_auc": 0.9957951724897706, "mean_logreg_auc": 0.8318962492242159, "mean_jaccard": 0.7047619047619047, "mean_fsds_cmean_jaccard": 0.037037037037037035, "mean_delta_Y": -0.007833333333333333, "reading": "strong transfer + stable top fea`

### `RC_PROBE_DISAGREE` — HGB 与 LogReg transfer 分歧
- severity: **watch** · family: probe
- pack: `tencent_gr` · chunk: `2000`
- copy: `tencent_gr` @N=2000: HGB=1.00 vs LogReg=0.66 分歧≥0.15。别把非线性探针优势误写成业务必然；对照 content/线性面板。
- next: 并列报告双探针, 优先信任两探针都同意的特征族
- evidence: `{"mean_auc": 0.9987730084570958, "mean_logreg_auc": 0.6634610688713213, "mean_jaccard": 0.4583333333333333, "mean_fsds_cmean_jaccard": 0.037037037037037035, "mean_delta_Y": -0.017333333333333336, "reading": "transfer holds; feature set part`

### `RC_SHIFTING_DRIVERS` — 高 transfer 但驱动名单在换
- severity: **watch** · family: stability
- pack: `metro_interstate` · chunk: `1000`
- copy: `metro_interstate` @N=1000: AUC 高但驱动名单 Jaccard≈0.23，预测还在、名单在换——按 regime/composition 读，勿锁死单一 tip。
- next: 对比 N=1000 vs 2000 粒度, 人审看当窗 top，不复用上窗话术
- evidence: `{"mean_auc": 0.9718807034328297, "mean_logreg_auc": 0.6090001912223876, "mean_jaccard": 0.2334656084656085, "mean_fsds_cmean_jaccard": 0.2556689342403628, "mean_delta_Y": -20.366857142857175, "reading": "strong transfer but shifting drivers`

### `RC_AD_BUY_INTENSITY` — 投放/曝光强度基线（买量侧）
- severity: **info** · family: ad_funnel
- pack: `tencent_gr` · chunk: `1000`
- copy: 买量/曝光强度在相邻窗稳传：这是投放量基线。可解释「量在不在」，不能单独解释「创意/落地好不好」。
- next: ops 行只报曝光强度, 创意问题看 content 面板
- evidence: `{"tops": ["e_log1p_exp", "e_n_exp", "i_log1p_n_exp", "i_credit_last", "i_share_last"]}`

### `RC_AD_FUNNEL_CONTEXT` — 广告漏斗语境：曝→点→转化
- severity: **info** · family: ad_funnel
- pack: `tencent_gr` · chunk: `1000`
- copy: 广告场景：边=曝光/点击流，Y=y_convert，商户键≈广告主。看板按每 N 条边切窗做 transfer 探针，服务审出分流而非出价模型。
- next: 用工单话术读 S1/S2/S3，不把 AUC 当 CTR 提升证明, 强度与末跳分两行写
- evidence: `{"funnel": "exp→clk→y_convert", "mean_delta_Y": -0.007833333333333333}`

### `RC_INTENSITY_BASELINE` — 强度基线（曝光量在传）
- severity: **info** · family: tencent_ops
- pack: `tencent_gr` · chunk: `1000`
- copy: Tencent ops：transfer 几乎由曝光/点击强度贡献。这是强度基线，不是 content tip；勿写成商户内容问题。
- next: 对照 tencent_gr_content 面板看缺口, 强度特征只做 ops 分流，不进定罪桶
- evidence: `{"mean_auc": 0.9957951724897706, "mean_logreg_auc": 0.8318962492242159, "mean_jaccard": 0.7047619047619047, "mean_fsds_cmean_jaccard": 0.037037037037037035, "mean_delta_Y": -0.007833333333333333, "reading": "strong transfer + stable top fea`

### `RC_INTENSITY_BASELINE` — 强度基线（曝光量在传）
- severity: **info** · family: tencent_ops
- pack: `tencent_gr` · chunk: `2000`
- copy: Tencent ops：transfer 几乎由曝光/点击强度贡献。这是强度基线，不是 content tip；勿写成商户内容问题。
- next: 对照 tencent_gr_content 面板看缺口, 强度特征只做 ops 分流，不进定罪桶
- evidence: `{"mean_auc": 0.9987730084570958, "mean_logreg_auc": 0.6634610688713213, "mean_jaccard": 0.4583333333333333, "mean_fsds_cmean_jaccard": 0.037037037037037035, "mean_delta_Y": -0.017333333333333336, "reading": "transfer holds; feature set part`

### `RC_PERSISTENT_ASSOCIATION` — 关联稳传且 top 稳定
- severity: **info** · family: stability
- pack: `beijing_pm25` · chunk: `1000`
- copy: `beijing_pm25` @N=1000: 高 transfer 且 top Jaccard≈0.56，关联稳传——先当 persistent correlate，勿直接当因果驱动。
- next: 记入 shortlist, 上 PO/ablation 前保持 L1
- evidence: `{"mean_auc": 0.9780904023588342, "mean_logreg_auc": 0.978103102615826, "mean_jaccard": 0.5634920634920634, "mean_fsds_cmean_jaccard": 0.5391156462585033, "mean_delta_Y": 7.190857142857142, "reading": "strong transfer + stable top feats (per`

### `RC_PERSISTENT_ASSOCIATION` — 关联稳传且 top 稳定
- severity: **info** · family: stability
- pack: `beijing_pm25` · chunk: `2000`
- copy: `beijing_pm25` @N=2000: 高 transfer 且 top Jaccard≈0.55，关联稳传——先当 persistent correlate，勿直接当因果驱动。
- next: 记入 shortlist, 上 PO/ablation 前保持 L1
- evidence: `{"mean_auc": 0.9797058426155137, "mean_logreg_auc": 0.97320706540759, "mean_jaccard": 0.5476190476190476, "mean_fsds_cmean_jaccard": 0.42857142857142855, "mean_delta_Y": 9.979999999999999, "reading": "strong transfer + stable top feats (per`

### `RC_PERSISTENT_ASSOCIATION` — 关联稳传且 top 稳定
- severity: **info** · family: stability
- pack: `metro_interstate` · chunk: `2000`
- copy: `metro_interstate` @N=2000: 高 transfer 且 top Jaccard≈0.55，关联稳传——先当 persistent correlate，勿直接当因果驱动。
- next: 记入 shortlist, 上 PO/ablation 前保持 L1
- evidence: `{"mean_auc": 0.9753765691624277, "mean_logreg_auc": 0.6396613442498914, "mean_jaccard": 0.5476190476190476, "mean_fsds_cmean_jaccard": 0.30952380952380953, "mean_delta_Y": 15.445333333333414, "reading": "strong transfer + stable top feats (`

### `RC_PERSISTENT_ASSOCIATION` — 关联稳传且 top 稳定
- severity: **info** · family: stability
- pack: `tencent_gr` · chunk: `1000`
- copy: `tencent_gr` @N=1000: 高 transfer 且 top Jaccard≈0.70，关联稳传——先当 persistent correlate，勿直接当因果驱动。
- next: 记入 shortlist, 上 PO/ablation 前保持 L1
- evidence: `{"mean_auc": 0.9957951724897706, "mean_logreg_auc": 0.8318962492242159, "mean_jaccard": 0.7047619047619047, "mean_fsds_cmean_jaccard": 0.037037037037037035, "mean_delta_Y": -0.007833333333333333, "reading": "strong transfer + stable top fea`

### `RC_PERSISTENT_ASSOCIATION` — 关联稳传且 top 稳定
- severity: **info** · family: stability
- pack: `waymo_proxy` · chunk: `1000`
- copy: `waymo_proxy` @N=1000: 高 transfer 且 top Jaccard≈0.72，关联稳传——先当 persistent correlate，勿直接当因果驱动。
- next: 记入 shortlist, 上 PO/ablation 前保持 L1
- evidence: `{"mean_auc": 0.9399096272322387, "mean_logreg_auc": 0.9457475548663155, "mean_jaccard": 0.7222222222222222, "mean_fsds_cmean_jaccard": 0.6802721088435374, "mean_delta_Y": 0.2307242796562466, "reading": "strong transfer + stable top feats (p`

### `RC_PERSISTENT_ASSOCIATION` — 关联稳传且 top 稳定
- severity: **info** · family: stability
- pack: `waymo_proxy` · chunk: `2000`
- copy: `waymo_proxy` @N=2000: 高 transfer 且 top Jaccard≈1.00，关联稳传——先当 persistent correlate，勿直接当因果驱动。
- next: 记入 shortlist, 上 PO/ablation 前保持 L1
- evidence: `{"mean_auc": 0.9692848500223602, "mean_logreg_auc": 0.9776463157951998, "mean_jaccard": 1.0, "mean_fsds_cmean_jaccard": 1.0, "mean_delta_Y": 0.46176069549992277, "reading": "strong transfer + stable top feats (persistent drivers)", "top_fsd`

### `RC_PLANTED_SANITY_OK` — 种植漂移包探针正常
- severity: **info** · family: sanity
- pack: `waymo_proxy` · chunk: `1000`
- copy: Waymo proxy 种植漂移上探针高 transfer——看板 sanity 通过；若此处失败，先修探针再读业务包。
- next: 保留为回归对照包
- evidence: `{"mean_auc": 0.9399096272322387, "mean_logreg_auc": 0.9457475548663155, "mean_jaccard": 0.7222222222222222, "mean_fsds_cmean_jaccard": 0.6802721088435374, "mean_delta_Y": 0.2307242796562466, "reading": "strong transfer + stable top feats (p`

### `RC_PLANTED_SANITY_OK` — 种植漂移包探针正常
- severity: **info** · family: sanity
- pack: `waymo_proxy` · chunk: `2000`
- copy: Waymo proxy 种植漂移上探针高 transfer——看板 sanity 通过；若此处失败，先修探针再读业务包。
- next: 保留为回归对照包
- evidence: `{"mean_auc": 0.9692848500223602, "mean_logreg_auc": 0.9776463157951998, "mean_jaccard": 1.0, "mean_fsds_cmean_jaccard": 1.0, "mean_delta_Y": 0.46176069549992277, "reading": "strong transfer + stable top feats (persistent drivers)", "top_fsd`

### `RC_RANK_MEAN_SPLIT` — FSDS 排名与 cmean 漂移不一致
- severity: **info** · family: selection
- pack: `diffusiondb` · chunk: `1000`
- copy: `diffusiondb` @N=1000: FSDS∩cmean≈0.00，排名选中的特征和均值漂移特征不是一回事——查法分开写。
- next: 卡上同时贴 top_fsds 与 top_cmean, 禁止只用 cmean 讲「预测力」
- evidence: `{"mean_auc": 0.6048427572700964, "mean_logreg_auc": 0.6161109153514835, "mean_jaccard": 0.2400793650793651, "mean_fsds_cmean_jaccard": 0.0, "mean_delta_Y": 0.014843780122963442, "reading": "weak transfer (association does not travel)", "top`

### `RC_RANK_MEAN_SPLIT` — FSDS 排名与 cmean 漂移不一致
- severity: **info** · family: selection
- pack: `diffusiondb` · chunk: `2000`
- copy: `diffusiondb` @N=2000: FSDS∩cmean≈0.00，排名选中的特征和均值漂移特征不是一回事——查法分开写。
- next: 卡上同时贴 top_fsds 与 top_cmean, 禁止只用 cmean 讲「预测力」
- evidence: `{"mean_auc": 0.6091376939024363, "mean_logreg_auc": 0.6187648813337822, "mean_jaccard": 0.3392857142857143, "mean_fsds_cmean_jaccard": 0.0, "mean_delta_Y": 0.029074908243802683, "reading": "weak transfer (association does not travel)", "top`

### `RC_RANK_MEAN_SPLIT` — FSDS 排名与 cmean 漂移不一致
- severity: **info** · family: selection
- pack: `tencent_gr` · chunk: `1000`
- copy: `tencent_gr` @N=1000: FSDS∩cmean≈0.04，排名选中的特征和均值漂移特征不是一回事——查法分开写。
- next: 卡上同时贴 top_fsds 与 top_cmean, 禁止只用 cmean 讲「预测力」
- evidence: `{"mean_auc": 0.9957951724897706, "mean_logreg_auc": 0.8318962492242159, "mean_jaccard": 0.7047619047619047, "mean_fsds_cmean_jaccard": 0.037037037037037035, "mean_delta_Y": -0.007833333333333333, "reading": "strong transfer + stable top fea`

### `RC_RANK_MEAN_SPLIT` — FSDS 排名与 cmean 漂移不一致
- severity: **info** · family: selection
- pack: `tencent_gr` · chunk: `2000`
- copy: `tencent_gr` @N=2000: FSDS∩cmean≈0.04，排名选中的特征和均值漂移特征不是一回事——查法分开写。
- next: 卡上同时贴 top_fsds 与 top_cmean, 禁止只用 cmean 讲「预测力」
- evidence: `{"mean_auc": 0.9987730084570958, "mean_logreg_auc": 0.6634610688713213, "mean_jaccard": 0.4583333333333333, "mean_fsds_cmean_jaccard": 0.037037037037037035, "mean_delta_Y": -0.017333333333333336, "reading": "transfer holds; feature set part`

### `RC_WEAK_TRANSFER` — 关联下一段传不过去
- severity: **info** · family: transfer
- pack: `diffusiondb` · chunk: `1000`
- copy: `diffusiondb` @N=1000: 下一段 AUC≈0.60，关联基本不传。不要把本窗 top 特征当成稳定过滤器。
- next: 降权该特征族的 tip 叙事, 若业务仍关心 Y，换 X（如 DiffDB 加超参/CLIP）再探针
- evidence: `{"mean_auc": 0.6048427572700964, "mean_logreg_auc": 0.6161109153514835, "mean_jaccard": 0.2400793650793651, "mean_fsds_cmean_jaccard": 0.0, "mean_delta_Y": 0.014843780122963442, "reading": "weak transfer (association does not travel)", "top`

### `RC_WEAK_TRANSFER` — 关联下一段传不过去
- severity: **info** · family: transfer
- pack: `diffusiondb` · chunk: `2000`
- copy: `diffusiondb` @N=2000: 下一段 AUC≈0.61，关联基本不传。不要把本窗 top 特征当成稳定过滤器。
- next: 降权该特征族的 tip 叙事, 若业务仍关心 Y，换 X（如 DiffDB 加超参/CLIP）再探针
- evidence: `{"mean_auc": 0.6091376939024363, "mean_logreg_auc": 0.6187648813337822, "mean_jaccard": 0.3392857142857143, "mean_fsds_cmean_jaccard": 0.0, "mean_delta_Y": 0.029074908243802683, "reading": "weak transfer (association does not travel)", "top`
