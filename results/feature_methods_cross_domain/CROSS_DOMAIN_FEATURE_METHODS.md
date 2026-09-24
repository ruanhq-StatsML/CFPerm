# 跨域特征方法面板 · smoke justify

同一套特征统计顺序（cmean → cov/PO-VIMP → FSDS → 共识/张力）可从电商图边迁到 prompt→image 与 tabular stream，无需改 Drill。

## Protocol（固定顺序）

- 1. shift tip: |μ_cur−μ_ref| (cmean)
- 2. shift importance: PO-VIMP or X→T cov VIMP
- 3. supervised tip: FSDS / SelectKBest F on y|support
- 4. optional: MMD-LOCO when pairwise feats exist
- 5. deliver consensus + FSDS-only + shift-only tension

## Domains

### tencent_gr_w1w2 (`ecommerce_graph_edges`)
- X/Y/T: user/item/path graph feats on K* → click/convert on localized edges · W1 early vs W2 late (gap≥30d)
- methods: cmean_abs, mmd_loco, po_vimp, fsds_f
- available: {'cmean': True, 'mmd_loco': True, 'po_vimp': True, 'fsds': True}
- consensus: i_credit_linear, i_share_linear, i_n_covisit_neighbors, i_n_exp, i_log1p_n_exp
- FSDS-only: i_share_first
- shift-only: e_log1p_exp, i_n_exp, i_n_users, i_share_last
- read: 多法共识(≥3): i_credit_linear, i_log1p_n_exp, i_n_covisit_neighbors, i_share_linear；仅FSDS(y判别、未必shift): i_share_first；仅shift描述(cmean/MMD/PO): e_log1p_exp, i_n_exp, i_n_users, i_share_last；仅PO/cov-VIMP: u_log1p_n_events, u_span_sec, ui_pop_mismatch

### diffusiondb_temporal (`prompt_tokens_to_image_nsfw`)
- X/Y/T: TF-IDF prompt tokens → image_nsfw · early vs late timestamp windows
- methods: abs_delta(cmean), vimp_cov(X→T), fsds_f
- available: {'cmean': True, 'mmd_loco': False, 'po_vimp': True, 'fsds': True}
- consensus: intricate, sharp, focus, greg, sharp focus
- FSDS-only: alphonse, alphonse mucha, artgerm, body
- shift-only: 3d, 4k, 8k, by
- read: 仅FSDS(y判别、未必shift): alphonse, alphonse mucha, artgerm, body；仅shift描述(cmean/MMD/PO): 3d, 4k, 8k, by；仅PO/cov-VIMP: intricate, mm

### waymo_proxy_stream (`tabular_stream_xy`)
- X/Y/T: proxy features d=9 → proxy label · index early/late half
- methods: cmean_abs, vimp_cov(|corr X,T|), fsds_f
- available: {'cmean': True, 'mmd_loco': False, 'po_vimp': True, 'fsds': True}
- consensus: x8, x5, x0, x4, x1
- FSDS-only: x7
- shift-only: x8
- read: 多法共识(≥3): x0, x1, x2, x3, x4；仅FSDS(y判别、未必shift): x7；仅shift描述(cmean/MMD/PO): x8

## Evidence

- TencentGR: 四法齐全，共识落在路径份额/归因/共现
- DiffusionDB: abs_delta + vimp_cov + FSDS 别名归一后同样出面板；ΔY>0 正向
- Waymo proxy: 无图也能跑 cmean + cov-proxy + F；证明骨架不绑 graph schema

## Not claimed

- 跨域 AUC/ATE 可比
- 自动定罪 / 改 Drill 门
- DiffusionDB CLIP 全量 token 网格（可用轻量 TF-IDF 协议）
