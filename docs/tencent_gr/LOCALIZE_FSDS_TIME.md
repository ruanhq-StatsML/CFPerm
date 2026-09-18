# TencentGR prototype: time-window localization → FSDS

两步：**Localization → FSDS**。FE 带时间范围；两窗间隔 ≥1 个月；无 network。

> TencentGR-1M labels = exposure(0)/click(1)；terminal success = **click**。

## Protocol (no leakage)
1. `feature_engineer(root, t_start, t_end)` — 只用窗内事件
2. W1 early / W2 late，gap = **30.0** days (≥ 30)
3. Localization subset = W1 `share_linear` top-k（三段启发式 aggregate through）
4. FSDS select+train 只在 W1；W2 仅作 temporal holdout

## Windows
- timeline span: **231.3** days
- W1: `W1_early` [1737673172, 1741993172) ≈ 50.0d — edges=129641 pos=0.0001
- W2: `W2_late` [1744585172, 1748905172) ≈ 50.0d — edges=4310736 pos=0.0003
- localize-k: **500** (W1-train shares) → W1-train=8435, W1-holdout=2714, W2=11963

## FSDS
- W1→W1 user-holdout: HGB AUC=0.642 AP=0.001 | LR AUC=0.147 AP=0.000 | n=8435/2714 sel=24
- W1→W2 temporal: HGB AUC=0.473 AP=0.001 | LR AUC=0.677 AP=0.002 | n=8435/11963 sel=24

## Selected features (W1 fit)
e_n_exp, u_n_events, u_n_exp, u_n_uniq_items, u_span_sec, u_log1p_n_events, u_log1p_n_uniq, u_user_activity_rank, i_n_exp, i_n_users, i_n_covisit_neighbors, i_log1p_n_exp, i_log1p_n_users, i_log1p_n_covisit, i_item_pop_rank, e_log1p_exp, ui_pop_mismatch, i_share_first, i_share_last, i_share_linear, i_credit_first, i_credit_last, i_credit_linear, i_item_credit_rank

## Guards
- no events outside window in FE
- drop `y_convert` / `*_n_cnv` / `*cvr*` from X
- subset + SelectKBest never see W2 labels
