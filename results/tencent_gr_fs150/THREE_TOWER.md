# 三塔 architecture

v1 是支路 concat 后接 HGB/LogReg/MLP。三塔把支路收成三路表征再融合。  
SKU 塔受 **daily mix 闸**：这批 `sku_clk_share=0.11%<1%`，闸关，SKU 塔输出 **0、不更新**。不是没做，是构成率不允许重训。

```
                    cnv_ts 已知（满窗 Y=y_post_clk_1d，按用户切）

         ┌──────────────────┬──────────────────┬──────────────────┐
         │   User 塔        │   Ctx 塔         │   SKU 塔         │
         │   这个人         │   这一单怎么来    │   这件漏斗        │
         │                  │   / 当场         │                  │
         │ n_clk_before_7d  │ empty_any        │ wo_prior_clk     │
         │ n_clk_before_1d  │ dt_any (+miss)   │ dt_item (+miss)  │
         │ lag_post_clk_1d  │ sess_clk_before  │ item_within_5m   │
         │ lag_empty_any    │ sess_pos         │ n_clk_same_before│
         │ n_prior_cnv      │ log1p_price      │                  │
         └────────┬─────────┴────────┬─────────┴────────┬─────────┘
                  │ MLP_u            │ MLP_c            │ MLP_s × GATE
                  │ d=16 ReLU        │ d=16 ReLU        │ d=16 ReLU
                  │                  │                  │ GATE=0 → 0向量
                  │                  │                  │ GATE=0 → freeze W_s
                  └────────┬─────────┴────────┬─────────┘
                           │ concat [u | c | s]  ∈ R^48
                           │ Fusion MLP  32 → 1
                           │ σ
                           ▼
              Ŷ = P(满窗后 1d 任意点击 | cnv_ts)
              ≠ P(哪次点击导致购买)
```

**GATE（每天一版）**

```
sku_clk_share = P(同品 last-clk 非空)     # 这批 0.0011
GATE = 1{share ≥ 1%}
GATE=0: s=0, 不反传 SKU 塔，fusion 只用 [u|c]
GATE=1: 打开 SKU 塔，整网（或只 fusion+SKU）重训
```

不要在 14/12866 上硬开闸。同品列继续进 daily mix 监控，占比跳了再开。

**Fusion 的显式版 = User × Ctx 交叉**（`scripts/tencent_gr/cross_feats.py`）。  
SKU 交叉跟 GATE。不要 `heat × empty`（empty ⇒ 7d 热度=0）。

**和 v1 的关系**

| | v1（已训） | 三塔 |
|---|---|---|
| 热度+滞后 | concat 特征 | User 塔 |
| 空路径+任意+当场 | concat 特征 | Ctx 塔 |
| SKU | 闸关，没进模型 | SKU 塔输出 0 |
| 头 | HGB / LogReg / MLP(16→8→1) | Fusion 32→1 |

v1 满窗 1d、user 70/30、HGB：heat 0.678 → +lag **0.746**（AP 0.374）。当场/任意路径几乎不动。SKU 支路 **没有重训**（闸关）。这就是预期效果。
