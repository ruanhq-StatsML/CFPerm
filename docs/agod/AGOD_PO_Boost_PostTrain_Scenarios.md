# PO-Boosted Post-Training — scenarios & gains (summary)

> LaTeX：[`AGOD_PO_Boost_PostTrain_Scenarios.tex`](AGOD_PO_Boost_PostTrain_Scenarios.tex)  
> PO 足够；增益 = 更新 FLOPs↓ / \(T(\mathrm{Acc}^\star)\)↓（Acc 约束）+ reject 后 √PO 行权。

## 场景 → 增益

| 场景 | PO 做什么 | 增益 |
|---|---|---|
| **S1** 不对称多模态 (Affec \(M\ge3\)) | Softmax/fuse → freeze 低残差塔 | FLOPs≈**0.70**（gated），Acc≈平 |
| **S2** 概念尖峰 | 短轨 ΔPO 倾倒 step；长轨稳集合 | \(T(\mathrm{Acc}^\star)\)↓ |
| **S3** 协变量糊 | PO−λMMD / gate | 不给纯 X 漂买 step |
| **S4** Reject 批 | \(w\propto\sqrt{\mathrm{PO}}\) | next-MSE；calm 仍 w=1 |
| **S5** 饱和 2-mod (COCO…) | 不冻，最多软 LR | **近 0 增益——别报** |
| **S6** 噪 vs 漂不清 | po_gated | 主 FLOPs 砍法 |

## 一句话

PO 不换 loss，只重分后训练预算；真增益在 S1/S2/S4/S6，S5 诚实报 noop。
