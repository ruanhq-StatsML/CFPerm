# Localization：输入、输出、颗粒度

不管映射死没死。这里只做 loc。

**颗粒度 = 实体**（不是列，不是整张图，不是 GNN 节点）。  
一块实体 = 一组列。一次 loc 只问这一块。

实体五块：`user` / `item` / `merchant` / `video` / `audio`。

---

## 输入

一张表，一行一个样本（这里是一笔 post-click 订单）。**不要 Y。**

| 列 | 含义 |
|---|---|
| `t_end`（或 `batch`） | 时钟。W = 1{t > 中位}，或相邻两窗 |
| `g_u_*` | **user** 实体 |
| `g_i_*` | **item** 实体 |
| `g_m_*` | **merchant** 实体（当前这一单的店） |
| `vid_*` | **video** 实体（店上向量 pool 到人） |
| `aud_*` | **audio** 实体（同上） |

列从哪来（手工特征，做完就能切 batch 对）：

- 左窗 clk/cnv 建人—货、人—店边（曝光不成边）
- user：该人度数；item：该货度数/共点；merchant：该店度数/PR/聚类
- video/audio：`(merchant_id → 向量)` attach 到店，沿人—店边 mean-pool 到 user  
  向量现在是假 DGP；换成编码器同一张表

---

## 输出

一张对照表，一行一个实体：

| entity | n | AUC | moved |
|---|---:|---:|---|
| audio | 64 | 0.853 | yes |
| video | 64 | 0.829 | yes |
| user | 2 | 0.605 | yes |
| item | 2 | 0.523 | no |
| merchant | 4 | 0.501 | no |

- `moved` = 这块 X 能分开 W（这里 AUC ≥ 0.55，门是日志）
- `quiet` 的实体停住，不再打开
- 输出 **不是** 谁导致 Y，也不是 fire

---

## 怎么做

对每个实体 `e`：

1. 取出该实体的列，得到 `X_e`（只有这些列）
2. `W` 来自时钟列（没有 Y）
3. 拟合 P(W | X_e)，hold-out AUC
4. AUC ≥ 0.55 → moved，否则 quiet

五块并排，同一批行、同一个 W。不要把五块拼成一块再事后分解。

```
for e in {user, item, merchant, video, audio}:
    X_e = 该实体列
    AUC_e = rf_domain(X_e, W)   # 无 Y
    moved_e = AUC_e >= 0.55
```

---

## 颗粒度不要混

| 实体 | 列 | 接到行上的键 | 问的是 |
|---|---|---|---|
| user | `g_u_item_deg`, `g_u_merch_deg` | `user_id` | 这个人连了多少货/店 |
| item | `g_i_user_deg`, `g_i_coclick_deg` | 当前单 `item_id` | 这件货的连通 |
| merchant | `g_m_*` | 当前单 `merchant_id` | 这家店的 4 个图标量 |
| video | `vid_00`…`vid_63` | `user_id` | 人走过的店的视频向量 |
| audio | `aud_00`…`aud_63` | `user_id` | 同上，音频 |

video/audio 接在 user 上，不是 merchant。所以 merchant 可以 quiet、video 可以 moved。

---

`scripts/tencent_gr/graph_loc_fsds_drill.py` 里 `localize_blocks` 就是上面这个循环。
