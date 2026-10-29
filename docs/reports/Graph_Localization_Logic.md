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

## 长名单 → 短名单 → 再拼起来

Loc **先不拼**。五块各打各的，输出五行，这是 **长名单**（五个实体都在，moved 或 quiet）。

**短名单** = 长名单里 `moved=yes` 的行。quiet 留在长名单上当记录，**不进短名单**。

**拼起来** 发生在短名单之后：同一批行，只把短名单实体的列水平拼成 `X_short`。  
不是把五行折成一行，也不是先把五块全拼再 loc。

```
长名单 L = {user, item, merchant, video, audio}     # 五行，并排 loc
短名单 S = {e ∈ L : moved_e}                        # 这次 {audio, video, user}
X_short = concat({ X_e : e ∈ S })                   # 只拼 S，行数不变
quiet 的 X_item、X_merchant 不进入 X_short
```

这次：长名单 5 行；短名单 3 行；`X_short` 是 64+64+2=130 列，n 仍是 12564。item/merchant 的列不在里面。

| | 长名单 | 短名单 | 拼起来之后 |
|---|---|---|---|
| 是什么 | 五实体对照表 | moved 的那几行 | `X_short` 一张宽表 |
| 列 | 各实体自己的 X_e | 还是各实体的 X_e | S 里实体的列水平拼接 |
| 行 | 五个实体 | 少于等于五 | 和输入样本同一批行 |
| 干什么 | loc（无 Y） | 决定谁进下一步 | 只在这张宽表上做下一步 |

长名单上 quiet 的实体：停。不拼、不打开、不进下一步。  
短名单若只有 1 个实体：`X_short = X_e`，实体之间没有可丢的块（LOGO 空），只能在列上走一步。  
短名单若是 5 个：S=L，拼起来等于全表，但那是 **loc 之后碰巧全动了**，不是一开始就拼。  
短名单若是空的：没有 `X_short`。不要把 quiet 硬拉进来凑名单（脚本里「全 quiet 则取 AUC 最高一个」只是防空，不是口径）。

拼起来之后还是同一批行、同一个 W。变的只是列集合：从「一次一块」变成「短名单那几块并排躺在一张表里」。  
下一步如果要丢块，丢的是短名单里的 **一个实体**（整块列），不是长名单里的 quiet 实体，也不是 64 个 loc。

---

## 还有其他的吗（大概）

就这三层：长名单、短名单、`X_short`。没有第四张名单。

顺带几条边：

- 两张表别混：长名单 5 行（实体），`X_short` 是 n 行（样本）。五行不是拼成一行样本。
- 短名单不按 AUC 排序当重要性；yes/no 收名单而已。
- video/audio 进短名单，仍是 user 键上的塔，不是把 merchant 从 quiet 改成 moved。
- 没有 Y 时，拼起来可以不做；短名单已经够指「去看哪几块输入」。
- 有 Y 时，评估只发生在 `X_short` 上，最多两步（丢实体、丢列）。Δ 不改长名单上的 moved。
- 下一窗再 loc：重新出一张长名单，短名单可以换人；不要沿用上一窗的 S 去拼这一窗。

`scripts/tencent_gr/graph_loc_fsds_drill.py`：`localize_blocks` 出长名单；`loc_cols` / `X2` 是短名单拼起来的 `X_short`。
