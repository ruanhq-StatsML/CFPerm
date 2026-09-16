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

## 五块并排，不是事后分解

Loc 的输出是 **五行方表**，不是一个 100% 切成五份的饼。

**并排：** 同一批行、同一个 W，五次独立的 P(W | X_e)。  
每次只用该实体的列。读法是五句 yes/no：这块像不像。  
user moved、merchant quiet 可以同时成立，不必加起来等于 1。

**分解（不要）：** 先拼 `X = [user | item | merchant | video | audio]`，训一个 P(W | X_all)，再用 VIMP / SHAP / LOGO-share 把可分性「分」回五个实体。  
那是联合分类器的贡献账，不是 loc。

| | 并排 | 分解 |
|---|---|---|
| 模型 | 五个，各吃自己的 X_e | 一个，吃拼起来的 X |
| 数字 | 五个 AUC，各管各的 | 一份份额，凑成 1 |
| quiet 实体 | 自己分不开 W 就停 | 仍可能分到质量（和 moved 块相关） |
| 问的问题 | 光看这块，W 分得开吗 | 在大杂烩里谁占权 |

拼起来再分解会把颗粒度毁掉：video 和 user 都接在 `user_id` 上，联合模型分不清是度数动了还是塔动了，却仍会吐出两个份额。  
并排则是两次单独问：只拿度数能不能分开 W；只拿塔能不能分开 W。这次就是 merchant 0.501 quiet、video 0.829 moved——并排才看得见。

---

## 每实体一行之后，下一步怎么评估

先看这五位：`moved` 还是 `quiet`。不要用 AUC 大小排名当成绩。

**quiet 的行：** 停。不进下一步，不打开列，不送 FSDS。  
这次：item、merchant。

**moved 的行：** 才是短名单。下一步只开这几块。  
这次：audio、video、user。

下一步仍是 **评估**，不是再 loc 一层（颗粒度已经是实体，不要拆成 64 个 loc）。

有 Y 时，最多两步，而且只在短名单上：

1. **LOGO（实体之间）：** 只用 moved 那几块拼起来的 X。丢掉某一个 moved 实体，R 降不降。  
   评估的是：这块 **和 Y 缺口** 有没有关系。Δ≤0 = 画像动了，但不是这个缺口来源。不要把 Δ 拿回去改 loc 的 moved。
2. **LOCO（实体内部）：** 还在短名单里，丢掉一列（64 维塔只 top-k）。  
   评估的是列，不是再发明一个实体。没有第三步。

没有 Y：下一步不是 FSDS。工单按实体走——video moved 就去看 `(merchant_id → 向量)` 那条 inference 和 pool；user moved 就去看度数怎么算的。quiet 的店标量不用动。

不要用的评估：

- 五块 AUC 排序当「谁最重要」
- 用全表 VIMP 检验 loc 对不对（那是分解）
- loc moved 却 LOGO Δ≤0 就说 loc 错了（两个问题）

这次短名单三块 LOGO 全是 Δ≤0：并排 loc 说画像动了；两步评估说和 Y 缺口没关系。两张表都留着。

`scripts/tencent_gr/graph_loc_fsds_drill.py` 里 `localize_blocks` 是并排循环；短名单上的 LOGO/LOCO 是下一步评估。
