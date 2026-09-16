# Graph localization：实体口径

Prototype。五粒够了。不是 GNN，不是因果图，不是归因。

**用起来就两句话：** 先按实体做 graph-localization，再只对 moved 的那几个实体做两步下钻。

---

## 0. 为什么必须先 loc

两件 FSDS 自己干不了，所以 loc 在前、下钻在后。

**第一，FSDS 没法做多步下钻。**  
LOGO 一次丢掉一整块，问 R 降不降。它不会先告诉你「人动了、店没动」，再往人里面走一步。没有实体名单，就只能对着全表或七个 hop 做一次 drop，那不是下钻，是一锅粥。  
所以先 loc：实体切成少数几块，quiet 的关掉。两步下钻只开 moved 的那几个——这正是 FSDS **能**做的那点事。

**第二，下一步 online training + LOGO 也要吃同一套实体口径。**  
Online 要有稳定的组名（user / item / merchant / video / audio），每组一块 X，来一个 batch 打一遍 moved/quiet，再对 moved 组做 LOGO。没有实体切，online 没有组，LOGO 没有可丢的块。  
本片只做到 **loc → 两步下钻**。online 不在这里实现；口径已经按实体钉死，下一步直接用。

---

## 1. 大体口径（钉死）

**实体 = 左窗图上的一种节点，或挂在节点上再沿边聚合的一座塔。**  
localization 的单位是实体，不是列，也不是「整张图」。

一张图、三种节点实体、两座挂在店上再聚到人的塔：

| 实体 | 图上是谁 | 这块 X | 接到订单行的键 | 这次 |
|---|---|---|---|---|
| user | `u` | 人—货度数、人—店度数 | `user_id` | 动 0.605 |
| item | `i` | 货—人度数、同场共点度数 | `item_id` | 静 0.523 |
| merchant | `m` | 店—人度数、投影 PR / 聚类 / 投影度 | **当前单** `merchant_id` | 静 0.501 |
| video | 属性在 `m` | 沿人—店边 mean-pool 到 `u`，64 维 | `user_id` | 动 0.829 |
| audio | 属性在 `m` | 同上 | `user_id` | 动 0.853 |

四条规矩，少一条就不是这个口径：

1. **先实体、后列。** loc 不问 `g_u_item_deg` 重不重要，问 user 这块画像动没动。
2. **一块 X 只含该实体。** 人度数和店 PR 不要拼成「graph」再事后分解。视频塔不要和音频塔捆成「多模态」。
3. **没有 Y。** 每实体只打 P(W \| X_实体)。moved = 这块 early vs late 不像了。
4. **同一时钟、同一批行。** 并排才叫对照。AUC 是日志门，不是成绩。

曝光不成边。clk/cnv 在左窗连边。行是 post-click 订单。W = `t_end` 中位。  
问的是：后来的人，在这张冻住的左窗图上，**哪个实体的画像**不像早的人。

---

## 2. 不同实体，口径不要混

三种节点看起来都是「度数」，不是同一口径。

**user**  
主体是人。X 是这个人在左窗图上连了多少货、多少店。  
merge 键是 `user_id`。一个用户多笔订单，共用同一份人侧画像。

**item**  
主体是货。X 是这件货被多少人点、同场还和谁一起被点。  
merge 键是这一单的 `item_id`。问的是 **当前这件货** 的连通，不是「用户爱点的货平均长什么样」。

**merchant**  
主体是店。X 是 **当前这一单挂的那家店** 的图结构（度数、投影 PageRank、聚类）。  
merge 键是这一单的 `merchant_id`。不是「这个人走过的店的平均结构」。

**video / audio**  
主体仍是人，但这块 X 不是度数。向量挂在店上，沿人—店边 visist-weighted mean-pool 到人。  
所以它们是 **user 实体上的 1-hop 属性塔**，不是 merchant 实体。  
口径后果：`merchant` 可以静（当前店 4 个图标量分不开 W），`video`/`audio` 可以动（人走过哪些店的 64 维味道能分开 W）。同一张图，实体不同，不要收成一句「店变了」或「视频是根因」。

塔是 DGP：1e4 个 id，64 维噪声，商户上乱七八糟 attach。协议和真塔相同——attach 在哪个实体、聚合到哪个实体，必须写在口径里。

---

## 3. 逻辑（按实体走一遍）

```
左窗 clk/cnv → 人—货、人—店、共点、店投影
     ↓
按实体抽出 X_user, X_item, X_merchant, X_video, X_audio
     ↓
localization（无 Y）
  每个实体：P(W | X_实体) → moved / quiet
  quiet 关掉
     ↓
两步下钻（有 Y，只开 moved 实体）
  ① LOGO：丢掉某一个 moved 实体，R 降不降
  ② LOCO：该实体内部的列（塔只 top-k）
```

FSDS 只出现在虚线以下。上面没有 Y，下面没有 quiet 实体。  
这就是「FSDS 做不了多步、所以 loc 先按实体收名单」。

---

## 4. 这次并排怎么读

| 实体 | 动？ | 读法 |
|---|---|---|
| user | 是 | 后来的人，点的货/店 **数量** 不像 |
| item | 否 | 当前货的连通还像 |
| merchant | 否 | 当前这一单的店，4 个图标量还像 |
| video | 是 | 人走过的店，视频 DGP 味道不像 |
| audio | 是 | 同上，音频 |

下钻只开 user / video / audio。item、merchant 停住。  
第一步 LOGO 三块全是 Δ≤0：这三块画像动了，不是这个时钟上的 Y 缺口来源。RF-mass 可以很响。不要把 share 写成贡献。

---

## 5. 和下一步 online + LOGO 的口径怎么接

不在本片实现。组名已经是实体，下一步不用再切一刀。

| 现在（离线 prototype） | 下一步（online，还没做） |
|---|---|
| 一个 cut，W = t_end 中位 | batch 到达，每实体 last-two / 域可分 |
| loc：P(W \| X_实体)，无 Y | 同一句，X_实体 随 batch 更新 |
| 两步下钻第一步 = LOGO 丢实体 | 只对 **本窗 moved 的实体** 做 LOGO |
| 第二步 LOCO | 仍只在名单内部；online 可先不做 |

Online training 若要做，训的是每实体各自的 P(W \| X_实体)（localization 那条），不是用 Y 去训「哪个实体重要」。  
映射 fire 另列，要有 Y，走 OnlineRFPerm。实体 loc **不替代** fire。

---

## 6. 不要做的

| 做 | 不做 |
|---|---|
| 五种实体并排 loc | 把人/货/店拼成一块「图特征」 |
| video 写成 user 的 1-hop 塔 | 写成 merchant 实体，或「视频导致转化」 |
| quiet 实体停住 | 用 FSDS 代替 loc 去做多步 |
| LOGO 组名 = 实体名 | 136 列当 136 个实体 |
| 口径留给 online | 本片里实现 online training |

`scripts/tencent_gr/graph_loc_fsds_drill.py`  
结果：`docs/reports/Recsys_Graph_Loc_FSDS.md`
