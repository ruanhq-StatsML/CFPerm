# Graph-localization package

FSDS 没法多步下钻，所以 **只讲 graph-localization**。  
Package 很小，直观：做 data-manipulation，再搭一条 inference pipeline，按实体打 moved / quiet。够了。

不是 GNN。不是因果。不是 FSDS 包装。

---

## Package 里有哪几块

就两块半。没有第三块业务。

| # | 块 | 干什么 | 现在 prototype | 换成真的 |
|---|---|---|---|---|
| 1 | **Data-manipulation** | 左窗切图；按实体抽画像；店上 attach 的向量沿边聚到人 | `split_lr` + networkx 四张图 + mean-pool | 同一套；换边表即可 |
| 2 | **Inference pipeline** | 给 item/merchant 打向量（视频/音频/文本） | 假 DGP：`randint` id + 64-d t/U | **LLM / 编码器 inference 一条 pipeline 就够** |
| ½ | **Detector** | 每实体 P(W \| X_实体)，无 Y → moved / quiet | RF-domain AUC ≥ 0.55 | last-two / 同一句分类器 |

没有：GNN、全图 embedding 当归因、FSDS 多步、online 训练器、社区发现。  
Detector 半块是因为没有新数据、没有新模型，只是每实体一块 X 打时钟。

**是的：data-manipulation + 一条 LLM-inference pipeline 就足够。**  
假 DGP 占的就是 pipeline 那个槽。换成真推理，loc 口径不动。

---

## 为什么不讲 FSDS

LOGO 一次丢掉一整块，不会自己按实体往下走。  
所以 FSDS 不是 localization package 的一部分。它最多是 loc 出名单之后、有 Y 时的日志。本 package 不靠它讲故事。

---

## Data-manipulation 具体做什么

左窗 clk/cnv（曝光不成边）→ 三种节点实体：

| 实体 | 抽出的 X | 接到行上的键 |
|---|---|---|
| user | 人—货度、人—店度 | `user_id` |
| item | 货—人度、共点度 | 当前单 `item_id` |
| merchant | 店—人度、投影 PR/聚类 | 当前单 `merchant_id` |

塔：inference 打在 **店** 上 → 沿人—店边 pool 到 **人**。  
video / audio 是 user 的 1-hop 属性，不是 merchant 实体。  
聚合键写进口径，不要事后猜。

时钟 W = `t_end` 中位。同一批行、五块 X 并排。这就是 manipulation 的全部。

---

## Inference pipeline 具体要什么

只要一张表：`(merchant_id 或 item_id) → 向量`。

现在：1e4 个假 id，64 维噪声，乱七八糟 attach 到店。  
换成 LLM-inference：对店的视频/音频/文案跑编码器，写出同一张表，后面 mean-pool 不用改。

一条 pipeline 够：encode → attach 到图上的实体 → loc。  
不要第二套训练、不要用 Y 去训这个 encoder 再回头说定位到了图。

---

## Detector（半块）

每实体：只有 `X_实体`，问 early vs late 像不像。  
moved / quiet。AUC 是日志。  
quiet 的实体停住——不是再丢给 FSDS 打开。

这次：user / video / audio 动；item / merchant 静。

---

## 还有其他的吗

没有必须的第三块。

| 可以后做 | 不进这个 package |
|---|---|
| 真 LLM pipeline 换掉假 DGP | GNN |
| Detector 换成 last-two（online 窗） | 用 FSDS 当多步 loc |
| 映射 fire 另列（要 Y） | 把 loc AUC 当成绩 |

Online training / LOGO 若做，组名仍是这五个实体；那是下一片，不是这个 package 里缺的零件。

`scripts/tencent_gr/graph_loc_fsds_drill.py`  
实体口径仍按 user / item / merchant / video / audio。
