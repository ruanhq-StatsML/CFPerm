# Graph-localization package

FSDS **最多两步下钻**（LOGO 丢实体 → LOCO 丢列），不能再往下走。所以 **只讲 graph-localization**。

手工做完那几步特征之后，同一套 X 按到达切窗，就是 **连续 batch 对**。loc 还是每实体 P(W | X_实体)，只是 W 换成相邻 batch。

Package 很小：data-manipulation + inference pipeline，按实体 moved / quiet。够了。

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

## FSDS 的硬限制（很明确）

最多两步，没有第三步：

1. LOGO：丢掉某一个 **实体**
2. LOCO：丢掉该实体里的 **列**（塔只 top-k）

不会按实体再 hop、再社区、再 2-hop。没有名单就连这两步也开不好。  
所以 FSDS 不是 loc package 的零件，也不能拿它当多步定位器。有 Y 时的日志，到此为止。

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

手工这几步做完，每行有五块实体 X。离线 W = `t_end` 中位；上线把行按到达切成 batch，相邻两窗就是连续 batch 对。特征口径不变。

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

**没有。** loc package 到 batch 对就齐了。

隔壁另列（不是这个包缺的零件）：映射 fire，要 Y，OnlineRFPerm。实体 loc 不替代 fire。

不进包：GNN、FSDS 第三步、更多实体、用 Y 训 encoder。

`scripts/tencent_gr/graph_loc_fsds_drill.py`  
实体口径仍按 user / item / merchant / video / audio。
