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

---

## 这里为什么不用 GNN

图在这里是 **切实体的坐标**，不是要训一个过图的预测器。GNN 的 hop 是消息传递半径，不是定位，也不是 FSDS 下钻。

别人在推荐图上掏 GNN，通常就这几条：

| 想用 GNN 干什么 | 这里为什么不对 |
|---|---|
| 消息传递自动吃 k-hop，少手写特征 | 手工那几步度数 / pool 已经够出 batch 对；再传会把实体揉进一个 hidden |
| 用图结构预测点击 / 转化 | 那是 P(Y\|X)，loc 问的是无 Y 的 P(W\|X_实体) |
| embedding 漂了 = 图变了 | 漂的是混在一起的表征，说不出 user 动、merchant 静 |
| GNNExplainer / 邻域重要性 | 读成归因；本包不做归因 |
| 换一版 GNN checkpoint | 那是模型 hop，不是实体画像 hop |

还有几条独立于「FSDS 只能两步」的原因：

1. **实体被消息传递混掉。** loc 的刀是 user / item / merchant / video / audio 并排。GNN 一层聚合就把人、货、店折进同一个向量，对照表没了。
2. **监督会把映射灌进 X。** GNN 几乎总是带着 Y 训。训完的向量不再是纯画像，不能当「无 Y 的 loc」。
3. **向量槽已经被 inference pipeline 占了。** 店/货上的塔来自 LLM/编码器（现在是假 DGP）。GNN 是第二套图编码器，loc 不需要。
4. **GNN 的 hop ≠ 变动。** 2-hop 邻居只是半径。变动是相邻 batch 上这块 X 像不像。半径写在特征口径里（度数是 0-hop，塔 pool 是 1-hop），不必用 GNN 再传一层。
5. **Online 时更糟。** 每个 batch 重训/微调 GNN，表征自己先 hop，实体 loc 和模型 hop 分不开。

所以：有图 ≠ 上 GNN。networkx 抽实体 X 就够。

---

## GNN 和 FSDS 的变动逻辑（不要混）

三种「变了」不是一件事。

| | **Graph-loc（本包）** | **FSDS** | **GNN** |
|---|---|---|---|
| 问 | 哪一 **实体画像** early vs late 不像 | 这块 X 和 **Y 缺口** 有没有关系 | 过图的 **表征/预测器** 换没换 |
| Y | **不要** | **要** | 通常要（监督） |
| 变动单位 | 实体（五块并排） | 事先切好的组：LOGO 实体、LOCO 列 | 节点 hidden / 一层权重 |
| 能走几步 | 并排一次，quiet 停 | **最多两步**，没有第三步 | k-hop 是半径，不是下钻 |
| hop 一词 | 不用 | 不是图 hop | 邻居半径 |
| 读法 | moved / quiet | Δ 是日志，不是贡献 | 不能读成「哪座塔动了」 |

**Loc 的变动：** 同一套手工实体 X，相邻 batch 对。P(W \| X_user) 分开了、P(W \| X_merchant) 没分开 → 人动店静。没有 Y。

**FSDS 的变动：** 已经有 Y、已经有实体名单之后。φ=(Y−μ)(W−e)。LOGO 丢掉某一个实体看 R 降不降；LOCO 再丢掉一列。这是「Y 缺口绑没绑在这块 X 上」，**不是**「这块画像动了」。最多两步，所以不能拿 FSDS 当多步定位器，也不能替代 loc。

**GNN 的变动：** 图一变或 checkpoint 一换，hidden 就变。那是表征漂了，还是预测器 hop 了，和「user 实体画像动了」不是同一句话。k-hop 只是聚合半径：2-hop GNN 不会告诉你第二步该下钻谁。若用 Y 训，变动里已经混进映射，loc 和 fire 分不开。

对照（同一时钟）：

```
GNN embedding 漂了     ≠  user 实体 moved
FSDS 某实体 LOGO Δ>0   ≠  该实体画像动了
某实体 loc moved       ≠  该实体是 Y 缺口来源
```

三句都要留着。本包只做第一列 loc。FSDS 两步是有 Y 时的日志上限。GNN 不进包。

`scripts/tencent_gr/graph_loc_fsds_drill.py`  
实体口径仍按 user / item / merchant / video / audio。
