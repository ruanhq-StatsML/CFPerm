# 图实体对照：接到现有口径，不是新方法

不 call 一个 package，也不发明方法论。  
推荐图这一片只是把 **已经有的** 对照 / loc / FSDS / fire 接到 Tencent-GR 的实体切上。

现有口径不动：

- 一个 Y，一堆 X；不是归因
- 对照 / loc：每块 X 对时钟 W，无 Y，moved / quiet
- fire：有 Y，OnlineRFPerm last-two，P(Y|X) hop
- FSDS：有 Y 之后最多两步（LOGO 组 → LOCO 列），Δ 是日志
- AUC / ratio 不是成绩

图只是 **切 X 的刀**（人 / 货 / 店 / 挂在店上再聚到人的塔）。不是 GNN，不是新 loc 算法。

---

## 怎么接（同一套话，换切刀）

| 现有 | 审核落地 | 推荐图落地（这一片） |
|---|---|---|
| 切 X 的刀 | 长度 / 语气 / 拒绝，或视频塔 / 文本塔 | user / item / merchant / video / audio |
| loc，无 Y | 这块文风 early vs late 像不像 | 这块实体画像像不像 |
| 手工特征 → batch 对 | `xy_*.csv` 按到达切窗 | 度数 + pool 做完，按到达切窗 |
| inference 槽 | style_vector / 审核器特征 | 店上的塔：假 DGP 或一条编码器 pipeline |
| fire，有 Y | OnlineRFPerm 审核决定 | 同一句，Y=`y_post_clk_1d`（另列，本片没当主叙事） |
| FSDS 两步 | 组 = 特征族 | 组 = loc 标成 moved 的实体；不能多步 |

Detector 还是 P(W | X_块)。没有新统计量。

---

## 数据侧只要两步半（落地，不是零件清单）

1. **手工特征：** 左窗图抽实体 X；塔 attach 在店、pool 到人。
2. **Inference 槽：** `(店或货) → 向量`。现在假 DGP；换成 LLM/编码器 pipeline 即可。
3. **同一句检测：** 每实体 moved / quiet。上线则相邻 batch 对。

FSDS 最多两步，所以 loc 故事不靠它。GNN 不接：会把实体揉进 hidden，还会把 Y 灌进 X。

---

## 结合时的建议（仍是原方法，不新造）

1. **对外不要起新名字。** 不写 Graph-Loc-FSDS、不写 graph localization package。就写：实体切的塔/族对照，接到 OnlineRFPerm / FSDS。
2. **一张表并排，不要合成一个分。** 和审核对照同一格式：

   | 窗 | user | item | merchant | video | audio | 映射 fire |
   |---|---|---|---|---|---|---|
   | t | moved/quiet | … | … | … | … | yes/no |

   loc 列无 Y。fire 列有 Y。缺 Y 就空着 fire，不要用 loc AUC 填。
3. **LOGO 的组名 = loc 的实体名。** 不要另切一套 family。quiet 实体不要送进 FSDS。两步封顶写进工单，免得有人要第三步。
4. **假 DGP 只占 inference 槽。** 真接法：encoder 写出 `(merchant_id, vec)`，后面 pool / loc 不动。不要用 Y 训这个 encoder，否则 loc 不再无 Y。
5. **上线只换时钟，不换切刀。** 离线 W=`t_end` 中位；上线相邻 batch。实体五块、无 Y loc、有 Y fire，三句话原样。
6. **审核和推荐用同一张对照语法。** 审核族、推荐实体，都是「一块 X、同一时钟、动/不动」。两份落地互证口径，不是两个方法。
7. **日志分开放。** loc AUC、FSDS Δ、fire bit、拒绝率，四列。不加权合成「图谱健康度」。
8. **GNN 若别人问起：** 图在这里是切刀；消息传递会毁掉实体对照。不接。不是「GNN 不 SOTA」，是和现有 loc 口径冲突。

---

## FSDS 两步（原限制，原样接）

1. LOGO：丢掉某一个 **moved 实体**
2. LOCO：丢掉该实体里的 **列**（塔 top-k）

没有第三步。Δ 不是贡献。画像 moved ≠ Δ>0 ≠ 根因。

---

## 实体切（刀，不是新 loc）

| 实体 | X | 接到行上 |
|---|---|---|
| user | 人—货度、人—店度 | `user_id` |
| item | 货—人度、共点度 | 当前单 `item_id` |
| merchant | 当前店度数 / PR / 聚类 | 当前单 `merchant_id` |
| video / audio | 店上向量 pool 到人 | `user_id` |

video/audio 是 user 的 1-hop 属性，不是 merchant。所以店标量可以 quiet、人侧塔可以 moved。这是切刀不同，不是新发现。

这次：user / video / audio moved；item / merchant quiet。quiet 停。

---

## 三种变动（原方法里本来就分开）

| | loc / 对照 | FSDS | GNN（不接） |
|---|---|---|---|
| 问 | 这块 X 像不像 | 和 Y 缺口有没有关系 | 表征/预测器换没换 |
| Y | 不要 | 要 | 通常要 |
| 几步 | 并排一次 | 最多两步 | 半径 ≠ 下钻 |

```
某实体 loc moved    ≠  该实体是 Y 缺口来源
FSDS LOGO Δ>0       ≠  该实体画像动了
GNN hidden 漂了     ≠  实体对照
```

映射 fire 另列，OnlineRFPerm。loc 不替代 fire。

`scripts/tencent_gr/graph_loc_fsds_drill.py` 是落地脚本，不是新方法实现。
