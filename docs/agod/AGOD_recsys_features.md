# 推荐系统的特征做法（因果时间 OOF）

看的是 Amazon 多模态推荐那条流水线：[`JingxiangQU/mmoe-multimodal-rec`](https://github.com/JingxiangQU/mmoe-multimodal-rec) 的 `data4moe_beam.py` / `data4model.py`，也就是 AGOD 吃的 WebDataset [`jingxiang11111/amazon_reviews_for_rec`](https://huggingface.co/datasets/jingxiang11111/amazon_reviews_for_rec)。

特征合同是 **时间**，不是折号。和上一张 Super Learner OOF 是同一件事：t 时刻的 X 不能由 y_t 算出来。模型头上的 MMoE 才是 `π(x)`——特征和门不要混。

---

## 他们实际写了什么

`CausalPosNegByUser`（按 user 排序）是整条协议：

```
for each user, sort reviews by timestamp:
    user_feat ← 只看 *过去* 的评论
    emit (user_feat, item_meta, image, label_good/best)
    *然后* 才把本条评论折进 running state
```

| 字段 | 从哪来 | 能不能看当前条 |
|---|---|---|
| `cat_hist` / `review_cnt` / `price_mean,std` | 过去评论的 Welford / 计数 | 否（先 emit 再 update） |
| `history`（最多 3 条 title+text） | 过去评论 | 否 |
| item 文本 | 类目 / 标题 / 价格 / features / description（目录元数据） | 目录，不是评分目标编码 |
| 图像 | 主图 → 14×14 patch | 目录 |
| `label_good` (≥4★) / `label_best` (5★) | **当前** 评分 | 标签，不当特征 |
| 当前评论 title/text | Enrich 里有，写 WebDataset 时丢掉 | **不准进 X** |

负样本：同一份因果 `user_feat`，随机抽没见过的 `parent_asin`，标签置 0。切分是 `SplitByDate`：train ≤ **2023-06-30** < valid ≤ **2023-09-30**，不是随机行切。

`build_user_text` / `build_item_text` 只是把结构化画像 **展成一段话**，给句级交叉注意力用。这是序列化，不是第二种泄漏。

---

## 和 OOF stacking 的同一张表

| | Super Learner / 上一张 | 这条推荐流水线 |
|---|---|---|
| 禁止 | `Z_{i,m}=f_m(x_i)` in-sample | 当前评论文本、未来的 user 统计 |
| 允许 | `f_m^{(-i)}(x_i)` | t 之前的评论 + 当时的目录 |
| 切分 | 折 / probe vs holdout | **事件时间** |
| 典型泄漏 | 噪声专家背 y | 用本条 review 预测本条星级（几乎是情感分类） |
| 不是这件事 | — | MMoE 门 `softmax(W·mean(experts))` 看的是 x |

本例玩具流：因果路径 current-review 泄漏率 **0.00**；若把当前条折进 `user_feat`，泄漏率 **1.00**。当前评论文本 vs `label_good` 相关 **1.00**，因果 `review_cnt` 只有 **0.03**。留下本条评论，任务就塌成读这句话。

---

## 他们还留的洞（特征侧，不是 MMoE）

- **负样本 item 池**是全局抽 10k `parent_asin`，不按 t 过滤上架时间。
- **item meta 是快照价格**，不是事件时刻的价。
- **5★ 降采样在 GroupByUser 之前**，`cat_hist` 不是真实历史，是抽过的历史。
- 过去评论的 *文本* 带情感（「did not work」），这是故意的用户画像，只要它是 **过去** 的就不是 y_t 泄漏。
- 没有 user_id / item_id embedding：冷启动友好，协同信号丢掉。

---

## AGOD 现在把合同扔了

`scripts/run_agod_amazon_modality_lr.py` 做的是 `text = item + "\n" + user`，再 `HashingVectorizer`。user 画像和 item 目录变成 **一个 ngram 袋**。因果切分还在样本行上，但模型看不到「这是过去的用户、那是当前的商品」。PO/VIMP 还在当前窗上 in-sample 拟合——相对这条推荐协议，那是 leaky stacking。

MMoE 头：`π=softmax(W mean(expert_vecs))`，专家是 user 文本 / item 文本 / 图 / 交叉。**那是 MoE。** 特征流水线不是。上一张说的 stacking `π(t)` 仍然不是这个门。

---

## 可复现

```bash
PYTHONPATH=. python3 -m tests.test_agod_recsys_features
PYTHONPATH=. python3 scripts/run_agod_recsys_features.py
```

源码对照：`CausalPosNegByUser` 先 yield 再 `hist.append`；`SplitByDate`；`build_user_text` 只用 `user_feat`。
