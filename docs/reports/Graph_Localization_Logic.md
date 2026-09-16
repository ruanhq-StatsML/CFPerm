# 图上 localization → 两步下钻

Prototype。这些维度够了，不需要再铺别的粒。

**用起来就两句话：** 先做 graph-localization，再只对 moved 的那几粒做两步下钻。

不是 GNN，不是因果图，不是归因。

---

## 用法

```
五粒（够了）
  user_connectivity
  item_connectivity
  merchant_structure
  video_tower      ← 商户上乱七八糟 attach，再聚到用户
  audio_tower      ← 同上

1. Graph localization（没有 Y）
   每一粒单独打 P(W | X_粒)
   moved / quiet
   quiet 的粒到此为止，不再打开

2. 两步下钻（有 Y，只打 moved 的那几粒）
   第一步  LOGO：这几粒之间谁和 Y 缺口有关
   第二步  LOCO：那几粒里面的列（64 维塔只 top-k）
```

漏斗、社区、2-hop、clk vs cnv、图对图，这个 prototype **都不加**。  
localization 出名单；下钻只开名单。只下钻几个维度，就是这样。

---

## 图不是判决，是坐标

左窗 clk/cnv 建成图（曝光不成边），五粒是图上切开的 X。每一粒问的是：

> 光看这块画像，early vs late 还像不像？

读法只有 **动了 / 没动**。AUC 是日志门，不是成绩。  
定位的是哪一块图画像换了，不是哪一个点导致了转化。

然后才把 moved 的几粒交给两步下钻。用 Y 去决定「哪块图重要」，就不是 localization 了。

---

## 五粒是什么

| 粒 | 这块 X | 这次 |
|---|---|---|
| user_connectivity | 人连了多少货、多少店 | 动 |
| item_connectivity | 货被多少人点、同场共点 | 静 |
| merchant_structure | **当前这一单的店** 的度数 / PR / 聚类 | 静 |
| video_tower | 人走过的店，64 维视频 DGP 均值 | 动 |
| audio_tower | 同上，音频 | 动 |

视频/音频不是生产塔。挂在商户上，沿进店边聚到人，再当一粒打 W。协议到此为止。

店侧 4 个数静、人侧 64 维塔动：同一张图，粒不同。对照写在 localization 表上，不要收成「所以是视频塔」。

---

## 两步下钻在干什么

只在 moved 的那几粒里：

| 步 | 问 | 不要读成 |
|---|---|---|
| 1 LOGO | 丢掉这一粒，Y 缺口的 R 会不会降 | 「这一粒贡献了 46%」 |
| 2 LOCO | 丢掉这一列，R 会不会降 | 「这一维是根因」 |

Δ≤0：这粒不是这个时钟上的 Y 缺口来源。RF-mass 仍可说明 P(X) 这块很响。两句话都留着。  
PO-risk 是早/晚对 Y 的距离，不是 treatment。量级 ~1e-4 就是日志。

画像动了 ≠ 要为 Y 缺口负责。localization 可以响，两步下钻可以全是 Δ≤0。这正是拆开的意义。

---

## 不要做的

| 做 | 不做 |
|---|---|
| 五粒并排 localization | 再加社区 / 边类型 / 漏斗进这个 prototype |
| quiet 的粒停住 | 把静的粒也拿去 FSDS |
| 两步下钻只开 moved | 三层 hop 板、全表 LOGO |
| PageRank、塔向量当画像坐标 | 当对 Y 的重要性 |

`scripts/tencent_gr/graph_loc_fsds_drill.py`  
结果：`docs/reports/Recsys_Graph_Loc_FSDS.md`
