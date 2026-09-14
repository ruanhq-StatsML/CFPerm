# 连续归因逻辑链

时钟是事件 unix 秒，不是自然日。相邻 `Δt > 30min` 切断一场。
归因是 asof 近邻，没有 GT，不是 CATE。任意点击 vs 同品是两问。

## 时间顺序（不能倒）

```
上架  item → merchant
  → 曝光
  → 点击（任意 / 同品）
  → 转化 此刻发 order  order → item → merchant
  → 买后点击  当场续逛 (dt≤30min) 或 跨场回访
```

货先挂在店上，行为才走得进店。单只在 `cnv_ts` 存在。不要先画图谱再补时间。

## 向后：这笔成交怎么点过来

对每笔 CNV：

- 任意 last-click：`dt_any`，空 = 没点过就买，**不填 0**
- 同品 last-click：`dt_item`，空 = `wo_prior_clk`
- 同秒先 CLK 后 CNV

同品漏斗闭合在这批数上几乎不发生（约 0.1%）。当日报 mix，闸在 1%，不要当模型主路径。

## 向前：买完还逛不逛（≠ 转化）

同一把尺子 W：

```
n_before = #{clk in (t−W, t]}
n_after  = #{clk in (t, t+W]}     # 跟不满 W → NaN，不当 0
next_clk = min {clk.ts | ts > t}  # 严格晚于
y_W      = 1{dt≤W}；follow≥W 没点到 → 0；否则 NaN
```

`dt≤30min` 是当场续逛，否则跨场。续逛不是转化。

## 特征（one-pass / rolling / Spark UDF 同构）

一条用户序列按 ts 扫一遍，同时积：漏斗窗、半衰期、场、last-touch、买后 pending、Markov。
过完再算比率。user/item 计数按进店人聚到 merchant。不做 DFS，不把店名 TF-IDF 进 X。

## 分解（两层，不是预报榜）

1. **RF-domain**：P(W|X)，谁来了（人 / tenure）。
2. **PO-risk** `φ=(Y−μ)(W−e)`：早/晚对这个 Y 的差。R 不当检验。
3. **LOGO**：整跳拿掉（user / order / merchant）。
4. **LOCO**：一列拿掉。impurity ≠ 漂移源。
5. **单位**：店 mix vs 店内。SKU 份额碎 → 停在店/人，不要点名 item。

这批订单粒：Y 缺口几乎全是店内（晚来的人 `n_prior_cnv` 更浅），不是换了一批店，更不是某件货。

## 实现

状态机 = `onepass_feats` + `onepass_post`。表上等价于 sort `(user, ts)` 后 rolling / merge_asof。Spark 同一套 UDF，不必上图。
