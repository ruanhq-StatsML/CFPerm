# 多层漏斗归因

不是 DFS，不是图谱。一层漏斗一问，时间不能倒。

```
曝光  →  点击  →  转化  →  买后点击
 CTR      CVR     路径      续逛 ≠ 转化
```

时钟 unix 秒。相邻 `Δt > 30min` 切断一场。同一 W：`(t−W, t]` 是 X，`(t, t+W]` 是 Y。跟不满不当 0。空 asof 不填 0。

## 1. 曝光 → 点击

Y = 这次曝光之后有没有点（或窗内 CTR = n_clk/n_exp，分母太小丢弃）。
X = 此前漏斗量 / 点率惯性 / 点击半衰期 / 场碎。
任意点击。还没有货、没有单。

## 2. 点击 → 转化

Y = 点了之后有没有买（CVR），或看见会不会买（CTCVR）。
X = 点率、买占轨迹、成交半衰期。
`pay_cnt` 是量，`cnv_share` 是结构，不要混。

## 3. 转化路径（这笔单怎么点过来）

粒 = 一笔 CNV。向后 asof，同秒 CLK 在 CNV 前。

- 任意 last-click：`dt_any`。空 = 没点过就买。
- 同品 last-click：`dt_item`。空 = `wo_prior_clk`。
两问，不是「多模态」。同品闭合在这批数上约 0.1%，当 mix 闸（1%），不当主路径。
转化当下才有 order；货此时已在店上。

## 4. 转化 → 买后点击

粒仍是这笔 CNV。向前 asof。续逛，不是第二笔转化。

```
next_clk = min {clk.ts | ts > t}
y_W      = 1{dt≤W}；follow≥W 没点 → 0；否则 NaN
当场     = dt ≤ 30min
跨场     = dt > 30min
```

`n_before` / `n_after` 同长窗。跟不满 W 的 after / lift = NaN。

## 扫一遍

一条序列按 ts 走：漏斗窗、半衰期、场、向后 last-touch、向前 pending。
Spark UDF 或 sort 后 rolling 同构。过完再算比率，不再扫，更不要 DFS 笛卡尔。

## 分解（每层自己的 Y）

RF-domain = 这层样本谁来了。  
PO-risk `φ=(Y−μ)(W−e)` + 列 LOCO = 这层早/晚差钉在哪列。R 不当检验。  
层与层不要共用一个 Y 去「总归因」。
