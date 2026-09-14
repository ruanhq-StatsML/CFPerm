# FSDS 怎么读这批特征（不是把名字翻译成故事）

OOD VIMP 只定位 **哪块特征在扛两批的差**，不是因果。
- **RF-domain** P(W|X)：后来的人长什么样（队列 / P(X)）
- **PO-risk** φ=(Y-μ)(W-e)：哪块和「早/晚对 Y 的差」绑在一起
- 看完 VIMP 再看 **E[X|Q1] vs E[X|Q4]**（README 里的 post-hoc 分组）

这根 `_seq_t_end` 时钟：Q1 pos=0.31，Q4=0。LOGO 在中位切上全是负 Δ——**没有家族是概念漂移源**。能讲的是：后来的人画像变了（RF-AUC 0.75/0.88），成交热度被抽干。

---

## pay_cnt ≠ 常态偶尔

| | 公式 | FSDS | Q1 → Q4 |
|---|---|---|---|
| **pay_cnt** | n_cnv 次数 | PO-mass **0.001**，几乎不是 OOD 驱动 | 3.76 → 0.47 |
| **life_cnv_share** | n_cnv / 总事件 | PO-top1 | 0.054 → 0.006 |

次数 = 买了几单（量）。share = 轨迹里买占多大（结构）。  
「常态还是偶尔」是 **share**。pay_cnt 进 cross（`x_pay_cnt__active_days`）只是密度，单列不被 PO 认。

---

## 衰减：为什么不是硬截 7 天

`dec = Σ 2^{-age/hl}`。硬窗把第 8 天砍死；半衰期让「昨天的点」比「上周的点」重。

FSDS：decay 的 **PO-mass 0.14（中位）/ 0.23（Q1Q4）**，高于 RF-mass 0.05/0.08——它更像绑在 Y 差上，不像纯谁来了。Q1Q4 的 PO-top1 是 `dec_hl7d_dec_cnv`。

Q1→Q4：`dec_hl3d_dec_clk` 0.29→0.01，`dec_hl7d_dec_cnv` 0.57→0.04。  
业务：后来的人 **最近既不点也不买**（热度掉光），不是「历史单量还在、只是窗口截错了」。

hl1d/3d/7d 同时进 PO-32：短窗看还在不在逛，长窗看成交热度掉没掉。

---

## 场次：30min 关场在定位什么

一场 = 连续逛；>30min 没动就关。bounce=看一条走；depth_cnv=一场成交几下；cnv_sess_rate=多少场真正成交。

FSDS：session mass 两边都 ~0.07。PO-top 有 `sess_depth_cnv_mean`、`sess_cnv_sess_rate`。  
Q1→Q4：depth 0.067→0.007，cnv_sess_rate 0.063→0.007。

后来的人不是「场次切错了」，是 **场里买不出来**。30min 是产品会话惯例；FSDS 认的是「场内成交深度」这个聚合，不是 gap 秒数本身。

（150 板把 `sess_n`/`bounce` 原子丢掉、只留 cross，这是 F 选的蠢，不是业务定义。）

---

## 归因：asof 刻画路径，FSDS 说它不是漂移源

三张 asof：同商品 last-click、任意 last-click、同商品 first-click。空 = 没点过就买。

FSDS：**attr RF-mass 0.01，PO-mass 0.03，LOGO Δ 为负**。  
不要讲「决策时长变了所以后来不买」。时长不是这根时钟的 shift driver。

PO-32 里的 `attr_anyclk2cnv_min_std`：Q1 均值 3036min vs Q4 174min——早期有转化的人决策时长又长又散；Q4 几乎没转化，std 被空值/少转化拖下去。这是 **有没有单可归因**，不是「后来下单更快」。

`x_life_ctcvr__attr_anyclk2cnv_min_p50` 只在 Q1Q4 的 PO-top 出现：会买 × 犹豫，是成交还在时的路径交互，LOGO 仍不认 attr 家族。

---

## Cross：业务就是「两件事同时成立」

原子乘：会买 × 场碎 / 爱点 × 来得勤 / 会买 × 还在点。

FSDS：cross **PO-mass 最大**（0.42 / Q1Q4 **0.51**），Q1Q4 LOGO share 也是 cross 头名。  
单看 pay_cnt 或 bounce 都不扛两批的差；**焊在一起**才进 OOD VIMP。

读法（方向用 Q1→Q4 均值，不是因果）：

| cross | Q1→Q4 | 后来的人 |
|---|---|---|
| ctcvr × bounce | 0.056→0.005 | 不再是「会买但场碎」，是两边都没了 |
| ctr × active_days | 1.61→0.31 | 不爱点、来的天也少 |
| pay_cnt × active_days | 137→16 | 单量密度掉 |

RF-top 另讲 P(X)：`ui_only_exp_share` 0.91→0.98，`trans_exp_to_exp` 0.83→0.97——后来的人 **只看不点、连着曝光**。这是队列/删失画像，不是转化机制 hop。

---

## 还有什么

- **Markov**：RF 爱 `exp→exp`（后来的人连着看）；PO 中位里 `cnv→exp` 较低。买完又逛是早期成交用户的痕迹。
- **ui_only_exp / exp_per_item**：RF-top1/2。后来 = 曝光堆着、点不出来。
- **money 单列**：PO 几乎 0。别用 ARPU 讲这根时钟。

一句话：FSDS 说这根时间轴上，**后来的人是「只曝不点、近期热度为 0」的队列**；share/衰减点击/场内成交/cross 是在定位这个差，不是在证明 last-click 或 pay_cnt 变了。
