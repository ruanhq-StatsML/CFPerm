# 路径、指针与 merge_asof

转化粒特征：怎么刻画、怎么评估。自己落地时按这个做即可。  
事件 `(user, item, act, ts)`，`act∈{0=曝光, 1=点击, 2=转化}`。先按 `(ts, act)` 排：同秒点击在转化前。`t_end` = 该用户序列最后一条 ts。

本文不是因果、不是 CATE。asof 只是「以某条事件为锚，找另一类事件的最近一次」。

---

## 1. 原语

```
merge_asof(left=锚, right=候选, by=组, left_on=ts, right_on=ts,
           direction=backward|forward, allow_exact_matches=exact)
```

| 参数 | 含义 |
|---|---|
| `by=user` | 同一个人 |
| `by=(user, item)` | 同一个人同一商品 |
| `backward` | 锚时刻之前（含/不含同秒由 exact 定）最近一条候选 |
| `forward` | 锚时刻之后最近一条候选 |
| `exact=True` | 同秒算匹配 |
| `exact=False` | 同秒不算；forward 时「下一次」必须严格晚于锚 |

pandas 3：时间键必须**全局单调**，不能只按 `(user, item)` 排。`by` 仍按组匹配。

指针扫一遍和 asof 同构：

- backward：正向扫，候选更新指针，锚出现时读指针。
- forward：锚推进 `pending`，候选出现时把还没闭上的锚闭上。

窗内计数用累计：

```
C(t) = #{ 候选 | ts ≤ t }
n_(t-W, t] = C(t) − C(t−W)
n_(t, t+W] = C(t+W) − C(t)
```

`C(t)` 本身也是 `merge_asof(锚, 候选.cumcount, backward, exact=True)`。

---

## 2. 已经成立的一块：最近一次点击 → 这一次转化

这是路径。锚 = 转化，候选 = 点击，`direction=backward`。

```
last_item = merge_asof(cnv, clk, by=(user, item), backward, exact=True)
last_any  = merge_asof(cnv, clk, by=user,         backward, exact=True)
first_item: 每个 (user, item) 只留最早一次点击，再 backward asof
```

落地字段（转化一行）：

| 字段 | 计算 | 在问什么 |
|---|---|---|
| `dt_item` | `(cnv_ts − last_item.ts) / 60`，空则空 | 这件商品点完多久下单 |
| `dt_any` | `(cnv_ts − last_any.ts) / 60` | 全局上次点完多久下单 |
| `wo_prior_clk` | `last_item` 空 | 这件从未点过就买 |
| `item_within_B` | `dt_item ≤ B` 且非空 | 5m / 30m / 1h / 1d 桶 |

指针（正向扫）：

```
last_any, last_item = None, {}
if CLK: last_any = ts; last_item[item] = ts
if CNV: dt_any = ts − last_any; dt_item = ts − last_item[item]
```

同一对 `(last_clk, cnv)` 从点击侧看就是「这次点击之后的下一次转化」：

```
next_cnv = merge_asof(clk, cnv, by=(user, item), forward, exact=False)
```

和 `last_item` 是同一批配对，只是粒从转化换成点击。路径分析停在转化粒；若做「点了会不会买」才用点击粒。

---

## 3. 对偶：这一次转化 → 下一次点击

锚 = 转化，候选 = 点击，`direction=forward`，`exact=False`。

```
next     = merge_asof(cnv, clk, by=user,         forward, exact=False)
next_i   = merge_asof(cnv, clk, by=(user, item), forward, exact=False)
dt       = next.ts − cnv_ts          # 没有 next → 空
follow   = t_end − cnv_ts
```

**只用跟满窗的单。** 窗 W∈{5m, 1h, 1d, 7d}：

```
y_post_clk_W  = 1{ dt ≤ W }
              ; follow ≥ W 且没点到 → 0
              ; follow < W 且没点到 → 丢掉（NaN，不当 0）
y_post_same_W = 同上，用 next_i
```

窗内已点到即使 `follow < W` 也算 1（事件已经发生）。没点到且没跟满，不知道，不能当负例。

`trans_cnv_to_exp` 是邻接 Markov（下一条事件），不是窗内点击，不要替代这一块。

---

## 4. 还有哪些 asof 值得做

同一原语，换锚 / 候选 / by / 方向。建议都做成转化粒或对应锚粒中间表，不要先堆用户一行。

| 锚 | 候选 | 方向 | by | 字段直觉 |
|---|---|---|---|---|
| cnv | clk | backward | item / user | 路径（§2，必做） |
| cnv | clk | backward | item，只留 first | 从第一次点到下单 |
| **cnv** | **clk** | **forward** | **user / item** | **买后下一次点（§3，必做）** |
| cnv | cnv | forward | user / item | 下一次买 / 复购同品 |
| cnv | exp | forward | user / item | 买后又看没看 |
| cnv | exp | backward | item | 最后一次曝光到下单（无点击路径） |
| clk | cnv | forward | item | 点完会不会买（点击粒） |
| clk | clk | forward | user | 点完还会不会点（较少单独做） |

场内版：匹配后再加 `dt ≤ 30min`（或先按 session id 切开再 asof）。用来区分当场续点 vs 隔场回访。

同长窗 before/after 计数（任意 / 同品）是 asof 的累计形式，不是新原语。

用户一行只对未删失转化做 `mean` / `median`，再 join。交叉最后乘。

---

## 5. 路径怎么刻画、怎么评估

路径 = 转化时刻已经知道的 backward asof，**不准含任何买后信息**。

刻画（转化粒，先不要模型）：

- `P(wo_prior_clk)`、`dt_item` / `dt_any` 的分位数、各 `within_B` 占比
- 空值单独一档，不要填 0（0 分钟 = 刚点完就买）
- 和当场一起切：`within_5m` × `sess_clk_before>0` = 热场短路径；`wo_prior_clk` × `sess_clk_before=0` = 冷场直达

评估（预测买后点击时）：

- 样本：`follow ≥ W` 的转化
- `Y = y_post_clk_W` 或 `y_post_same_W`（更稀、更干净）
- 对照：只拿买前点击量 `n_clk_before_W`；再加路径；看 AUC/AP 增量
- 短窗 Y（5m/1h）路径和当场才应该有用；1d/7d 任意点击多半是「人还在」，路径增量会很小
- 不要用 F-score 对着 `future_cnv` 选这些列

---

## 6. 当场：`sess_clk_before`

不是 asof。一场 = 相邻事件间隔 ≤ 30min；超过就关场。

```
sess_pos         = 本场第几条（含本条转化）
sess_clk_before  = 本场、本条之前的点击数
```

指针：

```
若 ts − prev_ts > 30min: sess_pos = 0; sess_clk = 0
sess_pos += 1
CLK 时 sess_clk += 1
CNV 时读出 sess_pos, sess_clk_before = sess_clk
```

**刻画什么：** 下单这一刻场热不热。`sess_clk_before=0`：这场还没点过就买（曝光即买 / 新场第一条附近）。`≥1`：买之前这场已经动过手。它是「买后立刻又点」的主混杂——场没关，下一击几乎是续逛。

**怎么评估：**

1. 分布：转化上 `sess_clk_before` 的直方图；对比随机事件上的同一量（转化是否更爱落在热场）。
2. 和 Y 的关系要按窗拆：  
   - 与 `y_post_clk_5m` / `next_clk_same_sess` 强相关 → 当场续点，不是购买带动。  
   - 控住 `n_clk_before_1d` 之后仍预测 `y_post_clk_1d` → 泄漏的是用户热度，不是场。
3. 消融：volume → +sess → 看短窗 Y 和长窗 Y 分别涨多少。短窗涨、长窗不涨，当场就做对了。
4. 不要和用户级 `sess_n` / `sess_bounce_rate` 混。那两个是整段历史的场统计；`sess_clk_before` 是**这一单所在场**。

---

## 7. 滞后

```
lag_post_clk_W_rate = mean( 此前各单的 y_post_clk_W | 非 NaN )
n_prior_cnv         = 此前转化条数
```

必须 `shift(1)`：本单的 Y 不能进本单的 X。此前单若当时没跟满窗，那一单不进均值。

刻画：用户「买完还点」的倾向，给下一单用。  
评估：volume 之上再加 lag，1d/7d Y 通常这是最大增量；短窗 Y 不该主要靠 lag（那是当场的事）。

---

## 8. 评估协议（自己做特征时按这个收）

**粒：** 转化一行。用户任务（如 `future_cnv`）只并历史倾向，不要拿本单买后量去预测用户未来转化。

**删失：** 预测窗 W 的模型，只用 `follow ≥ W` 的单；或窗内已观察到正例的单。负例必须跟满。

**X：** 只有 `cnv_ts` 已知的量——路径、买前窗计数、当场、价格、滞后。  
不准进：`n_after`、`lift`、`dt_next`、`y_post_*`。

**切分：** 按 `user_id` 切，不要随机切转化（lag 和用户热度会漏）。

**先刻画后模型：**

1. 未删失集合上的 `P(y_post_clk_W)`、`P(y_post_same_W)`  
2. `n_after` vs `n_before`（任意 / 同品）；`P(dt ≤ 30min | 有 next)`  
3. `E[X | Y=1]` vs `E[X | Y=0]`：路径、当场、lag 分列  
4. 再消融 AUC/AP：`volume` → `+path` → `+sess` → `+lag`

**读增量：**

| 消融涨在哪 | 说明 |
|---|---|
| 只有 volume | 爱点的人继续点 |
| +path 在短窗 | 决策时长 / 漏点在管当场余热 |
| +sess 在短窗、长窗不动 | 当场续逛，不是买完回访 |
| +lag 在长窗 | 用户稳定倾向 |
| 同品 Y 正例极少 | 「买完还看这件」事件本身稀，先报基数再建模 |

不要：全样本 F-score、删失当 0、用邻接转移代替窗、把用户级 bounce 当成当场。

---

## 9. 建议落地顺序

1. 转化中间表：§2 路径 asof（item / any / first）  
2. 同一张表：§3 forward asof + `follow` + `y_W`（先 1h/1d）  
3. 累计 `n_before_W` / `n_after_W`（after 只在跟满时留）  
4. 扫场：`sess_pos`、`sess_clk_before`  
5. 按转化时间 `lag`  
6. 按 §8 出一张刻画表 + 一组消融，再决定哪些进用户表

中间表 polish 后再 `user_id` merge。交叉最后做。
