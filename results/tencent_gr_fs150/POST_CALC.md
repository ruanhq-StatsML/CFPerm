# 买后点击：怎么算

代码：`scripts/tencent_gr/onepass_post.py`（状态机）≡ `block_tables.tab_post_cnv_events`（asof）。  
事件先按 `(ts, act)` 排：同秒 CLK=1 在 CNV=2 前。`t_end = 序列最后一条.ts`。

## 记号

```
clk  = { e | act=1 }
cnv  = { e | act=2 }          # 一行一笔转化
C(t)     = #{ clk.ts | ts <= t }                 # 含同秒
C_i(t)   = #{ clk.ts | item=i, ts <= t }
next     = min { clk.ts | ts > t }               # 严格晚于，同秒不算
next_i   = min { clk.ts | item=i, ts > t }
dt       = next - t          # 没有 next → 空
follow   = t_end - t
W ∈ {5m=300, 1h=3600, 1d=86400, 7d=7*86400}
```

## 每笔转化

```
n_before_W     = C(t) - C(t-W)        # (t-W, t]
n_after_W      = C(t+W) - C(t)        # (t, t+W]
n_same_before_W = C_i(t) - C_i(t-W)
n_same_after_W  = C_i(t+W) - C_i(t)
n_clk_same_before = C_i(t)

lift_W  = n_after_W / (n_before_W + 1)
delta_W = n_after_W - n_before_W

y_post_clk_W  = 1{ dt<=W } ; follow>=W 且没点到 → 0 ; follow<W 且没点到 → NaN
y_post_same_W = 同上，用 next_i

if follow < W:
    n_after_W = n_same_after_W = lift_W = delta_W = NaN

next_clk_same_sess = 1{ dt <= 30*60 }   # 没有 next → 空
```

last-touch（同一趟扫，CLK 记指针）：

```
wo_prior_clk = 该 item 从未 CLK（含同秒、在本条 CNV 前）
dt_item_min  = (t - last_clk[item]) / 60     空 → NaN
dt_any_min   = (t - last_any_clk) / 60
item_within_B = 1{ dt_item_min <= B } and not wo_prior_clk
```

场（>30min 没动作关场；本条计入 pos）：

```
sess_pos         = 本场第几条（含本条 CNV）
sess_clk_before  = 本场、本条之前的 CLK 数
```

滞后（按 cnv_ts 升序，只用已经收口的 y）：

```
lag_post_clk_1d_rate = mean( y_post_clk_1d of 此前各单 | 非 NaN )
n_prior_cnv = 此前转化条数
log1p_price = log1p(max(price, 0))     price 空 → 0
```

## 状态机（正向一遍）

```
clk_any=[], pending=[]
last_any, last_item = None, {}

for (iid, act, ts, price) in evs:          # 已按 (ts, act) 排
    if 距上场 > 30min: 关场

    if act==CLK:
        for p in pending:
            dt = ts - p.t
            if p.next is None: p.next = ts          # 第一次点，任意 item
            if iid==p.item and p.next_i is None: p.next_i = ts
            if 0 < dt <= W: p.n_after[W] += 1       # 同品再 + n_same_after
        clk_any.append(ts); last_any=ts; last_item[iid]=ts

    if act==CNV:
        n_before[W] = #{ clk in clk_any | t-W < ts <= t }   # bisect
        pending.append(新单)                                 # n_after 先 0

# 扫完用 t_end 收 y / 删失；再扫一遍转化算 lag
```

`pending` 不弹出：7d 窗内每来一个点都要给还在窗里的单 +1。只闭 `next`（第一次点）。

## asof（全表，等价）

pandas 3：`on` 键（时间）必须**全局**单调，不能只按 (user,item) 排。`by` 仍按组匹配。

```
next     = merge_asof(cnv, clk, by=user,       on=ts, forward, exact=False)
next_i   = merge_asof(cnv, clk, by=(user,item), on=ts, forward, exact=False)
C(t)     = merge_asof(cnv, clk.assign(cum=cumcount), by=user, on=ts, backward, exact=True).cum
```

## 用户一行

只平均非 NaN：

```
post_clk_W_rate  = mean(y_post_clk_W)
post_same_W_rate = mean(y_post_same_W)
post_lift_1d_p50 = median(lift_1d)
post_clk_dt_p50  = median(dt_next_clk_min)    全空 → -1
```

然后 `user_id` outer join 进 funnel/decay/sess/attr。cross 最后乘。
