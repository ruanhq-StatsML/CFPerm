# 连续时间团伙发现：我们的方法怎么接

> Drill 门不变：**rank-average(MMD, PO, cmean)** + Tail/CV/\(R\) 停钻。  
> 团伙 = **图上共现/共变的实体块**（检索支撑），**不是**因果效应子群 / ATE。  
> 图算法（Leiden 等）只产 **upstream key**，**不**用模块度 \(Q\) 当 Drill 门。

相关：[`Graph_CT_AD_Interface.md`](Graph_CT_AD_Interface.md) ·  
[`Graph_Localize_Algorithm_Elaboration.md`](Graph_Localize_Algorithm_Elaboration.md) ·  
数据：[`../../data/tencent_subset/README.md`](../../data/tencent_subset/README.md)

---

## 1. 一句话

| 问题 | 我们的答案 |
|---|---|
| **什么时候**出现可疑团伙活动？ | 连续 \(\Delta t\) 图流上 OnlineRFPerm：\(T_t^{\mathrm{MMD/PO/cmean}}\) → \(R_t\) |
| **哪一伙**在漂？ | Reject 后 Drill：实体键 = **社区 / 共现块**（或 merchant），分数仍是 MMD+PO+cmean |
| **团伙里 tip 是什么？** | \(K^\star\) 上一次 FSDS + 正向/负向 |
| Leiden / 大规模图算法？ | 只生成 `community_id` 钥匙；锁仍是我们的 Drill |

```text
(u,i,t) 边流
  → G_t（每 Δt 一张图）
  → [可选] Leiden_γ(G_t) → community_id
  → OnlineRFPerm(MMD, PO, cmean)        # WHEN
  → if R_t: Drill on community|merchant→user→order  # WHERE (团伙)
  → one FSDS on K* + sign(Δȳ), sign(δ_j)  # WHAT + 方向
```

---

## 2. 「团伙」在协议里的定义（讲武德）

在广告/推荐日志上，团伙 ≈ **行为上绑在一起的用户–物品–商户块**：

- 同社区内边密、跨社区疏（结构）  
- 或 W1→W2 / \(t\!-\!1\!\to\!t\) 上 **联合分布一起漂**（我们真正测的）

我们交付的是：

\[
K^\star_t=\{\text{entities whose joint shift explains } \delta(S_t)\}
\]

读成「这一时刻异常活动的载体块」。**禁止**读成「发现了因果作弊群体 / ATE 子群」。

报告条：

```text
[Gang-support] key=community|merchant; K*=...; n_entities=...; sign(Δȳ)=+/-/flat; tips={...}
```

---

## 3. 连续时间怎么跑（WHEN）

与 CT-AD 相同，三路并行：

\[
T_t^{\mathrm{MMD}}=\widehat{\mathrm{MMD}}^2(X_{t-1},X_t)-e_{\mathrm{MMD}},\quad
T_t^{\mathrm{PO}}=\mathrm{PO\text{-}gap}(t\!-\!1\to t)-e_{\mathrm{PO}},\quad
T_t^{\mathrm{cmean}}=\|\mu_t-\mu_{\mathrm{ref}}\|_2-e_{\mathrm{cm}}.
\]

- \(X_t\) = 窗 \([t,t+\Delta t)\) 内 `feature_engineer` 的 `grid`（或团伙级聚合特征）  
- \(R_t=1\) 才进入团伙 Drill；平时只记 \(T_t\) / Lead / AR  
- \(\Delta t\)：日级盯「团伙活动波次」；小时级盯突发

**团伙专用监测可选第四标量（不进门，只诊断）：**

\[
T_t^{\mathrm{comm}}=\widehat{\mathrm{MMD}}^2\big(X_{C^\star}^{t-1},X_{C^\star}^{t}\big)
\]

即「当前最大嫌疑社区」内部的分布差——仍用 OnlineRFPerm 骨架，**不**用 \(Q_t\) 报警。

---

## 4. 团伙键从哪来（WHERE 的 entity_key）

### 4.1 固定业务键（已有）

`merchant → user → order`：商户代理团伙 / 店铺簇。  
映射率低时大量 singleton shop → 键弱。

### 4.2 Leiden 多分辨率（推荐接法）

在切片图 \(G_t\)（或滚动 \(G_{t-h:t}\)）上：

1. 构图：二部 `user–item` 或投影 `user–user`（共点 item）/ `item–item`（共点 user）；边权 = 共现或 click  
2. Leiden 分辨率 \(\gamma_1>\gamma_2>\cdots\) → 粗→细 `community_id`  
3. 把 `community_id` 当作与 merchant **并行**的 `drill_key`  
4. 对每个社区 \(c\) 算  
   \(r_c^{\mathrm{blend}}=\mathrm{rank\text{-}avg}(\mathrm{cmean},\mathrm{MMD},\mathrm{PO})\)  
5. Drill 产品式照旧：\(\|\delta\|\ge\varepsilon\)、Tail/CV、\(R(K)\le R^\star\)、mass/\(\pi\)  
6. 停层 \(K^\star\) = 尖端社区（或社区∩商户）→ 一次 FSDS

**换钥匙，不换锁：** Leiden 产出层；是否下钻仍看 cmean/MMD/PO 谱。

### 4.3 与「社区发现」论文的差别

| 经典团伙/社区发现 | 我们 |
|---|---|
| 优化 \(Q\) / densest subgraph 本身就是目标 | \(Q\) **只**帮切块；目标是 **时段分布漂移定位** |
| 静态一张图 | **连续** \(G_t\) + reject 时刻再切 |
| 输出社区列表 | 输出 \(K^\star\) + tip 特征 + 正向/负向 |

---

## 5. 连续时间团伙生命周期（可落地的几件事）

1. **出现**：某 \(\gamma\) 下新社区在 \(t\) 首次进入 Top-\(r\) 且 \(R_t=1\)  
2. **加剧**：同一 `community_id`（需跨时刻对齐，见下）的 MMD/PO/cmean 持续升高  
3. **分裂/合并**：多分辨率下粗社区残差 \(R\) 变大 → Drill 到更细 \(\gamma\)  
4. **消退**：\(R_t=0\) 且该块 \(r_c\) 跌出 Tail → 停止交付  
5. **跨键核对**：同一 reject，比较 `merchant` vs `community` 的 \(\mathrm{eff}=\gamma/\alpha\)——哪个更尖用哪个讲故事（可并列报）

**社区跨时刻对齐（工程）：**  
Jaccard / Hungarian 匹配 \(C_t\) 与 \(C_{t-1}\)；或固定滚动图上跑一次 Leiden，切片只继承标签。对齐失败就当「新团伙候选」，勿硬追 ID。

---

## 6. 数据接口（接 TencentGR 子集）

```text
https://github.com/ruanhq-StatsML/CFPerm/tree/cursor/graph-localize-elab-abce/data/tencent_subset
```

```python
# 1) 时间片边
panel = feature_engineer(ROOT, t, t + delta, max_users=20000, terminal_action=1)
edges = panel["edge_df"]   # user_id, item_id, last_ts, ...
grid  = panel["grid"]

# 2) 可选：构图 → Leiden → community_id（伪）
#    A_ij = covisit / click co-occurrence on [t-h, t)
#    community_id = leiden(A, resolution=gamma)

# 3) 挂键后 Drill（分数仍是 MMD+PO+cmean）
#    score entities by rank_average(cmean_l2, mmd2, po_tau2)
#    nest: community → user → order   (or merchant → ...)

# 4) R_t from OnlineRFPerm; only then run step 3
```

边契约给图库：`(src=user, dst=item, ts, action, y=click)` + 可选 `community_id`。

---

## 7. 决策表（团伙版）

| \(\|\delta\|\) | 社区谱 \(r_c\) | 动作 |
|---|---|---|
| 小 | 平 | 无团伙报警 |
| 大 | 平 | 整层漂；FSDS 讲特征，**不**点名团伙 |
| 大 | 尖 | **团伙 Drill** → \(K^\star\) → FSDS |
| 尖但 \(R\approx 1\) | — | 假尖端；退回层 / 换 \(\gamma\) |

业务 cap：可规定最细交付到 `community` 或 `user`，禁止点到单订单。

---

## 8. 不要做的事

- 用 \(Q\) / densest / 「更显著」当 Drill 或 Fire 门  
- 把 \(K^\star\) 写成 ATE / 「作弊因果群」  
- 每个 \(\Delta t\) 无 reject 也全量 Leiden+FSDS（算力炸、故事碎）  
- 忽略 \(\mathrm{sign}(\Delta\bar y)\) / tip 符号  

---

## 9. 建议落地顺序

1. 日级 \(R_t\)（已有三路标量）+ reject 时 **merchant** Drill（已有脚本）  
2. 同 reject 加 **Leiden \(\gamma\)-scan** 社区键，比 \(\mathrm{eff}/\pi\)  
3. 社区跨日对齐 + 「团伙出现/加剧/消退」时间线  
4. 可选 tip-cmean×confirm 做团伙活动尖峰 timing  

---

## 10. Takeaway

连续时间团伙发现 = **CT 监测（何时）** + **结构键（哪一伙）** + **MMD+PO+cmean Drill（是否尖、停哪）** + **一次 FSDS（tip）**。  
Leiden 换的是「团伙候选钥匙」；锁始终是我们的分布漂移协议。
