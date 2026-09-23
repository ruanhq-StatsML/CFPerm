# PO-Boost：连续时间逻辑 + Post-training elaboration roadmap

> LaTeX：[`AGOD_PO_Boost_ContinuousTime_Roadmap.tex`](AGOD_PO_Boost_ContinuousTime_Roadmap.tex)  
> 场景增益：[`AGOD_PO_Boost_PostTrain_Scenarios.tex`](AGOD_PO_Boost_PostTrain_Scenarios.tex)

---

## 连续时间为什么显著

后训练是 **流** \((X_t,Y_t)\)，不是一次性重加权：

```text
PO_t  →  滤波(L=EMA, S=ΔPO, α)  →  执行器 u_{t+1}  →  θ_{t+1}
         低通 / 导数 / 因果保持（只作用下一窗）
```

| 原语 | 离散 | 连续读法 | 帮助 |
|---|---|---|---|
| EMA \(L\) | 慢性残差 | 低通 | freeze 不抖 |
| \(\Delta\mathrm{PO}\) \(S\) | 尖峰 | 短时导数 | Acc 掉之前倾倒 step |
| \(u_{t+1}=g(t)\) | 下一窗才动 | 零阶保持 | 无同窗泄漏 |
| RFPerm reject | 事件时刻 | 事件时间 | 只在报警抬 \(w\propto\sqrt{\mathrm{PO}}\) |

**增益故事：** 预期（更小 \(T(\mathrm{Acc}^\star)\)）+ 稳定（少 thrash）+ 选择性花费（calm 不抬权）。饱和流（COCO）连续控制器无可重分 → 诚实 noop。

### 不是 ε-greedy：传感器触发的 OOD

- 默认 \(\alpha=\mathrm{Softmax}(\mathrm{score}/\tau)\)（Boltzmann），**不是**以 ε 掷硬币换臂。
- 唯一 ε 味：`po_budget` / `structured_epsilon_alpha` —— \(f=\varepsilon/|M|\) 地板保冷塔温感，**禁止**随机改 α。
- OOD 响应链：Sense(PO,MMD,ΔPO,…) → Judge(drift-vs-noise) → Filter(L/S) → Actuate(freeze/dump/LR；reject 才 √PO) → 只作用 \(t{+}1\)。
- Justify：平静窗随机探会打穿 FLOPs+Acc 门；真漂移才该偏斜。

---

## Post-training elaboration roadmap

| Phase | 做什么 | 退出标准 |
|---|---|---|
| **R0** 已锁 | S1–S6、`po_fuse` 接线、ROI SQL 骨架 | — |
| **R1** 进行中 | `continuous_gain_metrics`：Jaccard / \(T(\mathrm{Acc}^\star)\) / cumFLOPs@★；已接 compare | Affec 跑满 fuse vs soft vs gated 填表 |
| **R2** 代码落地 | `realize_step_alloc` + block schedule；compare 默认 `step_mode=per_mod`；`--step-mode shared` 作 LR-only 对照 | Affec：per_mod vs shared 的 \(T(\mathrm{Acc}^\star)\) |
| **R3** 代码落地 | `po_iptw_weights` + `stream_reject_proxy`；compare 因果 \(w_{t+1}\)；calm \(w=1\) | Affec：reject 窗 next-MSE；换真 RFPerm flag |
| **R4** | 扫 \(\rho,\omega,\theta_{\mathrm{fr}}\) Acc@FLOPs Pareto | \(M\ge3\) / \(M=2\) 默认配置卡 |
| **R5** | 填 ROI SQL + ship gate | 业务只看 SQL，不报推理延迟 |

**不做：** α 反传当主路由；动 Drill；审出 NOT_READY 时堆 feature。

```text
R1 measure → R2 real steps → R3 event-time weights → R4 schedule → R5 ROI
```
