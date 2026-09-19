# 渐进 concept · 连续时间 · latency（可接管）

假设：**没有 drastic 跳变**，\(P(Y\mid X)\) 沿连续时钟 \(\tau\) 慢慢转。LOGO 代码好写；难的是 latency 怎么刻画、怎么评估。

Hop 和连续时间是同一条 \(S(\tau)\)（比 vs 水平）。其余好说。要跟 latency 绑在一起的是服务误差 vs 份额误差：`docs/serving_vs_share.md`。

---

## 表 1 · 为什么 last-two hop 在渐进里是错的 WHEN

窗口 \((\tau-W,\tau]\) 上 FrozenRF \(T(\tau)=\mathrm{MSE}(\tau)-E_{\mathrm{ref}}\)。

| 变动 | \(T(\tau)/T(\tau-W)\) | 该看什么 |
|---|---|---|
| 跳变（drastic） | 可以 ≥ 1.5 | hop 标 WHEN |
| 慢走 \(\alpha(\tau)\) 线性 0→1 | \(\to 1\)（相邻两窗几乎一样） | \(T(\tau)\) 相对 **参考水平** / pre-onset MA，不是相对上一窗 |

玩具：video \(\beta\) 对每一行连续旋转，\(P(X)\) 不动。n_new=40 时 hop **从未响**；PO 的 2× 门也从未响；看板一直 `keep`。MMD 全程 ~0，这是 concept 该有的。

---

## 表 2 · Latency 不是一个数（所以 tricky）

时钟：\(\tau\) = 观测序号（有墙钟就换成墙钟）。\(\alpha(\tau)\) 只在 DGP 里知道，评估用，生产没有。

| 钟 | 记号 | 定义 | 渐进里会发生什么 |
|---|---|---|---|
| 组批 | \(L_{\mathrm{batch}}=W\) | 必须攒满一窗才能打分 | \(W\) 大：分数稳、组批慢。\(W\) 小：更像连续、更吵、LOGO 更勤 |
| 检测 | \(L_{\mathrm{det}}\) | 真 onset \(\tau_0\) → WHEN 过程过线 | hop 经常 \(\infty\)。T-MA / Brier-MA 会在 \(\alpha\approx 0.3\) 过线 |
| 定位 | \(L_{\mathrm{loc}}\) | → \(\pi^{\mathrm{PO}}(\mathrm{video})\) 连续两次 ≥ 0.5 | 比检测稍晚或同时；早期份额会乱跳，必须 **连着看** |
| 动作 | \(L_{\mathrm{act}}\) | → 看板离开 keep / 开 train_top | 2× PO×MSE **可以永远不离开 keep**。不能把动作延迟当成检测延迟 |
| 恢复 | \(L_{\mathrm{rec}}\) | 开始更新 → online Brier 回到 oracle 附近 | 太早全量更新可能比冻住更差（见 area 为负） |
| 计算 | \(L_{\mathrm{cpu}}\) | 一次 LOGO/PO 的墙钟 | 本次 ~567 ms / 40 行。到达快于计算 → 永远欠账 |

没有唯一的「onset delay」。跳变没有发生，点延迟不适定。

---

## 表 3 · 该怎么评估（不要只报 hat − τ0）

| 指标 | 为什么要它 | 本次主窗 n_new=40 |
|---|---|---|
| 是否 never | hop / 2× 门会 never；never 不是实现 bug | hop ∞；PO 2× ∞；board keep ∞ |
| delay in obs | \(\times n_{\mathrm{new}}\)，才能跨窗宽比 | T-MA：160 行（α=0.35）；Brier+10%：120 行（α=0.27） |
| α at hat | 发现时机制已经转了多少 | ~0.27–0.35，不是 0，也还没到 1 |
| FA pre-onset | 过灵敏的代价 | 这几条检测器 FA=0 |
| **excess area** \(\sum_{\tau_0}^{\hat\tau}(L^{\mathrm{frozen}}-L^{\mathrm{oracle}})\) | **可评估的 latency**：等的这段时间多付的服务误差 | hop/PO2×/keep：0.207（等到结束）；T-MA / 份额：~0 甚至略负 |
| compute ms / 窗 | 统计延迟之外的管道延迟 | LOGO 中位 567 ms |

area 为负：过早 full update 的 oracle 在走的前半段 **不如冻住**。所以「更早检测 + 立刻全量」不是自动更优。渐进里默认是 **watch、攒份额、Brier 走起来再 train_top**，不是 hop 一响就冻层。

---

## 表 4 · 连续时间里的决策（逐步、无跳）

分数过程（同一 \(D_{\mathrm{ref}}\)）：

- \(T(\tau)\)：冻住 RF 的窗口 MSE − \(E_{\mathrm{ref}}\)
- \(\mathrm{Brier}(\tau)\)：现模型在窗口上的服务误差
- \(\mathrm{PO}(\tau)\)、\(\mathrm{MMD}(\tau)\)：机制 vs \(P(X)\)
- \(\pi^{\mathrm{PO}}_g(\tau)\)：LOGO 份额（连着两窗才算响）

| 过程 | 渐进 concept 该看到 | 下一窗做什么 |
|---|---|---|
| MMD vs ref | 安静 | 不是 X-shift，别 train_stem |
| T(\(\tau\)) 相对 pre-onset 水平 | 慢慢抬 | WHEN：开始 watch，不是 hop |
| Brier MA | 跟着抬，但不到 2× | 确认「房租」在走；仍可先不冻 |
| PO MA 2× | 经常不破 | **不要等这一格才动作** |
| \(\pi^{\mathrm{PO}}(\mathrm{video})\) 连着高 | 定位到 video 塔 | 记下 train_top 候选 |
| 看板 keep | 正常 | 子集/份额可以已经响了 |

动作句（渐进、无 drastic）：

1. hop 不响 → 不解释成「没有 drift」。
2. T-MA 或 Brier-MA 过 pre-onset 水平 → **watch**。只推理，全模型先不冻。
3. \(\pi^{\mathrm{PO}}\) 在某塔连着高、MMD 安静 → 这是 concept 定位，候选 `train_top` 那座塔。
4. 开 `train_top` 的确认仍是服务误差在走，不是 PO 2×。过早全量可能 area 为负。
5. 计算跟不上到达 → 先丢掉 LOGO 降频（每 k 窗打一次份额），WHEN 只用 FrozenRF T + Brier，这两条便宜。

---

## 表 5 · 窗宽：组批延迟 vs 检测延迟

同一条 640 行流，onset 后 α 线性走。

| n_new | 组批（行） | hop delay（行） | T-MA delay（行） | Brier+10%（行） | PO 2× | hop FA |
|---|---|---|---|---|---|---|
| 20 | 20 | 220 | 160 | 160 | ∞ | 0 |
| 40 | 40 | ∞ | 160 | 120 | ∞ | 0 |
| 80 | 80 | ∞ | 80 | 80 | ∞ | 0 |

窗越大，hop 越死（相邻两窗更像）；T-MA 的 **观测延迟** 可以变短，因为一窗里已经积了足够的 α。但你 **先付了 W 行组批**。总延迟 ≈ 组批 + 检测，不是其中一项。

n_new=20 时 hop 终于响了（delay 220），那是窗太碎、相邻比偶尔爆，不是更接近连续时间的胜利。连续时间应该让 hop 更安静，不是更吵。

---

## 表 6 · 和冻层看板、LOGO 怎么接

| 已有 | 渐进连续时间里改读 |
|---|---|
| last-two hop = WHEN | 只保留给 drastic。渐进 WHEN = T(\(\tau\)) / Brier 相对参考水平 |
| PO×MSE 2× → freeze | 渐进经常永远 keep。freeze 不是 WHEN |
| LOGO 两层比例 | 定位塔；要连着两窗。π_loss 仍可能全 0 |
| watch | **渐进的默认动作** |
| train_top / train_stem | 份额稳定 + 服务误差在走 之后才动那一座塔 |
| regret / online Brier | 用 excess area 评估「等了多久」，不要只报 hat |

图：`results/gradual_concept_latency/gradual_clocks.png`（α 直线走，T 抬，hop 没有点，π_PO 后期锁在 video）。数字：`results/gradual_concept_latency/TABLES.md`。

```bash
PYTHONPATH=Python/src:. python3 -m unittest tests.test_gradual_latency
PYTHONPATH=Python/src:. python3 scripts/run_gradual_concept_latency.py
```
