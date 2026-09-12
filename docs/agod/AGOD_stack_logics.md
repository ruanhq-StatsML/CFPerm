# 不是 MoE：有限步长和其他 stacking 逻辑

结论先说：**这不是 Mixture of Experts。** MoE（Jacobs et al. 1991；Switch Transformer 也是）是 **输入条件门** `π(x)=softmax(g(x))`，expert 按 x 的区域分工、和门一起训练。这里的 π 是 **全局（或慢变）的票组合权**，用 holdout / one-step-ahead 损失更新——Super Learner / forecast combination，门不看 x。

`softmax(s/τ)` 长得像 MoE 的门，但它的 logit 是 **expert 分数 s**，不是 `Wx`。没有 x，就不是 MoE。

---

## 有限步长（线性 gain 的下一项）

沿堆叠方向走真实一步：

```
L(θ − η Gπ) = L − η πᵀs + (η²/2) πᵀ H π + O(η³)
```

Gauss–Newton 下 `H ≈ GᵀG = R`。丢掉余项，单纯形上的 stacking 是

```
min_{π∈Δ}  −πᵀs + (η/2) πᵀ R π
```

| η | 行为 |
|---|---|
| → 0 | 线性项主导 → **顶点**（离散 OSL / Frank–Wolfe） |
| 中等 | 精度–方差权衡，内点组合 |
| → ∞ | 二次项主导 → **min-var** `min πᵀRπ`（趋近 MGDA） |
| 训练 LR ~ 3×10⁻³ | `quad/lin = η/2 ≪ 1`，真实 SGD 步几乎就是线性 gain |

所以：**不能指望实际学习率带来的曲率帮你组合。** 虚步长 η 是 meta-loss 上「敢不敢离开 winner-take-all」的旋钮，和训练 LR 不是同一个东西。

### 顺序 vs 联合：有限步长 *就是* GLS

把 (η, π) 一起放进二阶 Taylor：

1. **顺序**：先钉死训练 LR `η≈3e-3`，再选 π → 线性项主导 → **顶点**。
2. **线搜索（π 固定）**：Cauchy 步 `η* = (πᵀs)/(πᵀRπ)`。
3. **联合**：把 η* 代回去得到 `min_π −(πᵀs)² / (2 πᵀRπ)`，即最大化 Rayleigh / SNR。无约束解 `π ∝ R⁻¹ s`，就是方向匹配 / 静态 GLS。

有限步长不是另一种估计器，它是 **离开顶点的机制**。联合优化步长和权重，回到已经写过的 GLS；只用训练 LR 走一步，永远走不出离散 OSL。

Taylor 在 η=3e-3：linear=3.00e-03，quad/lin=1.50e-03（线性区）；η=1 才离开线性区。本例联合 Cauchy `η*=0.87`（图中红点线）。

η-path（s_a=0.85, s_b=0.45, s_c=0.10，a–b 相关 0.8）：

| η | π | N_eff | top |
|---|---|---:|---|
| 0.01 | a=1.00, b=0.00, c=0.00 | 1.00 | a |
| 0.0417 | a=1.00, b=0.00, c=0.00 | 1.00 | a |
| 0.174 | a=1.00, b=0.00, c=0.00 | 1.00 | a |
| 0.724 | a=1.00, b=0.00, c=0.00 | 1.00 | a |
| 3.02 | a=0.64, b=0.00, c=0.36 | 1.92 | a |
| 12.6 | a=0.36, b=0.20, c=0.44 | 2.85 | c |

互补单位梯度（holdout ∥ (e₀+e₁)/√2，R=I）：线性 / 训练 LR / Frank–Wolfe 都塌到一个轴；η=1 和联合 GLS 才是 ½–½。

| logic | π |
|---|---|
| vertex / train / FW | a=1.00, b=0.00, c=0.00 |
| finite-step η=1 | a=0.50, b=0.50, c=0.00 |
| joint GLS | a=0.50, b=0.50, c=0.00 |

---

## 其他逻辑（同样不是 MoE）

| 逻辑 | 公式 | 用 holdout 靶？ | 用 R？ | π 的路径 |
|---|---|---|---|---|
| 线性 gain | `max πᵀs` | 是（s） | 否 | 只有顶点 |
| Frank–Wolfe | 线性目标的 FW 步 | 是 | 否 | 每步仍是顶点 |
| **有限步长二次（顺序）** | `−πᵀs+(η/2)πᵀRπ`，η 给定 | 是 | 是 | 顶点 → min-var |
| **联合 (η,π)** | `max (πᵀs)²/(πᵀRπ)` | 是 | 是 | = GLS / 方向匹配 |
| 方向匹配 | `‖Gπ−ĝ_hold‖²` | 是 | 是 | 同上，η 无关的内点 |
| **熵正则** | `softmax(s/τ)` | 是 | **否** | 顶点 → **equal**（不拆 clone） |
| **MGDA** | `min ‖Gπ‖²` | **否** | 是 | 冲突几何的 min-norm，不是 stacking |
| mixloss / 聚合 | `−log Σ π_m e^{-L_m}` | 损失 | 否 | Vovk aggregating；proper scoring of the *mixture* |
| pseudo-BMA | `π ∝ e^{-n L}` | 边缘似然/损失 | 否 | `n=1` 温和；`n=|window|` 在 M-open 塌缩 |
| Hedge+share | 乘性 × 均匀 | 是 | 否 | tracking，不是门 |
| GLS-EWMA | `Σ̂^{-1}1` | 误差 | 是 | 在线 Bates–Granger |

同一组 (s, R) 上的快照：

| logic | π |
|---|---|
| vertex η=0 | a=1.00, b=0.00, c=0.00 |
| train η=3e-3 | a=1.00, b=0.00, c=0.00 |
| finite-step η=1 | a=0.92, b=0.00, c=0.08 |
| finite-step η=8 | a=0.41, b=0.16, c=0.43 |
| joint GLS (η*=0.87) | a=0.98, b=0.00, c=0.02 |
| entropy τ=1 | a=0.47, b=0.31, c=0.22 |
| MGDA（无靶） | a=0.26, b=0.26, c=0.47 |
| pseudo-BMA n=1 | a=0.47, b=0.31, c=0.22 |
| pseudo-BMA n=50 | a=1.00, b=0.00, c=0.00 |

读图：

- **训练 LR ≈ 顶点**：`train` 和 `vertex` 几乎一样。真实一步不会组合。
- **η=8** 把质量从共线的 a/b 分给独立的 c（二次项惩罚相关）。
- **联合 GLS** 在虚步长 η* 处取 Rayleigh 最优，不是「再把训练 LR 调大一点」。
- **熵正则不管 R**：τ 变大只是摊成 1/3，clone trap 解不了。
- **MGDA 没有 s**：更像「别打架」而不是「跟 holdout 对齐」。Sener & Koltun NeurIPS 2018 是多任务 Pareto，不是 stacking。
- **pseudo-BMA**：n=1 时损失差 0.4 nats 几乎还在混；n=50（ELPD 尺度）塌向 a。Yao, Vehtari, Simpson, Gelman (*Bayesian Analysis* 2018) 才主张 stacking 而不是 BMA——因为 stacking 优化的是预测风险，softmax **不乘 n**。

Armijo / 线搜索是有限步长的自适应 η：在 holdout 上沿 `d=Gπ` 找使 `L(θ−η d)` 下降的 η。训练 LR 太小，线搜索给出的虚步长才会进入组合区；把这条虚步长代回 π 的二次型，就是联合 GLS。

### 还相邻、但不要混进来的

| 名字 | 为什么不是这套 stacking |
|---|---|
| MoE / Switch Transformer | `π(x)=softmax(Wx)`，门看样本 |
| Gradient Blending (Wang et al. CVPR 2020) | 用过拟合间隙调模态权重，不是 OOF 组合 |
| PCGrad / OGM-GE | 梯度投影 / 冲突手术，不在单纯形上估 π |
| NCL (Liu & Yao) | 训练专家时罚相关，改的是专家本身不是 meta π |
| Mixup / 输入混合 | 混合的是 x，不是专家票 |

---

## 和 MoE 的对照（避免混）

| | MoE | 这里的 online stacking |
|---|---|---|
| π | `π_m(x)`，看样本 | `π_m(t)`，看 expert 票 / 损失 |
| 训练 | 门和 expert 联合（EM 或联合梯度） | expert 当黑盒，meta 只用 OOF |
| 目标 | 混合似然 `log Σ π_m(x) p_m(y\|x)` | 预测 / 方向风险 `L(π·v, y)` 或 `‖Gπ−g_hold‖²` |
| 多样性 | specialization by region of x | error–ambiguity / `N_eff` of votes |
| 文献 | Jacobs 1991；Switch Transformer | Wolpert；Super Learner；OSL；Bates–Granger |

熵正则的 softmax **不是** MoE 门，只是把分数温度化。若把 s 换成 `Wx`，那才滑向 MoE——不要滑。

```bash
PYTHONPATH=. python3 scripts/run_agod_stack_logics.py
PYTHONPATH=. python3 -m tests.test_agod_stack_logics
```
