# Online stacking of modality experts（怎么做、怎么评、文献在哪）

静态 `π ∝ R̃^{-1} α` 只是 stacking 的**一张快照**。代码好写；难的是 **online-stacking 的协议**：权重在单纯形上、只用 one-step-ahead 损失更新，而且评估必须对 oracle。这件事是 **modality-agnostic** 的——expert 只是带名字的投票器。

## 文献地图（这不是从零发明）

| 传统 | 文献 | 我们把它用在哪 |
|---|---|---|
| Stacked generalization | Wolpert 1992; Breiman 1996 stacked regressions | 模态 = base expert，π = meta-learner |
| Super Learner | van der Laan, Polley, Hubbard 2007; oracle inequality van der Laan & Dudoit 2003 | CV 组合渐近等价于 library 里最好的凸组合 |
| **Online Super Learner** | Benkeser, Ju, Lendle, van der Laan, *Stat Med* 2018 ([PMC5671383](https://pmc.ncbi.nlm.nih.gov/articles/PMC5671383/)) | 流上的 one-step-ahead CV；离散 OSL = 累加损失最小的 expert；凸 OSL = 单纯形上 SGD |
| Prequential | Dawid 1984 | **先用冻住的 π_t 计分，再 update** — 流式 OOF |
| Hedge / weighted majority + Fixed-Share | Littlestone & Warmuth 1994; Freund & Schapire 1997; **Herbster & Warmuth 1998 tracking** | `π ← (1-α) π⊙exp(−ηL) + α/|E|`，否则半程切换会被前半程锁死 |
| Tracking the best expert | Herbster & Warmuth 1998 | 非平稳窗（概念漂移）下的动态遗憾 |
| Forecast combination | Bates & Granger 1969; Newbold & Granger 1974; Timmermann 2006 | 误差协方差的最小方差权；**实践中常丢掉 off-diagonal**（Σ 估不稳） |
| Error–ambiguity / BVC | Krogh & Vedelsby 1995; Ueda & Nakano 1996 | `N_eff`、clone trap：相关 expert 不是新票 |
| Negative correlation learning | Liu & Yao 1999 | 负相关才是多样性；与 shared head 冲突不是 |
| Mixture of experts | Jacobs et al. 1991 | gating 网是另一类 meta；我们用损失驱动的单纯形权 |
| Gradient Blending | Wang, Tran, Feiszli CVPR 2020, arXiv:1905.12681 | **不是 stacking**：按 OGR 混 loss。我们混的是 *expert 的票* |
| PCGrad / OGM-GE | Yu et al. 2020; Peng et al. 2022 | 冲突投影 / 梯度调制，没有 online CV oracle inequality |

缺口很干净：Gradient Blending / PCGrad 不懂 prequential honesty，OSL / Hedge / Bates–Granger 不懂「票是梯度签名」。把模态塔当成 **online expert**，vote 的构造见下一节，π 走 OSL/Hedge/BG，再 `π → LR`。video/text/audio/tabular **同一段代码**。

## 梯度 expert 的票怎么造（这才是「怎么做」）

预测 stacking 里 `v_m = f_m(x)`、`y` 可观测。梯度 stacking 里 `g_m ∈ R^d`，真更新 `μ = −∇L_pop` **看不见**。holdout 梯度 `g_hold` 是 μ 的 noisy proxy。三种 meta-loss，几何完全不同：

| meta-loss | 公式 | π 的形状 | 文献对应 |
|---|---|---|---|
| 一阶 holdout gain | `s_m = ⟨ĝ_m, ĝ_hold⟩`，`max πᵀs` | **顶点**（离散 OSL / winner-take-all） | 线性 Super Learner 退化 |
| **方向匹配** | `min ||Gπ − ĝ_hold||²` | 内点：`π ∝ R^{−1}s`，`R=GᵀG` | Breiman stacked regressions in `R^d`；就是上一张 GLS 快照 |
| 均值–方差 | `max πᵀs − (λ/2) πᵀRπ` | λ=0 顶点，λ↑ 缩向 min-var | Markowitz / Bates–Granger |
| 有限步长 | `L_hold(θ − η Gπ)` | 非线性，η 大才离开线性 | 真 virtual step，贵 |

所以：**不能**把 directional cosine 直接塞进凸 stacking 还指望得到组合权——线性目标在单纯形上必崩到一个 expert。要组合，必须用方向匹配（平方损失）或方差惩罚。`direction_match_weights` 就是这件事：互补的两路单位梯度（轴 0 和轴 1）去拟合 `ĝ_hold = (e0+e1)/√2` 时，π 会拆成约 (0.5, 0.5, 0)，而 `linear_gain_is_vertex` 只会点其中一个。

诚实协议在梯度上还多一条：`g_m` 和 `g_hold` 必须来自 **不同 batch**（probe vs holdout）。同一窗的 holdout 既当票又当靶，就是 leaky stacking。

```
P_t probe → g_m          # 票（OOF 特征）
H_t holdout → g_hold     # 靶（one-step-ahead y）
score ||G π_t − ĝ_hold||²
train on A_t
update π_{t+1}
```

## 协议（这才是 online-stacking，不是 `R^{-1}α` 一行）

```
for window t:
  1. freeze π_t
  2. each expert emits vote v_{m,t}  *before* being trained on window t
  3. score ŷ_t = π_t · v_t on the holdout of window t     # prequential loss
  4. train bases on the adapt split                       # FWD still on
  5. update π_{t+1} from (v_t, y_t, L_{m,t})
  6. actuator: LR_{t+1,m} = β + (1-β)|E| π_{t+1,m}
```

漏写第 2–3 步就是 **leaky stacking**（用 in-sample 票训 meta），等价于 Wolpert 警告过的用训练集预测做 stacking 特征。流上的正确类比是 Dawid 的 prequential / OSL 的 one-step-ahead CV。

Meta-learners：

| method | 更新 | 何时该赢 |
|---|---|---|
| `equal` | 不动 | 相关高、Σ 估不稳（Bates–Granger 自己的告诫） |
| `hedge` | 乘性权重 | expert 会切换（漂移） |
| `osl_disc` | argmin 累加损失 | 有一个稳定的最好 expert |
| `osl_sgd` | 投影 SGD on `(π·v − y)²` | 要凸组合，K 小 |
| `bg` | `π ∝ 1/MSE_m`（忽略相关） | 异方差、Σ 噪 |
| `gls_ewma` | `π ∝ Σ̂^{-1} 1` + 特征值收缩 | clone / 共线票，且 T 够估 Σ |

## 怎么评估（对 oracle，不只看 Acc lift）

1. **Prequential risk** `R_n = n^{-1} Σ_t (π_t·v_t − y_t)²` — 唯一诚实的流风险  
2. **vs best expert**（离散 oracle）和 **vs best fixed convex combo**（网格上的 Super Learner oracle）  
3. **Tracking**：半程切换最好 expert，Hedge/OSL 必须把质量搬过去  
4. **Clone trap**：两个几乎相同的好 expert + 一个独立 expert；等权把 2/3 砸在同一个因子上，GLS 必须丢掉冗余票  
5. **Honesty gap**：leaky − honest；leaky 应系统性低估风险  
6. **N_eff(π)**：stacking 是否还剩多样性  
7. 模态无关：expert 名叫 `a,b,c` 或 `video,text,audio` 行为应一样  

### Switch stream（最好 expert 在 t=T/2 切换）

| method | preq MSE | regret vs best | regret vs combo | N_eff |
|---|---:|---:|---:|---:|
| equal | 0.226 | -0.302 | +0.000 | 3.00 |
| hedge | 0.075 | -0.453 | -0.151 | 1.76 |
| osl_disc | 0.546 | +0.017 | +0.320 | 1.00 |
| osl_sgd | 0.058 | -0.470 | -0.168 | 1.41 |
| bg | 0.076 | -0.452 | -0.150 | 1.15 |
| gls_ewma | 0.097 | -0.431 | -0.129 | 1.91 |

负 regret vs best expert 是正常的：全时段最好的**单个** expert 往往是那个「一直中等」的 c，组合/跟踪能赢它。`osl_disc` 锁在前半程赢家 a 上，后半程崩掉（preq 0.55）——离散 Super Learner 在非平稳流上必须加 Fixed-Share / 窗，不能只用累加损失。`hedge`+share 和 `osl_sgd` 把质量从 a 搬到 b。

### Clone trap（e0≈e1 共线，e2 独立）

| method | preq MSE | regret vs best | regret vs combo | N_eff |
|---|---:|---:|---:|---:|
| equal | 0.136 | +0.029 | +0.055 | 3.00 |
| hedge | 0.092 | -0.015 | +0.011 | 2.81 |
| osl_disc | 0.109 | +0.002 | +0.028 | 1.00 |
| osl_sgd | 0.092 | -0.015 | +0.010 | 1.31 |
| bg | 0.102 | -0.005 | +0.021 | 2.74 |
| gls_ewma | 0.092 | -0.015 | +0.011 | 2.26 |

Clone 终盘 π：equal e0+e1 = 0.67；GLS e0+e1 = 0.31（冗余票被拆掉）。

### Honesty gap（Hedge）

| stream | honest | leaky | gap |
|---|---:|---:|---:|
| switch | 0.087 | 0.070 | +0.017 |
| clones | 0.093 | 0.079 | +0.014 |

## Amazon / MSR-VTT replay

把 `align_gain = ½(1+cos(g_m,g_shared))` 当 vote、目标=1，跑诚实 Hedge。T=6 且 align_gain 常挤在 0.5 附近，**不够当 OSL 渐近，也不该指望看到 Amazon 共线 → N_eff 下降**；replay 只证明同一套 stacker 吃 `video/text/audio` 或 `text/image` 这些名字。真正的评估是上面的 switch / clone / honesty 三条。

| traj | top expert | π_top | N_eff | preq |
|---|---:|---:|---:|---:|
| amazon:equal | image | 0.54 | 2.00 | 0.243 |
| amazon:soft | image | 0.53 | 2.00 | 0.244 |
| amazon:soft_gradcos | image | 0.53 | 2.00 | 0.242 |
| msrvtt:equal | audio | 0.34 | 3.00 | 0.254 |
| msrvtt:soft | audio | 0.34 | 3.00 | 0.253 |
| msrvtt:soft_gradcos | audio | 0.34 | 3.00 | 0.253 |

下一步 LR 仍是 `pi_to_lr(π)`，FWD 不关。T 短时不要用 replay 的 N_eff 下结论。

## 和上一张 GLS 快照的关系

`soft_stat` 的 `π* ∝ R^{-1}α` 是 **窗口内闭式 stacking**（Bates–Granger / 匹配滤波）。Online stacking 多了三件闭式没有的东西：

- **时间**：π 是状态，能 track 切换（Hedge 遗憾界）  
- **诚实**：one-step-ahead，不拿本窗训练票更新本窗权  
- **oracle inequality**：OSL 对 library 里最好凸组合渐近等价（Benkeser et al. 2018）  

相关矩阵仍在：`gls_ewma` 用误差协方差的 EWMA，clone 时它就是 `N_eff` 故事的在线版。

```bash
PYTHONPATH=. python3 scripts/run_agod_online_stacking.py
PYTHONPATH=. python3 -m tests.test_agod_online_stacking
```
