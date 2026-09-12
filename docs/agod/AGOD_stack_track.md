# 半程换人、OOF 特征、π 的变动

上一张讲的是 **选哪一种 meta-loss**（线性顶点 / 有限步长 / GLS）。这一张是 **协议里还没讲透的三件事**：最好的 expert 会换人、meta 只能看 OOF 特征、π 在单纯形上怎么走。还是不是 MoE：门不看 x。

---

## 半程换人（tracking，不是静态遗憾）

最好 expert 在 `t=T/2` 从 a 换成 b。对 **全程最好的固定 expert** 的遗憾是错的 oracle——那个固定赢家往往是一直中等的 c。对的 oracle 是 **最多换 k 次的比较器**（Herbster & Warmuth 1998 switching regret）。

三种锁死 / 三种花钱方式：

| 机制 | 会怎样 |
|---|---|
| 离散 OSL（累加损失） | 选出全程最好的**固定**专家（本例是中等的 c），看不见 a→b |
| Hedge `share=0` | 乘性更新从 ~0 质量爬升是对数慢的，半程经常来不及 |
| Hedge **Fixed-Share** `π ← (1-α)π⊙e^{-ηL} + α/\|E\|` | 每个死专家每步至少 `α/\|E\|`，能复活 |
| 滑窗 OSL（只看最近 W 步） | 忘掉前半程，W 步后跳到 b |
| 睡眠专家 / specialists | 前半程没投票的人不当成「很差」，只当没出场 |

share 不是学习率。它是 **允许 π 花掉的复活预算**：α=0 跟踪不动；α 太大则全程抖动，定常流上也在烧 TV。

本例 `T=80`、a→b 在 t=40：

| method | preq | path TV | TV after switch | delay | leader | final π |
|---|---:|---:|---:|---:|---|---|
| `hedge_share0` | 0.399 | 1.75 | 0.55 | 35 | b | a=0.04, b=0.96, c=0.00 |
| `hedge_share` | 0.062 | 3.04 | 0.55 | 5 | b | a=0.05, b=0.84, c=0.11 |
| `osl_disc` | 0.546 | 1.67 | 0.60 | ∞ | c | a=0.00, b=0.00, c=1.00 |
| `osl_window` | 0.157 | 2.67 | 0.75 | 10 | b | a=0.00, b=1.00, c=0.00 |
| `osl_sgd` | 0.060 | 5.61 | 0.58 | 2 | b | a=0.12, b=0.88, c=0.00 |

读表：离散 OSL 的累加损失选出的是全程最好的**固定**专家 c（一直中等），π 停在 c 上，a→b 的切换它看不见。Fixed-Share delay=5 把质量搬到 b，preq 从 0.40 掉到 0.06。`share=0` 最后也到了 b，但晚了 35 步，整段风险已经付过了。滑窗 OSL 是一次迟到的跳跃（delay=10）。Hedge `share=0` 从 ~0 质量爬升是对数慢的。

**负 regret vs best fixed expert 是预期，不是 bug。** 跟踪赢的是换人的比较器；离散 OSL 对齐的是那个「一直中等」的固定赢家。

---

## OOF 特征（Wolpert 的 level-1 矩阵）

Super Learner 的 meta 特征是 `Z_{i,m} = f_m^{(-i)}(x_i)`：第 m 个专家在 **没见过 i** 时对 i 的预测。用 `f_m(x_i)`（in-sample）当 stacking 特征，就是 Wolpert 警告过的 leaky stacking。流上的翻译：

```
P_t probe  → v_{m,t}     # level-1 票（OOF）
H_t holdout → y_t / g_hold # 靶，和 P_t 不相交
score with frozen π_t
train bases on A_t
update π_{t+1}
```

三种泄漏，严重程度不一样：

| 泄漏 | 发生了什么 | 看起来 |
|---|---|---|
| in-sample 票 | 专家在同一点上训过再投票 | meta 以为噪声专家也好 |
| 同一窗既当票又当靶 | `g` 同时是 vote 和 `g_hold` | 方向余弦被抬成 1 |
| 先 update 再 score | `leaky_step` | 低估 prequential risk |

本例：信号专家看 1 维 x，噪声专家看 60 维纯噪声（in-sample 能背下来）。

| | π_good | π_noise | corr(Z_noise, y) | in-sample MSE | holdout MSE | optimism |
|---|---:|---:|---:|---:|---:|---:|
| OOF | 1.00 | 0.00 | -0.09 | 0.035 | 0.030 | -0.005 |
| leaky | 0.86 | 0.14 | 0.87 | 0.030 | 0.107 | 0.077 |

诊断是 **corr(Z_noise, y)**：leaky 票是噪声专家在同一点上背下来的，会跟 y 相关；OOF 票几乎不相关。π 跟着走——OOF 把质量放在 good 上，leaky 给噪声专家更多票。holdout MSE / optimism 是后果，小样本会抖；相关才是特征有没有漏。

梯度协议里同一 batch 的对齐余弦 = 1.00，probe/holdout 不相交才是 0.77。票和靶必须拆开。

能进 level-1 的列：OOF 预测、OOF 残差、OOF 方向分数 `⟨ĝ_m, ĝ_hold⟩`。不能塞：本窗 in-sample 损失、x 本身、时间下标 t——把 x 或 t 塞进门，就滑向 MoE / 时变门，不是 stacking。

---

## π 的变动（要花在刀刃上）

`TV(π_t, π_{t+1}) = ½‖π_{t+1}−π_t‖_1`。`path TV = Σ_t TV_t` 是这条轨迹的复杂度（和 switching 次数是一家）。

| 流 | 想要的 π |
|---|---|
| 半程换人 | 切换后 **局部** TV 尖峰，把质量从 a 搬到 b |
| 定常 | path TV 小；大 TV = 在拟合上一窗噪声 |
| clone | TV 可以有，但应拆冗余票，不是来回抖 |

本例切换流 Hedge+share path TV = 3.04。同样方法在**定常**流上是 8.28——三个差不多好的专家时，Fixed-Share 每步掺 `α/K`，会在没事的时候烧 TV（复活税）。sticky λ=0.12 把执行器 TV 压下来。跟踪算法的 TV 应该花在切换点，而不是当成默认抖动。

离散 OSL：TV 是 **一次跳到顶点**（不平滑）。Hedge：连续、乘性。滑窗：窗口一满就跳。Fixed-Share：每步至少掺 `α/K`，底噪 TV 换复活。

**Two-timescale（给执行器）：** meta 可以用快的 π（跟踪），LR 用慢的 `π_slow ← (1-λ)π_slow + λ π_fast`。FWD 一直开着的时候，π 一抖 LR 就抖。λ=0.12 的 sticky 路径 path TV 更小、延迟稍长——这是执行器和跟踪之间的预算。

share α、窗宽 W、sticky λ 是 **三个花 TV 的旋钮**，不是三个新模型。α/W 管「半程换人复不复活」，λ 管「下一步 LR 跟不跟得上」。

---

## 和前面几张的关系

| 已经说清 | 这一张补的 |
|---|---|
| 线性 gain → 顶点 | 半程换人时顶点会 **锁死在旧人** |
| 有限步长 / 联合 GLS | 那是一张窗里怎么组合；这里是 π **跨窗怎么走** |
| honesty gap | 把 gap 拆成 OOF 特征泄漏 vs 先 update 再 score |
| Hedge+share 能跟踪 | share、窗、sticky 各自花掉多少 TV、换来多少 delay |

```bash
PYTHONPATH=. python3 scripts/run_agod_stack_track.py
PYTHONPATH=. python3 -m tests.test_agod_stack_track
```
