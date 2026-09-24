# 服务误差 vs 份额误差（要跟 latency 一起看）

Hop 和连续时间是 **同一条** \(S(\tau)\)：hop 读 \(S(\tau)/S(\tau-W)\)，水平读 \(S(\tau)\) vs \(S_{\mathrm{ref}}\)。其余都好说。

真正 tricky 的是两条 **不同** 的误差过程，共用一个钟，过线时刻不一样：

| 过程 | 是什么 | 错了意味着 |
|---|---|---|
| 服务误差 | 现模型在窗口上的 Brier/MSE 相对 pre-onset | 房租在走；现在的策略付得出去付不出去 |
| 份额误差 | 定位指针有多歪：onset 后 \(1-\pi_{\mathrm{true}}\)；onset 前 \(\pi_{\mathrm{true}}\) 本身就是误指 | 塔指错了；还没 drift 就指 |

份额不是真分解。份额误差是指针噪声，不是 Shapley 残差。

---

## 表 1 · 四格（不是 WHEN 的两种理论）

同一窗、两个 loud：

- serve loud：Brier MA ≥ 1.1 × pre-onset
- share loud：\(\pi^{\mathrm{PO}}(\mathrm{video})\) **连着两窗** ≥ 0.5（单窗 π 会在 t=0 就到 0.55，不能当钟）

| 格 | 服务 | 份额 | 跟 latency 绑在一起时该干什么 |
|---|---|---|---|
| quiet | 安静 | 安静 | 全量或接着看。单窗 π 乱跳也算安静 |
| share_only | 安静 | 响 | **watch**。指针先动、房租还在。这时更新/冻层，付的是服务误差（过早适应） |
| serve_only | 响 | 安静 | 房租在走，塔还没锁。可以准备更新，**不要按份额去冻某一塔** |
| both | 响 | 响 | 才把 `train_top` 落到那座塔。确认是服务误差，指向是份额 |

Latency 不是 hat−τ0。是：你待在 share_only / serve_only 里 **多久**，以及这段时间的 online Brier。

---

## 表 2 · 这次慢走（同一条渐进 concept）

hat_serve=6，hat_share=7，hat_both=7。份额钟比服务钟 **晚 1 窗**。t=0 的 π=0.55 被连着两窗挡掉了，否则会在 onset 前误指。

| t | α | 服务超额 | 份额误差 | 格 |
|---|---|---|---|---|
| 0–2 | 0 | ~0 | 最高 0.55 | quiet（π 已经在晃） |
| 3–5 | 0.04–0.19 | ~0 | **0.58–1.0** | quiet：机制开始转，指针更歪，房租还没走 |
| 6 | 0.27 | 0.059 | 0.095 | **serve_only** 一窗 |
| 7–15 | 0.35–0.96 | 0.09→0.25 | 0–0.41 | **both** |

onset 后格计数：quiet 3 · share_only **0** · serve_only 1 · both 9。

读法：渐进里份额误差在服务误差前面先变脏，再变干净。如果你用单窗份额当 WHEN，你会在房租还没动、指针最脏的时候动手。连着两窗之后，份额钟几乎贴着服务钟（晚 1 窗），share_only 可以为空——那是过滤器的代价，也是 latency。

---

## 表 3 · 用哪口钟去更新，评估必须是服务误差

各策略：到 hat 之前冻住，hat 起才 refit。面积 = onset 之后 online Brier 相对「每窗都更新」。

| 策略 | hat | area vs always |
|---|---|---|
| never | ∞ | **+0.207**（一直不更新，房租全付） |
| from_serve | 6 | −0.008 |
| from_share | 7 | −0.014 |
| from_both | 7 | −0.014 |
| always | 0 | 0 |

负数：等服务/份额都响了再更新，比从 t=0 就全量 **更省服务误差**。Always 不是 oracle。Latency 评估是这条面积，不是检测器谁先响。

所以综合 latency 的规则就一句：

**动塔看份额，动手看服务；两口钟的错开用面积计，不用 hop。**

- 只有份额响 → watch，把 latency 花在等服务误差
- 只有服务响 → 可以开始更新，但还不能按塔冻
- 两口都响 → `train_top` 那座塔
- 用连着两窗降份额误差，付的是 1 窗定位延迟；这次换来 onset 前的误指被挡住

图：`results/serving_vs_share/serving_vs_share.png`。上：两条误差过程不同步。下：never 后期抬起来，from_both 和 from_share 贴在一起。
