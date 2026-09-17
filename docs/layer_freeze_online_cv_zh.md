# PO × MSE 对照：评估逻辑

看板不是假设检验，是两套**不同问题**的对照。崩 = 因果 MA ≥ 2× 各自的 baseline。不做 online-bootstrap。

## 两条序列在量什么

**PO-risk** 问的是：新 batch（\(T=1\)）相对 \(D_{\mathrm{ref}}\)（\(T=0\)），\(P(Y\mid X)\) 有没有 hop。

单独拟合 \(\mu(Y\mid X)\) 和 \(e(T\mid X)\)，不是 serving MLP：

\[
\varphi=(Y-\mu)(T-e),\quad \widehat\tau(X)\approx\varphi,\quad \mathrm{risk}=\mathrm{mean}(\widehat\tau^2)
\]

baseline 是把 \(D_{\mathrm{ref}}\) 对半切、假装一段 \(T=1\)。那是「没有 hop」时这个数该有的量级。stream MA ≥ 2× 这条线，只说明残差里 \(T\) 和 \(Y\) 的关联起来了，**不说明当前 MLP 已经不能用**。

**serving MSE** 问的是：当前 MLP 在新 batch 上 \(\mathrm{mean}((Y-\hat y)^2)\) 是否还停在预训练时 \(D_{\mathrm{ref}}\) 的误差量级。这是「全量每一层 backprop」这套策略还付不付得起房租。baseline 是预训练模型在 \(D_{\mathrm{ref}}\) 上的 MSE。同样用因果 MA、同样 2×。评估发生在 **update 之前**：先看当前策略在新 batch 上的数，再决定训还是冻。

两条可以分开走，这是对照存在的理由：

| PO-risk | serving MSE | 读法 | 动作 |
|---|---|---|---|
| 安静 | 安静 | 没有 hop，模型还在拟合 | **接着 train** |
| 崩了 | 没崩 | 机制可能动了，但现模型误差还在线内 | **再观察，不冻** |
| 崩了 | 崩了 | hop 在残差里可见，而且现模型也崩了 | **冻住**（策略不行） |
| 没崩 | 崩了 | 没有 hop 证据，只是这一段更难或更噪 | **接着 train**（不当成冻） |

## 为什么这样 justify

1. **冻要两个都崩。** PO 单独过线只是 hop 探测器。electricity 上 PO MA 过线而 MSE MA 不过（0.21 vs 2×0.121），全量 backprop 的误差并没有跟着翻倍，冻下层会把还在工作的更新掐掉。必须 MSE 也过线，才说「这套全量更新策略不太行」。

2. **只崩 PO 就观察。** hop 可以真实存在，同时 \(\hat y\) 仍够用（或正在跟上）。看板的职责是标黄，不是立刻改 freeze-depth。业务要不要切域、要不要停更，不是这个数能 justify 的。

3. **两个都正常就接着训。** covertype 的 raw PO 有过单个点偏高，但 PO MA 和 MSE MA 都在 2× 下面。读 MA 不是读单个点。安静段全开每一层，符合「没有大 hop 就不要冻」。

4. **MSE 单独偏高不冻。** 没有 PO hop 就没有「\(P(Y\mid X)\) 换了」的证据。单段 MSE 升高可以是同一机制下更难的 batch。对照表里这一格归 keep training，避免把噪声当策略失败。

5. **2× 和 MA 不是 p-value。** raw 点估计在小 `n_new` 上会抖；因果 MA（窗口约 1000 条流）是稳定性读数。2× 相对 baseline 只是看板上「明显偏离」。重复推 MLP 做 online-bootstrap 顶不住，也不该用来 justify 何时 update。

6. **冻住以后才归因层。** `PO_Dict` / `MSE_Dict` 只在 freeze hop 上填。那是「从哪一层开始动」的第二步，不是第一步的门。

何时 **update** 仍是业务逻辑。看板只给这三档对照。

```bash
PYTHONPATH=Python/src:. python3 scripts/run_layer_freeze_online_cv.py --dataset all --replay-json
```
