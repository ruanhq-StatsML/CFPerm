# Use-case 4 — LLM moderation (车道 A)

Drop-in manuscript section. One prediction table: **one Y, many X**. Not attribution.

## English

**Business.** A reviewer (human, policy model, or external LLM judge) writes a decision as traffic arrives. The object to monitor is whether \(P(Y\mid X)\) is still the same auditor, not the refusal rate and not HH preference.

**Table.** \(Y\in\{0,1\}\) is pass/fail. \(X=(x_{\mathrm{n\_toks}},\ldots,x_{\mathrm{refuse}},\ldots,x_{\mathrm{thank}})\) are reply features. Arrival windows are `batch`. For CFPerm, two reviewer queues are encoded as \(T\in\{0,1\}\).

Anthropic HH-RLHF `chosen`/`rejected` is **not** \(Y\). Those fields only label which queue the row came from (helpfulness vs harmlessness). The deployed auditor’s pass/fail is \(Y\).

**Already on disk (prototype A).**

- Helpful / harmless streams: `results/manuscript/llm_audit/xy_hh_*.csv` (\(n=1200\) each).
- Two-stream CFPerm table: `xy_hh_two_stream.csv` (\(T=0\) helpful, \(T=1\) harmless, \(n=2400\)).

**More realistic labels (still the same table shape).**

| Source | \(Y\) | File | \(n\) | pass rate |
|---|---|---|---:|---:|
| BeaverTails | `is_safe` | `xy_beavertails.csv` | 1200 | 0.43 |
| WildGuard | `response_harm_label = unharmful` | `xy_wildguard.csv` | 1200 | 0.92 |
| ToxicChat | human `toxicity = 0` | `xy_toxicchat.csv` | 1200 | 0.78 |
| BeaverTails vs ToxicChat | real labels, two queues | `xy_real_two_stream.csv` | 2400 | — |

A hop of \(P(Y\mid X)\) is a policy-pack / judge swap, not “this feature caused the fail.” After a fire: Top-\(k\) re-review, do not treat the current judge as gold.

**Two streaming gates (same table, not interchangeable).**

1. **Frozen-ref online AR-bootstrap** (Palm & Nagler; same protocol as the concept-drift board). Fit lstsq once on \(D_{\mathrm{ref}}\) (first `n_ref` batches). \(s_t\) is that probe’s MSE on trail batch \(t\), \(\mu_{\mathrm{ref}}\) is the mean MSE on the \(D_{\mathrm{ref}}\) mini-batches, \(\Delta_t=s_t-\mu_{\mathrm{ref}}\). Fire when the online AR-bootstrap \(\mathrm{CI}_{\mathrm{lo}}(\Delta)>0\). Needs two trail updates before a CI exists. The last-two gate still uses the shallow RF.
2. **Last-two `hop_fires`** (OnlineRFPerm). Refit on \(B_{t-1}\), score \(B_t\). Fire iff \(e_{\mathrm{prev}}\ge e_{\mathrm{floor}}\) and \(e_{\mathrm{now}}/e_{\mathrm{prev}}\ge\gamma\). First hop is always quiet.

Prototype (HH consistent/hop + BeaverTails / WildGuard / ToxicChat native and Y-flip overlay): `PYTHONPATH=. python3 scripts/llm_audit_online_bootstrap_prototype.py`. Numbers and CI plots live in `results/manuscript/llm_audit_online_bootstrap/`.

```r
tab <- read.csv("results/manuscript/llm_audit/xy_hh_two_stream.csv")
X <- as.matrix(tab[, grep("^x_", names(tab))])
Y <- tab$Y
T <- tab$T
set.seed(2026)
cfperm(X = X, Y = Y, T = T, n_perm = 30, vimp = "permuCATE",
       seed = 2026, level_feature = 0.05, level_across_feature = 0.05, top_k = 1)
```

Rebuild: `python3 scripts/build_llm_audit_xy.py`. Raw prompts are not stored; only \((Y,X,\mathrm{batch})\).

## 中文

**业务。** 审核员（人 / 政策模型 / 外挂 LLM）按到达时间写决定。要盯的是 \(P(Y\mid X)\) 还是不是同一套审核逻辑，不是拒绝率，也不是 HH 偏好对。

**一张表。** \(Y\in\{0,1\}\) 过/不过。\(X=x_{\mathrm{n\_toks}}\ldots x_{\mathrm{refuse}}\) 是待审回复特征。`batch` 是到达窗。CFPerm 里两路审核队列写成 \(T\in\{0,1\}\)。

HH 的 chosen/rejected **不当 \(Y\)**。那只说明流量从 helpful 还是 harmless 队列来。\(Y\) 是部署审核器当场写的过/不过。

**现成 prototype A** 在 `xy_hh_*.csv`。**更像真审核的标签** 用 BeaverTails / WildGuard / ToxicChat，仍然是同一张预测表，见 `xy_beavertails.csv`、`xy_wildguard.csv`、`xy_toxicchat.csv`，以及两路真标签 `xy_real_two_stream.csv`。

补真 \(Y\) 之后，probe 仍然只学 \(P(Y\mid X)\)。hop = 政策包 / judge 换代。不要把 13 个 \(x\) 当成审核逻辑的归因维。

**两个闸，同一张表，不要混。** 冻参考窗的 Palm–Nagler online AR-bootstrap：\(\Delta_t=s_t-\mu_{\mathrm{ref}}\)，\(\mathrm{CI}_{\mathrm{lo}}>0\) 才叫显著变差。Last-two `hop_fires`：相邻窗 OOS 比 \(\ge\gamma\) 才叫 map hop。CI 至少要两个 trail 点；last-two 的第一跳永远 quiet。原型脚本 `scripts/llm_audit_online_bootstrap_prototype.py`。

连续时间 serving 的同一套表（Graph-RAG / 混合检索 / Agent 下一步）见 `docs/manuscript/use_case_05_rag_agent_stream.md`。
