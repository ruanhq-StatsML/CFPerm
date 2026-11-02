# Adaptation opportunity for novel-category cold-start

## One sentence
**归因 + 分离可以做，但真正要刻画的是：这一批有没有「适应到新类」的机会 \(O_t\)；下一批长什么样未知，所以机会是 batch-local 的，持久复用要打折。**

## Track benchmarks（这个赛道对标什么）

| 赛道 | 标准板 | 主指标 | 和咱们的关系 |
| --- | --- | --- | --- |
| **C-GCD** | CIFAR-100, ImageNet-100, Tiny-ImageNet, CUB（MetaGCD / Happy 协议：Stage-0 已知类，之后 session 混 known+novel） | 每 session **Old / New / All** acc | 主赛道：连续发现新类 + 防遗忘 |
| **GCD / NCD** | CIFAR-10/100, ImageNet-100, CUB, Cars, Herbarium | All / Old / New（单次 unlabeled mix） | 非连续版；表示分离的邻居 |
| **PDA / OSDA / UniDA** | Office-Home, VisDA, DomainNet | shared acc / **H-score**（known↔unknown） | 标签支撑不一致；我们的 Jaccard 边界同精神 |
| **我们的代理** | synth C-GCD、Office-Home PDA proxy、tabular HAR/mice/… | novel_acc、forgetting、feat_hit、\(s(w)\)、**\(O_t\)** | 工程可跑的 cold-start + 机会审计 |

外部方法多在表示对比 / meta / 蒸馏上卷 New acc；我们多一个可运维量：**本批适应机会 \(O_t\)**，用来决定要不要 sep / expand / 是否保留 \(S\)-offset。

## 算法刻画：适应新类的机会

对时刻 \(t\) 的 LIVE 批，定义

\[
O^{\mathrm{local}}_t = (D_t\, L_t\, R_t\, M_t)^{1/4}, \qquad
O^{\mathrm{persist}}_t = O^{\mathrm{local}}_t \cdot \mathrm{Stab}_t.
\]

| 因子 | 含义 | 估计 |
| --- | --- | --- |
| \(D\) Detect | 有没有必要动 | LIVE/REF error ratio（RFPerm gate） |
| \(L\) Localize | 新类差在哪 | FSDS VIMP 在 \(\mathcal{S}\) 上的质量 \(v(\mathcal{S})\) |
| \(R\) Separate | 还能不能推离旧中心 | membership 浓度 / gap + 相对 known prototype 的 residual |
| \(M\) Mass | 有没有行可花 | high-mem / private-live 计数 vs \(N^\star\) |
| \(\mathrm{Stab}\) | 下一批未知下的复用折扣 | 与上批 VIMP-top / feature 支撑 overlap |

**读法**

- \(O^{\mathrm{local}}\) 高 → **本批**值得归因+分离（feature-level 聚焦 + 几何推开）。
- \(O^{\mathrm{local}}\) 高且 \(M\) 够 → 才 **expand**（sample-level 权重在这里赋能 SGD 份额）。
- \(O^{\mathrm{persist}}\) 低 → **不要**赌 \(\mathcal{S}_t\)/offset 对 \(t+1\) 仍好使；短记忆、下批重归因。

这就是「连续时间归因链路」的门控：两层权重都能赋能新类冷启动，但只有机会够时才动手；机会不够就 hold。

## 政策（锁定）

```
if incomparable:          skip + drop offset
elif O_local ≥ τ_sep:     representation separation on S
   if also expand-ready:  head expand with sample_w / fused
   if O_persist ≥ τ_keep: EMA-keep S-offset
   else:                  replace offset (batch-local only)
else:                     hold; drop stale offset
```

默认：每批重算 VIMP（feature-level）；sample-level 只在 expand 点火；几何状态可作废。

## 和「不知下一批」的关系

机会函数把「下一批未知」写进 \(\mathrm{Stab}\) / \(O^{\mathrm{persist}}\)，而不是假装 \(P_{t+1}=P_t\)。  
**类节点可以累积；\(\mathcal{S}\)-几何默认不累积。**

## Code

- `Python/src/adaptation_opportunity.py` — `adaptation_opportunity(...)`
- `repr_sep_cold_start.adapt_batch_repr_sep` — 用 \(O_t\) 门控 sep / expand / keep_offset
- `scripts/run_adaptation_opportunity_probe.py`

```bash
python3 scripts/run_adaptation_opportunity_probe.py --quick
python3 scripts/run_repr_sep_cold_start.py --quick
```
