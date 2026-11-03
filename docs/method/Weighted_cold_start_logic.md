# Cold-start logic (v4): opportunity-gated attribution + separation

## One sentence
**每批重估适应新类的机会 \(O_t\)：机会够 → 归因聚焦 \(S\) → 分离 →（量够再）扩头加权；机会不够或足迹不稳 → hold / 作废几何。下一批未知，故 \(O^{\mathrm{persist}}\) 默认打折。**

Feature-level（RFPerm/FSDS）和 sample-level 权重都能赋能新类冷启动：一层管 where，一层管 whom；**\(O_t\) 决定这一批动不动。**

## Adaptation opportunity（核心）

\[
O^{\mathrm{local}}_t=(D_t L_t R_t M_t)^{1/4},\quad
O^{\mathrm{persist}}_t=O^{\mathrm{local}}_t\cdot\mathrm{Stab}_t.
\]

| 因子 | 角色 |
| --- | --- |
| \(D\) | detect：error rise |
| \(L\) | localize：VIMP mass on \(S\) |
| \(R\) | separate：还能推离旧原型 |
| \(M\) | mass：有行可花 |
| \(\mathrm{Stab}\) | 下一批未知下的复用折扣 |

政策：`O_local≥τ_sep` → sep；再够 → expand（sample weights）；`O_persist<τ_keep` → 不盲留 offset。

详见 `docs/method/Adaptation_opportunity_novel_cold_start.md`。

## Track benchmarks

- **C-GCD**：CIFAR-100 / ImageNet-100 / Tiny-ImageNet / CUB — Old/New/All per session（MetaGCD, Happy）
- **PDA/OSDA/UniDA**：Office-Home / VisDA / DomainNet — H-score
- **我们**：synth C-GCD + PDA proxy + tabular；额外审计 **\(O_t\)**, novel_acc, forgetting, gap_lift, \(s(w)\)

## Four steps（仍成立，但由 \(O_t\) 门控）

1. **Detect + attribute** — 每批重算；role=focus  
2. **Representation separation** — 仅当 \(O^{\mathrm{local}}\) 过线；只推 private-live/high-mem  
3. **Dynamic head expand** — buffer / new-private / \(O\) 过 expand 线；sample_w×subset 在此点火  
4. **Optional TTA on \(S\)** — 谨慎；offset 仅当 \(O^{\mathrm{persist}}\) 高才 EMA

## 两层权重如何赋能

| 层 | 赋能什么 | 何时 |
| --- | --- | --- |
| Feature-level VIMP | 找到 \(S\)，抬高 \(L\)（从而抬 \(O\)） | 每批归因 |
| Sample-level \(w_i\) | 抬高 novel SGD 份额 | expand 时 |

连续时间：归因链路连续重跑；不是一次 VIMP 吃到永远。

## Code

- `adaptation_opportunity.py` — \(O_t\)
- `repr_sep_cold_start.py` — opportunity-gated pipeline
- `scripts/run_adaptation_opportunity_probe.py`

```bash
python3 scripts/run_adaptation_opportunity_probe.py --quick
python3 scripts/run_repr_sep_cold_start.py --quick
```
