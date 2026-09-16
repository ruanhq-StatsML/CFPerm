# HH-RLHF → tidy 表 → OnlineRFPerm 连续实时检测

## 数据是什么

| 项 | 值 |
|----|-----|
| 源 | **Anthropic/hh-rlhf** `helpful-base`（对齐 / 偏好） |
| 缓存 | `data/hf_cache/hh_rlhf_helpful_base_2500.jsonl` |
| 整理后 | `results/agod/hh_online_stream/stream_table.parquet` |
| 行数 | 2000（chosen+rejected 展开） |
| featurizer | `style10+hash128` |
| n_per / cut_batch | 80 / 4 |

### 你要的 form（样例）

| t_idx | batch | question | answer | y | feature |
|------:|------:|----------|--------|--:|---------|
| 0 | 0 | How can you learn to be polite… | I have an idea about a system I’d like t… | 1 | dim=138 |
| 1 | 0 | How are LEDs made?… | Maybe!  But it also depends on the type … | 1 | … |

真实线上：`t_idx` 换成 serving / 标注时间戳；`y` 换成人工偏好或裁判分；**不要**做 demo 里的标签翻转。

## 怎么连续实时检测（核心）

```text
新 batch t 到齐
    → μ0 = fit_online_probe(batch t-1)     # 上一窗当制度探针
    → e_now = err(μ0, batch t)
    → fire ⇔ e_now / e_prev ≥ γ 且 e_prev≥floor
    → po_i  = po_risk0_rows(μ0, batch t)
    → fire? w=√po , audit Top-k(po)
      quiet? w=1 退火
    → e_prev ← e_now                       # 连续：下一跳接着比
```

包调用（就这几行）::

```python
from agod.online_rfperm import fit_online_probe, probe_err, hop_fires, po_risk0_rows
from agod.po_iptw import po_iptw_weights

probe = fit_online_probe(X_prev, y_prev, task="acc")
e_now = probe_err(probe, X_cur, y_cur, task="acc")
fire  = hop_fires(e_now, e_prev, gate=1.25)
po    = po_risk0_rows(probe, X_cur, y_cur, task="acc")
w     = po_iptw_weights(po, mode="sqrt") if fire else np.ones_like(po)
```

## 本跑 hop 表（连续）

| batch | t_idx | fired | ratio | e_now | P@k | action |
|------:|-------|:-----:|------:|------:|----:|--------|
| 1 | 80–159 | False | 1.000 | 0.475 | 0.60 | `ACCEPT_merge + w=1…` |
| 2 | 160–239 | False | 1.079 | 0.512 | 1.00 | `ACCEPT_merge + w=1…` |
| 3 | 240–319 | False | 0.878 | 0.450 | 0.60 | `ACCEPT_merge + w=1…` |
| 4 | 320–399 | True | 1.222 | 0.550 | 0.10 | `REJECT_or_LIMIT_merge + audi…` |
| 5 | 400–479 | False | 1.000 | 0.550 | 0.90 | `ACCEPT_merge + w=1…` |
| 6 | 480–559 | False | 0.773 | 0.425 | 0.60 | `ACCEPT_merge + w=1…` |
| 7 | 560–639 | True | 1.265 | 0.537 | 0.40 | `REJECT_or_LIMIT_merge + audi…` |
| 8 | 640–719 | False | 0.907 | 0.488 | 0.50 | `ACCEPT_merge + w=1…` |
| 9 | 720–799 | False | 1.077 | 0.525 | 0.40 | `ACCEPT_merge + w=1…` |
| 10 | 800–879 | False | 0.952 | 0.500 | 0.70 | `ACCEPT_merge + w=1…` |
| 11 | 880–959 | False | 1.000 | 0.500 | 0.10 | `ACCEPT_merge + w=1…` |
| 12 | 960–1039 | False | 0.950 | 0.475 | 0.70 | `ACCEPT_merge + w=1…` |
| 13 | 1040–1119 | False | 1.105 | 0.525 | 0.20 | `ACCEPT_merge + w=1…` |
| 14 | 1120–1199 | False | 1.071 | 0.562 | 0.80 | `ACCEPT_merge + w=1…` |
| 15 | 1200–1279 | False | 0.822 | 0.463 | 0.00 | `ACCEPT_merge + w=1…` |
| 16 | 1280–1359 | False | 0.838 | 0.387 | 0.70 | `ACCEPT_merge + w=1…` |
| 17 | 1360–1439 | True | 1.290 | 0.500 | 0.30 | `REJECT_or_LIMIT_merge + audi…` |
| 18 | 1440–1519 | False | 0.875 | 0.438 | 0.50 | `ACCEPT_merge + w=1…` |
| 19 | 1520–1599 | False | 1.086 | 0.475 | 0.90 | `ACCEPT_merge + w=1…` |
| 20 | 1600–1679 | False | 1.079 | 0.512 | 0.10 | `ACCEPT_merge + w=1…` |
| 21 | 1680–1759 | False | 0.927 | 0.475 | 0.80 | `ACCEPT_merge + w=1…` |
| 22 | 1760–1839 | False | 0.737 | 0.350 | 1.00 | `ACCEPT_merge + w=1…` |
| 23 | 1840–1919 | True | 1.536 | 0.537 | 0.10 | `REJECT_or_LIMIT_merge + audi…` |
| 24 | 1920–1999 | False | 0.884 | 0.475 | 0.90 | `ACCEPT_merge + w=1…` |

- 首次 fire：batch=`4`，ratio≈`1.2222222222222225`
- cut=4 处：fired=`True`，ratio≈`1.2222222222222225`，P@k=`0.1`
- 总 fires：4 / 24 hops

## HashingVectorizer 靠谱吗？Intuition + 替代

| 方案 | Intuition | 何时用 | 风险 |
|------|-----------|--------|------|
| **style 10 维** | 语气/长度/正式度 = P(X) 画像，便宜可审计 | 文风轴、快速探针 | 不解语义偏好 |
| **HashingVectorizer** | 固定维 n-gram 投影，streaming、无词表 | **仅 demo / 无 embedding 时** | 无语义、碰撞噪声；**生产不推荐当唯一特征** |
| **style + hash**（本 demo 默认） | 可解释轴 + 弱字面信号，保证 hop 可复现 | 把 OnlineRFPerm 闭环跑通 | 仍非语义 |
| **sentence-transformers / 业务 embedding** | 语义几何接近裁判/模型内部表征 | **生产默认** | 要模型与版本钉扎 |
| **LLM hidden / reward head** | 与对齐目标同空间 | 已有 reward model 时 | 贵；要防泄漏 |

**结论：** HashingVectorizer **不够靠谱当最终特征**；它只是「没有 GPU/没有 embedding 服务时，仍能把时间序 + OnlineRFPerm 闭环演示出来」的弱代理。你接生产时把 `feature` 列换成你们的 embedding 即可，**检测逻辑一行都不用改**。

## 下一步评估 / 刻画

1. **制度**：fire 率、ratio 分布、cut 对齐（本表）
2. **排序**：po_risk0 P@k / AUROC（相对 y）
3. **动作**：fire→拒合并/审计；quiet→放行
4. **接到¥**：同一 fire 窗 → 客服 act-vs-ignore 账（另见 biz SQL）

## 怎么跑

```bash
PYTHONPATH=. python3 scripts/agod/hh_online_rfperm_stream.py
PYTHONPATH=. python3 scripts/agod/hh_online_rfperm_stream.py --featurizer style_only
```
