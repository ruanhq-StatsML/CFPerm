# 画风/文风漂移 — 迭代 Prototype

> 数据：`Anthropic/hh-rlhf` helpful-base（本地 cache）。  
> **另账**：CTR/品牌路径；禁止并进客服主账 / 偏好合并门禁。

## 对外一句

画风漂移：Before AUC(pre/post)=1.000（fired=True，主因 length/bang）→ 修模板/decoding 后 After=0.535（fired=False）；迭代有效；另账，不动偏好头/客服账。

## 迭代流程

```text
HH answer
  → style_feature(10)          # tidy，只打在 answer
  → detect domain AUC          # pre/post + portrait tercile
  → attribute                  # 文本指标 + length/punct/register
  → style-only ticket          # decoding / 模板 / 创意配比
  → apply fix → re-detect      # 一轮迭代
```

| 步 | 做什么 | 本跑 |
|----|--------|------|
| 1 tidy | answer → 10 维文风 | n=`3200`，cut_batch=`20` |
| 2 detect | pre/post AUC / portrait AUC | Before `1.000` / `0.988` |
| 3 attribute | Top 指标 + 特征维 | `length` / `bang` |
| 4 action | style-only 工单 | `style_only_creative_decoding` |
| 5 iterate | 修后重检 | After AUC `0.535`，fired `True`→`False` |

## Before → After

| 项 | Before | After |
|----|--------|-------|
| style AUC (pre/post) | **1.000** | **0.535** |
| style AUC (portrait) | 0.988 | 0.997 |
| fired | True | False |
| top family / metric | length / bang | register / bang |
| 迭代有效 | — | **True** |

## Round-0 文本指标漂移（Top）

| 指标 | mean_pre | mean_post | Δ |
|------|----------|-----------|---|
| bang | 0.041 | 0.471 | +0.430 |
| hedge | 0.071 | 0.394 | +0.323 |
| tok_len | 0.090 | 0.411 | +0.320 |
| char_len | 0.116 | 0.424 | +0.308 |
| first_person | 0.034 | 0.298 | +0.264 |

## Round-0 特征维份额

| 维 | LOGO share | RF mass |
|----|------------|---------|
| length | 1.000 | 0.136 |
| punct | 0.000 | 0.379 |
| register | 0.000 | 0.485 |

## 动作（另账）

**Do：** 控长度：max_tokens / 模板短答 / 去啰嗦句, 盯 Top 指标 `bang` 回落, 改 decoding / 模板 / 创意配比后重跑 detect  
**Don't：** 不要把 style_auc 当成对齐失败去动偏好头, 不要并进客服工单/退款账, 不要用 PSI(长度) 解释偏好掉点

## 口径

- 画风轴 = P(X) 画像，不是偏好对错 P(Y|X)。  
- 归因先收成特征维 `length/punct/register`（与图谱并维同句式）。  
- 迭代成功判据：pre/post AUC 明显下降，或 fired 熄灭。

## 怎么跑

```bash
PYTHONPATH=. python3 scripts/agod/style_drift_iter_proto.py
PYTHONPATH=. python3 scripts/agod/style_drift_iter_proto.py --synth
# → results/agod/style_drift_iter/
```

相关：`HH_DATA_MANIP_STYLE_PROCEDURES.md` · `FEATURE_DIM_UNIFIED_ATTR.md` · `LANDING_BULLETIN_POLISH.md`
