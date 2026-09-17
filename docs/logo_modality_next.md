# LOGO 多模态：初步探索（可接管）

口径与冻层看板相同：**定位，不是唯一分解**。T 是 batch 标签，不是要识别的处理。PO-risk 问的是 \(P(Y\mid X)\) 有没有 hop。MMD 只对 \(D_{\mathrm{ref}}\)。份额是 ReLU(Δ) 的归一化，不是 Shapley，也不是 CATE。

继续逻辑（这一页的 shape）：

```
新 batch
  → 0. 全局 PO × serving MSE/Brier × MMD²(X_new, X_ref)     （已有看板）
  → 1. 子集 excess（哪一段数据在拖误差）
  → 2. LOGO：每个模态对 Brier/MSE、对 PO-risk、对 MMD 的 Δ → 两层比例
  → 3. 按比例决定下一 batch：全量 / 只推理 / 哪一座塔冻哪一层
  → 4. 记 online Brier、regret、loss curve
```

全局 2× 门安静，**不等于**没有 drift。玩具里只污染一个 slice 时，全 batch 一直 `keep_training`，但子集表已经能标出 slice 3。下一步先子集，再在子集上做 LOGO。

---

## 表 1 · 三个量，两种 LOGO 口径

| 量 | 问什么 | 口径 | 模型 |
|---|---|---|---|
| serving MSE / Brier | 现模型在新 batch 还付不付得起房租 | 回归 MSE；二分类 Brier（概率的 MSE） | 服务模型，可 refit 少一模态 |
| PO-risk | \(P(Y\mid X)\) 有没有 hop | \(\varphi=(Y-\mu)(T-e)\)，T=batch 标签 | 独立 RF nuisance，不是 ViT |
| MMD² | \(P(X)\) 相对参考窗 | \(\mathrm{MMD}^2(X_{\mathrm{new}},X_{\mathrm{ref}})\)，σ 钉在该空间的 \(D_{\mathrm{ref}}\) | 无模型 |

LOGO 对模态 \(g\)：

\[
\Delta L_g = L(\text{full}) - L(-g),\quad
\Delta \mathrm{PO}_g = \mathrm{PO}(\text{full}) - \mathrm{PO}(-g),\quad
\Delta \mathrm{MMD}_g = \mathrm{MMD}(\text{full}) - \mathrm{MMD}(-g).
\]

丢掉 \(g\) 之后度量下降 → \(\Delta>0\) → \(g\) 在往这个度量上送。负的 Δ 记下来，但 **不进份额**（不当成唯一分解的交叉项）。

---

## 表 2 · 两层比例

| 层 | 符号 | 定义 | 什么时候才响 |
|---|---|---|---|
| 1 · 服务误差 | \(\pi^{\mathrm{loss}}_g = \mathrm{ReLU}(\Delta L_g)/\sum_h \mathrm{ReLU}(\Delta L_h)\) | 衰减被定位到哪座塔 | serving 真的变差时。安静时经常全 0。和看板「MSE 先崩」同一句 |
| 2 · 机制 vs 协变量 | \(\pi^{\mathrm{PO}}_g,\ \pi^{\mathrm{MMD}}_g\) 同上归一化；再 \(\mathrm{mix}^{\mathrm{PO}}_g=\pi^{\mathrm{PO}}_g/(\pi^{\mathrm{PO}}_g+\pi^{\mathrm{MMD}}_g)\) | 这座塔更像 \(P(Y\mid X)\) hop 还是 \(P(X)\) 在走 | MMD 对 covariate 最干净；PO 对 concept 较慢、较吵 |

读法：层 2 选塔、选动作；层 1 只做确认。不要用 \(\pi^{\mathrm{loss}}\) 单独开冻层。

---

## 表 3 · 全局动作 × 模态份额 → 下一 batch

全局五格仍是冻层看板（PO × MSE × MMD）。LOGO 只在格子内部回答「哪座塔」。

| 全局 | 该塔 mix | 下一 batch（这座塔） | fusion | 更新 |
|---|---|---|---|---|
| keep | — | `full_train` | `fusion_full` | 全量 |
| watch | — | `infer_only` | `fusion_infer` | 只推理，不冻 |
| tricky | — | `infer_only` | `fusion_infer` | 只推理 |
| x_shift | mix_MMD ≥ 0.6 | `train_stem`（只训底一组） | `fusion_infer` | 不当成 concept 去冻顶 |
| x_shift | quiet | `freeze_tower` | 同上 | 别在安静塔上花梯度 |
| freeze | mix_PO ≥ 0.6 | `train_top`（冻底、训顶） | `fusion_head` | 选择冻层 |
| freeze | quiet | `freeze_tower` | 同上 | 梯度留给 PO-loud 塔 |

`mix` 阈值 0.6 是启发式。未过线的塔 `infer_only`。何时真正上线更新仍是业务。

---

## 表 4 · 塔内层（ViT / Transformer / fusion）

对应现有 `apply_train_top_i`（concept）和这次加的 `apply_train_stem`（covariate）。每座塔自己一组 k。

| 塔 | 层组（建议 k=4） | `train_stem` | `train_top` | `freeze_tower` | `full_train` |
|---|---|---|---|---|---|
| 视频 ViT | patch+pos · early blocks · late blocks · head | 只训 patch/early | 冻 stem，训 late+head | 全冻 | 全训 |
| 音频 (AST / wav2vec) | conv stem · early · late · head | 只训 stem | 冻 stem，训 late+head | 全冻 | 全训 |
| 文本 Transformer | embed · early · late · head | 只训 embed | 冻 embed，训 late+head | 全冻 | 全训 |
| fusion | concat / cross-attn · classifier | 不动（X-shift 时 \(Y\mid X\) 还在） | 训 head | — | 全训 |

FineTuneWrapper 已有 `freeze_backbone` + `unfreeze_last_n`（单塔 train_top）。多塔用 `ModuleDict` 对每个 key 分别 `train_top` / `train_stem`。代码里还没有三塔 net，这是接管时要接的形状。

---

## 表 5 · 玩具：三块特征 ≈ video / audio / text

DGP：12 维，4+4+4。只污染 **slice 3**。onset = batch 2。n_ref=640，n_new=160。Y 二分类，层 1 是 Brier。

| kind | 种下的事 | 全 batch 2× 门 | onset 后该看到 |
|---|---|---|---|
| concept_video | 只在 slice 3 把 video 系数翻号 | 一直 `keep`（污染太稀） | 子集 3 冒头；π_PO 往 video 走 |
| covariate_audio | 只在 slice 3 走 audio 均值 | 一直 `keep` | π_MMD = 0/1/0（audio 吃满） |
| both | 两件同时 | 一直 `keep` | MMD 仍在 audio；PO 后期偏向 video |

实测（节选，完整表 `results/logo_modality/TABLES.md`）：

| kind | t | onset | 全局 | π_loss v/a/t | π_PO v/a/t | π_MMD v/a/t | top subset |
|---|---|---|---|---|---|---|---|
| concept_video | 3 | yes | keep | **1/0/0** | **0.80**/0.04/0.16 | **1/0/0** | **3** |
| concept_video | 4 | yes | keep | 0/0/0 | **0.47**/0.23/0.30 | **1/0/0** | **3** |
| covariate_audio | 2–4 | yes | keep | 0/0/0 | 吵 | **0/1/0** | 不稳 |
| both | 3 | yes | keep | 0/0/0 | **0.83**/0/0.17 | **0/1/0** | **3** |

初步结论：

1. **π_MMD 是最干净的塔定位**（covariate / both 的 audio）。
2. **π_PO 在 concept 上会指向 video，但有滞后、有噪声**（t=2 还在 text）。
3. **π_loss 大部分时候是 0**——服务模型还没崩。只有 concept t=3 给了 video=1，当作确认，不要当开关。
4. **全 batch 门会漏掉 slice-local drift**。concept/both 的 top subset 从 t=2/3 起就是 3；covariate 的 Brier 子集排序不稳，改看子集 MMD（单测已覆盖）。
5. 全局一直 keep → LOGO 策略 = 全量更新，regret 相对 oracle 约 0。策略表要等门被打破、或改成「子集响了也可以 stem/top」才会分叉。

---

## 表 6 · 评估（已接线）

| 指标 | 定义 | 玩具里看到什么 |
|---|---|---|
| online Brier / MSE | 更新**前**在新 batch 上的服务误差 | concept onset 后 frozen 升到 ~0.23；full update 略低 |
| regret | \(\sum_s(L_s^{\mathrm{policy}}-L_s^{\mathrm{oracle}})\) | concept oracle=full update → 0；covariate oracle=frozen → 后期略负（更新有时更好） |
| loss curve | `results/logo_modality/online_brier.png` | 三列 DGP；灰点线 = labeled onset |
| 份额曲线 | `results/logo_modality/logo_shares.png` | π_loss 稀疏；π_PO 在 concept t=3 才集中到 video |

Oracle 必须事先指定（concept → 更新；covariate → 冻住），**不是**从数据里唯一分解出来的。

---

## 表 7 · MSR-VTT 落地形状（未下载视频）

Hub：`friedrichor/MSR-VTT`（train_7k / train_9k / test_1k），字段 `video_id, video, caption, category, url, start time, end time`。音频在 mp4 里，没有单独列。category ∈ {0…19}，正好当子集键。

| 步骤 | 建议 | 不要做的 |
|---|---|---|
| Y | 第一刀：20 类 `category` 的 Brier。检索任务第二刀（对比损失当 serving loss） | 不要把 caption 当 Y 又当 text 模态 |
| 三模态 | 冻住抽 embedding：ViT/CLIP 帧、CLAP/AST 音频、文本 encoder。LOGO 打在三块 embedding 上 | 第一轮不要从像素训 ViT |
| \(D_{\mathrm{ref}}\) | train_7k 抽一次，钉死 | 不要用 rolling pairwise MMD |
| 流 | test_1k 按 category 或时间切 batch，n_new 先 50–100 | 不要 n_new=20 再 bootstrap |
| 子集 | `category`；可选 caption 长度、时长 | 子集 n 太小就跳过（`min_n`） |
| 塔 | embedding 上 LOGO 选塔 → 再对该塔 `train_top` / `train_stem` | 不要三塔共用一个 i* |

---

## 表 8 · 接管时剩下的 shape

| 优先级 | 做什么 | 现成接口 |
|---|---|---|
| 1 | 全局 keep 时也跑子集；子集 excess / 子集 MMD 响了 → 在该子集上再 LOGO | `subset_excess` + `logo_batch` |
| 2 | 子集响、全局还 keep：记下塔，**先不冻全模型**；等 serving Brier 崩再 `train_top`/`train_stem` | 表 3；与「PO 崩 MSE 没崩 → watch」同句 |
| 3 | 真三塔 `ModuleDict`：每塔自己的 `apply_train_top_i` / `apply_train_stem` | `dl_model_registry.FineTuneWrapper`、`_default_multi` |
| 4 | MSR-VTT embedding 缓存 + category 流 | 表 7 |
| 5 | regret 的 oracle 按 DGP/业务写死：concept→更新，X-shift→stem 或冻住 | `cumulative_regret` |
| 6 | 把 LOGO 表接到老板 HTML 作为冻层看板的下一页 | `render_flow_html.py` 还没接 |

代码：

- `Python/src/logo_modality.py` — LOGO、两层比例、`plan_next_batch`
- `Python/src/stream_dgps.py` — `make_trimodal_stream`
- `scripts/run_logo_modality.py` — 出表和曲线
- `tests/test_logo_modality.py`

```bash
PYTHONPATH=Python/src:. python3 -m unittest tests.test_logo_modality
PYTHONPATH=Python/src:. python3 scripts/run_logo_modality.py
```
