# RSI 10 分钟自迭代：模型 / 统计效率（不锁场景）

> 节奏：**每 10 分钟一轮** — 自我更新 · 自我修复 · 自我迭代。  
> 视角：**模型与统计方法**（null、校准、FLOPs、流式探针、多模态 OOD、PO-risk），  
> **不**强制挂广告 S1/S2/S3。发散找机会；持续轮次后再 polish。

账本：[`RSI_Iteration_Log.md`](./RSI_Iteration_Log.md)

---

## 1. 单轮 Protocol（掐表）

| 分 | 动作 | 产出 |
|---:|---|---|
| 0–1 | 抽一张「方法约束卡」 | 本轮不许重复上一轮点 |
| 1–4 | **发散**：写 1 个跨 pack 的模型/统计点子 | 1 句 idea + 假设 |
| 4–7 | **Justify**：用可检的统计/ML 语言说明为什么值 | 定义量 / H0 / 效率 |
| 7–9 | **落地**：改 ≤1 模块 + 测试；自修失败则先修 | diff / green tests |
| 9–10 | 写账本一行 + commit/push；等下一 tick | log 行 |

**方法约束卡（轮流抽）：**
1. 必须引入或使用 **null / 置换 / 重采样** 之一  
2. 必须谈 **校准**（Brier / ECE），不许只报 AUC  
3. 必须谈 **相对 FLOPs / 探针效率**，跨数据集可比  
4. 必须跨 **非广告** pack（DiffDB / metro / PM25 / multimodal）  
5. 必须触及 **PO / IPTW / rfperm** 之一的效率或稳健  
6. 自我修复：优先修红测 / 破文档，不堆新名词  
7. 只 polish 读法/门控文案，不改算法核（本轮）  

---

## 2. 核心统计主张（本栈）

### 2.1 原始 AUC 不可识别「技能」

相邻窗 transfer AUC 回答的是：*在 chunk t 上选出的关联，在 t+1 是否还能排 Y*。  
它在下列情况下会**虚高**：

- 类别极稀 + 分数与先验同向  
- 特征几乎是标签的单调变换（买量强度）  
- 探针过拟合小窗  

**纠正：** label-permutation null —— 固定训练好的分数，打乱 **测试标签**，得 \(AUC_{null}\)。  

\[
\text{excess\_auc} = AUC_{obs} - \mathbb{E}[AUC_{null}]
\]

`excess≈0` ⇒ 看起来很漂亮的 AUC 几乎是机会水平；跨 pack 用同一尺子（广告 / 污染 / token）。

### 2.2 效率 = 技能 / 相对计算

\[
\text{probe\_eff} = \frac{\text{excess\_auc}}{\text{relative\_flops}/10^6}
\]

`relative_flops ≈ n_train · k · (HGB_iters + LR_iters)` —— **适应代价代理**，不是 serving 延迟。  
跨 pack 比的是「每单位探针算力换来多少超出偶然的 transfer」，不是比绝对 AUC。

### 2.3 校准 ≠ 排序（Iter2）

`excess_auc` 说的是排序技能；**ECE**（等宽 bin，transfer 概率）说的是：探针报 \(p\) 时，下一窗频率是否 ≈ \(p\)。

| 读法 | 含义 |
|---|---|
| 高 excess · 低 ECE | 能排且概率可信 |
| 高 excess · 高 ECE | 能排但别当概率用（稀标签/量特征常见） |
| 低 excess | 先别谈校准 |

与 Brier 并存：Brier 是整体 proper score；ECE 局部化偏差。

### 2.4 不确定性：块 bootstrap（Iter3）

相邻对 `(t,t+1),(t+1,t+2)` **共享端点** → 对上的 `excess_auc` 序列相关。  
IID bootstrap 会**低估** mean excess 的方差。

**做法：** moving-block bootstrap（默认 `block_size=2`），对有序 pair 序列重采样，报 90% CI：`excess_auc_ci90` / `hgb_ece_ci90`。

读法：CI 盖住 0 → 别把 pack 均值写成「稳定技能」；跨 pack 比区间宽度，不比点估计。

### 2.5 与业务案由正交

`sign_Dy` 定 S1/S2/S3；excess/probe_eff/ECE/CI 定「故事有没有超出偶然的统计含量 / 概率能不能信 / 均值稳不稳」。  
高 excess 的买量基线仍 **≠ 限投**（业务门另写）。

---

## 3. 发散机会池（不锁广告）

| ID | 点子 | Justify | 状态 |
|---|---|---|---|
| R1 | Transfer null + excess + probe_eff | §2 | **Iter1 落地** |
| R2 | Blocked bootstrap CI on mean excess / ECE | 相邻对共享端点 → 块 bootstrap 保守区间 | **Iter3 落地** |
| R3 | ECE bins on LogReg/HGB transfer | 校准 vs 排序；高 excess+高 ECE=能排不能信 | **Iter2 落地** |
| R4 | DiffusionDB token excess vs Tencent | 弱特征族的 null 对照 | backlog |
| R5 | Streaming reservoir subsample of pairs | 线性扫描 vs 全对；偏差-方差 | backlog |
| R6 | AGOD gate FLOPs vs probe_eff 对照表 | 适应 FLOPs 双尺子 | backlog |
| R7 | PO refit vs ref IPTW MSE / FLOP | 已有脚本串联 scorecard | backlog |
| R8 | Multimodal image-OOD bench 抽 1 指标进账本 | OOD 效率 | backlog |

---

## 4. 代码挂载

- `agod/transfer_null.py` — null / excess / probe_eff  
- `scripts/run_sample_chunk_adjacent_board.py` — 每对自动 enrich  
- `tests/test_transfer_null.py`

---

## 5. 自修复优先级

1. 红测 / import 破  
2. 文档与字段不一致  
3. 新点子（仅当 1–2 绿）
