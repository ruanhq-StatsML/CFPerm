# 广告漏斗场景：看板 × reason-code 怎么讲

> 具体 pack = **TencentGR**（曝→点→转化边；商户键 ≈ 广告主 `item_feat.122`）。  
> 看板仍是 transfer probe；本页把读数翻成 **广告审出话术**，对齐 S1/S2/S3。

相关：[`Board_Reason_Codes.md`](./Board_Reason_Codes.md) ·
[`Business_Scenarios_Brush_vs_Inject.md`](./Business_Scenarios_Brush_vs_Inject.md) ·
[`Graph_CT_AD_Interface.md`](./Graph_CT_AD_Interface.md)

---

## 1. 场景一张图

```text
投放账户 / 创意 item
        │ 曝光 e_n_exp …
        ▼
     点击 e_ctr …
        ▼
   转化 Y = y_convert          ← 看板监督标签
        │
每 N=1000/2000 条边（按 e_last_ts）切窗
        │
   ops 面板（留量） ──► 买量/曝光强度基线
   content 面板（去量）──► 末跳 credit/share 候选
        │
   sign_Dy(board) = sign(mean Δȳ) ──► S1 刷量 / S2 灌入 / S3 漂移
```

**Justify**

| 层 | 取什么 | 为什么 |
|---|---|---|
| Y | `y_convert` | 广告业务对立轴是转化变多/变少，不是 tip 本身 |
| X_ops | 曝光/点击强度 | 回答「量在不在」——买量基线 |
| X_content | credit/share/covisit | 回答「路径/末跳像不像被操控」——审出候选 |
| 切窗 | 每 N 条边 | 业务粒=样本量，不是日历小时（投放日志不均） |

---

## 2. 和三族怎么对齐

看板给出 `mean_Δȳ`（相邻窗转化率差），压成 `sign_Dy_board`：

| sign | 族 | 广告体感 | 默认队 |
|---|---|---|---|
| **pos** | S1 刷量 | 「这账户/创意突然转化好得离谱」 | 末跳/刷量队（有 credit tip）或弱刷量·爆款辨 |
| **neg** | S2 灌入 | 「曝光还在，转化被差流冲了」 | 差流/落地劣化/打压残留 |
| **flat** | S3 漂移 | 「结构漂了，转化没动」 | 策略/定向/创意轮换，慎升强动作 |

**本机 smoke（N=1000）**：mean_Δȳ ≈ −0.008 → **neg → S2_inject**；  
ops tops = `e_log1p_exp…`；content tops = `i_credit_last / i_share_last`。

读法：**转化在掉**，但 content 仍甩出末跳特征 ——  
先当「差流 + 路径结构异常」盯梢，**不要**在 neg 上讲刷量，也**不要**把强度 AUC≈1 写成「模型可上线」。

---

## 3. reason-code 在广告里的分工

| code | 广告一句话 |
|---|---|
| `RC_AD_FUNNEL_CONTEXT` | 本卡是广告漏斗探针，不是 CTR 出价证明 |
| `RC_AD_BUY_INTENSITY` / `RC_INTENSITY_BASELINE` | 买量强度基线：量在传 |
| `RC_AD_LAST_TOUCH_CANDIDATE` / `RC_CONTENT_CREDIT_CANDIDATE` | 末跳/份额候选 → 映射审出 tip 桶 |
| `RC_OPS_CONTENT_GAP` | 强度解释大半 transfer；创意/落地看残差 |
| `RC_AD_CONVERT_DIP` | 转化下行 → S2 |
| `RC_AD_CONVERT_SPIKE` | 转化上行 → S1（须辨爆款） |
| `RC_AD_STRUCTURE_DRIFT` | 转化平 → S3 |
| `RC_BOARD_NOT_SHIP` | 禁止因探针 AUC 上线 HGB/限投模型 |

---

## 4. 人审读卡顺序（广告）

1. **sign_Dy(board)**：升 / 降 / 平？  
2. **ops 一行**：强度基线（买量）——只解释量。  
3. **content 一行**：credit/share ——末跳候选，进 L1。  
4. **gap**：大 → 先报「大半是量」；小 → 非量结构更值得查。  
5. **永远**：`allows_ship_model=false`；真驱动走 PO/ablation。

---

## 5. 产物

```bash
PYTHONPATH=. python3 scripts/generate_board_reason_codes.py \
  --summary results/sample_chunk_adjacent_board/summary.json \
  --out results/sample_chunk_adjacent_board/reason_codes
```

- `ad_scenario.json` / `ad_scenario_paste.txt`  
- `paste_for_agent.txt` 顶部已嵌广告场景段  

---

## 6. 边界（广告特有）

| 不要宣称 | 为什么 |
|---|---|
| 「AUC 高 = 投放模型可上」 | 探针 ≠ 在线打分；缺延迟标签与校准 |
| 「曝光 tip = 作弊」 | 强度是买量内生，政策/预算都推它 |
| 「credit tip = 定罪末跳」 | shortlist only；需路径与转化极性同向才进强队 |
| 「日历小时切窗」 | 投放日志突发；本场景坚持样本数切窗 |
