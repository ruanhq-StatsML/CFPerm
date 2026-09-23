# X / Y / 特征→业务（concise）

> 审出加速包：效果、商业门禁、打标/预测——先把 **Y、X、特征能不能直接上业务** 说死。

---

## 1. 两套 (X,Y)——别混

| | **算法层**（出卡原料） | **商业层**（门禁 + 打标） |
|---|---|---|
| **Y** | 块上成功：click/convert 均值差 \(\Delta\bar y\) → `sign_Dy`∈{pos,neg,flat} | **主 KPI** 件均审核人分钟↓；**门禁代理** `label`∈{useful, not_useful} |
| **X** | \(K^\star\) 上的图/商户/用户特征（FSDS tip） | 审核看到的**整张卡**（极性+队列+tip 话术） |
| 现在 | tip 有；`sign_Dy` 样例常缺→当 flat | 卡能贴；**n_human=0** → gate NOT_READY |

一句话：算法 Y 回答「成功变多/变少/没动」；商业 Y 回答「人有没有少花时间 / 卡有没有用」。

---

## 2. 选出的特征 → 业务 map（能直接上）

| 特征（X 列） | 业务读法（审核一句话） | 直接动作？ |
|---|---|---|
| `sign_Dy=pos` | 刷量族：成功变多 | **能**：进刷量/末跳队（L1） |
| `sign_Dy=neg` | 灌入族：成功变少 | **能**：进差流/短会话队（L1） |
| `sign_Dy=flat` / 缺 | 漂移族：结构漂、成功没动 | **能**：进供给/活动队；**慎升** L2 |
| `i_share_last` / `i_credit_last` + | 成功压在末跳 | **能**：末跳盯梢话术 |
| `i_share_linear` / `i_credit_linear` | 路径结构份额变了 | **能**：漂移队（当前样例多这个） |
| `u_span_sec` − | 会话变短 | **能**：劣质短会话（配 neg 才像灌入） |
| `i_n_covisit_*` | 共现变密 | **半直接**：只提高抽检，**不定罪团伙** |
| `gap_days` 大 | 波次偏老 | **能**：同队先清（SLA） |

**不能直接上业务的：** 裸 tip 名当「刷量」；无 `sign_Dy` 就 L2/封号；PO/MMD 分数本身不当案由。

---

## 3. 商业置信门禁 = 对商业 Y 的门

```text
Y_biz_proxy = useful|not_useful（真人）
X_biz       = 卡（由算法 X+Y 渲染）

n_human=0  →  不敢信 UR  →  feature_freeze（现在）
n≥1        →  READY_THIN（只改话术）
n≥5        →  才用 UR 预测调灰度
```

刻画：`n_human`（不是 Cover，不是 tip AUC）。  
预测（等人中）：只有「催打标」期望 >0；堆 feature / 拉 r / 解冻 ETA 对人分钟 **不可预测为正**。

---

## 4. 现在效果（一行）

| 层 | 效果 |
|---|---|
| 工程 X→卡 | **有**：出卡/灰度/handoff |
| 算法 Y `sign_Dy` | **弱**：样例常 flat/缺 |
| 特征→业务词典 | **有表**：上表可直接说队列话术 |
| 商业 Y（人分钟/有用） | **无实测**：门禁挡住瞎报 ROI |

---

## 5. 打标在公式里是什么

\[
\hat{\mathrm{UR}}
=\frac{\#\{Y_{\mathrm{biz}}=\mathrm{useful}\}}{\#\{\mathrm{useful}\}+\#\{\mathrm{not\_useful}\}}
\quad\text{（仅真人；smoke 剔除）}
\]

预测调灰度：\(\widehat{\mathrm{harm}}\propto r(1-\mathrm{UR})\)——**Y 必须是真人标**，否则 r 的业务预测无定义。

Handoff：`artifacts/review_accel_human_handoff.md` · 详态：[`Human_Label_Wait_Characterize_Predict.md`](Human_Label_Wait_Characterize_Predict.md)
