# 审出加速：角色痛点 elaborations + 真人反馈闭环

> 补全 [`Scope_Review_Agent_Card.md`](Scope_Review_Agent_Card.md) §0.2：还有哪些角色、痛点怎么拆、下一步怎么走；  
> 以及 **为什么必须闭环真人反馈**、具体怎么转一圈、什么算 READY。  
> 仍停在 L1；不改 Drill；不解冻 ETA/规则包。

---

## 1. 角色图谱：比「五人表」再全一点

| # | 角色 | 是否本包主用户 | 今天痛点（具体） | 本包触点 | 下一步（按优先级） |
|---|---|---|---|---|---|
| R1 | **一线审核同学** | **主** | 打开工单后要自己猜「这波像刷量还是供给漂」；翻图/规则/历史耗时 | 卡 paste、队列桶、读法 | ① 贴 1 张真卡审 ② `useful/not_useful` ③ 记 note「少翻了哪一页」 |
| R2 | **审核 Agent** | **主** | 上下文只有 ID/规则命中，没有图谱变动叙事 | `paste_for_agent`、工单字段 | ① 工单 POST dry-run→真 URL ② Agent 提示词强制读 `sign_Dy`+队列 ③ 回流同 R1 |
| R3 | **审核组长 / 值班长** | 主（管理） | 不知道哪类队列有用率高、老波次谁清 | useful_rate 按队列、SLA | ① 周看 `summarize_review_feedback` ② 低有用率队列降 `queue_sample_rates` ③ 定 OWNER |
| R4 | **风控策略 / 灰度负责人** | 辅 | 怕一刀切误伤、缺回滚把手 | kill_switch、sample_rate 建议 | ① 人工确认 suggested flags ② 单队列缩流试点 ③ **禁止**自动覆盖现网 |
| R5 | **监测 / 算法同学** | 辅 | `summary.json` 没人看；direction 偶发缺失 | export / smoke / batch | ① 保证 direction 写入 ② 波次批出卡进值班目录 ③ 用反馈反查 tip 词典是否胡说 |
| R6 | **内容安审**（复用） | 旁路 | 另建一套 UI 成本高 | tip overlay | READY 后再扩；现阶段只验证「同一卡换词典可读」 |
| R7 | **调查 / 案件（L3）** | **非主** | 要证据链与资金动作 | 卡只作时间线附件 | **Out**：不自动开案；人要拷卡自行升级 |
| R8 | **ETA / 调度** | **非主** | 异常日客诉 | `eta_soft_hint` 只读 | **冻结**至 human-gate READY；本包不改 ETA 服务 |
| R9 | **规则引擎同学** | **非主** | 规则∩结构异常难合流 | 无接通 | **冻结**；解冻后另包「规则∩K\*」 |
| R10 | **商户运营 / 客服** | 旁听 | 被限频时要可解释话术 | 暂无 | Out；L2 解冻后用 direction 模板（B10） |
| R11 | **平台 / 工单 SOAR** | 辅 | 要稳定字段契约 | `ticket_custom_fields` + schema 习惯 | 字段名冻结；dry-run→真 webhook；改字段要 version |

**主路径只服务 R1–R5**；R7–R10 知道存在，但不在本 scope 交付。

---

## 2. 痛点怎么 elaborate →「下一步」模板

每个痛点拆成四格，避免空喊：

```text
痛点: <一句话现象>
证据: <现在缺什么字段/动作>
本包已有: <脚本/字段>
下一步 1 周内: <可点的一件事 + 验收人>
不做: <明确 Out>
```

### 2.1 R1 一线审核（核心）
| | |
|---|---|
| 痛点 | 不知道先走 S1 末跳还是 S3 供给漂移，来回切系统 |
| 证据 | 无卡时工单只有 ID；有卡但没人说「有没有用」 |
| 已有 | S1–S7 读法、handoff 粘贴段 |
| 下一步 | 组长指定 1 人审 **当前 flat 样例卡**；打标+note；验收= gate 变 READY |
| 不做 | 为「好用」再加 Web UI / 新模型 |

### 2.2 R2 审核 Agent
| | |
|---|---|
| 痛点 | Agent 胡编队列或忽略 flat「慎升」 |
| 证据 | 无强制字段进 prompt |
| 已有 | `paste_for_agent`、`ticket_custom_fields` |
| 下一步 | Agent 系统提示：必须引用 `graph_shift_sign_dy` + `queue_bucket`；缺字段则拒答升权 |
| 不做 | 让 Agent 直接 L2 |

### 2.3 R3 值班长
| | |
|---|---|
| 痛点 | 无法按队列调人力 |
| 证据 | useful_rate 无真人切片 |
| 已有 | `summarize_review_feedback.py` |
| 下一步 | 凑满 ≥5 条真人反馈后出首周报；低有用队列写进 `queue_sample_rates` |
| 不做 | 没样本就改 tip 排序算法 |

### 2.4 R4 灰度
| | |
|---|---|
| 痛点 | 全开怕炸、全关没价值 |
| 证据 | 只有全局 sample_rate 时粒度粗 |
| 已有 | `queue_sample_rates`、suggest CLI |
| 下一步 | 对「路径线性份额」等样例队列先 100% 开给试点组，其它 0 |
| 不做 | 脚本自动 `cp` suggested→现网 |

### 2.5 R5 算法
| | |
|---|---|
| 痛点 | direction 空导致全 flat，业务误读 |
| 证据 | 当前 `summary.json` 可无 `direction` 键 |
| 已有 | `direction_report`、export 缺省 flat |
| 下一步 | **工程债**：localize 跑完强制 attach direction；验收=卡上 Dy 非空或显式 `dy_missing=true` |
| 不做 | 用 AUC 证明 direction「对」 |

---

## 3. 真人反馈闭环（具体谈）

### 3.1 为什么必须有闭环
商业计分（见 Motivate 文）：

\[
\mathrm{Score}=1_{\text{现网可点}}\cdot(1_{\text{接通}}+1_{\text{有人用}}+1_{\text{ROI粗数}})
\]

| 只有工程接通 | Score≤1 |
|---|---|
| + 真人用过并打标 | → 有人用 |
| + 有用率/耗时粗数 | → ROI |

**没有真人反馈 = 永远停在「可演示」，进不了「可卖」。**  
所以 human-gate 把「再堆 B35」卡死，不是偷懒，是逼 Score 从 1→2。

### 3.2 一圈里有什么（现成脚本）
```text
① 出卡   export_review_agent_card.py
② 贴卡   工单 / Agent / handoff §1
③ 人审   按 S1–S7 读法分流（可改队）
④ 打标   log_review_card_feedback.py
         --label useful|not_useful
         --reviewer <真名>          ← 关键：非 smoke/ood_timer
         --note "少翻了X / 队列错了"
⑤ 汇总   summarize_review_feedback.py
         → overall + by queue_bucket useful_rate
⑥ 建议   suggest_review_card_sample_rate.py
         → flags.suggested.json（人工确认）
⑦ 回写灰度  人改 configs/review_agent_card_flags.json
         或 queue_sample_rates 单队列缩/扩
⑧ 再出卡  下一波次 → ①
```

Gate：`check_review_accel_human_gate.py` 读的是 ④ 里 **真人** 行数。

### 3.3 一条反馈长什么样（契约）
```json
{
  "ts": "...",
  "label": "useful | not_useful",
  "note": "人话：哪里有用/坑",
  "reviewer": "真名",
  "sign_Dy": "flat|pos|neg",
  "queue_bucket": "路径线性份额变动",
  "tips_top3": ["i_share_linear", "..."],
  "ticket_fields_echo": { "...": "..." }
}
```

| 字段 | 闭环用途 |
|---|---|
| `label` | 有用率分子分母 |
| `reviewer` | 过滤 bot；READY 判定 |
| `queue_bucket` | 按队切片；喂 `queue_sample_rates` |
| `sign_Dy` / tips | 以后改词典、查「flat 是否总被标 not_useful」 |
| `note` | 质性：少翻页？队列错？话术看不懂？ |

### 3.4 useful vs not_useful 分别驱动什么
| 标签 | 含义（约定） | 系统侧下一步 | 产品侧下一步 |
|---|---|---|---|
| **useful** | 卡帮你少翻页或分流对了 | 该队列 sample_rate 可维持/略升 | 记「哪句读法有用」进词典 |
| **not_useful** | 干扰、队列错、或废话 | 该队列降 sample_rate；连续差可 kill_switch | 改桶话术 / 查 direction 是否空 |
| 不打标 | 等于没用过 | gate 永不 READY | **禁止**宣称人分钟↓ |

**注意**：useful ≠「真是欺诈」。只评 **线索对审出是否省事**。

### 3.5 闭环状态机（给值班长看）
```text
[工程可点] ──贴卡──► [等待真人]
                │
                ▼
         打标 n_human≥1 ──► [READY / Score≥2]
                │                │
                │                ├─ 周报 useful_rate
                │                ├─ 调 queue_sample_rates
                │                └─ 才议 ETA/规则解冻
                ▼
         一直 0 真人 ──► feature_freeze + 6h 催打标
                         （现在就在这）
```

### 3.6 最小真人实验设计（1 周内可做完）
| 步 | 谁 | 做什么 | 完成定义 |
|---|---|---|---|
| 1 | OWNER（填 handoff） | 指定 1 名审核 | 名字写进 handoff |
| 2 | 算法 | 导出当前样例卡（S3 flat） | md 可打开 |
| 3 | 审核 | 假装审 1 单：只许用卡+必要系统 | 口头说多花/少花几分钟 |
| 4 | 审核 | CLI 打标 + note | `reviewer=真名` 进 jsonl |
| 5 | OWNER | 跑 gate | exit 0 = READY |
| 6 | 值班长 | 若 not_useful：改话术或查 direction 债 | 再请第 2 人，不堆 feature |

样本量：READY 门槛 **1**；调灰度建议 **≥5**；报「人分钟↓」至少 **一小队一周** 前后对比（工单时长或自评分钟）。

### 3.7 反模式（闭环常见翻车）
| 反模式 | 为何坏 | 正确做法 |
|---|---|---|
| `reviewer=ood_timer/smoke` 算有人用 | Score 造假 | gate 黑名单过滤 |
| 自动覆盖 flags | 不可回滚 | 只写 suggested |
| 用 AUC 代替 useful | 与人分钟脱节 | 只认打标与耗时 |
| 没打标继续加 overlay/UI | 扩面逃避验证 | freeze |
| useful 当定罪 | 法律/误伤 | 文案保持 clue_not_conviction |

---

## 4. 「还有其他的吗」——场景/角色增量清单

### 4.1 可加但 **仍属审出加速**（解冻后或 READY 后小步）
- 审核自评「本单多/少花了几分钟」数值字段（比二元 label 更贴 KPI）  
- 队列错标：`label=not_useful` + note 含「应去 S1」→ 以后校准词典  
- 双人抽检：同卡两人打标看一致率（降话术歧义）  

### 4.2 像相关但 **不是本包**（保持 Kill/Backlog）
- ETA 降置信真接通、规则∩K\*、名单 TTL、商户 App 推送、支付失败侧通道、退货边…  

### 4.3 现在唯一缺口（比「再想场景」更要紧）
1. **真人打标 = 0** → 闭环断在 ④  
2. **direction 可能空** → S3 过多「假 flat」，反馈会被噪声污染  

---

## 5. Takeaway

1. 角色不止五个：主路径 R1–R5；L3/ETA/规则/商户是旁路或冻结。  
2. 每个痛点都要落到「1 周内可点的下一步 + 验收人」。  
3. 真人反馈闭环 = 出卡→贴卡→打标→汇总→**人工**调灰度→再出卡；gate 卡的是「有人用」。  
4. useful 评的是 **省审出时间**，不是案情真伪。  
5. 下一步不是新场景脑暴，是 **跑通 §3.6 最小真人实验**。
