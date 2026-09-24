# 审核 Agent 上下文卡

> 图谱分布变动线索，非定罪结论；供审核分流与上下文，不自动封禁。

- source: `results/tencent_gr_w1w2_mmd_po_fsds/summary.json`
- sign_Dy: **flat** (Dy=None, dy_missing)
- 场景: **漂移族** / 结构漂移
- 建议队列: **路径线性份额变动**
- 读法: Dy 缺失按 flat：X 结构漂了，成功率未证实联动 → 漂移族，禁自称刷量
- 建议动作级: `L1_watch`

## Tips

| feature | sign | bucket |
|---|:---:|---|
| `i_share_linear` | - | 路径线性份额变动 |
| `i_credit_linear` | + | 路径线性归因变动 |
| `i_share_first` | - | 首次份额偏高 |
| `i_credit_first` | + | 首次归因偏高 |
| `i_log1p_n_users` | + | 触达用户规模(log) |
| `i_log1p_n_covisit` | + | 共现规模(log) |
| `i_n_covisit_neighbors` | + | 共现邻域变动（团伙共点） |
| `i_log1p_n_exp` | + | 曝光规模(log) |
| `i_n_users` | + | 触达用户数变动 |
| `i_n_exp` | + | 曝光规模变动 |
| `ui_pop_mismatch` | - | 热度-活跃错配 |
| `u_log1p_n_events` | + | 其它图特征 tip |

## 工单自定义字段（可直接 POST）

```json
{
  "graph_shift_sign_dy": "flat",
  "graph_shift_dy": null,
  "graph_shift_dy_missing": true,
  "graph_shift_scenario_family": "S3_drift",
  "graph_shift_scenario_sub": "structure_shift",
  "graph_shift_queue_bucket": "路径线性份额变动",
  "graph_shift_action_level": "L1_watch",
  "graph_shift_tip_top3": "i_share_linear,i_credit_linear,i_share_first",
  "graph_shift_tip_signs_top3": "i_share_linear:-,i_credit_linear:+,i_share_first:-",
  "graph_shift_localize_k": 200,
  "graph_shift_disclaimer": "clue_not_conviction",
  "graph_shift_gray_allow": true,
  "graph_shift_gray_reason": "ok",
  "graph_shift_tip_industry": "default",
  "graph_shift_sla_level": "urgent",
  "graph_shift_sla_gap_days": 30.0,
  "graph_shift_eta_soft_hint": "lower_confidence"
}
```

## 粘贴给审核 Agent

```
【审核上下文·图谱变动线索】
灰度: allow=True reason=ok
词典: industry=default
方向: sign_Dy=flat Dy=None (dy_missing)
场景: 漂移族 / 结构漂移
支撑: localize_k=200 edges={'W1_train': 962, 'W1_holdout': 314, 'W2': 8815}
建议队列: 路径线性份额变动
SLA: urgent (波次偏老：优先清本队列)
共享: eta_soft_hint=lower_confidence | 图谱在动：优先本队列；ETA 侧建议降置信（只读提示）
读法: Dy 缺失按 flat：X 结构漂了，成功率未证实联动 → 漂移族，禁自称刷量
Tips:
  - i_share_linear (-): 路径线性份额变动
  - i_credit_linear (+): 路径线性归因变动
  - i_share_first (-): 首次份额偏高
  - i_credit_first (+): 首次归因偏高
  - i_log1p_n_users (+): 触达用户规模(log)
  - i_log1p_n_covisit (+): 共现规模(log)
  - i_n_covisit_neighbors (+): 共现邻域变动（团伙共点）
  - i_log1p_n_exp (+): 曝光规模(log)
声明: 图谱分布变动线索，非定罪结论；供审核分流与上下文，不自动封禁。
```
