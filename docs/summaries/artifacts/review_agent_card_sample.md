# 审核 Agent 上下文卡

> 图谱分布变动线索，非定罪结论；供审核分流与上下文，不自动封禁。

- source: `results/tencent_gr_w1w2_mmd_po_fsds/summary.json`
- sign_Dy: **flat** (Dy=None)
- 建议队列: **路径线性份额变动**
- 读法: X 漂了但 click 率平：先当供给/分布漂移，慎升强动作
- 建议动作级: `L1_watch`

## Tips

| feature | sign | bucket |
|---|:---:|---|
| `i_share_linear` | 0 | 路径线性份额变动 |
| `i_credit_linear` | 0 | 路径线性归因变动 |
| `i_share_first` | 0 | 首次份额偏高 |
| `i_credit_first` | 0 | 首次归因偏高 |
| `i_log1p_n_users` | 0 | 触达用户规模(log) |
| `i_log1p_n_covisit` | 0 | 共现规模(log) |
| `i_n_covisit_neighbors` | 0 | 共现邻域变动（团伙共点） |
| `i_log1p_n_exp` | 0 | 曝光规模(log) |
| `i_n_users` | 0 | 触达用户数变动 |
| `i_n_exp` | 0 | 曝光规模变动 |
| `ui_pop_mismatch` | 0 | 热度-活跃错配 |
| `u_log1p_n_events` | 0 | 其它图特征 tip |

## 粘贴给审核 Agent

```
【审核上下文·图谱变动线索】
方向: sign_Dy=flat Dy=None
支撑: localize_k=200 edges={'W1_train': 962, 'W1_holdout': 314, 'W2': 8815}
建议队列: 路径线性份额变动
读法: X 漂了但 click 率平：先当供给/分布漂移，慎升强动作
Tips:
  - i_share_linear (0): 路径线性份额变动
  - i_credit_linear (0): 路径线性归因变动
  - i_share_first (0): 首次份额偏高
  - i_credit_first (0): 首次归因偏高
  - i_log1p_n_users (0): 触达用户规模(log)
  - i_log1p_n_covisit (0): 共现规模(log)
  - i_n_covisit_neighbors (0): 共现邻域变动（团伙共点）
  - i_log1p_n_exp (0): 曝光规模(log)
声明: 图谱分布变动线索，非定罪结论；供审核分流与上下文，不自动封禁。
```
