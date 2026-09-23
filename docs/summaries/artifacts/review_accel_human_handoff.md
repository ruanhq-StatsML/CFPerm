# 审出加速 · 真人打标交接纸

> Gate: **NOT_READY** — 需要 1 位真人审核打标后才解冻 ETA/规则 / 新 feature。


## 0. IM 一句话（发给审核同学）

> 帮看一眼图谱变动审出卡（非定罪），有没有帮你少翻页？有空回 useful / not_useful。卡：`docs/summaries/artifacts/review_accel_human_handoff.md` §1；打标命令在 §2，`--reviewer` 用你真名。

## 1. 把下面整段贴进工单 / 发给审核同学

```
【审核上下文·图谱变动线索】
灰度: allow=True reason=ok
词典: industry=default
方向: sign_Dy=flat Dy=None
支撑: localize_k=200 edges={'W1_train': 962, 'W1_holdout': 314, 'W2': 8815}
建议队列: 路径线性份额变动
SLA: urgent (波次偏老：优先清本队列)
共享: eta_soft_hint=lower_confidence | 图谱在动：优先本队列；ETA 侧建议降置信（只读提示）
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

## 2. 打标（把 YOUR_NAME 换成真名）

```bash
PYTHONPATH=. python3 scripts/tencent_gr/log_review_card_feedback.py \
  --card results/tencent_gr_review_agent_card/review_agent_card.json \
  --label useful \
  --reviewer YOUR_NAME \
  --note "一句话：线索有没有帮你少翻页"
```

`not_useful` 同样合法；`reviewer` 不要用 smoke/ood_timer/agent。

## 3. 验收 READY

```bash
PYTHONPATH=. python3 scripts/tencent_gr/check_review_accel_human_gate.py
# exit 0 = READY
```

## 责任人（填空）

- OWNER（催打标）: `________________`
- 审核同学: `________________`
- 目标: 本周内拿到 1 条 `--reviewer` 真名的 useful/not_useful
- 禁止: 自动 `cp` `flags.suggested.json` 覆盖现网；禁止开 ETA/规则包


