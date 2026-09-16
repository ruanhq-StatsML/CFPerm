-- =============================================================================
-- 06) 客服助手经营总看板 + 动作推荐
-- 对外只讲：before→after 费率、少工单/退款/多承接、毛/净¥、优先动作。
-- =============================================================================

-- Before / After 对照（三臂并排，经营会第一屏）
CREATE OR REPLACE VIEW vw_cs_assist_before_after AS
SELECT
  r.surface_id,
  ROUND(100.0 * r.ticket_rate_quiet, 2) AS ticket_rate_quiet_pct,
  ROUND(100.0 * r.ticket_rate_ignored, 2) AS ticket_rate_ignored_pct,
  ROUND(100.0 * r.ticket_rate_acted, 2) AS ticket_rate_acted_pct,
  ROUND(100.0 * (r.ticket_rate_ignored - r.ticket_rate_acted), 2) AS ticket_rate_saved_pp,
  ROUND(100.0 * r.refund_rate_quiet, 2) AS refund_rate_quiet_pct,
  ROUND(100.0 * r.refund_rate_ignored, 2) AS refund_rate_ignored_pct,
  ROUND(100.0 * r.refund_rate_acted, 2) AS refund_rate_acted_pct,
  ROUND(100.0 * (r.refund_rate_ignored - r.refund_rate_acted), 2) AS refund_rate_saved_pp,
  ROUND(100.0 * r.contain_rate_quiet, 2) AS contain_rate_quiet_pct,
  ROUND(100.0 * r.contain_rate_ignored, 2) AS contain_rate_ignored_pct,
  ROUND(100.0 * r.contain_rate_acted, 2) AS contain_rate_acted_pct,
  ROUND(100.0 * (r.contain_rate_acted - r.contain_rate_ignored), 2) AS contain_rate_lift_pp,
  e.tickets_avoided,
  e.refunds_avoided,
  e.extra_sessions_contained,
  e.incremental_yen AS gross_yen,
  e.yen_from_tickets,
  e.yen_from_refunds,
  e.yen_from_containment,
  e.monthly_runrate_yen,
  e.days_acted,
  e.days_ignored,
  e.sessions_acted,
  n.net_incremental_yen AS net_yen,
  n.audit_cost_acted,
  n.net_yen_per_1k_sessions
FROM vw_cs_assist_rate_compare r
LEFT JOIN vw_cs_assist_exec_summary e ON e.surface_id = r.surface_id
LEFT JOIN vw_cs_assist_net_increment n ON n.surface_id = r.surface_id;

-- 动作推荐：按「日均净贡献」排序，给 on-call 直接执行
CREATE OR REPLACE VIEW vw_cs_assist_action_recommend AS
SELECT
  surface_id,
  action_type,
  n_days,
  sessions,
  tickets_avoided,
  refunds_avoided,
  extra_contained,
  gross_yen,
  action_day_cost_yen,
  net_yen_after_action_cost AS net_yen,
  ROUND(net_yen_after_action_cost * 1.0 / NULLIF(n_days, 0), 0) AS net_yen_per_day,
  ROUND(gross_yen * 1.0 / NULLIF(n_days, 0), 0) AS gross_yen_per_day,
  avg_rag_hit,
  CASE
    WHEN action_type = 'retrieval_refresh' THEN 'avg_rag_hit 低 → 先刷检索/索引'
    WHEN action_type = 'model_rollback' THEN '生成制度跳变 → 回滚/切稳定版'
    WHEN action_type = 'audit_topk' THEN 'po_risk Top-k 人工审计清差例'
    ELSE '按费率监控'
  END AS playbook,
  ROW_NUMBER() OVER (
    PARTITION BY surface_id
    ORDER BY net_yen_after_action_cost * 1.0 / NULLIF(n_days, 0) DESC
  ) AS recommend_rank
FROM vw_cs_assist_action_net;

-- 经营总看板：一行对外
CREATE OR REPLACE VIEW vw_cs_assist_exec_dashboard AS
SELECT
  b.surface_id,
  b.ticket_rate_quiet_pct,
  b.ticket_rate_ignored_pct,
  b.ticket_rate_acted_pct,
  b.ticket_rate_saved_pp,
  b.refund_rate_saved_pp,
  b.contain_rate_lift_pp,
  b.tickets_avoided,
  b.refunds_avoided,
  b.extra_sessions_contained,
  b.gross_yen,
  b.net_yen,
  b.audit_cost_acted,
  b.monthly_runrate_yen,
  b.days_acted,
  b.sessions_acted,
  a.action_type AS top_action,
  a.net_yen_per_day AS top_action_net_yen_per_day,
  a.playbook AS top_action_playbook,
  ROUND(b.net_yen * 1000.0 / NULLIF(b.sessions_acted, 0), 0) AS net_yen_per_1k,
  ROUND(b.net_yen * 1000.0 / NULLIF(b.sessions_acted, 0) * 10000 * 30 / 1000.0, 0)
    AS prod_mid_monthly_net_yen
FROM vw_cs_assist_before_after b
LEFT JOIN vw_cs_assist_action_recommend a
  ON a.surface_id = b.surface_id AND a.recommend_rank = 1;
