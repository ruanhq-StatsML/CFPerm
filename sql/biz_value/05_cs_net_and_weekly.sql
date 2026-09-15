-- =============================================================================
-- 05) 客服助手净贡献（扣审计人力）+ 周经营看板
-- =============================================================================

-- 审计人力：从 fct_audit_queue 按日计件，单价来自 dim_value_assumption
CREATE OR REPLACE VIEW vw_cs_assist_audit_daily AS
SELECT
  a.dt,
  a.surface_id,
  COUNT(*) AS n_audits,
  SUM(CASE WHEN a.audit_status = 'confirmed_bad' THEN 1 ELSE 0 END) AS n_confirmed_bad,
  SUM(CASE WHEN a.audit_status = 'false_alarm' THEN 1 ELSE 0 END) AS n_false_alarm,
  SUM(CASE WHEN a.audit_status = 'confirmed_bad' THEN 1 ELSE 0 END) * 1.0
    / NULLIF(COUNT(*), 0) AS precision_proxy
FROM fct_audit_queue a
WHERE a.surface_id IN ('shop_assistant', 'cs_bot', 'rag_qa')
GROUP BY 1, 2;

CREATE OR REPLACE VIEW vw_cs_assist_net_daily AS
SELECT
  l.dt,
  l.surface_id,
  l.arm,
  l.n_sessions,
  l.n_tickets,
  l.n_refunds,
  l.n_contained,
  l.quality_cost_yen,
  l.containment_value_yen,
  l.net_ops_value_yen,
  COALESCE(a.n_audits, 0) AS n_audits,
  COALESCE(a.n_confirmed_bad, 0) AS n_confirmed_bad,
  COALESCE(a.precision_proxy, 0) AS audit_precision,
  COALESCE(a.n_audits, 0) * COALESCE(v.audit_unit_cost, 0) AS audit_cost_yen,
  l.net_ops_value_yen
    - COALESCE(a.n_audits, 0) * COALESCE(v.audit_unit_cost, 0) AS net_contrib_yen
FROM vw_cs_assist_daily_ledger l
LEFT JOIN vw_cs_assist_audit_daily a
  ON a.dt = l.dt AND a.surface_id = l.surface_id
LEFT JOIN (
  SELECT
    surface_id,
    MAX(CASE WHEN metric = 'audit_unit_cost' THEN unit_value END) AS audit_unit_cost
  FROM dim_value_assumption
  GROUP BY surface_id
) v ON v.surface_id = l.surface_id;

CREATE OR REPLACE VIEW vw_cs_assist_net_by_arm AS
SELECT
  surface_id,
  arm,
  COUNT(*) AS n_days,
  SUM(n_sessions) AS sessions,
  SUM(n_tickets) AS tickets,
  SUM(n_refunds) AS refunds,
  SUM(n_audits) AS audits,
  ROUND(AVG(audit_precision), 3) AS avg_audit_precision,
  ROUND(SUM(quality_cost_yen), 0) AS quality_cost_yen,
  ROUND(SUM(audit_cost_yen), 0) AS audit_cost_yen,
  ROUND(SUM(net_ops_value_yen), 0) AS net_ops_value_yen,
  ROUND(SUM(net_contrib_yen), 0) AS net_contrib_yen
FROM vw_cs_assist_net_daily
GROUP BY surface_id, arm;

-- 相对 ignored：acted 臂的净增量（已扣审计人力）
-- 口径：运营毛增量（费率差×会话×单价）− acted 侧审计人力成本
-- 不把 ignored 臂多烧的审计当成“贡献”（避免净>毛的反直觉）。
CREATE OR REPLACE VIEW vw_cs_assist_net_increment AS
WITH a AS (
  SELECT * FROM vw_cs_assist_net_by_arm WHERE arm = 'fire_acted'
),
i AS (
  SELECT * FROM vw_cs_assist_net_by_arm WHERE arm = 'fire_ignored'
),
g AS (
  SELECT surface_id, incremental_yen_realized, sessions_acted
  FROM vw_cs_assist_increment
)
SELECT
  COALESCE(a.surface_id, i.surface_id) AS surface_id,
  a.n_days AS days_acted,
  i.n_days AS days_ignored,
  a.sessions AS sessions_acted,
  a.audits AS audits_acted,
  a.avg_audit_precision,
  a.audit_cost_yen AS audit_cost_acted,
  i.quality_cost_yen - a.quality_cost_yen AS quality_cost_saved,
  a.net_ops_value_yen - i.net_ops_value_yen AS ops_value_lift,
  ROUND(g.incremental_yen_realized - COALESCE(a.audit_cost_yen, 0), 0) AS net_incremental_yen,
  ROUND(
    (g.incremental_yen_realized - COALESCE(a.audit_cost_yen, 0))
      * 1000.0 / NULLIF(a.sessions, 0),
    0
  ) AS net_yen_per_1k_sessions
FROM a
FULL OUTER JOIN i ON a.surface_id = i.surface_id
LEFT JOIN g ON COALESCE(a.surface_id, i.surface_id) = g.surface_id;

-- 周经营看板（ISO 周）
CREATE OR REPLACE VIEW vw_cs_assist_weekly AS
SELECT
  date_trunc('week', CAST(dt AS TIMESTAMP))::DATE AS week_start,
  surface_id,
  SUM(n_sessions) AS sessions,
  SUM(n_tickets) AS tickets,
  SUM(n_refunds) AS refunds,
  SUM(n_contained) AS contained,
  SUM(n_audits) AS audits,
  ROUND(SUM(n_tickets) * 1.0 / NULLIF(SUM(n_sessions), 0), 4) AS ticket_rate,
  ROUND(SUM(n_refunds) * 1.0 / NULLIF(SUM(n_sessions), 0), 4) AS refund_rate,
  ROUND(SUM(n_contained) * 1.0 / NULLIF(SUM(n_sessions), 0), 4) AS containment_rate,
  ROUND(SUM(quality_cost_yen), 0) AS quality_cost_yen,
  ROUND(SUM(audit_cost_yen), 0) AS audit_cost_yen,
  ROUND(SUM(net_contrib_yen), 0) AS net_contrib_yen,
  SUM(CASE WHEN arm = 'fire_acted' THEN 1 ELSE 0 END) AS days_acted,
  SUM(CASE WHEN arm = 'fire_ignored' THEN 1 ELSE 0 END) AS days_ignored
FROM vw_cs_assist_net_daily
GROUP BY 1, 2
ORDER BY 1, 2;

-- 火情周环比（WoW）：只看有 fire 的周，给经营简报用
CREATE OR REPLACE VIEW vw_cs_assist_weekly_wow AS
SELECT
  week_start,
  surface_id,
  sessions,
  tickets,
  refunds,
  contained,
  ticket_rate,
  refund_rate,
  containment_rate,
  net_contrib_yen,
  days_acted,
  days_ignored,
  LAG(ticket_rate) OVER (PARTITION BY surface_id ORDER BY week_start) AS prev_ticket_rate,
  LAG(containment_rate) OVER (PARTITION BY surface_id ORDER BY week_start) AS prev_containment_rate,
  LAG(net_contrib_yen) OVER (PARTITION BY surface_id ORDER BY week_start) AS prev_net_contrib_yen,
  ROUND(
    ticket_rate - LAG(ticket_rate) OVER (PARTITION BY surface_id ORDER BY week_start),
    4
  ) AS delta_ticket_rate,
  ROUND(
    containment_rate - LAG(containment_rate) OVER (PARTITION BY surface_id ORDER BY week_start),
    4
  ) AS delta_containment_rate,
  ROUND(
    net_contrib_yen - LAG(net_contrib_yen) OVER (PARTITION BY surface_id ORDER BY week_start),
    0
  ) AS delta_net_contrib_yen
FROM vw_cs_assist_weekly
WHERE days_acted + days_ignored > 0;

-- 动作毛增量 − 动作日成本 − 审计件成本分摊（仅 audit_topk）
CREATE OR REPLACE VIEW vw_cs_assist_action_net AS
WITH costs AS (
  SELECT
    surface_id,
    MAX(CASE WHEN metric = 'action_cost_retrieval_refresh_day' THEN unit_value END) AS cost_retrieval_day,
    MAX(CASE WHEN metric = 'action_cost_model_rollback_day' THEN unit_value END) AS cost_rollback_day,
    MAX(CASE WHEN metric = 'action_cost_audit_topk_day' THEN unit_value END) AS cost_audit_day,
    MAX(CASE WHEN metric = 'audit_unit_cost' THEN unit_value END) AS audit_unit_cost
  FROM dim_value_assumption
  GROUP BY surface_id
),
aud AS (
  SELECT surface_id, COUNT(*) AS n_audits
  FROM fct_audit_queue
  WHERE surface_id IN ('shop_assistant', 'cs_bot', 'rag_qa')
  GROUP BY surface_id
)
SELECT
  a.surface_id,
  a.action_type,
  a.n_days,
  a.sessions,
  a.tickets_avoided,
  a.refunds_avoided,
  a.extra_contained,
  a.incremental_yen AS gross_yen,
  ROUND(
    a.n_days * CASE a.action_type
      WHEN 'retrieval_refresh' THEN COALESCE(c.cost_retrieval_day, 0)
      WHEN 'model_rollback' THEN COALESCE(c.cost_rollback_day, 0)
      WHEN 'audit_topk' THEN COALESCE(c.cost_audit_day, 0)
      ELSE 0
    END
  , 0) AS action_day_cost_yen,
  ROUND(
    a.incremental_yen
    - a.n_days * CASE a.action_type
        WHEN 'retrieval_refresh' THEN COALESCE(c.cost_retrieval_day, 0)
        WHEN 'model_rollback' THEN COALESCE(c.cost_rollback_day, 0)
        WHEN 'audit_topk' THEN COALESCE(c.cost_audit_day, 0)
        ELSE 0
      END
    - CASE
        WHEN a.action_type = 'audit_topk'
          THEN COALESCE(aud.n_audits, 0) * COALESCE(c.audit_unit_cost, 0)
             * a.sessions / NULLIF(SUM(a.sessions) OVER (PARTITION BY a.surface_id), 0)
        ELSE 0
      END
  , 0) AS net_yen_after_action_cost,
  a.avg_rag_hit
FROM vw_cs_assist_action_increment a
LEFT JOIN costs c ON c.surface_id = a.surface_id
LEFT JOIN aud ON aud.surface_id = a.surface_id;
