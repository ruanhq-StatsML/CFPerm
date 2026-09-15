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
