-- =============================================================================
-- 04) 客服助手增量贡献账 (CS Assistant incremental contribution)
-- 产出只讲：少多少工单、少多少退款、多承接多少会话、贡献多少元。
-- 对照臂：fire_acted vs fire_ignored（同一制度跳变下有没有动作）。
-- =============================================================================

CREATE OR REPLACE VIEW vw_cs_assist_daily AS
SELECT
  e.dt,
  e.surface_id,
  e.model_version,
  COUNT(*) AS n_sessions,
  SUM(
    CASE
      WHEN COALESCE(e.cs_ticketed, 0) = 0 AND COALESCE(e.refunded, 0) = 0 THEN 1
      ELSE 0
    END
  ) AS n_contained,
  SUM(COALESCE(e.cs_ticketed, 0)) AS n_tickets,
  SUM(COALESCE(e.refunded, 0)) AS n_refunds,
  SUM(COALESCE(e.converted, 0)) AS n_conversions,
  SUM(COALESCE(e.gmv, 0)) AS gmv,
  AVG(COALESCE(e.y_halluc, 0)) AS halluc_rate,
  AVG(COALESCE(e.rag_hit, 0)) AS avg_rag_hit
FROM fct_serve_event e
WHERE e.surface_id IN ('shop_assistant', 'cs_bot', 'rag_qa')
GROUP BY 1, 2, 3;

CREATE OR REPLACE VIEW vw_cs_assist_with_signal AS
SELECT
  d.*,
  MAX(CASE WHEN s.signal_type = 'rfperm_fire' AND s.axis = 'concept'
           THEN s.fired ELSE 0 END) AS regime_fired,
  MAX(CASE WHEN s.notes LIKE 'acted%' THEN 1 ELSE 0 END) AS acted,
  MAX(CASE WHEN s.signal_type = 'rfperm_fire' AND s.axis = 'concept'
           THEN s.score END) AS fire_ratio,
  MAX(CASE
        WHEN s.notes LIKE 'acted:%' THEN regexp_extract(s.notes, 'acted:(.*)', 1)
        WHEN s.notes LIKE 'acted%' THEN 'acted'
        ELSE NULL
      END) AS action_type
FROM vw_cs_assist_daily d
LEFT JOIN fct_shift_signal s
  ON s.dt = d.dt
 AND s.surface_id = d.surface_id
 AND (s.model_version = d.model_version OR s.model_version IS NULL)
GROUP BY
  d.dt, d.surface_id, d.model_version, d.n_sessions, d.n_contained,
  d.n_tickets, d.n_refunds, d.n_conversions, d.gmv, d.halluc_rate, d.avg_rag_hit;

CREATE OR REPLACE VIEW vw_cs_assist_rates AS
SELECT
  *,
  n_contained * 1.0 / NULLIF(n_sessions, 0) AS containment_rate,
  n_tickets * 1.0 / NULLIF(n_sessions, 0) AS ticket_rate,
  n_refunds * 1.0 / NULLIF(n_sessions, 0) AS refund_rate,
  n_conversions * 1.0 / NULLIF(n_sessions, 0) AS cvr
FROM vw_cs_assist_with_signal;

CREATE OR REPLACE VIEW vw_cs_assist_arms AS
SELECT
  *,
  CASE
    WHEN regime_fired = 1 AND acted = 1 THEN 'fire_acted'
    WHEN regime_fired = 1 AND acted = 0 THEN 'fire_ignored'
    ELSE 'quiet'
  END AS arm
FROM vw_cs_assist_rates;

-- 按动作类型拆分：检索刷新 / 模型回滚 / Top-k 审计 各自贡献多少
CREATE OR REPLACE VIEW vw_cs_assist_arm_stats AS
SELECT
  surface_id,
  arm,
  COUNT(*) AS n_days,
  SUM(n_sessions) AS sessions,
  SUM(n_contained) AS contained,
  SUM(n_tickets) AS tickets,
  SUM(n_refunds) AS refunds,
  SUM(n_conversions) AS conversions,
  SUM(gmv) AS gmv,
  AVG(containment_rate) AS avg_containment_rate,
  AVG(ticket_rate) AS avg_ticket_rate,
  AVG(refund_rate) AS avg_refund_rate,
  AVG(cvr) AS avg_cvr,
  SUM(n_tickets) * 1.0 / NULLIF(SUM(n_sessions), 0) AS ticket_per_session,
  SUM(n_refunds) * 1.0 / NULLIF(SUM(n_sessions), 0) AS refund_per_session,
  SUM(n_contained) * 1.0 / NULLIF(SUM(n_sessions), 0) AS contained_per_session
FROM vw_cs_assist_arms
GROUP BY surface_id, arm;

CREATE OR REPLACE VIEW vw_cs_assist_action_stats AS
SELECT
  surface_id,
  COALESCE(action_type, 'none') AS action_type,
  COUNT(*) AS n_days,
  SUM(n_sessions) AS sessions,
  SUM(n_tickets) AS tickets,
  SUM(n_refunds) AS refunds,
  SUM(n_contained) AS contained,
  SUM(n_tickets) * 1.0 / NULLIF(SUM(n_sessions), 0) AS ticket_per_session,
  SUM(n_refunds) * 1.0 / NULLIF(SUM(n_sessions), 0) AS refund_per_session,
  SUM(n_contained) * 1.0 / NULLIF(SUM(n_sessions), 0) AS contained_per_session,
  AVG(avg_rag_hit) AS avg_rag_hit
FROM vw_cs_assist_arms
WHERE arm = 'fire_acted'
GROUP BY surface_id, COALESCE(action_type, 'none');

CREATE OR REPLACE VIEW vw_cs_assist_action_increment AS
WITH baseline AS (
  SELECT * FROM vw_cs_assist_arm_stats WHERE arm = 'fire_ignored'
),
v AS (
  SELECT
    surface_id,
    MAX(CASE WHEN metric = 'cs_ticket_cost' THEN unit_value END) AS cs_ticket_cost,
    MAX(CASE WHEN metric = 'refund_unit_cost' THEN unit_value END) AS refund_unit_cost,
    MAX(CASE WHEN metric = 'contained_session_value' THEN unit_value END) AS contained_session_value
  FROM dim_value_assumption
  GROUP BY surface_id
)
SELECT
  a.surface_id,
  a.action_type,
  a.n_days,
  a.sessions,
  a.avg_rag_hit,
  ROUND((b.ticket_per_session - a.ticket_per_session) * a.sessions, 1) AS tickets_avoided,
  ROUND((b.refund_per_session - a.refund_per_session) * a.sessions, 1) AS refunds_avoided,
  ROUND((a.contained_per_session - b.contained_per_session) * a.sessions, 1) AS extra_contained,
  ROUND(
    (b.ticket_per_session - a.ticket_per_session) * a.sessions * COALESCE(v.cs_ticket_cost, 0)
    + (b.refund_per_session - a.refund_per_session) * a.sessions * COALESCE(v.refund_unit_cost, 0)
    + (a.contained_per_session - b.contained_per_session) * a.sessions * COALESCE(v.contained_session_value, 0)
  , 0) AS incremental_yen,
  ROUND(a.ticket_per_session, 4) AS ticket_rate,
  ROUND(b.ticket_per_session, 4) AS ticket_rate_ignored
FROM vw_cs_assist_action_stats a
JOIN baseline b ON a.surface_id = b.surface_id
LEFT JOIN v ON a.surface_id = v.surface_id;

-- 流量情景：按每千会话单价 × 日会话量外推
CREATE OR REPLACE VIEW vw_cs_assist_increment AS
WITH a AS (
  SELECT * FROM vw_cs_assist_arm_stats WHERE arm = 'fire_acted'
),
i AS (
  SELECT * FROM vw_cs_assist_arm_stats WHERE arm = 'fire_ignored'
),
q AS (
  SELECT * FROM vw_cs_assist_arm_stats WHERE arm = 'quiet'
),
v AS (
  SELECT
    surface_id,
    MAX(CASE WHEN metric = 'cs_ticket_cost' THEN unit_value END) AS cs_ticket_cost,
    MAX(CASE WHEN metric = 'refund_unit_cost' THEN unit_value END) AS refund_unit_cost,
    MAX(CASE WHEN metric = 'contained_session_value' THEN unit_value END) AS contained_session_value
  FROM dim_value_assumption
  GROUP BY surface_id
)
SELECT
  COALESCE(a.surface_id, i.surface_id) AS surface_id,
  a.n_days AS days_acted,
  i.n_days AS days_ignored,
  a.sessions AS sessions_acted,
  i.sessions AS sessions_ignored,
  a.ticket_per_session - i.ticket_per_session AS delta_ticket_per_session,
  a.refund_per_session - i.refund_per_session AS delta_refund_per_session,
  a.contained_per_session - i.contained_per_session AS delta_contained_per_session,
  a.avg_containment_rate - i.avg_containment_rate AS delta_containment_rate,
  a.avg_cvr - i.avg_cvr AS delta_cvr,
  (i.ticket_per_session - a.ticket_per_session) * 1000 AS tickets_avoided_per_1k_sessions,
  (i.refund_per_session - a.refund_per_session) * 1000 AS refunds_avoided_per_1k_sessions,
  (a.contained_per_session - i.contained_per_session) * 1000 AS extra_contained_per_1k_sessions,
  (i.ticket_per_session - a.ticket_per_session) * 1000 * COALESCE(v.cs_ticket_cost, 0)
    AS yen_from_tickets_avoided_per_1k,
  (i.refund_per_session - a.refund_per_session) * 1000 * COALESCE(v.refund_unit_cost, 0)
    AS yen_from_refunds_avoided_per_1k,
  (a.contained_per_session - i.contained_per_session) * 1000 * COALESCE(v.contained_session_value, 0)
    AS yen_from_extra_containment_per_1k,
  (i.ticket_per_session - a.ticket_per_session) * 1000 * COALESCE(v.cs_ticket_cost, 0)
    + (i.refund_per_session - a.refund_per_session) * 1000 * COALESCE(v.refund_unit_cost, 0)
    + (a.contained_per_session - i.contained_per_session) * 1000 * COALESCE(v.contained_session_value, 0)
    AS incremental_yen_per_1k_sessions,
  (i.ticket_per_session - a.ticket_per_session) * a.sessions AS tickets_avoided_realized,
  (i.refund_per_session - a.refund_per_session) * a.sessions AS refunds_avoided_realized,
  (a.contained_per_session - i.contained_per_session) * a.sessions AS extra_contained_realized,
  (i.ticket_per_session - a.ticket_per_session) * a.sessions * COALESCE(v.cs_ticket_cost, 0)
    + (i.refund_per_session - a.refund_per_session) * a.sessions * COALESCE(v.refund_unit_cost, 0)
    + (a.contained_per_session - i.contained_per_session) * a.sessions * COALESCE(v.contained_session_value, 0)
    AS incremental_yen_realized,
  a.ticket_per_session - q.ticket_per_session AS acted_vs_quiet_ticket_delta,
  i.ticket_per_session - q.ticket_per_session AS ignored_vs_quiet_ticket_delta,
  v.cs_ticket_cost,
  v.refund_unit_cost,
  v.contained_session_value
FROM a
FULL OUTER JOIN i ON a.surface_id = i.surface_id
LEFT JOIN q ON COALESCE(a.surface_id, i.surface_id) = q.surface_id
LEFT JOIN v ON COALESCE(a.surface_id, i.surface_id) = v.surface_id;

CREATE OR REPLACE VIEW vw_cs_assist_traffic_scenarios AS
SELECT
  i.surface_id,
  s.daily_sessions,
  s.scenario,
  ROUND(i.incremental_yen_per_1k_sessions * s.daily_sessions * 30 / 1000.0, 0) AS monthly_yen,
  ROUND(i.tickets_avoided_per_1k_sessions * s.daily_sessions * 30 / 1000.0, 1) AS monthly_tickets,
  ROUND(i.refunds_avoided_per_1k_sessions * s.daily_sessions * 30 / 1000.0, 1) AS monthly_refunds,
  ROUND(i.extra_contained_per_1k_sessions * s.daily_sessions * 30 / 1000.0, 1) AS monthly_extra_contained
FROM vw_cs_assist_increment i
CROSS JOIN (
  SELECT 200 AS daily_sessions, 'demo_seed' AS scenario
  UNION ALL SELECT 500, 'pilot_small'
  UNION ALL SELECT 2000, 'pilot_mid'
  UNION ALL SELECT 10000, 'prod_mid'
  UNION ALL SELECT 50000, 'prod_large'
) s;

CREATE OR REPLACE VIEW vw_cs_assist_daily_ledger AS
SELECT
  r.dt,
  r.surface_id,
  r.arm,
  r.n_sessions,
  r.n_contained,
  r.n_tickets,
  r.n_refunds,
  r.containment_rate,
  r.ticket_rate,
  r.refund_rate,
  r.n_tickets * COALESCE(v.cs_ticket_cost, 0) AS cs_cost_yen,
  r.n_refunds * COALESCE(v.refund_unit_cost, 0) AS refund_cost_yen,
  r.n_contained * COALESCE(v.contained_session_value, 0) AS containment_value_yen,
  r.n_tickets * COALESCE(v.cs_ticket_cost, 0)
    + r.n_refunds * COALESCE(v.refund_unit_cost, 0) AS quality_cost_yen,
  r.n_contained * COALESCE(v.contained_session_value, 0)
    - r.n_tickets * COALESCE(v.cs_ticket_cost, 0)
    - r.n_refunds * COALESCE(v.refund_unit_cost, 0) AS net_ops_value_yen
FROM vw_cs_assist_arms r
LEFT JOIN (
  SELECT
    surface_id,
    MAX(CASE WHEN metric = 'cs_ticket_cost' THEN unit_value END) AS cs_ticket_cost,
    MAX(CASE WHEN metric = 'refund_unit_cost' THEN unit_value END) AS refund_unit_cost,
    MAX(CASE WHEN metric = 'contained_session_value' THEN unit_value END) AS contained_session_value
  FROM dim_value_assumption
  GROUP BY surface_id
) v ON v.surface_id = r.surface_id;

CREATE OR REPLACE VIEW vw_cs_assist_exec_summary AS
SELECT
  surface_id,
  ROUND(tickets_avoided_realized, 1) AS tickets_avoided,
  ROUND(refunds_avoided_realized, 1) AS refunds_avoided,
  ROUND(extra_contained_realized, 1) AS extra_sessions_contained,
  ROUND(incremental_yen_realized, 0) AS incremental_yen,
  ROUND(incremental_yen_per_1k_sessions, 0) AS incremental_yen_per_1k_sessions,
  ROUND(tickets_avoided_per_1k_sessions, 1) AS tickets_avoided_per_1k,
  ROUND(refunds_avoided_per_1k_sessions, 1) AS refunds_avoided_per_1k,
  ROUND(100.0 * delta_containment_rate, 2) AS containment_rate_lift_pp,
  ROUND(yen_from_tickets_avoided_per_1k * sessions_acted / 1000.0, 0)
    AS yen_from_tickets,
  ROUND(yen_from_refunds_avoided_per_1k * sessions_acted / 1000.0, 0)
    AS yen_from_refunds,
  ROUND(yen_from_extra_containment_per_1k * sessions_acted / 1000.0, 0)
    AS yen_from_containment,
  days_acted,
  days_ignored,
  sessions_acted,
  -- 按 acted 日均会话外推 30 天跑率（给业务周会用）
  ROUND(incremental_yen_per_1k_sessions * (sessions_acted / NULLIF(days_acted, 0)) * 30 / 1000.0, 0)
    AS monthly_runrate_yen,
  ROUND(tickets_avoided_per_1k_sessions * (sessions_acted / NULLIF(days_acted, 0)) * 30 / 1000.0, 1)
    AS monthly_tickets_avoided,
  ROUND(refunds_avoided_per_1k_sessions * (sessions_acted / NULLIF(days_acted, 0)) * 30 / 1000.0, 1)
    AS monthly_refunds_avoided
FROM vw_cs_assist_increment;

-- 费率对照表：quiet / fire_ignored / fire_acted 三臂并排（给业务直接看 before→after）
CREATE OR REPLACE VIEW vw_cs_assist_rate_compare AS
SELECT
  surface_id,
  MAX(CASE WHEN arm = 'quiet' THEN ticket_per_session END) AS ticket_rate_quiet,
  MAX(CASE WHEN arm = 'fire_ignored' THEN ticket_per_session END) AS ticket_rate_ignored,
  MAX(CASE WHEN arm = 'fire_acted' THEN ticket_per_session END) AS ticket_rate_acted,
  MAX(CASE WHEN arm = 'quiet' THEN refund_per_session END) AS refund_rate_quiet,
  MAX(CASE WHEN arm = 'fire_ignored' THEN refund_per_session END) AS refund_rate_ignored,
  MAX(CASE WHEN arm = 'fire_acted' THEN refund_per_session END) AS refund_rate_acted,
  MAX(CASE WHEN arm = 'quiet' THEN contained_per_session END) AS contain_rate_quiet,
  MAX(CASE WHEN arm = 'fire_ignored' THEN contained_per_session END) AS contain_rate_ignored,
  MAX(CASE WHEN arm = 'fire_acted' THEN contained_per_session END) AS contain_rate_acted
FROM vw_cs_assist_arm_stats
GROUP BY surface_id;
