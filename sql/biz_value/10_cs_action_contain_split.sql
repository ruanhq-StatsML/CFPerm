-- =============================================================================
-- 10) 按动作拆多承接贡献（可对账）
-- 业务问题：316 次多承接 / ¥1106 里，rollback / retrieval / audit 各贡献多少？
-- 口径与 vw_cs_assist_action_increment 同一套费率差，只把毛¥拆成三项。
-- =============================================================================

CREATE OR REPLACE VIEW vw_cs_assist_action_value_split AS
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
),
raw AS (
  SELECT
    a.surface_id,
    a.action_type,
    a.n_days,
    a.sessions,
    a.avg_rag_hit,
    (b.ticket_per_session - a.ticket_per_session) * a.sessions AS tickets_avoided,
    (b.refund_per_session - a.refund_per_session) * a.sessions AS refunds_avoided,
    (a.contained_per_session - b.contained_per_session) * a.sessions AS extra_contained,
    (b.ticket_per_session - a.ticket_per_session) * a.sessions * COALESCE(v.cs_ticket_cost, 0)
      AS yen_from_tickets,
    (b.refund_per_session - a.refund_per_session) * a.sessions * COALESCE(v.refund_unit_cost, 0)
      AS yen_from_refunds,
    (a.contained_per_session - b.contained_per_session) * a.sessions
      * COALESCE(v.contained_session_value, 0) AS yen_from_containment,
    COALESCE(v.contained_session_value, 0) AS contain_unit_price
  FROM vw_cs_assist_action_stats a
  JOIN baseline b ON a.surface_id = b.surface_id
  LEFT JOIN v ON a.surface_id = v.surface_id
)
SELECT
  surface_id,
  action_type,
  n_days,
  sessions,
  ROUND(avg_rag_hit, 3) AS avg_rag_hit,
  ROUND(tickets_avoided, 1) AS tickets_avoided,
  ROUND(refunds_avoided, 1) AS refunds_avoided,
  ROUND(extra_contained, 1) AS extra_contained,
  ROUND(yen_from_tickets, 0) AS yen_from_tickets,
  ROUND(yen_from_refunds, 0) AS yen_from_refunds,
  ROUND(yen_from_containment, 0) AS yen_from_containment,
  ROUND(yen_from_tickets + yen_from_refunds + yen_from_containment, 0) AS gross_yen,
  ROUND(
    100.0 * yen_from_containment
      / NULLIF(yen_from_tickets + yen_from_refunds + yen_from_containment, 0),
    1
  ) AS contain_share_pct_of_action_gross,
  ROUND(extra_contained * 1.0 / NULLIF(n_days, 0), 1) AS extra_contained_per_day,
  ROUND(yen_from_containment * 1.0 / NULLIF(n_days, 0), 0) AS contain_yen_per_day,
  contain_unit_price,
  ROUND(
    100.0 * extra_contained
      / NULLIF(SUM(extra_contained) OVER (PARTITION BY surface_id), 0),
    1
  ) AS pct_of_total_extra_contained,
  ROUND(
    100.0 * yen_from_containment
      / NULLIF(SUM(yen_from_containment) OVER (PARTITION BY surface_id), 0),
    1
  ) AS pct_of_total_contain_yen
FROM raw
ORDER BY surface_id, yen_from_containment DESC;

-- 对账：动作承接合计 ≈ exec_summary 多承接 / 承接¥
CREATE OR REPLACE VIEW vw_cs_assist_contain_attribution_check AS
SELECT
  e.surface_id,
  e.extra_sessions_contained AS ledger_extra_contained,
  e.yen_from_containment AS ledger_contain_yen,
  ROUND(SUM(a.extra_contained), 1) AS action_sum_extra_contained,
  ROUND(SUM(a.yen_from_containment), 0) AS action_sum_contain_yen,
  ROUND(SUM(a.extra_contained) - e.extra_sessions_contained, 1) AS contain_count_gap,
  ROUND(SUM(a.yen_from_containment) - e.yen_from_containment, 0) AS contain_yen_gap
FROM vw_cs_assist_exec_summary e
LEFT JOIN vw_cs_assist_action_value_split a ON a.surface_id = e.surface_id
GROUP BY
  e.surface_id,
  e.extra_sessions_contained,
  e.yen_from_containment;
