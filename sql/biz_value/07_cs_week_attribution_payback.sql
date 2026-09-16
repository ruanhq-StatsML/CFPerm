-- =============================================================================
-- 07) 周归因贡献账 + 动作回本天数
-- 财务对账用：哪一周贡献了多少¥；哪个动作几天回本。
-- =============================================================================

-- 不动作臂的基准费率（全样本），作为 acted 日的反事实
CREATE OR REPLACE VIEW vw_cs_assist_ignored_baseline AS
SELECT
  surface_id,
  ticket_per_session AS ignored_ticket_rate,
  refund_per_session AS ignored_refund_rate,
  contained_per_session AS ignored_contain_rate
FROM vw_cs_assist_arm_stats
WHERE arm = 'fire_ignored';

-- 每个 acted 日的反事实贡献（相对 ignored 费率）
CREATE OR REPLACE VIEW vw_cs_assist_acted_day_value AS
SELECT
  d.dt,
  d.surface_id,
  d.n_sessions,
  d.n_tickets,
  d.n_refunds,
  d.n_contained,
  d.ticket_rate,
  d.refund_rate,
  d.containment_rate,
  b.ignored_ticket_rate,
  b.ignored_refund_rate,
  b.ignored_contain_rate,
  (b.ignored_ticket_rate - d.ticket_rate) * d.n_sessions AS tickets_avoided_day,
  (b.ignored_refund_rate - d.refund_rate) * d.n_sessions AS refunds_avoided_day,
  (d.containment_rate - b.ignored_contain_rate) * d.n_sessions AS extra_contained_day,
  (b.ignored_ticket_rate - d.ticket_rate) * d.n_sessions * COALESCE(v.cs_ticket_cost, 0)
    + (b.ignored_refund_rate - d.refund_rate) * d.n_sessions * COALESCE(v.refund_unit_cost, 0)
    + (d.containment_rate - b.ignored_contain_rate) * d.n_sessions * COALESCE(v.contained_session_value, 0)
    AS gross_yen_day
FROM vw_cs_assist_arms d
JOIN vw_cs_assist_ignored_baseline b ON b.surface_id = d.surface_id
LEFT JOIN (
  SELECT
    surface_id,
    MAX(CASE WHEN metric = 'cs_ticket_cost' THEN unit_value END) AS cs_ticket_cost,
    MAX(CASE WHEN metric = 'refund_unit_cost' THEN unit_value END) AS refund_unit_cost,
    MAX(CASE WHEN metric = 'contained_session_value' THEN unit_value END) AS contained_session_value
  FROM dim_value_assumption
  GROUP BY surface_id
) v ON v.surface_id = d.surface_id
WHERE d.arm = 'fire_acted';

-- 按周汇总归因（财务周账）
CREATE OR REPLACE VIEW vw_cs_assist_week_attribution AS
SELECT
  date_trunc('week', CAST(dt AS TIMESTAMP))::DATE AS week_start,
  surface_id,
  COUNT(*) AS days_acted,
  SUM(n_sessions) AS sessions_acted,
  ROUND(SUM(tickets_avoided_day), 1) AS tickets_avoided,
  ROUND(SUM(refunds_avoided_day), 1) AS refunds_avoided,
  ROUND(SUM(extra_contained_day), 1) AS extra_contained,
  ROUND(SUM(gross_yen_day), 0) AS gross_yen,
  ROUND(
    SUM(gross_yen_day) * 100.0 / NULLIF(SUM(SUM(gross_yen_day)) OVER (PARTITION BY surface_id), 0),
    1
  ) AS pct_of_total_gross
FROM vw_cs_assist_acted_day_value
GROUP BY 1, 2
ORDER BY 1, 2;

-- 动作回本：动作日成本 / 日均毛贡献 = 回本天数（<1 表示当天回本）
CREATE OR REPLACE VIEW vw_cs_assist_action_payback AS
SELECT
  a.surface_id,
  a.action_type,
  a.n_days,
  a.gross_yen,
  a.action_day_cost_yen,
  a.net_yen_after_action_cost AS net_yen,
  ROUND(a.gross_yen * 1.0 / NULLIF(a.n_days, 0), 0) AS gross_yen_per_day,
  ROUND(a.action_day_cost_yen * 1.0 / NULLIF(a.n_days, 0), 0) AS cost_yen_per_day,
  ROUND(
    (a.action_day_cost_yen * 1.0 / NULLIF(a.n_days, 0))
      / NULLIF(a.gross_yen * 1.0 / NULLIF(a.n_days, 0), 0),
    3
  ) AS payback_days,
  CASE
    WHEN a.gross_yen * 1.0 / NULLIF(a.n_days, 0) <= 0 THEN 'no_payback'
    WHEN (a.action_day_cost_yen * 1.0 / NULLIF(a.n_days, 0))
           / NULLIF(a.gross_yen * 1.0 / NULLIF(a.n_days, 0), 0) <= 1.0
      THEN 'same_day_payback'
    WHEN (a.action_day_cost_yen * 1.0 / NULLIF(a.n_days, 0))
           / NULLIF(a.gross_yen * 1.0 / NULLIF(a.n_days, 0), 0) <= 3.0
      THEN 'within_3_days'
    ELSE 'slow_payback'
  END AS payback_bucket,
  ROUND(
    a.net_yen_after_action_cost * 1.0 / NULLIF(a.action_day_cost_yen, 0),
    1
  ) AS net_roi_multiple
FROM vw_cs_assist_action_net a;

-- 对账摘要：周归因合计应≈总毛增量
CREATE OR REPLACE VIEW vw_cs_assist_attribution_check AS
SELECT
  e.surface_id,
  e.incremental_yen AS ledger_gross_yen,
  ROUND(SUM(w.gross_yen), 0) AS week_attributed_gross_yen,
  ROUND(SUM(w.gross_yen) - e.incremental_yen, 0) AS attribution_gap_yen,
  ROUND(100.0 * ABS(SUM(w.gross_yen) - e.incremental_yen) / NULLIF(ABS(e.incremental_yen), 0), 2)
    AS attribution_gap_pct
FROM vw_cs_assist_exec_summary e
LEFT JOIN vw_cs_assist_week_attribution w ON w.surface_id = e.surface_id
GROUP BY e.surface_id, e.incremental_yen;
