-- =============================================================================
-- 17) 全成本净贡献 + delay/留白总捕获
-- 业务问题：扣审计人力 + 动作日成本后还净多少？桌上还留多少（ignored 留白）？
-- 口径：
--   fully_loaded_net = 毛增量 − 审计件成本 − 动作日成本合计
--   uncaptured_gross = ignored 留白毛（已知火情未动作）
--   total_addressable = realized_gross + uncaptured_gross
-- =============================================================================

CREATE OR REPLACE VIEW vw_cs_assist_fully_loaded_capture AS
WITH action_cost AS (
  SELECT
    surface_id,
    ROUND(SUM(action_day_cost_yen), 0) AS action_day_cost_yen_total,
    ROUND(SUM(gross_yen), 0) AS action_gross_yen_check,
    ROUND(SUM(net_yen_after_action_cost), 0) AS action_net_yen_sum
  FROM vw_cs_assist_action_net
  GROUP BY surface_id
),
delay1 AS (
  SELECT surface_id, lost_gross_yen AS delay1_lost_gross_yen, lost_net_yen AS delay1_lost_net_yen
  FROM vw_cs_assist_detection_delay_profit
  WHERE delay_days = 1
)
SELECT
  e.surface_id,
  ROUND(e.incremental_yen, 0) AS realized_gross_yen,
  ROUND(n.audit_cost_acted, 0) AS audit_cost_yen,
  ROUND(COALESCE(a.action_day_cost_yen_total, 0), 0) AS action_day_cost_yen,
  ROUND(n.net_incremental_yen, 0) AS net_after_audit_yen,
  ROUND(
    e.incremental_yen
      - n.audit_cost_acted
      - COALESCE(a.action_day_cost_yen_total, 0),
    0
  ) AS fully_loaded_net_yen,
  ROUND(o.opportunity_gross_yen, 0) AS uncaptured_ignored_gross_yen,
  ROUND(o.days_ignored, 0) AS days_ignored,
  ROUND(o.realized_gross_yen + o.opportunity_gross_yen, 0) AS total_addressable_gross_yen,
  ROUND(
    100.0 * e.incremental_yen
      / NULLIF(e.incremental_yen + o.opportunity_gross_yen, 0),
    1
  ) AS gross_capture_pct,
  ROUND(COALESCE(d.delay1_lost_gross_yen, 0), 0) AS delay1_lost_gross_yen,
  CONCAT(
    '全成本净¥',
    CAST(CAST(ROUND(
      e.incremental_yen - n.audit_cost_acted - COALESCE(a.action_day_cost_yen_total, 0),
      0
    ) AS BIGINT) AS VARCHAR),
    '（毛¥',
    CAST(CAST(ROUND(e.incremental_yen, 0) AS BIGINT) AS VARCHAR),
    '−审计¥',
    CAST(CAST(ROUND(n.audit_cost_acted, 0) AS BIGINT) AS VARCHAR),
    '−动作日¥',
    CAST(CAST(ROUND(COALESCE(a.action_day_cost_yen_total, 0), 0) AS BIGINT) AS VARCHAR),
    '）；已捕获 ',
    CAST(ROUND(
      100.0 * e.incremental_yen / NULLIF(e.incremental_yen + o.opportunity_gross_yen, 0),
      1
    ) AS VARCHAR),
    '%，ignored 留白毛¥',
    CAST(CAST(ROUND(o.opportunity_gross_yen, 0) AS BIGINT) AS VARCHAR),
    '；delay+1 天再少拿毛¥',
    CAST(CAST(ROUND(COALESCE(d.delay1_lost_gross_yen, 0), 0) AS BIGINT) AS VARCHAR),
    '。'
  ) AS external_one_liner_cn
FROM vw_cs_assist_exec_summary e
JOIN vw_cs_assist_net_increment n ON n.surface_id = e.surface_id
JOIN vw_cs_assist_ignored_opportunity o ON o.surface_id = e.surface_id
LEFT JOIN action_cost a ON a.surface_id = e.surface_id
LEFT JOIN delay1 d ON d.surface_id = e.surface_id;

-- 周维度：毛归因 − 当周动作日成本分摊（按 acted 天均摊总动作日成本）
CREATE OR REPLACE VIEW vw_cs_assist_weekly_action_cost_net AS
WITH tot AS (
  SELECT
    surface_id,
    SUM(action_day_cost_yen) AS action_day_cost_yen_total,
    SUM(n_days) AS acted_days_total
  FROM vw_cs_assist_action_net
  GROUP BY surface_id
),
w AS (
  SELECT
    week_start,
    surface_id,
    days_acted,
    yen_from_tickets,
    yen_from_refunds,
    yen_from_containment,
    gross_yen
  FROM vw_cs_assist_week_value_split
)
SELECT
  w.week_start,
  w.surface_id,
  w.days_acted,
  ROUND(w.gross_yen, 0) AS gross_yen,
  ROUND(
    COALESCE(t.action_day_cost_yen_total, 0)
      * w.days_acted
      / NULLIF(t.acted_days_total, 0),
    0
  ) AS action_day_cost_yen_alloc,
  ROUND(
    w.gross_yen
      - COALESCE(t.action_day_cost_yen_total, 0)
          * w.days_acted
          / NULLIF(t.acted_days_total, 0),
    0
  ) AS net_yen_after_action_cost_alloc,
  ROUND(w.yen_from_tickets, 0) AS yen_from_tickets,
  ROUND(w.yen_from_refunds, 0) AS yen_from_refunds,
  ROUND(w.yen_from_containment, 0) AS yen_from_containment
FROM w
LEFT JOIN tot t ON t.surface_id = w.surface_id
ORDER BY w.week_start;
