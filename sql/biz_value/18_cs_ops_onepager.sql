-- =============================================================================
-- 18) 值班一页纸：Before/After + 全成本净 + delay 损失 + 优先动作
-- 业务问题：开会/值班一眼看到可对账数字，不用翻 10 张表
-- =============================================================================

CREATE OR REPLACE VIEW vw_cs_assist_ops_onepager AS
WITH top_act AS (
  SELECT
    surface_id,
    action_type AS top_action,
    ROUND(net_yen_after_action_cost, 0) AS top_action_net_yen,
    ROUND(net_yen_after_action_cost * 1.0 / NULLIF(n_days, 0), 0) AS top_action_net_yen_per_day
  FROM (
    SELECT
      surface_id,
      action_type,
      net_yen_after_action_cost,
      n_days,
      ROW_NUMBER() OVER (
        PARTITION BY surface_id
        ORDER BY net_yen_after_action_cost * 1.0 / NULLIF(n_days, 0) DESC
      ) AS rn
    FROM vw_cs_assist_action_net
  ) z
  WHERE rn = 1
),
delay AS (
  SELECT
    surface_id,
    MAX(CASE WHEN delay_days = 0 THEN lost_gross_yen END) AS delay0_lost_gross,
    MAX(CASE WHEN delay_days = 1 THEN lost_gross_yen END) AS delay1_lost_gross,
    MAX(CASE WHEN delay_days = 3 THEN lost_gross_yen END) AS delay3_lost_gross,
    MAX(gross_yen_per_acted_day) AS gross_yen_per_acted_day
  FROM vw_cs_assist_detection_delay_profit
  GROUP BY surface_id
)
SELECT
  f.surface_id,
  ROUND(100.0 * r.ticket_rate_ignored, 2) AS ticket_rate_ignored_pct,
  ROUND(100.0 * r.ticket_rate_acted, 2) AS ticket_rate_acted_pct,
  ROUND(100.0 * r.refund_rate_ignored, 2) AS refund_rate_ignored_pct,
  ROUND(100.0 * r.refund_rate_acted, 2) AS refund_rate_acted_pct,
  ROUND(100.0 * r.contain_rate_ignored, 2) AS contain_rate_ignored_pct,
  ROUND(100.0 * r.contain_rate_acted, 2) AS contain_rate_acted_pct,
  ROUND(
    100.0 * (r.contain_rate_acted - r.contain_rate_ignored),
    2
  ) AS contain_lift_pp,
  f.realized_gross_yen,
  f.audit_cost_yen,
  f.action_day_cost_yen,
  f.net_after_audit_yen,
  f.fully_loaded_net_yen,
  f.uncaptured_ignored_gross_yen,
  f.gross_capture_pct,
  d.gross_yen_per_acted_day,
  d.delay0_lost_gross,
  d.delay1_lost_gross,
  d.delay3_lost_gross,
  t.top_action,
  t.top_action_net_yen,
  t.top_action_net_yen_per_day,
  CONCAT(
    'Before→After 工单 ',
    CAST(ROUND(100.0 * r.ticket_rate_ignored, 2) AS VARCHAR),
    '%→',
    CAST(ROUND(100.0 * r.ticket_rate_acted, 2) AS VARCHAR),
    '%，承接 ',
    CAST(ROUND(100.0 * r.contain_rate_ignored, 2) AS VARCHAR),
    '%→',
    CAST(ROUND(100.0 * r.contain_rate_acted, 2) AS VARCHAR),
    '%（+',
    CAST(ROUND(100.0 * (r.contain_rate_acted - r.contain_rate_ignored), 2) AS VARCHAR),
    'pp）；全成本净¥',
    CAST(CAST(f.fully_loaded_net_yen AS BIGINT) AS VARCHAR),
    '（毛¥',
    CAST(CAST(f.realized_gross_yen AS BIGINT) AS VARCHAR),
    '−审计¥',
    CAST(CAST(f.audit_cost_yen AS BIGINT) AS VARCHAR),
    '−动作¥',
    CAST(CAST(f.action_day_cost_yen AS BIGINT) AS VARCHAR),
    '）；捕获 ',
    CAST(f.gross_capture_pct AS VARCHAR),
    '%；delay0/1/3 少拿毛¥',
    CAST(CAST(d.delay0_lost_gross AS BIGINT) AS VARCHAR),
    '/',
    CAST(CAST(d.delay1_lost_gross AS BIGINT) AS VARCHAR),
    '/',
    CAST(CAST(d.delay3_lost_gross AS BIGINT) AS VARCHAR),
    '；优先 ',
    CAST(t.top_action AS VARCHAR),
    '（日净≈¥',
    CAST(CAST(t.top_action_net_yen_per_day AS BIGINT) AS VARCHAR),
    '）。'
  ) AS external_one_liner_cn
FROM vw_cs_assist_fully_loaded_capture f
JOIN vw_cs_assist_rate_compare r ON r.surface_id = f.surface_id
JOIN delay d ON d.surface_id = f.surface_id
JOIN top_act t ON t.surface_id = f.surface_id;
