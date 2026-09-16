-- =============================================================================
-- 16) Early detection → detection delay → 利润差
-- 业务问题：火情晚发现 delay 天，相对零延迟少拿多少毛/净？
-- 口径：delay × 每动作日贡献（与「火了仍 ignored」留白同构）
-- =============================================================================

CREATE OR REPLACE VIEW vw_cs_assist_detection_delay_profit AS
WITH base AS (
  SELECT
    m.surface_id,
    m.days_acted,
    m.days_ignored,
    m.sessions_per_acted_day,
    m.gross_yen_per_acted_day,
    m.net_yen_per_acted_day,
    m.contain_yen_per_acted_day,
    m.remaining_opportunity_gross_yen,
    e.incremental_yen AS realized_gross_yen,
    n.net_incremental_yen AS realized_net_yen,
    r.ticket_rate_ignored,
    r.ticket_rate_acted,
    r.refund_rate_ignored,
    r.refund_rate_acted,
    r.contain_rate_ignored,
    r.contain_rate_acted
  FROM vw_cs_assist_marginal_day m
  JOIN vw_cs_assist_exec_summary e ON e.surface_id = m.surface_id
  JOIN vw_cs_assist_net_increment n ON n.surface_id = m.surface_id
  JOIN vw_cs_assist_rate_compare r ON r.surface_id = m.surface_id
),
delays AS (
  SELECT * FROM (VALUES
    (0), (1), (2), (3)
  ) AS t(delay_days)
)
SELECT
  b.surface_id,
  d.delay_days,
  b.gross_yen_per_acted_day,
  b.net_yen_per_acted_day,
  ROUND(d.delay_days * b.gross_yen_per_acted_day, 0) AS lost_gross_yen,
  ROUND(d.delay_days * b.net_yen_per_acted_day, 0) AS lost_net_yen,
  ROUND(b.realized_gross_yen - d.delay_days * b.gross_yen_per_acted_day, 0)
    AS gross_if_delayed_vs_zero,
  ROUND(b.realized_net_yen - d.delay_days * b.net_yen_per_acted_day, 0)
    AS net_if_delayed_vs_zero,
  -- before/after 费率（火情窗）
  ROUND(100.0 * b.ticket_rate_ignored, 2) AS ticket_rate_ignored_pct,
  ROUND(100.0 * b.ticket_rate_acted, 2) AS ticket_rate_acted_pct,
  ROUND(100.0 * b.contain_rate_ignored, 2) AS contain_rate_ignored_pct,
  ROUND(100.0 * b.contain_rate_acted, 2) AS contain_rate_acted_pct,
  CASE
    WHEN d.delay_days = 0 THEN 'zero_delay_full_capture'
    WHEN d.delay_days = 1 THEN 'one_day_late'
    ELSE 'multi_day_late'
  END AS delay_bucket
FROM base b
CROSS JOIN delays d
ORDER BY d.delay_days;

CREATE OR REPLACE VIEW vw_cs_assist_detection_delay_summary AS
SELECT
  surface_id,
  MAX(CASE WHEN delay_days = 0 THEN lost_gross_yen END) AS lost_gross_at_delay0,
  MAX(CASE WHEN delay_days = 1 THEN lost_gross_yen END) AS lost_gross_at_delay1,
  MAX(CASE WHEN delay_days = 2 THEN lost_gross_yen END) AS lost_gross_at_delay2,
  MAX(CASE WHEN delay_days = 3 THEN lost_gross_yen END) AS lost_gross_at_delay3,
  MAX(CASE WHEN delay_days = 1 THEN lost_net_yen END) AS lost_net_at_delay1,
  MAX(gross_yen_per_acted_day) AS gross_yen_per_acted_day,
  MAX(net_yen_per_acted_day) AS net_yen_per_acted_day,
  MAX(ticket_rate_ignored_pct) AS ticket_rate_ignored_pct,
  MAX(ticket_rate_acted_pct) AS ticket_rate_acted_pct,
  MAX(contain_rate_ignored_pct) AS contain_rate_ignored_pct,
  MAX(contain_rate_acted_pct) AS contain_rate_acted_pct,
  CONCAT(
    'Early detection：delay 0→1→3 天少拿毛¥',
    CAST(CAST(MAX(CASE WHEN delay_days = 0 THEN lost_gross_yen END) AS BIGINT) AS VARCHAR),
    '/',
    CAST(CAST(MAX(CASE WHEN delay_days = 1 THEN lost_gross_yen END) AS BIGINT) AS VARCHAR),
    '/',
    CAST(CAST(MAX(CASE WHEN delay_days = 3 THEN lost_gross_yen END) AS BIGINT) AS VARCHAR),
    '（每动作日≈¥',
    CAST(CAST(MAX(gross_yen_per_acted_day) AS BIGINT) AS VARCHAR),
    '）；Before→After 工单 ',
    CAST(MAX(ticket_rate_ignored_pct) AS VARCHAR),
    '%→',
    CAST(MAX(ticket_rate_acted_pct) AS VARCHAR),
    '%，承接 ',
    CAST(MAX(contain_rate_ignored_pct) AS VARCHAR),
    '%→',
    CAST(MAX(contain_rate_acted_pct) AS VARCHAR),
    '%。'
  ) AS external_one_liner_cn
FROM vw_cs_assist_detection_delay_profit
GROUP BY surface_id;
