-- =============================================================================
-- 11) 周贡献三项拆分 + 火情覆盖率外推
-- 业务问题：
--   (1) 每周毛¥里，工单/退款/承接各多少？
--   (2) 覆盖率从 50%→75%→100% 时，多拿多少毛/净？
-- =============================================================================

-- 日级三项¥（在 acted_day_value 上拆）
CREATE OR REPLACE VIEW vw_cs_assist_acted_day_value_split AS
SELECT
  d.*,
  (d.ignored_ticket_rate - d.ticket_rate) * d.n_sessions * COALESCE(v.cs_ticket_cost, 0)
    AS yen_tickets_day,
  (d.ignored_refund_rate - d.refund_rate) * d.n_sessions * COALESCE(v.refund_unit_cost, 0)
    AS yen_refunds_day,
  (d.containment_rate - d.ignored_contain_rate) * d.n_sessions
    * COALESCE(v.contained_session_value, 0) AS yen_contain_day
FROM vw_cs_assist_acted_day_value d
LEFT JOIN (
  SELECT
    surface_id,
    MAX(CASE WHEN metric = 'cs_ticket_cost' THEN unit_value END) AS cs_ticket_cost,
    MAX(CASE WHEN metric = 'refund_unit_cost' THEN unit_value END) AS refund_unit_cost,
    MAX(CASE WHEN metric = 'contained_session_value' THEN unit_value END) AS contained_session_value
  FROM dim_value_assumption
  GROUP BY surface_id
) v ON v.surface_id = d.surface_id;

-- 周归因 + 三项¥
CREATE OR REPLACE VIEW vw_cs_assist_week_value_split AS
SELECT
  date_trunc('week', CAST(dt AS TIMESTAMP))::DATE AS week_start,
  surface_id,
  COUNT(*) AS days_acted,
  SUM(n_sessions) AS sessions_acted,
  ROUND(SUM(tickets_avoided_day), 1) AS tickets_avoided,
  ROUND(SUM(refunds_avoided_day), 1) AS refunds_avoided,
  ROUND(SUM(extra_contained_day), 1) AS extra_contained,
  ROUND(SUM(yen_tickets_day), 0) AS yen_from_tickets,
  ROUND(SUM(yen_refunds_day), 0) AS yen_from_refunds,
  ROUND(SUM(yen_contain_day), 0) AS yen_from_containment,
  ROUND(SUM(gross_yen_day), 0) AS gross_yen,
  ROUND(
    100.0 * SUM(yen_contain_day) / NULLIF(SUM(gross_yen_day), 0),
    1
  ) AS contain_share_pct,
  ROUND(
    SUM(gross_yen_day) * 100.0
      / NULLIF(SUM(SUM(gross_yen_day)) OVER (PARTITION BY surface_id), 0),
    1
  ) AS pct_of_total_gross
FROM vw_cs_assist_acted_day_value_split
GROUP BY 1, 2
ORDER BY 1, 2;

-- 周三项加总对账
CREATE OR REPLACE VIEW vw_cs_assist_week_split_check AS
SELECT
  e.surface_id,
  e.incremental_yen AS ledger_gross_yen,
  e.yen_from_tickets AS ledger_tickets_yen,
  e.yen_from_refunds AS ledger_refunds_yen,
  e.yen_from_containment AS ledger_contain_yen,
  ROUND(SUM(w.gross_yen), 0) AS week_sum_gross,
  ROUND(SUM(w.yen_from_tickets), 0) AS week_sum_tickets,
  ROUND(SUM(w.yen_from_refunds), 0) AS week_sum_refunds,
  ROUND(SUM(w.yen_from_containment), 0) AS week_sum_contain,
  ROUND(SUM(w.gross_yen) - e.incremental_yen, 0) AS gross_gap,
  ROUND(SUM(w.yen_from_containment) - e.yen_from_containment, 0) AS contain_yen_gap
FROM vw_cs_assist_exec_summary e
LEFT JOIN vw_cs_assist_week_value_split w ON w.surface_id = e.surface_id
GROUP BY
  e.surface_id,
  e.incremental_yen,
  e.yen_from_tickets,
  e.yen_from_refunds,
  e.yen_from_containment;

-- 火情覆盖率外推：当前已实现 + 按目标覆盖率吃掉留白
CREATE OR REPLACE VIEW vw_cs_assist_coverage_expansion AS
WITH base AS (
  SELECT
    o.surface_id,
    o.fire_day_coverage_pct AS coverage_now_pct,
    o.realized_gross_yen,
    o.opportunity_gross_yen,
    o.potential_full_coverage_gross_yen,
    n.net_incremental_yen AS realized_net_yen,
    n.audit_cost_acted,
    -- 审计成本按覆盖天数近似线性（当前 50% 覆盖）
    n.audit_cost_acted * 1.0 / NULLIF(o.fire_day_coverage_pct / 100.0, 0)
      AS audit_cost_at_full_coverage
  FROM vw_cs_assist_ignored_opportunity o
  JOIN vw_cs_assist_net_increment n ON n.surface_id = o.surface_id
),
scen AS (
  SELECT * FROM (VALUES
    (50.0, 'cover_50_current'),
    (75.0, 'cover_75'),
    (100.0, 'cover_100_full')
  ) AS t(target_coverage_pct, scenario)
)
SELECT
  b.surface_id,
  s.scenario,
  s.target_coverage_pct,
  b.coverage_now_pct,
  ROUND(
    b.realized_gross_yen
      + b.opportunity_gross_yen
        * GREATEST(s.target_coverage_pct - b.coverage_now_pct, 0)
        / NULLIF(100.0 - b.coverage_now_pct, 0),
    0
  ) AS projected_gross_yen,
  ROUND(
    (
      b.realized_gross_yen
        + b.opportunity_gross_yen
          * GREATEST(s.target_coverage_pct - b.coverage_now_pct, 0)
          / NULLIF(100.0 - b.coverage_now_pct, 0)
    )
    - b.audit_cost_at_full_coverage * (s.target_coverage_pct / 100.0),
    0
  ) AS projected_net_yen,
  ROUND(
    (
      b.realized_gross_yen
        + b.opportunity_gross_yen
          * GREATEST(s.target_coverage_pct - b.coverage_now_pct, 0)
          / NULLIF(100.0 - b.coverage_now_pct, 0)
    )
    - b.realized_gross_yen,
    0
  ) AS incremental_gross_vs_now,
  ROUND(b.potential_full_coverage_gross_yen, 0) AS full_coverage_gross_yen
FROM base b
CROSS JOIN scen s
ORDER BY s.target_coverage_pct;

-- 经营一页：before/after + 覆盖外推一句
CREATE OR REPLACE VIEW vw_cs_assist_coverage_decision AS
SELECT
  ba.surface_id,
  ba.ticket_rate_ignored_pct,
  ba.ticket_rate_acted_pct,
  ba.contain_rate_ignored_pct,
  ba.contain_rate_acted_pct,
  ba.contain_rate_lift_pp,
  ba.tickets_avoided,
  ba.refunds_avoided,
  ba.extra_sessions_contained,
  ba.gross_yen,
  ba.net_yen,
  o.fire_day_coverage_pct,
  o.opportunity_gross_yen,
  c75.projected_net_yen AS net_at_cover_75,
  c100.projected_net_yen AS net_at_cover_100,
  c100.incremental_gross_vs_now AS gross_uplift_if_full_coverage,
  d.top_action,
  CONCAT(
    'Before→After 工单 ',
    CAST(ba.ticket_rate_ignored_pct AS VARCHAR), '%→',
    CAST(ba.ticket_rate_acted_pct AS VARCHAR), '%，承接 ',
    CAST(ba.contain_rate_ignored_pct AS VARCHAR), '%→',
    CAST(ba.contain_rate_acted_pct AS VARCHAR), '%（+',
    CAST(ba.contain_rate_lift_pp AS VARCHAR), 'pp）；',
    '现净¥', CAST(CAST(ba.net_yen AS BIGINT) AS VARCHAR),
    '（覆盖 ', CAST(o.fire_day_coverage_pct AS VARCHAR), '%）；',
    '若覆盖拉到 75%/100%，净≈¥',
    CAST(CAST(c75.projected_net_yen AS BIGINT) AS VARCHAR), '/',
    CAST(CAST(c100.projected_net_yen AS BIGINT) AS VARCHAR),
    '；优先 ', CAST(d.top_action AS VARCHAR), '。'
  ) AS external_one_liner_cn
FROM vw_cs_assist_before_after ba
JOIN vw_cs_assist_ignored_opportunity o ON o.surface_id = ba.surface_id
JOIN vw_cs_assist_exec_dashboard d ON d.surface_id = ba.surface_id
JOIN vw_cs_assist_coverage_expansion c75
  ON c75.surface_id = ba.surface_id AND c75.scenario = 'cover_75'
JOIN vw_cs_assist_coverage_expansion c100
  ON c100.surface_id = ba.surface_id AND c100.scenario = 'cover_100_full';
