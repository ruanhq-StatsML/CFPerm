-- =============================================================================
-- 08) 累计贡献曲线 + 单位经济敏感度（财务对账）
-- 不改流量，只改单价 → 看贡献¥带；按日累计毛贡献给财务画曲线。
-- =============================================================================

-- 日累计毛贡献（acted 日）
CREATE OR REPLACE VIEW vw_cs_assist_cumulative_curve AS
SELECT
  dt,
  surface_id,
  n_sessions,
  ROUND(tickets_avoided_day, 2) AS tickets_avoided_day,
  ROUND(refunds_avoided_day, 2) AS refunds_avoided_day,
  ROUND(extra_contained_day, 2) AS extra_contained_day,
  ROUND(gross_yen_day, 0) AS gross_yen_day,
  ROUND(SUM(gross_yen_day) OVER (
    PARTITION BY surface_id ORDER BY dt
    ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW
  ), 0) AS cumulative_gross_yen,
  ROUND(
    100.0 * SUM(gross_yen_day) OVER (
      PARTITION BY surface_id ORDER BY dt
      ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW
    ) / NULLIF(SUM(gross_yen_day) OVER (PARTITION BY surface_id), 0),
    1
  ) AS cumulative_pct_of_total
FROM vw_cs_assist_acted_day_value
ORDER BY surface_id, dt;

-- 单位经济敏感度：固定量（少工单/退款/多承接），扫单价
CREATE OR REPLACE VIEW vw_cs_assist_unit_econ_sensitivity AS
WITH base AS (
  SELECT
    e.surface_id,
    e.tickets_avoided,
    e.refunds_avoided,
    e.extra_sessions_contained AS extra_contained,
    v.cs_ticket_cost AS base_ticket_cost,
    v.refund_unit_cost AS base_refund_cost,
    v.contained_session_value AS base_contain_value,
    e.incremental_yen AS base_gross_yen,
    n.net_incremental_yen AS base_net_yen,
    n.audit_cost_acted
  FROM vw_cs_assist_exec_summary e
  JOIN vw_cs_assist_net_increment n ON n.surface_id = e.surface_id
  JOIN (
    SELECT
      surface_id,
      MAX(CASE WHEN metric = 'cs_ticket_cost' THEN unit_value END) AS cs_ticket_cost,
      MAX(CASE WHEN metric = 'refund_unit_cost' THEN unit_value END) AS refund_unit_cost,
      MAX(CASE WHEN metric = 'contained_session_value' THEN unit_value END) AS contained_session_value
    FROM dim_value_assumption
    GROUP BY surface_id
  ) v ON v.surface_id = e.surface_id
),
scen AS (
  SELECT * FROM (VALUES
    ('base', 1.00),
    ('ticket_plus20', 1.00),
    ('ticket_minus20', 1.00),
    ('refund_plus20', 1.00),
    ('refund_minus20', 1.00),
    ('all_plus20', 1.20),
    ('all_minus20', 0.80)
  ) AS t(scenario, dummy)
)
SELECT
  b.surface_id,
  s.scenario,
  ROUND(CASE
    WHEN s.scenario = 'ticket_plus20' THEN b.base_ticket_cost * 1.20
    WHEN s.scenario = 'ticket_minus20' THEN b.base_ticket_cost * 0.80
    WHEN s.scenario = 'all_plus20' THEN b.base_ticket_cost * 1.20
    WHEN s.scenario = 'all_minus20' THEN b.base_ticket_cost * 0.80
    ELSE b.base_ticket_cost
  END, 2) AS ticket_cost,
  ROUND(CASE
    WHEN s.scenario = 'refund_plus20' THEN b.base_refund_cost * 1.20
    WHEN s.scenario = 'refund_minus20' THEN b.base_refund_cost * 0.80
    WHEN s.scenario = 'all_plus20' THEN b.base_refund_cost * 1.20
    WHEN s.scenario = 'all_minus20' THEN b.base_refund_cost * 0.80
    ELSE b.base_refund_cost
  END, 2) AS refund_cost,
  ROUND(CASE
    WHEN s.scenario = 'all_plus20' THEN b.base_contain_value * 1.20
    WHEN s.scenario = 'all_minus20' THEN b.base_contain_value * 0.80
    ELSE b.base_contain_value
  END, 2) AS contain_value,
  b.tickets_avoided,
  b.refunds_avoided,
  b.extra_contained,
  ROUND(
    b.tickets_avoided * CASE
      WHEN s.scenario = 'ticket_plus20' THEN b.base_ticket_cost * 1.20
      WHEN s.scenario = 'ticket_minus20' THEN b.base_ticket_cost * 0.80
      WHEN s.scenario = 'all_plus20' THEN b.base_ticket_cost * 1.20
      WHEN s.scenario = 'all_minus20' THEN b.base_ticket_cost * 0.80
      ELSE b.base_ticket_cost
    END
    + b.refunds_avoided * CASE
      WHEN s.scenario = 'refund_plus20' THEN b.base_refund_cost * 1.20
      WHEN s.scenario = 'refund_minus20' THEN b.base_refund_cost * 0.80
      WHEN s.scenario = 'all_plus20' THEN b.base_refund_cost * 1.20
      WHEN s.scenario = 'all_minus20' THEN b.base_refund_cost * 0.80
      ELSE b.base_refund_cost
    END
    + b.extra_contained * CASE
      WHEN s.scenario = 'all_plus20' THEN b.base_contain_value * 1.20
      WHEN s.scenario = 'all_minus20' THEN b.base_contain_value * 0.80
      ELSE b.base_contain_value
    END
  , 0) AS gross_yen,
  ROUND(
    b.tickets_avoided * CASE
      WHEN s.scenario = 'ticket_plus20' THEN b.base_ticket_cost * 1.20
      WHEN s.scenario = 'ticket_minus20' THEN b.base_ticket_cost * 0.80
      WHEN s.scenario = 'all_plus20' THEN b.base_ticket_cost * 1.20
      WHEN s.scenario = 'all_minus20' THEN b.base_ticket_cost * 0.80
      ELSE b.base_ticket_cost
    END
    + b.refunds_avoided * CASE
      WHEN s.scenario = 'refund_plus20' THEN b.base_refund_cost * 1.20
      WHEN s.scenario = 'refund_minus20' THEN b.base_refund_cost * 0.80
      WHEN s.scenario = 'all_plus20' THEN b.base_refund_cost * 1.20
      WHEN s.scenario = 'all_minus20' THEN b.base_refund_cost * 0.80
      ELSE b.base_refund_cost
    END
    + b.extra_contained * CASE
      WHEN s.scenario = 'all_plus20' THEN b.base_contain_value * 1.20
      WHEN s.scenario = 'all_minus20' THEN b.base_contain_value * 0.80
      ELSE b.base_contain_value
    END
    - b.audit_cost_acted
  , 0) AS net_yen,
  b.base_gross_yen,
  b.base_net_yen
FROM base b
CROSS JOIN scen s;

-- 敏感度摘要：base / 全+20 / 全-20 三档
CREATE OR REPLACE VIEW vw_cs_assist_unit_econ_band AS
SELECT
  surface_id,
  MAX(CASE WHEN scenario = 'base' THEN net_yen END) AS net_yen_base,
  MAX(CASE WHEN scenario = 'all_minus20' THEN net_yen END) AS net_yen_low,
  MAX(CASE WHEN scenario = 'all_plus20' THEN net_yen END) AS net_yen_high,
  MAX(CASE WHEN scenario = 'base' THEN gross_yen END) AS gross_yen_base,
  MAX(CASE WHEN scenario = 'all_minus20' THEN gross_yen END) AS gross_yen_low,
  MAX(CASE WHEN scenario = 'all_plus20' THEN gross_yen END) AS gross_yen_high,
  ROUND(
    (MAX(CASE WHEN scenario = 'all_plus20' THEN net_yen END)
      - MAX(CASE WHEN scenario = 'all_minus20' THEN net_yen END))
    / NULLIF(MAX(CASE WHEN scenario = 'base' THEN net_yen END), 0),
    2
  ) AS net_band_width_vs_base
FROM vw_cs_assist_unit_econ_sensitivity
GROUP BY surface_id;
