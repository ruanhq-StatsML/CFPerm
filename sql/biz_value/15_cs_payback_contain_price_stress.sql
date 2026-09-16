-- =============================================================================
-- 15) 承接单价承压下的回本（仅承接分子）
-- 业务问题：自助会话价值砍到 ¥1.75（-50%）时，各动作「仅靠承接」还能当天回本吗？
-- =============================================================================

CREATE OR REPLACE VIEW vw_cs_assist_payback_contain_price_stress AS
WITH base AS (
  SELECT
    p.surface_id,
    p.action_type,
    p.n_days,
    p.action_day_cost_yen,
    p.cost_yen_per_day,
    p.payback_days AS payback_days_full_gross,
    s.extra_contained,
    s.contain_unit_price AS contain_price_base,
    s.yen_from_containment AS contain_yen_base
  FROM vw_cs_assist_action_payback p
  JOIN vw_cs_assist_action_value_split s
    ON s.surface_id = p.surface_id AND s.action_type = p.action_type
),
scen AS (
  SELECT * FROM (VALUES
    ('contain_price_base', 1.00),
    ('contain_price_minus50', 0.50),
    ('contain_price_plus50', 1.50)
  ) AS t(scenario, price_scale)
)
SELECT
  b.surface_id,
  b.action_type,
  s.scenario,
  ROUND(b.contain_price_base * s.price_scale, 2) AS contain_unit_price,
  b.extra_contained,
  ROUND(b.contain_yen_base * s.price_scale, 0) AS contain_yen,
  ROUND(b.contain_yen_base * s.price_scale / NULLIF(b.n_days, 0), 0)
    AS contain_yen_per_day,
  b.cost_yen_per_day,
  b.payback_days_full_gross,
  ROUND(
    b.cost_yen_per_day
      / NULLIF(b.contain_yen_base * s.price_scale / NULLIF(b.n_days, 0), 0),
    3
  ) AS payback_days_contain_only,
  CASE
    WHEN b.contain_yen_base * s.price_scale / NULLIF(b.n_days, 0) <= 0
      THEN 'no_contain_payback'
    WHEN b.cost_yen_per_day
           / NULLIF(b.contain_yen_base * s.price_scale / NULLIF(b.n_days, 0), 0) <= 1.0
      THEN 'same_day_from_contain'
    WHEN b.cost_yen_per_day
           / NULLIF(b.contain_yen_base * s.price_scale / NULLIF(b.n_days, 0), 0) <= 3.0
      THEN 'within_3_days_from_contain'
    ELSE 'slow_from_contain'
  END AS contain_payback_bucket,
  ROUND(
    (b.contain_yen_base * s.price_scale) / NULLIF(b.action_day_cost_yen, 0),
    1
  ) AS contain_roi_vs_action_cost
FROM base b
CROSS JOIN scen s
ORDER BY b.action_type, s.price_scale;

-- 摘要：-50% 下是否仍全部 same_day
CREATE OR REPLACE VIEW vw_cs_assist_payback_contain_stress_summary AS
SELECT
  surface_id,
  COUNT(*) AS n_actions,
  SUM(CASE WHEN contain_payback_bucket = 'same_day_from_contain' THEN 1 ELSE 0 END)
    AS n_same_day_at_minus50,
  MIN(payback_days_contain_only) AS best_contain_payback_days,
  MAX(payback_days_contain_only) AS worst_contain_payback_days,
  CASE
    WHEN SUM(CASE WHEN contain_payback_bucket = 'same_day_from_contain' THEN 1 ELSE 0 END)
         = COUNT(*)
    THEN 'all_arms_same_day_at_contain_minus50'
    ELSE 'some_arms_need_ticket_refund_to_payback'
  END AS stress_status,
  CONCAT(
    '承接单价 -50% 时，仅靠多承接回本：最快 ',
    CAST(MIN(payback_days_contain_only) AS VARCHAR),
    ' 天 / 最慢 ',
    CAST(MAX(payback_days_contain_only) AS VARCHAR),
    ' 天；',
    CAST(SUM(CASE WHEN contain_payback_bucket = 'same_day_from_contain' THEN 1 ELSE 0 END) AS VARCHAR),
    '/',
    CAST(COUNT(*) AS VARCHAR),
    ' 臂仍当天回本。'
  ) AS external_one_liner_cn
FROM vw_cs_assist_payback_contain_price_stress
WHERE scenario = 'contain_price_minus50'
GROUP BY surface_id;
