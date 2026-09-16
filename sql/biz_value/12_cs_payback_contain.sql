-- =============================================================================
-- 12) 动作回本分子拆分（承接贡献进回本）
-- 业务问题：回本天数的分子毛¥里，有多少来自多承接？若只靠承接，几天回本？
-- =============================================================================

CREATE OR REPLACE VIEW vw_cs_assist_action_payback_contain AS
SELECT
  p.surface_id,
  p.action_type,
  p.n_days,
  p.gross_yen,
  s.yen_from_tickets,
  s.yen_from_refunds,
  s.yen_from_containment,
  s.extra_contained,
  p.action_day_cost_yen,
  p.net_yen,
  p.gross_yen_per_day,
  p.cost_yen_per_day,
  p.payback_days,
  p.payback_bucket,
  p.net_roi_multiple,
  ROUND(s.yen_from_containment * 1.0 / NULLIF(p.n_days, 0), 0) AS contain_yen_per_day,
  ROUND(
    100.0 * s.yen_from_containment / NULLIF(p.gross_yen, 0),
    1
  ) AS contain_share_pct_of_gross,
  -- 若分子只用承接¥：成本 / 日均承接¥
  ROUND(
    (p.action_day_cost_yen * 1.0 / NULLIF(p.n_days, 0))
      / NULLIF(s.yen_from_containment * 1.0 / NULLIF(p.n_days, 0), 0),
    3
  ) AS payback_days_contain_only,
  CASE
    WHEN s.yen_from_containment * 1.0 / NULLIF(p.n_days, 0) <= 0 THEN 'no_contain_payback'
    WHEN (p.action_day_cost_yen * 1.0 / NULLIF(p.n_days, 0))
           / NULLIF(s.yen_from_containment * 1.0 / NULLIF(p.n_days, 0), 0) <= 1.0
      THEN 'same_day_from_contain'
    WHEN (p.action_day_cost_yen * 1.0 / NULLIF(p.n_days, 0))
           / NULLIF(s.yen_from_containment * 1.0 / NULLIF(p.n_days, 0), 0) <= 3.0
      THEN 'within_3_days_from_contain'
    ELSE 'slow_from_contain'
  END AS contain_payback_bucket,
  ROUND(
    s.yen_from_containment * 1.0 / NULLIF(p.action_day_cost_yen, 0),
    1
  ) AS contain_roi_vs_action_cost
FROM vw_cs_assist_action_payback p
JOIN vw_cs_assist_action_value_split s
  ON s.surface_id = p.surface_id AND s.action_type = p.action_type
ORDER BY p.payback_days, s.yen_from_containment DESC;
