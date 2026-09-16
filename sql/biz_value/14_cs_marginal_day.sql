-- =============================================================================
-- 14) 边际动作日贡献（多动作一天值多少）
-- 业务问题：当前每个 acted 日平均贡献多少¥？再把 1 个 ignored 日改成动作，期望多拿多少？
-- =============================================================================

CREATE OR REPLACE VIEW vw_cs_assist_marginal_day AS
SELECT
  e.surface_id,
  e.days_acted,
  o.days_ignored,
  e.sessions_acted,
  ROUND(e.sessions_acted * 1.0 / NULLIF(e.days_acted, 0), 0) AS sessions_per_acted_day,
  e.incremental_yen AS realized_gross_yen,
  n.net_incremental_yen AS realized_net_yen,
  e.yen_from_containment AS realized_contain_yen,
  ROUND(e.incremental_yen * 1.0 / NULLIF(e.days_acted, 0), 0) AS gross_yen_per_acted_day,
  ROUND(n.net_incremental_yen * 1.0 / NULLIF(e.days_acted, 0), 0) AS net_yen_per_acted_day,
  ROUND(e.yen_from_containment * 1.0 / NULLIF(e.days_acted, 0), 0) AS contain_yen_per_acted_day,
  ROUND(e.tickets_avoided * 1.0 / NULLIF(e.days_acted, 0), 1) AS tickets_per_acted_day,
  ROUND(e.refunds_avoided * 1.0 / NULLIF(e.days_acted, 0), 1) AS refunds_per_acted_day,
  ROUND(e.extra_sessions_contained * 1.0 / NULLIF(e.days_acted, 0), 1) AS contain_per_acted_day,
  -- 把 1 个 ignored 日改成动作：按留白毛 / ignored 天数
  ROUND(o.opportunity_gross_yen * 1.0 / NULLIF(o.days_ignored, 0), 0)
    AS expected_gross_if_one_more_acted_day,
  ROUND(
    (o.opportunity_gross_yen * 1.0 / NULLIF(o.days_ignored, 0))
      - (n.audit_cost_acted * 1.0 / NULLIF(e.days_acted, 0)),
    0
  ) AS expected_net_if_one_more_acted_day,
  o.days_ignored AS remaining_ignored_days,
  ROUND(o.opportunity_gross_yen, 0) AS remaining_opportunity_gross_yen,
  CONCAT(
    '每动作日毛≈¥',
    CAST(CAST(ROUND(e.incremental_yen * 1.0 / NULLIF(e.days_acted, 0), 0) AS BIGINT) AS VARCHAR),
    ' / 净≈¥',
    CAST(CAST(ROUND(n.net_incremental_yen * 1.0 / NULLIF(e.days_acted, 0), 0) AS BIGINT) AS VARCHAR),
    '（含承接≈¥',
    CAST(CAST(ROUND(e.yen_from_containment * 1.0 / NULLIF(e.days_acted, 0), 0) AS BIGINT) AS VARCHAR),
    '）；再动作 1 天期望毛≈¥',
    CAST(CAST(ROUND(o.opportunity_gross_yen * 1.0 / NULLIF(o.days_ignored, 0), 0) AS BIGINT) AS VARCHAR),
    '，仍余 ',
    CAST(o.days_ignored AS VARCHAR),
    ' 个 ignored 日留白毛¥',
    CAST(CAST(ROUND(o.opportunity_gross_yen, 0) AS BIGINT) AS VARCHAR),
    '。'
  ) AS external_one_liner_cn
FROM vw_cs_assist_exec_summary e
JOIN vw_cs_assist_net_increment n ON n.surface_id = e.surface_id
JOIN vw_cs_assist_ignored_opportunity o ON o.surface_id = e.surface_id;
