-- =============================================================================
-- 09) 不动作留白（机会成本）+ 成本/单价盈亏平衡
-- 业务问题：
--   (1) fire 但没动作的日子，若按 acted 费率，桌上留了多少工单/退款/承接/¥？
--   (2) 审计单价高到多少、单价缩到多少，净贡献才会归零？
-- =============================================================================

-- acted 臂费率（用作 ignored 日的反事实「若当时动作」）
CREATE OR REPLACE VIEW vw_cs_assist_acted_baseline AS
SELECT
  surface_id,
  ticket_per_session AS acted_ticket_rate,
  refund_per_session AS acted_refund_rate,
  contained_per_session AS acted_contain_rate
FROM vw_cs_assist_arm_stats
WHERE arm = 'fire_acted';

-- 每个 ignored 日的留白（相对 acted 费率）
CREATE OR REPLACE VIEW vw_cs_assist_ignored_day_opportunity AS
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
  a.acted_ticket_rate,
  a.acted_refund_rate,
  a.acted_contain_rate,
  (d.ticket_rate - a.acted_ticket_rate) * d.n_sessions AS tickets_left_on_table,
  (d.refund_rate - a.acted_refund_rate) * d.n_sessions AS refunds_left_on_table,
  (a.acted_contain_rate - d.containment_rate) * d.n_sessions AS contain_left_on_table,
  (d.ticket_rate - a.acted_ticket_rate) * d.n_sessions * COALESCE(v.cs_ticket_cost, 0)
    + (d.refund_rate - a.acted_refund_rate) * d.n_sessions * COALESCE(v.refund_unit_cost, 0)
    + (a.acted_contain_rate - d.containment_rate) * d.n_sessions * COALESCE(v.contained_session_value, 0)
    AS opportunity_gross_yen_day
FROM vw_cs_assist_arms d
JOIN vw_cs_assist_acted_baseline a ON a.surface_id = d.surface_id
LEFT JOIN (
  SELECT
    surface_id,
    MAX(CASE WHEN metric = 'cs_ticket_cost' THEN unit_value END) AS cs_ticket_cost,
    MAX(CASE WHEN metric = 'refund_unit_cost' THEN unit_value END) AS refund_unit_cost,
    MAX(CASE WHEN metric = 'contained_session_value' THEN unit_value END) AS contained_session_value
  FROM dim_value_assumption
  GROUP BY surface_id
) v ON v.surface_id = d.surface_id
WHERE d.arm = 'fire_ignored';

-- 留白汇总 + 覆盖率（已实现 vs 潜在全覆盖）
CREATE OR REPLACE VIEW vw_cs_assist_ignored_opportunity AS
WITH opp AS (
  SELECT
    surface_id,
    COUNT(*) AS days_ignored,
    SUM(n_sessions) AS sessions_ignored,
    ROUND(SUM(tickets_left_on_table), 1) AS tickets_left_on_table,
    ROUND(SUM(refunds_left_on_table), 1) AS refunds_left_on_table,
    ROUND(SUM(contain_left_on_table), 1) AS contain_left_on_table,
    ROUND(SUM(opportunity_gross_yen_day), 0) AS opportunity_gross_yen
  FROM vw_cs_assist_ignored_day_opportunity
  GROUP BY surface_id
),
real AS (
  SELECT
    surface_id,
    days_acted,
    days_ignored AS days_ignored_ledger,
    sessions_acted,
    incremental_yen AS realized_gross_yen,
    tickets_avoided,
    refunds_avoided,
    extra_sessions_contained
  FROM vw_cs_assist_exec_summary
)
SELECT
  r.surface_id,
  r.days_acted,
  o.days_ignored,
  r.sessions_acted,
  o.sessions_ignored,
  ROUND(
    100.0 * r.days_acted / NULLIF(r.days_acted + o.days_ignored, 0),
    1
  ) AS fire_day_coverage_pct,
  r.tickets_avoided AS tickets_realized,
  o.tickets_left_on_table,
  r.refunds_avoided AS refunds_realized,
  o.refunds_left_on_table,
  r.extra_sessions_contained AS contain_realized,
  o.contain_left_on_table,
  r.realized_gross_yen,
  o.opportunity_gross_yen,
  ROUND(r.realized_gross_yen + o.opportunity_gross_yen, 0) AS potential_full_coverage_gross_yen,
  ROUND(
    100.0 * r.realized_gross_yen
      / NULLIF(r.realized_gross_yen + o.opportunity_gross_yen, 0),
    1
  ) AS yen_capture_pct
FROM real r
JOIN opp o ON o.surface_id = r.surface_id;

-- 盈亏平衡：审计单价上限；单价整体缩放到净=0 的倍数；动作日成本上限（相对毛）
CREATE OR REPLACE VIEW vw_cs_assist_cost_breakeven AS
SELECT
  e.surface_id,
  e.incremental_yen AS gross_yen,
  n.net_incremental_yen AS net_yen,
  n.audits_acted,
  n.audit_cost_acted,
  v.audit_unit_cost AS audit_unit_cost_now,
  v.cs_ticket_cost,
  v.refund_unit_cost,
  v.contained_session_value,
  -- 净=0 时审计单价上限（件）
  ROUND(e.incremental_yen / NULLIF(n.audits_acted, 0), 1) AS max_audit_unit_cost_at_net0,
  ROUND(
    (e.incremental_yen / NULLIF(n.audits_acted, 0)) / NULLIF(v.audit_unit_cost, 0),
    1
  ) AS audit_unit_headroom_multiple,
  -- 毛贡献可承受的总审计成本上限 = 毛¥（净=0）
  ROUND(e.incremental_yen, 0) AS max_audit_spend_at_net0,
  -- 单价整体缩放 s：s*毛 - 审计 = 0 → s = 审计/毛（再低净为负）
  ROUND(n.audit_cost_acted / NULLIF(e.incremental_yen, 0), 4) AS min_price_scale_at_net0,
  ROUND(1.0 - n.audit_cost_acted / NULLIF(e.incremental_yen, 0), 4) AS price_cut_headroom_frac,
  -- 解读：现价还可砍多少比例仍净>0
  ROUND(100.0 * (1.0 - n.audit_cost_acted / NULLIF(e.incremental_yen, 0)), 1)
    AS price_cut_headroom_pct
FROM vw_cs_assist_exec_summary e
JOIN vw_cs_assist_net_increment n ON n.surface_id = e.surface_id
LEFT JOIN (
  SELECT
    surface_id,
    MAX(CASE WHEN metric = 'audit_unit_cost' THEN unit_value END) AS audit_unit_cost,
    MAX(CASE WHEN metric = 'cs_ticket_cost' THEN unit_value END) AS cs_ticket_cost,
    MAX(CASE WHEN metric = 'refund_unit_cost' THEN unit_value END) AS refund_unit_cost,
    MAX(CASE WHEN metric = 'contained_session_value' THEN unit_value END) AS contained_session_value
  FROM dim_value_assumption
  GROUP BY surface_id
) v ON v.surface_id = e.surface_id;

-- 对外一页：before/after + 留白 + 盈亏平衡（经营粘贴）
CREATE OR REPLACE VIEW vw_cs_assist_business_oneliner AS
SELECT
  e.surface_id,
  ba.ticket_rate_quiet_pct,
  ba.ticket_rate_ignored_pct,
  ba.ticket_rate_acted_pct,
  ba.contain_rate_quiet_pct,
  ba.contain_rate_ignored_pct,
  ba.contain_rate_acted_pct,
  ba.contain_rate_lift_pp,
  e.tickets_avoided,
  e.refunds_avoided,
  e.extra_sessions_contained,
  e.incremental_yen AS gross_yen,
  n.net_incremental_yen AS net_yen,
  o.opportunity_gross_yen,
  o.yen_capture_pct,
  o.fire_day_coverage_pct,
  o.potential_full_coverage_gross_yen,
  b.max_audit_unit_cost_at_net0,
  b.price_cut_headroom_pct,
  d.top_action,
  CONCAT(
    'Before→After：工单率 ',
    CAST(ba.ticket_rate_ignored_pct AS VARCHAR), '%→',
    CAST(ba.ticket_rate_acted_pct AS VARCHAR), '%，承接率 ',
    CAST(ba.contain_rate_ignored_pct AS VARCHAR), '%→',
    CAST(ba.contain_rate_acted_pct AS VARCHAR), '%（+',
    CAST(ba.contain_rate_lift_pp AS VARCHAR), 'pp）；',
    '少工单 ', CAST(e.tickets_avoided AS VARCHAR),
    ' / 少退款 ', CAST(e.refunds_avoided AS VARCHAR),
    ' / 多承接 ', CAST(e.extra_sessions_contained AS VARCHAR),
    '；净¥', CAST(CAST(n.net_incremental_yen AS BIGINT) AS VARCHAR),
    '；不动作留白毛¥', CAST(CAST(o.opportunity_gross_yen AS BIGINT) AS VARCHAR),
    '（火情日覆盖 ', CAST(o.fire_day_coverage_pct AS VARCHAR), '%）；',
    '优先 ', CAST(d.top_action AS VARCHAR), '。'
  ) AS external_one_liner_cn
FROM vw_cs_assist_exec_summary e
JOIN vw_cs_assist_before_after ba ON ba.surface_id = e.surface_id
JOIN vw_cs_assist_exec_dashboard d ON d.surface_id = e.surface_id
JOIN vw_cs_assist_net_increment n ON n.surface_id = e.surface_id
JOIN vw_cs_assist_ignored_opportunity o ON o.surface_id = e.surface_id
JOIN vw_cs_assist_cost_breakeven b ON b.surface_id = e.surface_id;
