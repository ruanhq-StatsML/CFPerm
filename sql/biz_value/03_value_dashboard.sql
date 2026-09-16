-- =============================================================================
-- 3) Unified value dashboard (what leadership / on-call reads)
-- =============================================================================

-- Executive one-pager: quality cost + style monetization + axis hygiene
CREATE OR REPLACE VIEW vw_biz_value_dashboard AS
SELECT
  'hallucination' AS theme,
  h.surface_id,
  h.n_fire_acted_days,
  h.n_fire_ignored_days,
  h.avg_cost_fire_acted,
  h.avg_cost_fire_ignored,
  h.avg_cost_quiet,
  h.est_daily_save_act_vs_ignore AS primary_value_proxy,
  CAST(NULL AS DOUBLE) AS style_only_days,
  CAST(NULL AS DOUBLE) AS avg_ctr_style_only
FROM vw_halluc_roi_rollup h
UNION ALL
SELECT
  'style_creative' AS theme,
  s.surface_id,
  CAST(NULL AS BIGINT) AS n_fire_acted_days,
  CAST(NULL AS BIGINT) AS n_fire_ignored_days,
  CAST(NULL AS DOUBLE) AS avg_cost_fire_acted,
  CAST(NULL AS DOUBLE) AS avg_cost_fire_ignored,
  CAST(NULL AS DOUBLE) AS avg_cost_quiet,
  -- value proxy: brand-cost elevation on style-only days vs quiet (descriptive)
  (s.avg_brand_cost_style_only - s.avg_brand_cost_quiet) AS primary_value_proxy,
  s.n_style_only_days AS style_only_days,
  s.avg_ctr_style_only
FROM vw_style_roi_rollup s;

-- Merge-gate hygiene: refuse DPO/RLHF merge when judge_err_ratio high OR concept fire
CREATE OR REPLACE VIEW vw_alignment_merge_gate AS
SELECT
  s.dt,
  s.surface_id,
  s.model_version,
  MAX(CASE WHEN s.signal_type = 'judge_err_ratio' THEN s.score END) AS judge_err_ratio,
  MAX(CASE WHEN s.signal_type = 'style_domain_auc' THEN s.score END) AS style_domain_auc,
  MAX(CASE WHEN s.signal_type = 'rfperm_fire' AND s.axis = 'concept' THEN s.fired ELSE 0 END) AS concept_fired,
  CASE
    WHEN MAX(CASE WHEN s.signal_type = 'rfperm_fire' AND s.axis = 'concept' THEN s.fired ELSE 0 END) = 1
      THEN 'BLOCK_MERGE — concept/preference hop'
    WHEN MAX(CASE WHEN s.signal_type = 'judge_err_ratio' THEN s.score END) >= 1.5
      THEN 'BLOCK_MERGE — RFPerm-as-Judge degraded'
    WHEN MAX(CASE WHEN s.signal_type = 'style_domain_auc' THEN s.score END) >= 0.70
         AND MAX(CASE WHEN s.signal_type = 'rfperm_fire' AND s.axis = 'concept' THEN s.fired ELSE 0 END) = 0
      THEN 'ALLOW_WITH_STYLE_REBALANCE — portrait only'
    ELSE 'ALLOW'
  END AS merge_decision
FROM fct_shift_signal s
GROUP BY 1, 2, 3;

-- On-call combined queue
CREATE OR REPLACE VIEW vw_oncall_queue AS
SELECT
  'halluc' AS source,
  dt,
  surface_id,
  biz_bucket AS bucket,
  recommended_action,
  total_quality_cost AS severity_proxy
FROM vw_halluc_workorders
WHERE regime_fired = 1
UNION ALL
SELECT
  'style' AS source,
  dt,
  surface_id,
  biz_axis_bucket AS bucket,
  recommended_action,
  brand_flag_rate AS severity_proxy
FROM vw_style_workorders
WHERE biz_axis_bucket IN ('style_only_portrait', 'concept_only_quality', 'joint_style_and_concept');
