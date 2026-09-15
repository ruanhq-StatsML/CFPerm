-- =============================================================================
-- 2) Style / Creative (文风·画风) DRIFT → CTR / brand / fatigue value
-- Axis: covariate P(X). Do NOT treat style AUC as "model quality failed".
-- =============================================================================

CREATE OR REPLACE VIEW vw_style_daily AS
SELECT
  e.dt,
  e.surface_id,
  COALESCE(c.style_cluster, 'unknown') AS style_cluster,
  e.model_version,
  COUNT(*) AS n_impr,
  AVG(e.clicked) AS ctr,
  AVG(e.converted) AS cvr,
  SUM(e.gmv) AS gmv,
  AVG(e.brand_flagged) AS brand_flag_rate,
  AVG(e.style_formal) AS mean_formal,
  AVG(e.style_len) AS mean_len,
  AVG(e.style_emoji) AS mean_emoji
FROM fct_serve_event e
LEFT JOIN dim_creative c ON c.creative_id = e.creative_id
WHERE e.surface_id IN ('feed_caption', 'ad_creative', 'shop_image_gen', 'video_cover')
   OR e.modality IN ('image', 'video', 'multi')
GROUP BY 1, 2, 3, 4;

CREATE OR REPLACE VIEW vw_style_shift_day AS
SELECT
  d.dt,
  d.surface_id,
  d.model_version,
  SUM(d.n_impr) AS n_impr,
  SUM(d.n_impr * d.ctr) / NULLIF(SUM(d.n_impr), 0) AS ctr,
  SUM(d.gmv) AS gmv,
  SUM(d.n_impr * d.brand_flag_rate) / NULLIF(SUM(d.n_impr), 0) AS brand_flag_rate,
  MAX(CASE WHEN s.signal_type = 'style_domain_auc' THEN s.score END) AS style_domain_auc,
  MAX(CASE WHEN s.signal_type = 'rfperm_fire' AND s.axis = 'covariate' THEN s.fired ELSE 0 END) AS covar_fired,
  MAX(CASE WHEN s.signal_type = 'rfperm_fire' AND s.axis = 'concept' THEN s.fired ELSE 0 END) AS concept_fired,
  MAX(CASE WHEN s.notes LIKE 'acted%' THEN 1 ELSE 0 END) AS acted
FROM vw_style_daily d
LEFT JOIN fct_shift_signal s
  ON s.dt = d.dt AND s.surface_id = d.surface_id
 AND (s.model_version = d.model_version OR s.model_version IS NULL)
GROUP BY 1, 2, 3;

-- Split portrait vs preference/quality — prevents wrong retrain
CREATE OR REPLACE VIEW vw_style_vs_concept_split AS
SELECT
  s.*,
  CASE
    WHEN COALESCE(s.style_domain_auc, 0) >= 0.70 AND COALESCE(s.concept_fired, 0) = 0
      THEN 'style_only_portrait'
    WHEN COALESCE(s.concept_fired, 0) = 1 AND COALESCE(s.style_domain_auc, 0) < 0.60
      THEN 'concept_only_quality'
    WHEN COALESCE(s.style_domain_auc, 0) >= 0.70 AND COALESCE(s.concept_fired, 0) = 1
      THEN 'joint_style_and_concept'
    ELSE 'quiet_or_weak'
  END AS biz_axis_bucket,
  CASE
    WHEN COALESCE(s.style_domain_auc, 0) >= 0.70 AND COALESCE(s.concept_fired, 0) = 0
      THEN 'Rebalance creative mix / decoding; do NOT retrain preference head'
    WHEN COALESCE(s.concept_fired, 0) = 1 AND COALESCE(s.style_domain_auc, 0) < 0.60
      THEN 'Preference/quality hop: RFPerm-as-Judge gate + √PO; keep style sampler'
    WHEN COALESCE(s.style_domain_auc, 0) >= 0.70 AND COALESCE(s.concept_fired, 0) = 1
      THEN 'Two tickets: (1) style mix (2) preference map — do not conflate'
    ELSE 'No action'
  END AS recommended_action
FROM vw_style_shift_day s;

-- Descriptive cluster gap (not causal)
CREATE OR REPLACE VIEW vw_style_cluster_gap AS
WITH bounds AS (
  SELECT MAX(dt) AS max_dt FROM vw_style_daily
),
ref AS (
  SELECT d.surface_id, d.style_cluster,
         AVG(d.ctr) AS ctr_ref, SUM(d.gmv) AS gmv_ref, SUM(d.n_impr) AS impr_ref
  FROM vw_style_daily d, bounds b
  WHERE d.dt BETWEEN b.max_dt - INTERVAL 27 DAY AND b.max_dt - INTERVAL 14 DAY
  GROUP BY 1, 2
),
cur AS (
  SELECT d.surface_id, d.style_cluster,
         AVG(d.ctr) AS ctr_cur, SUM(d.gmv) AS gmv_cur, SUM(d.n_impr) AS impr_cur
  FROM vw_style_daily d, bounds b
  WHERE d.dt BETWEEN b.max_dt - INTERVAL 13 DAY AND b.max_dt
  GROUP BY 1, 2
)
SELECT
  COALESCE(c.surface_id, r.surface_id) AS surface_id,
  COALESCE(c.style_cluster, r.style_cluster) AS style_cluster,
  r.ctr_ref, c.ctr_cur,
  (c.ctr_cur - r.ctr_ref) AS ctr_delta,
  r.impr_ref, c.impr_cur,
  (COALESCE(c.impr_cur, 0) - COALESCE(r.impr_ref, 0)) AS impr_delta,
  r.gmv_ref, c.gmv_cur
FROM ref r
FULL OUTER JOIN cur c
  ON r.surface_id = c.surface_id AND r.style_cluster = c.style_cluster;

CREATE OR REPLACE VIEW vw_style_value_daily AS
SELECT
  b.*,
  v_ctr.unit_value AS value_per_ctr_point,
  v_brand.unit_value AS brand_incident_cost,
  (b.brand_flag_rate * b.n_impr * COALESCE(v_brand.unit_value, 0)) AS brand_cost,
  (b.n_impr * b.ctr * COALESCE(v_ctr.unit_value, 0)) AS ctr_monetized_proxy
FROM vw_style_vs_concept_split b
LEFT JOIN dim_value_assumption v_ctr
  ON v_ctr.surface_id = b.surface_id AND v_ctr.metric = 'value_per_ctr_point'
LEFT JOIN dim_value_assumption v_brand
  ON v_brand.surface_id = b.surface_id AND v_brand.metric = 'brand_incident_cost';

CREATE OR REPLACE VIEW vw_style_roi_rollup AS
SELECT
  surface_id,
  SUM(CASE WHEN biz_axis_bucket = 'style_only_portrait' THEN 1 ELSE 0 END) AS n_style_only_days,
  SUM(CASE WHEN biz_axis_bucket = 'concept_only_quality' THEN 1 ELSE 0 END) AS n_concept_only_days,
  SUM(CASE WHEN biz_axis_bucket = 'joint_style_and_concept' THEN 1 ELSE 0 END) AS n_joint_days,
  AVG(CASE WHEN biz_axis_bucket = 'style_only_portrait' THEN ctr END) AS avg_ctr_style_only,
  AVG(CASE WHEN biz_axis_bucket = 'quiet_or_weak' THEN ctr END) AS avg_ctr_quiet,
  AVG(CASE WHEN biz_axis_bucket = 'style_only_portrait' THEN brand_cost END) AS avg_brand_cost_style_only,
  AVG(CASE WHEN biz_axis_bucket = 'quiet_or_weak' THEN brand_cost END) AS avg_brand_cost_quiet
FROM vw_style_value_daily
GROUP BY surface_id;

CREATE OR REPLACE VIEW vw_style_workorders AS
SELECT
  dt,
  surface_id,
  biz_axis_bucket,
  style_domain_auc,
  concept_fired,
  covar_fired,
  ctr,
  brand_flag_rate,
  recommended_action
FROM vw_style_vs_concept_split
ORDER BY
  CASE biz_axis_bucket
    WHEN 'joint_style_and_concept' THEN 0
    WHEN 'concept_only_quality' THEN 1
    WHEN 'style_only_portrait' THEN 2
    ELSE 3
  END,
  brand_flag_rate DESC;
