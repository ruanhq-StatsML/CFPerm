-- =============================================================================
-- 1) Hallucination REGIME → CS / refund / trust value
-- Axis: concept P(Y|X). Separate retrieval_gap vs generation_regime.
-- =============================================================================

CREATE OR REPLACE VIEW vw_halluc_daily AS
SELECT
  e.dt,
  e.surface_id,
  e.model_version,
  COUNT(*) AS n_serve,
  AVG(CASE WHEN e.y_halluc IS NOT NULL THEN e.y_halluc END) AS halluc_rate_labeled,
  SUM(CASE WHEN e.y_halluc IS NOT NULL THEN 1 ELSE 0 END) AS n_labeled,
  AVG(e.rag_hit) AS avg_rag_hit,
  AVG(e.cs_ticketed) AS ticket_rate,
  AVG(e.refunded) AS refund_rate,
  AVG(e.converted) AS cvr,
  SUM(e.gmv) AS gmv,
  SUM(e.cs_ticketed) AS n_tickets,
  SUM(e.refunded) AS n_refunds
FROM fct_serve_event e
WHERE e.surface_id IN ('shop_assistant', 'cs_bot', 'rag_qa')
GROUP BY 1, 2, 3;

CREATE OR REPLACE VIEW vw_halluc_fire_day AS
SELECT
  d.dt,
  d.surface_id,
  d.model_version,
  d.n_serve,
  d.halluc_rate_labeled,
  d.n_labeled,
  d.avg_rag_hit,
  d.ticket_rate,
  d.refund_rate,
  d.cvr,
  d.gmv,
  d.n_tickets,
  d.n_refunds,
  MAX(CASE WHEN s.signal_type = 'rfperm_fire' AND s.axis = 'concept' THEN s.fired ELSE 0 END) AS regime_fired,
  MAX(CASE WHEN s.signal_type = 'rfperm_fire' AND s.axis = 'concept' THEN s.score END) AS fire_ratio,
  MAX(CASE WHEN s.signal_type = 'judge_err_ratio' THEN s.score END) AS judge_err_ratio,
  AVG(CASE WHEN s.signal_type = 'po_risk0' THEN s.score END) AS mean_po_risk0,
  MAX(CASE WHEN s.notes = 'acted' THEN 1 ELSE 0 END) AS acted
FROM vw_halluc_daily d
LEFT JOIN fct_shift_signal s
  ON s.dt = d.dt
 AND s.surface_id = d.surface_id
 AND (s.model_version = d.model_version OR s.model_version IS NULL)
GROUP BY 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13;

-- Root-cause BUCKET for ops (not causal proof — routing ticket type)
CREATE OR REPLACE VIEW vw_halluc_root_bucket AS
SELECT
  f.*,
  CASE
    WHEN f.regime_fired = 1 AND COALESCE(f.avg_rag_hit, 1) < 0.35 THEN 'retrieval_gap'
    WHEN f.regime_fired = 1 AND COALESCE(f.avg_rag_hit, 0) >= 0.35 THEN 'generation_regime'
    WHEN f.regime_fired = 0 AND COALESCE(f.ticket_rate, 0) > 0.02 THEN 'labeler_or_product_other'
    ELSE 'quiet'
  END AS biz_bucket
FROM vw_halluc_fire_day f;

CREATE OR REPLACE VIEW vw_halluc_audit_yield AS
SELECT
  a.dt,
  a.surface_id,
  COUNT(*) AS n_queued,
  SUM(CASE WHEN a.audit_status = 'confirmed_bad' THEN 1 ELSE 0 END) AS n_confirmed_bad,
  SUM(CASE WHEN a.audit_status = 'false_alarm' THEN 1 ELSE 0 END) AS n_false_alarm,
  AVG(CASE WHEN a.rank_in_batch <= 10 AND a.audit_status = 'confirmed_bad' THEN 1.0 ELSE 0.0 END)
    AS precision_at_10_proxy,
  AVG(a.po_risk0) AS mean_po_queued
FROM fct_audit_queue a
GROUP BY 1, 2;

-- Unit economics join
CREATE OR REPLACE VIEW vw_halluc_value_daily AS
SELECT
  b.*,
  v_cs.unit_value AS cs_ticket_cost,
  v_rf.unit_value AS refund_unit_cost,
  b.n_tickets * COALESCE(v_cs.unit_value, 0) AS cs_cost,
  b.n_refunds * COALESCE(v_rf.unit_value, 0) AS refund_cost,
  b.n_tickets * COALESCE(v_cs.unit_value, 0)
    + b.n_refunds * COALESCE(v_rf.unit_value, 0) AS total_quality_cost
FROM vw_halluc_root_bucket b
LEFT JOIN dim_value_assumption v_cs
  ON v_cs.surface_id = b.surface_id AND v_cs.metric = 'cs_ticket_cost'
LEFT JOIN dim_value_assumption v_rf
  ON v_rf.surface_id = b.surface_id AND v_rf.metric = 'refund_unit_cost';

-- ROI: acting on fire days vs ignoring them
CREATE OR REPLACE VIEW vw_halluc_roi_rollup AS
SELECT
  surface_id,
  SUM(CASE WHEN regime_fired = 1 AND acted = 1 THEN 1 ELSE 0 END) AS n_fire_acted_days,
  SUM(CASE WHEN regime_fired = 1 AND acted = 0 THEN 1 ELSE 0 END) AS n_fire_ignored_days,
  AVG(CASE WHEN regime_fired = 1 AND acted = 1 THEN total_quality_cost END) AS avg_cost_fire_acted,
  AVG(CASE WHEN regime_fired = 1 AND acted = 0 THEN total_quality_cost END) AS avg_cost_fire_ignored,
  AVG(CASE WHEN regime_fired = 0 THEN total_quality_cost END) AS avg_cost_quiet,
  AVG(CASE WHEN regime_fired = 1 AND acted = 0 THEN total_quality_cost END)
    - AVG(CASE WHEN regime_fired = 1 AND acted = 1 THEN total_quality_cost END)
    AS est_daily_save_act_vs_ignore
FROM vw_halluc_value_daily
GROUP BY surface_id;

CREATE OR REPLACE VIEW vw_halluc_workorders AS
SELECT
  dt,
  surface_id,
  biz_bucket,
  regime_fired,
  fire_ratio,
  avg_rag_hit,
  ticket_rate,
  refund_rate,
  total_quality_cost,
  CASE biz_bucket
    WHEN 'retrieval_gap' THEN 'Refresh retriever/index; do NOT full SFT yet'
    WHEN 'generation_regime' THEN 'Gate sqrt(PO) / rollback model_version; Top-k po_risk0 audit'
    WHEN 'labeler_or_product_other' THEN 'Check labeler drift or product policy; RFPerm quiet'
    ELSE 'No action — anneal uniform'
  END AS recommended_action
FROM vw_halluc_value_daily
ORDER BY regime_fired DESC, total_quality_cost DESC;
