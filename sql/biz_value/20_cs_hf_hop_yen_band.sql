-- =============================================================================
-- 20) HF hop / rag 情景贡献带（由 demo 重 seed 写入 fct_hop_scenario_yen）
-- 业务问题：hop 强弱、rag 好坏直接驱动 seed 后，毛/净/全成本净差多少？动作混比怎么变？
-- =============================================================================

CREATE TABLE IF NOT EXISTS fct_hop_scenario_yen (
  surface_id                 VARCHAR,
  scenario                   VARCHAR,
  hop_ratio                  DOUBLE,
  fire_halluc                DOUBLE,
  rag_low                    DOUBLE,
  rag_ok                     DOUBLE,
  rag_threshold              DOUBLE,
  tickets_avoided            DOUBLE,
  refunds_avoided            DOUBLE,
  extra_sessions_contained   DOUBLE,
  gross_yen                  DOUBLE,
  net_yen                    DOUBLE,
  fully_loaded_net_yen       DOUBLE,
  retrieval_day_share        DOUBLE,
  top_action                 VARCHAR,
  delta_gross_vs_base        DOUBLE,
  delta_net_vs_base          DOUBLE,
  delta_fully_loaded_vs_base DOUBLE,
  delta_retrieval_share_vs_base DOUBLE,
  PRIMARY KEY (surface_id, scenario)
);

CREATE OR REPLACE VIEW vw_cs_assist_hop_yen_band AS
SELECT
  surface_id,
  scenario,
  ROUND(hop_ratio, 2) AS hop_ratio,
  ROUND(fire_halluc, 3) AS fire_halluc,
  ROUND(rag_low, 3) AS rag_low,
  ROUND(rag_ok, 3) AS rag_ok,
  ROUND(rag_threshold, 3) AS rag_threshold,
  ROUND(tickets_avoided, 1) AS tickets_avoided,
  ROUND(refunds_avoided, 1) AS refunds_avoided,
  ROUND(extra_sessions_contained, 1) AS extra_sessions_contained,
  ROUND(gross_yen, 0) AS gross_yen,
  ROUND(net_yen, 0) AS net_yen,
  ROUND(fully_loaded_net_yen, 0) AS fully_loaded_net_yen,
  ROUND(100.0 * retrieval_day_share, 1) AS retrieval_day_share_pct,
  top_action,
  ROUND(delta_gross_vs_base, 0) AS delta_gross_vs_base,
  ROUND(delta_net_vs_base, 0) AS delta_net_vs_base,
  ROUND(delta_fully_loaded_vs_base, 0) AS delta_fully_loaded_vs_base,
  ROUND(100.0 * delta_retrieval_share_vs_base, 1) AS delta_retrieval_share_pp
FROM fct_hop_scenario_yen
ORDER BY
  CASE scenario
    WHEN 'weak_hop' THEN 1
    WHEN 'base_hop' THEN 2
    WHEN 'strong_hop' THEN 3
    WHEN 'rag_worse' THEN 4
    WHEN 'rag_better' THEN 5
    ELSE 9
  END;

CREATE OR REPLACE VIEW vw_cs_assist_hop_yen_band_summary AS
WITH p AS (
  SELECT
    surface_id,
    MAX(CASE WHEN scenario = 'weak_hop' THEN fully_loaded_net_yen END) AS weak_fl,
    MAX(CASE WHEN scenario = 'base_hop' THEN fully_loaded_net_yen END) AS base_fl,
    MAX(CASE WHEN scenario = 'strong_hop' THEN fully_loaded_net_yen END) AS strong_fl,
    MAX(CASE WHEN scenario = 'base_hop' THEN net_yen END) AS base_net,
    MAX(CASE WHEN scenario = 'base_hop' THEN retrieval_day_share END) AS base_ret,
    MAX(CASE WHEN scenario = 'rag_worse' THEN retrieval_day_share END) AS worse_ret,
    MAX(CASE WHEN scenario = 'rag_better' THEN retrieval_day_share END) AS better_ret,
    MAX(CASE WHEN scenario = 'rag_worse' THEN fully_loaded_net_yen END) AS worse_fl,
    MAX(CASE WHEN scenario = 'rag_better' THEN fully_loaded_net_yen END) AS better_fl
  FROM fct_hop_scenario_yen
  GROUP BY surface_id
)
SELECT
  surface_id,
  ROUND(weak_fl, 0) AS weak_fully_loaded_net_yen,
  ROUND(base_fl, 0) AS base_fully_loaded_net_yen,
  ROUND(strong_fl, 0) AS strong_fully_loaded_net_yen,
  ROUND(base_net, 0) AS base_net_yen,
  ROUND(100.0 * base_ret, 1) AS base_retrieval_share_pct,
  ROUND(100.0 * worse_ret, 1) AS rag_worse_retrieval_share_pct,
  ROUND(100.0 * better_ret, 1) AS rag_better_retrieval_share_pct,
  ROUND(worse_fl, 0) AS rag_worse_fully_loaded_net_yen,
  ROUND(better_fl, 0) AS rag_better_fully_loaded_net_yen,
  CONCAT(
    'HF hop 弱→强：全成本净¥',
    CAST(CAST(ROUND(weak_fl, 0) AS BIGINT) AS VARCHAR),
    '→',
    CAST(CAST(ROUND(base_fl, 0) AS BIGINT) AS VARCHAR),
    '→',
    CAST(CAST(ROUND(strong_fl, 0) AS BIGINT) AS VARCHAR),
    '；rag 变差时 retrieval 天占比 ',
    CAST(ROUND(100.0 * base_ret, 0) AS VARCHAR),
    '%→',
    CAST(ROUND(100.0 * worse_ret, 0) AS VARCHAR),
    '%，变好时→',
    CAST(ROUND(100.0 * better_ret, 0) AS VARCHAR),
    '%。'
  ) AS external_one_liner_cn
FROM p;
