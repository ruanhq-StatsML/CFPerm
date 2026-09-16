-- =============================================================================
-- Business-value SQL: Hallucination regimes + Style/Creative drift
-- Dialect: DuckDB (also mostly BigQuery/Snowflake ANSI)
--
-- Maps OnlineRFPerm / PO-risk / RF-domain signals → operational KPIs so the
-- toolkit can justify CS cost, refund, CTR, brand, and merge-gate spend.
--
-- Stance: fire / po_risk0 / style_domain_auc are MONITOR + PRIORITY signals.
-- They do NOT claim fact-check truth or causal creative attribution.
-- =============================================================================

CREATE TABLE IF NOT EXISTS dim_surface (
  surface_id   VARCHAR PRIMARY KEY,  -- shop_assistant | cs_bot | rag_qa | feed_caption | ad_creative | shop_image_gen
  biz_line     VARCHAR,              -- llm_assist | content | ads
  owner_team   VARCHAR
);

CREATE TABLE IF NOT EXISTS dim_creative (
  creative_id    VARCHAR PRIMARY KEY,
  campaign_id    VARCHAR,
  style_cluster  VARCHAR,
  artist_tone    VARCHAR,
  launch_dt      DATE
);

-- One row per generation / impression unit
CREATE TABLE IF NOT EXISTS fct_serve_event (
  event_id       VARCHAR PRIMARY KEY,
  dt             DATE,
  ts             TIMESTAMP,
  surface_id     VARCHAR,
  session_id     VARCHAR,
  user_id        VARCHAR,
  model_version  VARCHAR,
  creative_id    VARCHAR,
  modality       VARCHAR,   -- text | image | video | multi
  task_type      VARCHAR,
  -- labels (sparse OK)
  y_halluc       INTEGER,   -- 1 unfaithful / halluc, 0 ok, NULL unlabeled
  y_pref_win     INTEGER,   -- 1 preferred / chosen
  -- style register features
  style_len      DOUBLE,
  style_formal   DOUBLE,
  style_emoji    DOUBLE,
  style_lang     VARCHAR,
  -- RAG support
  rag_hit        DOUBLE,
  -- outcomes
  clicked        INTEGER,
  converted      INTEGER,
  gmv            DOUBLE,
  refunded       INTEGER,
  cs_ticketed    INTEGER,
  brand_flagged  INTEGER
);

-- Monitor outputs (OnlineRFPerm / PO / RF-domain)
CREATE TABLE IF NOT EXISTS fct_shift_signal (
  signal_id      VARCHAR PRIMARY KEY,
  dt             DATE,
  ts             TIMESTAMP,
  surface_id     VARCHAR,
  batch_id       VARCHAR,
  signal_type    VARCHAR,  -- rfperm_fire | po_risk0 | style_domain_auc | judge_err_ratio
  axis           VARCHAR,  -- concept | covariate | joint
  fired          INTEGER,
  score          DOUBLE,
  n_batch        INTEGER,
  model_version  VARCHAR,
  notes          VARCHAR   -- acted | ignored | NULL
);

CREATE TABLE IF NOT EXISTS fct_audit_queue (
  audit_id       VARCHAR PRIMARY KEY,
  dt             DATE,
  surface_id     VARCHAR,
  event_id       VARCHAR,
  batch_id       VARCHAR,
  po_risk0       DOUBLE,
  rank_in_batch  INTEGER,
  audit_status   VARCHAR,  -- queued | confirmed_bad | false_alarm | skipped
  auditor        VARCHAR,
  audited_at     TIMESTAMP
);

-- Controlled unit economics (do NOT hardcode in queries)
CREATE TABLE IF NOT EXISTS dim_value_assumption (
  as_of          DATE,
  surface_id     VARCHAR,
  metric         VARCHAR,
  unit_value     DOUBLE,
  currency       VARCHAR,
  source_note    VARCHAR,
  PRIMARY KEY (as_of, surface_id, metric)
);

-- HF landing hop knobs that drive CS seed fire/rag intensity
CREATE TABLE IF NOT EXISTS dim_hf_hop_knobs (
  surface_id           VARCHAR PRIMARY KEY,
  source               VARCHAR,
  dataset              VARCHAR,
  quiet_halluc         DOUBLE,
  fire_halluc          DOUBLE,
  raw_fire_halluc      DOUBLE,
  acted_halluc_scale   DOUBLE,
  hop_ratio            DOUBLE,
  rag_low              DOUBLE,
  rag_ok               DOUBLE,
  rag_threshold        DOUBLE,
  precision_at_10      DOUBLE,
  fired_at_cut         INTEGER,
  ignore_ticket_bump   DOUBLE,
  ignore_refund_bump   DOUBLE
);
