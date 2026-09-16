-- =============================================================================
-- 19) 周全成本净周报：周毛 − 审计分摊 − 动作日成本分摊 + WoW
-- 业务问题：这一周扣完审计人力和动作日成本，还净多少？周环比怎么变？
-- 口径：
--   audit_cost_alloc     = 当周 acted 日实际审计件成本（从 net_daily）
--   action_day_cost_alloc = 总动作日成本 × 当周 acted 天 / 总 acted 天
--   fully_loaded_net     = 周毛归因 − audit_cost_alloc − action_day_cost_alloc
-- =============================================================================

CREATE OR REPLACE VIEW vw_cs_assist_weekly_fully_loaded AS
WITH action_tot AS (
  SELECT
    surface_id,
    SUM(action_day_cost_yen) AS action_day_cost_yen_total,
    SUM(n_days) AS acted_days_total
  FROM vw_cs_assist_action_net
  GROUP BY surface_id
),
-- 当周 acted 臂真实审计件成本（不是按会话均摊，直接按日记账）
week_audit AS (
  SELECT
    date_trunc('week', CAST(dt AS TIMESTAMP))::DATE AS week_start,
    surface_id,
    ROUND(SUM(audit_cost_yen), 0) AS audit_cost_yen_week,
    ROUND(SUM(n_audits), 0) AS audits_week
  FROM vw_cs_assist_net_daily
  WHERE arm = 'fire_acted'
  GROUP BY 1, 2
),
w AS (
  SELECT
    week_start,
    surface_id,
    days_acted,
    yen_from_tickets,
    yen_from_refunds,
    yen_from_containment,
    gross_yen
  FROM vw_cs_assist_week_value_split
)
SELECT
  w.week_start,
  w.surface_id,
  w.days_acted,
  ROUND(w.gross_yen, 0) AS gross_yen,
  ROUND(w.yen_from_tickets, 0) AS yen_from_tickets,
  ROUND(w.yen_from_refunds, 0) AS yen_from_refunds,
  ROUND(w.yen_from_containment, 0) AS yen_from_containment,
  ROUND(COALESCE(a.audit_cost_yen_week, 0), 0) AS audit_cost_yen_alloc,
  ROUND(COALESCE(a.audits_week, 0), 0) AS audits_week,
  ROUND(
    COALESCE(t.action_day_cost_yen_total, 0)
      * w.days_acted
      / NULLIF(t.acted_days_total, 0),
    0
  ) AS action_day_cost_yen_alloc,
  ROUND(
    w.gross_yen
      - COALESCE(a.audit_cost_yen_week, 0)
      - COALESCE(t.action_day_cost_yen_total, 0)
          * w.days_acted
          / NULLIF(t.acted_days_total, 0),
    0
  ) AS fully_loaded_net_yen,
  ROUND(
    100.0 * (
      w.gross_yen
        - COALESCE(a.audit_cost_yen_week, 0)
        - COALESCE(t.action_day_cost_yen_total, 0)
            * w.days_acted
            / NULLIF(t.acted_days_total, 0)
    ) / NULLIF(w.gross_yen, 0),
    1
  ) AS fully_loaded_pct_of_gross
FROM w
LEFT JOIN week_audit a
  ON a.week_start = w.week_start AND a.surface_id = w.surface_id
LEFT JOIN action_tot t ON t.surface_id = w.surface_id
ORDER BY w.week_start;

-- 周全成本净环比（WoW）
CREATE OR REPLACE VIEW vw_cs_assist_weekly_fully_loaded_wow AS
SELECT
  week_start,
  surface_id,
  days_acted,
  gross_yen,
  audit_cost_yen_alloc,
  action_day_cost_yen_alloc,
  fully_loaded_net_yen,
  fully_loaded_pct_of_gross,
  LAG(fully_loaded_net_yen) OVER (
    PARTITION BY surface_id ORDER BY week_start
  ) AS prev_fully_loaded_net_yen,
  ROUND(
    fully_loaded_net_yen
      - LAG(fully_loaded_net_yen) OVER (
          PARTITION BY surface_id ORDER BY week_start
        ),
    0
  ) AS delta_fully_loaded_net_yen,
  ROUND(
    100.0 * (
      fully_loaded_net_yen
        - LAG(fully_loaded_net_yen) OVER (
            PARTITION BY surface_id ORDER BY week_start
          )
    ) / NULLIF(
      ABS(LAG(fully_loaded_net_yen) OVER (
        PARTITION BY surface_id ORDER BY week_start
      )),
      0
    ),
    1
  ) AS delta_fully_loaded_pct
FROM vw_cs_assist_weekly_fully_loaded;

-- 对账：周全成本净合计 ≈ 总账全成本净
CREATE OR REPLACE VIEW vw_cs_assist_weekly_fully_loaded_check AS
WITH week_sum AS (
  SELECT
    surface_id,
    ROUND(SUM(gross_yen), 0) AS week_gross_sum,
    ROUND(SUM(audit_cost_yen_alloc), 0) AS week_audit_sum,
    ROUND(SUM(action_day_cost_yen_alloc), 0) AS week_action_sum,
    ROUND(SUM(fully_loaded_net_yen), 0) AS week_fully_loaded_sum
  FROM vw_cs_assist_weekly_fully_loaded
  GROUP BY surface_id
)
SELECT
  f.surface_id,
  f.realized_gross_yen AS total_gross,
  f.audit_cost_yen AS total_audit,
  f.action_day_cost_yen AS total_action,
  f.fully_loaded_net_yen AS total_fully_loaded,
  w.week_gross_sum,
  w.week_audit_sum,
  w.week_action_sum,
  w.week_fully_loaded_sum,
  ROUND(w.week_fully_loaded_sum - f.fully_loaded_net_yen, 0) AS fully_loaded_gap_yen,
  CONCAT(
    '周全成本净合计¥',
    CAST(CAST(w.week_fully_loaded_sum AS BIGINT) AS VARCHAR),
    ' vs 总账¥',
    CAST(CAST(f.fully_loaded_net_yen AS BIGINT) AS VARCHAR),
    '（缺口¥',
    CAST(CAST(ROUND(w.week_fully_loaded_sum - f.fully_loaded_net_yen, 0) AS BIGINT) AS VARCHAR),
    '）；两周毛¥',
    CAST(CAST(w.week_gross_sum AS BIGINT) AS VARCHAR),
    '−审计¥',
    CAST(CAST(w.week_audit_sum AS BIGINT) AS VARCHAR),
    '−动作¥',
    CAST(CAST(w.week_action_sum AS BIGINT) AS VARCHAR),
    '。'
  ) AS external_one_liner_cn
FROM vw_cs_assist_fully_loaded_capture f
JOIN week_sum w ON w.surface_id = f.surface_id;

-- 对外一句（取最近一周 + 环比）
CREATE OR REPLACE VIEW vw_cs_assist_weekly_fully_loaded_summary AS
WITH latest AS (
  SELECT *
  FROM vw_cs_assist_weekly_fully_loaded_wow
  QUALIFY ROW_NUMBER() OVER (
    PARTITION BY surface_id ORDER BY week_start DESC
  ) = 1
),
chk AS (
  SELECT * FROM vw_cs_assist_weekly_fully_loaded_check
)
SELECT
  l.surface_id,
  l.week_start AS latest_week_start,
  l.days_acted AS latest_days_acted,
  l.gross_yen AS latest_gross_yen,
  l.audit_cost_yen_alloc AS latest_audit_yen,
  l.action_day_cost_yen_alloc AS latest_action_yen,
  l.fully_loaded_net_yen AS latest_fully_loaded_net_yen,
  l.prev_fully_loaded_net_yen,
  l.delta_fully_loaded_net_yen,
  c.week_fully_loaded_sum,
  c.total_fully_loaded,
  c.fully_loaded_gap_yen,
  CONCAT(
    '最近一周（',
    CAST(l.week_start AS VARCHAR),
    '）全成本净¥',
    CAST(CAST(l.fully_loaded_net_yen AS BIGINT) AS VARCHAR),
    '（毛¥',
    CAST(CAST(l.gross_yen AS BIGINT) AS VARCHAR),
    '−审计¥',
    CAST(CAST(l.audit_cost_yen_alloc AS BIGINT) AS VARCHAR),
    '−动作¥',
    CAST(CAST(l.action_day_cost_yen_alloc AS BIGINT) AS VARCHAR),
    '）',
    CASE
      WHEN l.delta_fully_loaded_net_yen IS NULL THEN ''
      WHEN l.delta_fully_loaded_net_yen >= 0 THEN CONCAT(
        '；较上周 +¥',
        CAST(CAST(l.delta_fully_loaded_net_yen AS BIGINT) AS VARCHAR)
      )
      ELSE CONCAT(
        '；较上周 ¥',
        CAST(CAST(l.delta_fully_loaded_net_yen AS BIGINT) AS VARCHAR)
      )
    END,
    '；周合计 vs 总账缺口¥',
    CAST(CAST(c.fully_loaded_gap_yen AS BIGINT) AS VARCHAR),
    '。'
  ) AS external_one_liner_cn
FROM latest l
JOIN chk c ON c.surface_id = l.surface_id;
