-- =============================================================================
-- 13) HF hop knobs → 客服费率桥接（可对账）
-- 业务问题：贡献账的 before/after 费率，到底是哪几个 HF knobs 驱动的？
-- =============================================================================

CREATE OR REPLACE VIEW vw_cs_assist_hf_knob_bridge AS
SELECT
  k.surface_id,
  k.source AS hf_source,
  k.dataset AS hf_dataset,
  k.quiet_halluc,
  k.fire_halluc AS fire_halluc_biz,
  k.raw_fire_halluc,
  k.acted_halluc_scale,
  k.hop_ratio,
  k.rag_low,
  k.rag_ok,
  k.rag_threshold,
  k.precision_at_10,
  k.fired_at_cut,
  k.ignore_ticket_bump,
  k.ignore_refund_bump,
  ROUND(100.0 * r.ticket_rate_quiet, 2) AS ticket_rate_quiet_pct,
  ROUND(100.0 * r.ticket_rate_ignored, 2) AS ticket_rate_ignored_pct,
  ROUND(100.0 * r.ticket_rate_acted, 2) AS ticket_rate_acted_pct,
  ROUND(100.0 * r.contain_rate_quiet, 2) AS contain_rate_quiet_pct,
  ROUND(100.0 * r.contain_rate_ignored, 2) AS contain_rate_ignored_pct,
  ROUND(100.0 * r.contain_rate_acted, 2) AS contain_rate_acted_pct,
  ROUND(100.0 * (r.contain_rate_acted - r.contain_rate_ignored), 2) AS contain_lift_pp,
  e.tickets_avoided,
  e.refunds_avoided,
  e.extra_sessions_contained,
  e.incremental_yen AS gross_yen,
  n.net_incremental_yen AS net_yen,
  -- 驱动链是否「站得住」：跳变后 ignored 工单率应明显高于 quiet；acted 应压回去
  CASE
    WHEN r.ticket_rate_ignored > r.ticket_rate_quiet * 2
     AND r.ticket_rate_acted < r.ticket_rate_ignored * 0.5
     AND k.hop_ratio >= 2.0
     AND k.fired_at_cut = 1
    THEN 'knobs_align_with_rates'
    ELSE 'check_seed_or_knobs'
  END AS bridge_status,
  CONCAT(
    'HF knobs（hop_ratio≈', CAST(ROUND(k.hop_ratio, 1) AS VARCHAR),
    ', fire_halluc=', CAST(k.fire_halluc AS VARCHAR),
    ', rag_thr=', CAST(k.rag_threshold AS VARCHAR),
    '）→ 工单率 ',
    CAST(ROUND(100.0 * r.ticket_rate_ignored, 2) AS VARCHAR), '%→',
    CAST(ROUND(100.0 * r.ticket_rate_acted, 2) AS VARCHAR),
    '%，承接 +',
    CAST(ROUND(100.0 * (r.contain_rate_acted - r.contain_rate_ignored), 2) AS VARCHAR),
    'pp；净¥', CAST(CAST(n.net_incremental_yen AS BIGINT) AS VARCHAR), '。'
  ) AS external_one_liner_cn
FROM dim_hf_hop_knobs k
JOIN vw_cs_assist_rate_compare r ON r.surface_id = k.surface_id
JOIN vw_cs_assist_exec_summary e ON e.surface_id = k.surface_id
JOIN vw_cs_assist_net_increment n ON n.surface_id = k.surface_id;
