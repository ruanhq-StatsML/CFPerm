-- PO-risk post-train: modality emphasis + hard up-weight → ROI mapping
-- Goal: accelerate post-training. Acc is a constraint, not the ROI numerator alone.
-- Maintain this SQL; other dashboards are optional.

-- ---------------------------------------------------------------------------
-- 1) Dimension: two logics (frozen vocabulary)
-- ---------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS po_roi_logic (
  logic_id        TEXT PRIMARY KEY,   -- modality_emphasis | hard_upweight
  resolution      TEXT NOT NULL,      -- modality | observation
  sensor          TEXT NOT NULL,      -- PO_m | PO_i
  gate            TEXT,               -- NULL = every window; reject = RFPerm/CFPerm
  actuator        TEXT NOT NULL,
  cost_proxy      TEXT NOT NULL,      -- what we save
  benefit_proxy   TEXT NOT NULL,      -- what we gain (under Acc constraint)
  notes           TEXT
);

INSERT OR REPLACE INTO po_roi_logic VALUES
('modality_emphasis', 'modality', 'PO_m (long L + short S fuse)',
 NULL,
 'freeze<-L; step_dump<-S; LR<-alpha',
 'flops_rel, wall_clock_update, steps_used',
 't_to_acc_star, delta_acc_vs_equal',
 'High-residual concept tower gets budget; low-L BWD frozen.'),
('hard_upweight', 'observation', 'PO_i = |Y-mu| (opt. blend batch PO)',
 'reject',
 'w_i = sqrt(PO_i) / mean(sqrt(PO)); sample_weight on Fit_t+1',
 'fit_wall_clock, n_weighted_rows',
 'next_mse_drop_vs_uniform, hard_rank_precision_at_k',
 'Only after gate reject; calm windows stay w=1 (soft sqrt, not raw PO).');

-- ---------------------------------------------------------------------------
-- 2) Per-window / per-reject logs (fill from training loop)
-- ---------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS po_roi_window_log (
  run_id          TEXT NOT NULL,
  window_t        INTEGER NOT NULL,
  pack            TEXT,
  -- modality emphasis (logic A)
  flops_rel       REAL,          -- active_mods / M
  wall_clock_s    REAL,          -- update only
  steps_used      INTEGER,
  t_to_acc_star   INTEGER,       -- nullable until bar hit
  acc             REAL,
  acc_equal       REAL,          -- baseline equal policy same window
  top_concept_mod TEXT,
  top_spike_mod   TEXT,
  -- hard up-weight (logic B)
  rejected        INTEGER NOT NULL DEFAULT 0,  -- 0/1
  next_mse        REAL,
  next_mse_uniform REAL,
  hard_p_at_20    REAL,
  PRIMARY KEY (run_id, window_t)
);

-- ---------------------------------------------------------------------------
-- 3) ROI views (the only mapping we need to maintain)
-- ---------------------------------------------------------------------------
-- A: modality emphasis ROI ≈ benefit / cost
--    benefit = max(0, acc - (acc_equal - eps)) is constraint gate;
--    speed benefit = (1 - flops_rel) and/or drop in t_to_acc_star
CREATE VIEW IF NOT EXISTS v_roi_modality_emphasis AS
SELECT
  run_id,
  pack,
  AVG(flops_rel)                         AS avg_flops_rel,
  AVG(wall_clock_s)                      AS avg_wall_clock_s,
  AVG(acc - acc_equal)                   AS avg_delta_acc,
  AVG(1.0 - flops_rel)                   AS avg_flops_saved,
  -- ROI proxy: FLOPs saved per unit Acc risk (Acc drop counts against)
  AVG(1.0 - flops_rel)
    / NULLIF(AVG(CASE WHEN (acc - acc_equal) < -0.005 THEN 1.0 ELSE 0.05 END), 0)
                                         AS roi_flops_per_acc_risk,
  MIN(t_to_acc_star)                     AS t_to_acc_star
FROM po_roi_window_log
GROUP BY run_id, pack;

-- B: hard up-weight ROI only on rejected windows
CREATE VIEW IF NOT EXISTS v_roi_hard_upweight AS
SELECT
  run_id,
  pack,
  COUNT(*)                               AS n_reject,
  AVG(next_mse_uniform - next_mse)       AS avg_next_mse_drop,
  AVG(hard_p_at_20)                      AS avg_hard_p20,
  AVG(wall_clock_s)                      AS avg_fit_wall_clock_s,
  -- ROI proxy: next-MSE drop per reject-fit second
  AVG(next_mse_uniform - next_mse)
    / NULLIF(AVG(wall_clock_s), 0)       AS roi_mse_drop_per_sec
FROM po_roi_window_log
WHERE rejected = 1
GROUP BY run_id, pack;

-- Combined ship gate: A saves FLOPs with Acc OK; B helps on reject
CREATE VIEW IF NOT EXISTS v_roi_ship_gate AS
SELECT
  a.run_id,
  a.pack,
  a.avg_flops_rel,
  a.avg_delta_acc,
  a.roi_flops_per_acc_risk,
  b.n_reject,
  b.avg_next_mse_drop,
  b.roi_mse_drop_per_sec,
  CASE
    WHEN a.avg_delta_acc >= -0.005
     AND a.avg_flops_rel < 1.0
    THEN 1 ELSE 0
  END AS pass_modality_emphasis,
  CASE
    WHEN b.n_reject IS NULL THEN NULL
    WHEN b.avg_next_mse_drop > 0 THEN 1 ELSE 0
  END AS pass_hard_upweight
FROM v_roi_modality_emphasis a
LEFT JOIN v_roi_hard_upweight b
  ON a.run_id = b.run_id AND a.pack = b.pack;
