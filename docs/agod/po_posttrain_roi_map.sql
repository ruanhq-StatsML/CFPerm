-- PO-risk post-train ROI map
-- Direct map: "next focus (which modality / which part)" → ROI dimension
-- Maintain THIS file only.

-- ===========================================================================
-- 0) Answer in one line
-- ===========================================================================
-- Yes: we *constrain update focus* to a subset —
--   A) modality subset (which towers get BWD/steps/LR)
--   B) hard-row subset (which observations get weight after reject)
-- Not: permanently discard the rest of the dataset; FWD/serving still sees all.

-- ===========================================================================
-- 1) Focus → ROI dimension (the vocabulary to maintain)
-- ===========================================================================
CREATE TABLE IF NOT EXISTS po_roi_focus_dim (
  focus_id        TEXT PRIMARY KEY,
  -- what we decide "next"
  asks            TEXT NOT NULL,      -- decision question
  subset_of       TEXT NOT NULL,      -- modality | row | (not dataset_partition)
  how_chosen      TEXT NOT NULL,      -- rule
  roi_dim         TEXT NOT NULL,      -- cost / benefit / constraint columns
  roi_formula     TEXT NOT NULL,      -- how ROI is read
  ship_rule       TEXT NOT NULL
);

INSERT OR REPLACE INTO po_roi_focus_dim VALUES
(
  'next_modality',
  '下一步应侧重哪一个模态（概念塔）？',
  'modality',
  'top_concept_mod = argmax L_m (chronic); top_spike_mod = argmax S_m (this window); '
  || 'active_set = {m: L_m >= q30(L)}; dump steps to top_spike among active; LR from fused alpha',
  'cost=flops_rel,wall_clock; benefit=t_to_acc_star↓; constraint=delta_acc',
  'ROI_A = flops_saved / acc_risk   with acc_risk=1 if delta_acc<-0.005 else 0.05',
  'pass if flops_rel<1 AND delta_acc>=-0.005'
),
(
  'next_hard_rows',
  'reject 后应强调数据集/batch 的哪一部分（难行）？',
  'row',
  'only if rejected=1; rank by PO_i=|Y-mu|; w_i∝sqrt(PO_i) (soft); calm => w=1',
  'cost=fit_wall_clock; benefit=next_mse_drop_vs_uniform, hard_p_at_20',
  'ROI_B = next_mse_drop / fit_wall_clock_s',
  'pass if rejected windows have next_mse_drop>0'
),
(
  'joint_focus',
  '同窗：塔侧重 +（若 reject）难行侧重如何一起读 ROI？',
  'modality+row',
  'every t: apply next_modality; if reject: also next_hard_rows; orthogonal actuators',
  'see ROI_A and ROI_B; ship_gate = pass_A AND (pass_B OR no_reject)',
  'report both; do not average away reject-only B into calm windows',
  'v_roi_ship_gate.pass_modality_emphasis=1; pass_hard_upweight in (1,NULL)'
);

-- ===========================================================================
-- 2) Logic registry (A/B)
-- ===========================================================================
CREATE TABLE IF NOT EXISTS po_roi_logic (
  logic_id        TEXT PRIMARY KEY,
  focus_id        TEXT NOT NULL REFERENCES po_roi_focus_dim(focus_id),
  resolution      TEXT NOT NULL,
  sensor          TEXT NOT NULL,
  gate            TEXT,
  actuator        TEXT NOT NULL,
  notes           TEXT
);

INSERT OR REPLACE INTO po_roi_logic VALUES
('modality_emphasis', 'next_modality', 'modality',
 'PO_m → L (EMA·proto) + S (ΔPO) → α',
 NULL,
 'freeze<-L; step_dump<-S; LR<-α',
 'Subset focus = active modality towers only for BWD/steps.'),
('hard_upweight', 'next_hard_rows', 'observation',
 'PO_i=|Y-μ|',
 'reject',
 'w_i∝√PO_i on Fit_{t+1}',
 'Subset focus = high-residual rows inside rejected batch (soft weights).');

-- ===========================================================================
-- 3) Window log (fill from train loop)
-- ===========================================================================
CREATE TABLE IF NOT EXISTS po_roi_window_log (
  run_id           TEXT NOT NULL,
  window_t         INTEGER NOT NULL,
  pack             TEXT,
  -- decision outputs (the "next focus" answer)
  top_concept_mod  TEXT,              -- next chronic modality
  top_spike_mod    TEXT,              -- next step-dump modality
  active_mods_csv  TEXT,              -- modality subset this window
  n_active_mods    INTEGER,
  n_mods           INTEGER,
  -- hard-row focus
  rejected         INTEGER NOT NULL DEFAULT 0,
  n_rows           INTEGER,
  n_hard_top20     INTEGER,           -- |{i: rank by PO in top 20%}|
  -- ROI numerics
  flops_rel        REAL,
  wall_clock_s     REAL,
  steps_used       INTEGER,
  t_to_acc_star    INTEGER,
  acc              REAL,
  acc_equal        REAL,
  next_mse         REAL,
  next_mse_uniform REAL,
  hard_p_at_20     REAL,
  PRIMARY KEY (run_id, window_t)
);

-- ===========================================================================
-- 4) Per-window "what to emphasize next" (direct ROI map)
-- ===========================================================================
CREATE VIEW IF NOT EXISTS v_roi_next_focus AS
SELECT
  run_id,
  window_t,
  pack,
  -- A: which modality next
  top_concept_mod                         AS next_modality_chronic,
  top_spike_mod                           AS next_modality_spike,
  active_mods_csv                         AS modality_subset,
  CAST(n_active_mods AS REAL) / NULLIF(n_mods, 0)
                                          AS modality_subset_frac,
  (1.0 - flops_rel)                       AS flops_saved,
  (acc - acc_equal)                       AS delta_acc,
  CASE
    WHEN flops_rel < 1.0 AND (acc - acc_equal) >= -0.005 THEN 1
    ELSE 0
  END                                     AS roi_pass_modality,
  -- B: which rows next (only meaningful if rejected)
  rejected,
  CASE WHEN rejected = 1
       THEN CAST(n_hard_top20 AS REAL) / NULLIF(n_rows, 0)
       ELSE NULL END                      AS hard_row_subset_frac,
  CASE WHEN rejected = 1
       THEN (next_mse_uniform - next_mse) ELSE NULL END
                                          AS next_mse_drop,
  CASE
    WHEN rejected = 0 THEN NULL
    WHEN (next_mse_uniform - next_mse) > 0 THEN 1
    ELSE 0
  END                                     AS roi_pass_hard_rows
FROM po_roi_window_log;

-- ===========================================================================
-- 5) Aggregate ROI by focus dimension
-- ===========================================================================
CREATE VIEW IF NOT EXISTS v_roi_modality_emphasis AS
SELECT
  run_id,
  pack,
  -- which mods were emphasized most often
  (
    SELECT top_concept_mod FROM po_roi_window_log w2
    WHERE w2.run_id = w.run_id AND IFNULL(w2.pack,'') = IFNULL(w.pack,'')
    GROUP BY top_concept_mod
    ORDER BY COUNT(*) DESC LIMIT 1
  )                                       AS dominant_next_modality,
  AVG(flops_rel)                          AS avg_flops_rel,
  AVG(1.0 - flops_rel)                    AS avg_flops_saved,
  AVG(acc - acc_equal)                    AS avg_delta_acc,
  AVG(CAST(n_active_mods AS REAL) / NULLIF(n_mods, 0))
                                          AS avg_modality_subset_frac,
  AVG(1.0 - flops_rel)
    / NULLIF(AVG(CASE WHEN (acc - acc_equal) < -0.005 THEN 1.0 ELSE 0.05 END), 0)
                                          AS roi_flops_per_acc_risk,
  MIN(t_to_acc_star)                      AS t_to_acc_star
FROM po_roi_window_log w
GROUP BY run_id, pack;

CREATE VIEW IF NOT EXISTS v_roi_hard_upweight AS
SELECT
  run_id,
  pack,
  COUNT(*)                                AS n_reject,
  AVG(CAST(n_hard_top20 AS REAL) / NULLIF(n_rows, 0))
                                          AS avg_hard_row_subset_frac,
  AVG(next_mse_uniform - next_mse)        AS avg_next_mse_drop,
  AVG(hard_p_at_20)                       AS avg_hard_p20,
  AVG(wall_clock_s)                       AS avg_fit_wall_clock_s,
  AVG(next_mse_uniform - next_mse)
    / NULLIF(AVG(wall_clock_s), 0)        AS roi_mse_drop_per_sec
FROM po_roi_window_log
WHERE rejected = 1
GROUP BY run_id, pack;

CREATE VIEW IF NOT EXISTS v_roi_ship_gate AS
SELECT
  a.run_id,
  a.pack,
  a.dominant_next_modality,
  a.avg_modality_subset_frac,
  a.avg_flops_rel,
  a.avg_delta_acc,
  a.roi_flops_per_acc_risk,
  b.n_reject,
  b.avg_hard_row_subset_frac,
  b.avg_next_mse_drop,
  b.roi_mse_drop_per_sec,
  CASE WHEN a.avg_delta_acc >= -0.005 AND a.avg_flops_rel < 1.0
       THEN 1 ELSE 0 END                  AS pass_modality_emphasis,
  CASE WHEN b.n_reject IS NULL THEN NULL
       WHEN b.avg_next_mse_drop > 0 THEN 1 ELSE 0 END
                                          AS pass_hard_upweight
FROM v_roi_modality_emphasis a
LEFT JOIN v_roi_hard_upweight b
  ON a.run_id = b.run_id AND IFNULL(a.pack,'') = IFNULL(b.pack,'');

-- ===========================================================================
-- 6) Convenience: answer "subset focus?" for operators
-- ===========================================================================
CREATE VIEW IF NOT EXISTS v_roi_subset_meaning AS
SELECT
  focus_id,
  subset_of,
  asks,
  'YES: constrain *update focus* to this subset; NOT: drop other data from the corpus forever'
    AS interpretation,
  roi_dim,
  ship_rule
FROM po_roi_focus_dim;
