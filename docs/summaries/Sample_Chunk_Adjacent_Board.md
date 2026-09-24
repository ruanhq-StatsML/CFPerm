# Sample-chunk adjacent board — what the HGB number actually means

## Short answer

**HGB AUC on this board is a next-chunk transfer probe**, not a model claim.

Protocol per adjacent pair \(t \to t{+}1\):

1. Sort rows (time / stamp / index).
2. Cut every \(N\) samples (\(N \in \{1000,2000\}\)) → chunk ids.
3. On chunk \(t\): binarize \(Y\) (native binary, else early quantile), run
   `StandardScaler → VarianceThreshold → SelectKBest(f_classif)` → fit a shallow HGB.
4. On chunk \(t{+}1\): **freeze** that feature map + HGB; score ROC-AUC.

So AUC asks: *associations that looked good for \(Y\) on this batch — do they still
rank \(Y\) on the next batch?*

It does **not** mean:

- causal tip / graph tip / “we found the driver”
- that HGB is the right production model
- that high AUC = “good OOD detector”
- that the window cut is an estimand for policy (it is a **business grain**)

## Read three numbers together

| signal | question |
|---|---|
| \(\Delta\bar Y\) | Did the outcome *level* move across the cut? |
| HGB AUC | Did the *selected association* transfer? |
| Jaccard(top-5) | Did *which features* stay the same? |

Typical readings (judgment, not a hard rule):

| pattern | reading |
|---|---|
| low AUC | association does not travel (e.g. DiffusionDB prompt tokens ≈ 0.60) |
| high AUC + high Jaccard | persistent drivers (often intensity / meteorology / planted kinematics) |
| high AUC + low Jaccard | still predictive, but the *which* set is shifting (regime / mix change) |

## Why HGB at all?

HGB is a **cheap nonlinear probe** after univariate SelectKBest — same role as in
the DiffusionDB FSDS smoke. Alternatives (LogReg, RF) would answer the same
transfer question; we keep HGB so boards stay comparable across packs.

If you only want “which tokens moved,” look at `top_fsds` / `top_cmean` and
Jaccard — the AUC is optional corroboration that the ranking was not pure noise.

## Judgment calls this iteration (not mechanical follow)

1. **Sample-count grain** stays primary (\(N=1000/2000\)) — calendar width is a
   different estimand; do not mix them in one cell.
2. **Tencent split panels**:
   - `tencent_gr` (ops): convert-leak dropped, volume kept → expect ~1.0 AUC
     (“more exposure transfers”).
   - `tencent_gr_content`: volume also dropped → AUC falls on early pairs
     (e.g. 0.99 → ~0.65) and tops flip to credit/share ranks. Residual mean
     AUC can stay mid-high if later chunks recover — still shows intensity
     was doing most of the ops-board work.
3. **Cross-pack scatter**: mean AUC vs \(|\Delta\bar Y|\), marker size ∝ Jaccard —
   one glance for “transfer vs level-shift vs driver stability.”
   Empirical contrast: DiffusionDB weak (~0.60); Metro@1000 high AUC + low
   Jaccard (shifting drivers); Waymo high AUC + high Jaccard (planted persist).

## Run

```bash
PYTHONPATH=. python3 scripts/run_sample_chunk_adjacent_board.py \
  --chunk-sizes 1000,2000 \
  --out results/sample_chunk_adjacent_board
```

Artifacts: per-pack `*_board.png`, `cross_pack_transfer.png`, `SAMPLE_CHUNK_BOARD.md`,
`summary.json`.

## Datasets

| pack | X | Y | sort | expected board reading |
|---|---|---|---|---|
| diffusiondb | prompt TF-IDF | `image_nsfw` | timestamp | weak transfer |
| tencent_gr | edge numerics (−convert leak) | `y_convert` | `e_last_ts` | strong transfer via intensity |
| tencent_gr_content | above − volume | `y_convert` | `e_last_ts` | AUC should drop if intensity dominated |
| waymo_proxy | kinematics proxy | next_disp | row | strong + stable (planted drift) |
| metro_interstate | weather + calendar | traffic_volume | date_time | strong transfer |
| beijing_pm25 | meteo + lag | pm2.5 | stamp | strong transfer |

## Related

- Calendar-width board: [`DiffusionDB_Temporal_Attribution.md`](./DiffusionDB_Temporal_Attribution.md)
- Cross-domain feature methods: [`Cross_Domain_Feature_Methods.md`](./Cross_Domain_Feature_Methods.md)
