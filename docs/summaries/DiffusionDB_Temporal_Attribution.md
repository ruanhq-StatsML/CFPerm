# DiffusionDB temporal attribution (FSDS)

## Definition (prototype)

| Symbol | Meaning | Implementation here |
|---|---|---|
| **X** | Prompt token embedding | Lightweight TF-IDF / hashing bag-of-tokens (`d≤512`) |
| **Y** | Image attribute score | `image_nsfw` from metadata (no PNG download) |
| **T** | Time window | Equal-count early vs late along `timestamp` |

Attribution goal: which prompt tokens drive temporal drift in **Y**.

Signed direction: \(\Delta\bar Y = \bar Y_{\mathrm{late}}-\bar Y_{\mathrm{early}}\) (正向/负向).

> Full CLIP token grid `(77, 768)` can replace the light embedding later; the protocol stays the same.

## Data

```bash
hf download poloclub/diffusiondb --type dataset \
  --include metadata.parquet --local-dir data/diffusiondb
```

`metadata.parquet` ≈ 186MB, 2M rows (gallery span ≈ 2022-08-06 → 2022-08-20).

## Run

```bash
PYTHONPATH=. python3 scripts/run_diffusiondb_temporal_fsds.py \
  --n-sample 4000 --max-features 512 --select-k 40 \
  --concat-hyperparams --window-scheme equal_count \
  --out results/diffusiondb_temporal_fsds_hp

PYTHONPATH=. python3 scripts/run_diffusiondb_window_sweep.py \
  --n-sample 2000 --out results/diffusiondb_window_sweep
```

## Methods (feature selection — not graph tip)

1. **Covariate drift**: RF Domain VIMP (`X → T`)
2. **cmean**: signed \(\mu_{\mathrm{late}}-\mu_{\mathrm{early}}\) per feature
3. **FSDS**: `StandardScaler → VarianceThreshold → SelectKBest → HGB/LogReg`  
   predicting high-`Y` (threshold = early-window quantile); train early / test late
4. **Blend rank**: mean of percentile ranks (`f_score`, `vimp_cov`, `|Δ|`)

### X = prompt ⊕ hyperparams

Default \(X\) = TF-IDF/hash of prompt. Optional `--concat-hyperparams` stacks
`cfg`, `step`, sampler one-hot onto \(X\) so generation knobs are not absorbed
into token VIMP. (Qwen/CLIP embedding can concat the same way later.)

### T window schemes (not interchangeable)

| scheme | balances | side effect |
|---|---|---|
| `equal_count` | \(n\) per window | calendar width differs |
| `equal_time` | calendar span | \(n\) imbalances |
| `width` | fixed hours/bin | early/late = first/last occupied bin |

Sweep artifact: `results/diffusiondb_window_sweep/` — \(\Delta\bar Y\), HGB AUC,
`n_by_T`, `span_hours_by_T`, top features all move with the cut.

### Adjacent-batch board (业务切窗看板)

Same FS, no new estimator: pick width (5/10/20/60 min …), `groupby T`,
adjacent pairs \(t\to t{+}1\) as train→test, dump joint table + PNG.

```bash
PYTHONPATH=. python3 scripts/run_diffusiondb_adjacent_board.py \
  --widths-min 5,10,20,60 --out results/diffusiondb_adjacent_board
```

Just a rolling visualization board for business time grains.

### Sample-chunk adjacent board (每 N 条切窗)

Same board, no calendar weirdness: sort → cut every **1000 / 2000** rows →
adjacent chunks \(t\to t{+}1\) as train→test. Ported across DiffusionDB,
Tencent-GR edges, Waymo proxy, Metro Interstate, Beijing PM2.5.

```bash
PYTHONPATH=. python3 scripts/run_sample_chunk_adjacent_board.py \
  --chunk-sizes 1000,2000 --out results/sample_chunk_adjacent_board
```

Artifact: `results/sample_chunk_adjacent_board/` — per-dataset PNG + `SAMPLE_CHUNK_BOARD.md`.

## Smoke notes

- \(Y\) = continuous `image_nsfw` (not JSON aggregate); FSDS uses early-quantile binary
- With `--concat-hyperparams`, `hp_cfg` / `hp_step` appear in blend/VIMP (`hp_vimp_share>0`)

## Related

See [`StableDiffusion_Eval_Logic.md`](./StableDiffusion_Eval_Logic.md) for why SD-side evaluation (CFG / sampler / which Y head) is tricky relative to this prototype.
