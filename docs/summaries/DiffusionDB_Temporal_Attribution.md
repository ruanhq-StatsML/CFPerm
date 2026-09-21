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
  --n-sample 8000 --max-features 512 --select-k 40 \
  --out results/diffusiondb_temporal_fsds
```

## Methods

1. **Covariate drift**: RF Domain VIMP (`X → T`)
2. **Tip-cmean**: signed \(\mu_{\mathrm{late}}-\mu_{\mathrm{early}}\) per token feature
3. **FSDS**: `StandardScaler → VarianceThreshold → SelectKBest → HGB/LogReg`  
   predicting high-`Y` (threshold = early-window quantile); train early / test late
4. **Blend rank**: mean of percentile ranks (`f_score`, `vimp_cov`, `|Δ|`)

## Smoke result (n=8000)

- \(\bar Y\): early 0.167 → late 0.218, \(\Delta\bar Y=+0.051\) → **正向**
- FSDS HGB AUC(late) ≈ 0.63
- Top blend tokens include style/quality phrases (`sharp focus`, `greg rutkowski`, `art by`, …)

Artifacts: `results/diffusiondb_temporal_fsds/`.

## Related

See [`StableDiffusion_Eval_Logic.md`](./StableDiffusion_Eval_Logic.md) for why SD-side evaluation (CFG / sampler / which Y head) is tricky relative to this prototype.
