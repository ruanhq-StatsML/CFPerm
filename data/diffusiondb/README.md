# DiffusionDB local data

Download (not committed):

```bash
hf download poloclub/diffusiondb --type dataset \
  --include metadata.parquet --local-dir data/diffusiondb
```

Then run:

```bash
PYTHONPATH=. python3 scripts/run_diffusiondb_temporal_fsds.py --n-sample 8000
```
