# Sample-chunk adjacent boards

Business grain = **sample count**, not calendar width.

1. Sort by time (or row order when no stamp).
2. Assign chunk id = `row_index // N` with \(N \in \{1000, 2000\}\).
3. Adjacent chunks \(t \to t{+}1\) = train → test.
4. Same FS smoke: SelectKBest → HGB; report \(\Delta\bar Y\), AUC, top tokens/feats.
5. Joint PNG board per dataset — visualization only, not a new estimator.

## Run

```bash
PYTHONPATH=. python3 scripts/run_sample_chunk_adjacent_board.py \
  --chunk-sizes 1000,2000 \
  --datasets diffusiondb,tencent_gr,waymo_proxy,metro_interstate,beijing_pm25 \
  --out results/sample_chunk_adjacent_board
```

## Datasets

| pack | X | Y | sort |
|---|---|---|---|
| diffusiondb | prompt TF-IDF | `image_nsfw` | timestamp |
| tencent_gr | edge `feature_grid` numerics | `y_convert` | `e_last_ts` |
| waymo_proxy | kinematics proxy | next_disp | row index |
| metro_interstate | weather + calendar OH | traffic_volume | date_time |
| beijing_pm25 | meteo + lag | pm2.5 | stamp |

## Related

- Calendar-width board: [`DiffusionDB_Temporal_Attribution.md`](./DiffusionDB_Temporal_Attribution.md) (adjacent-batch / width mins).
- Cross-domain feature methods: [`Cross_Domain_Feature_Methods.md`](./Cross_Domain_Feature_Methods.md).
