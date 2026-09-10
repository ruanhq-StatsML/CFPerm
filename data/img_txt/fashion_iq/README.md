# Fashion-IQ CLIP embeddings (FSDS layout)

```
data/img_txt/fashion_iq/
  clip_embedding.npy   # (n, 1025) = [image(512) | text(512) | label]
  df1.npy              # train batch
  df2.npy              # test batch
  img_feats.npy / txt_feats.npy / labels.npy / df_metadata.csv
```

Last column = label; first half of features = image; second half = text.

```bash
python3 scripts/run_fashion_iq_fsds_benchmark.py
```

Runs FSDS `benchmark_feature_selection` subset: **RF Domain**, **coord-wise MMD**, **PO-risk**.

Note: Drive file `1LEmgx…` (4285×1024) is corrupted (all rows identical); this board uses the
train/test CLIP zip (`X_{train,test}_{img,txt}.npy` + `y_*.npy`).
