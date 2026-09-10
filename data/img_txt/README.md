# Drop image/text feature files under dataset subfolders:

```
data/img_txt/indiana_cxr/
  img_feats.npy   # (7426, 512) — Indiana/OpenI CXR image embeddings
  txt_feats.npy   # (7426, 512) — report text embeddings
  labels.npy
  df_metadata.csv # batch = projection (Frontal vs Lateral)

data/img_txt/mm_train_test/
  img_feats.npy   # (4034, 512) — train∪test image embeddings
  txt_feats.npy   # (4034, 512)
  labels.npy
  df_metadata.csv # batch = split (train vs test)
```

Then run:

```bash
python3 scripts/run_img_txt_fsds.py
# or one dataset:
python3 scripts/run_img_txt_fsds.py --data-dir data/img_txt/indiana_cxr
```

Source Drive IDs (user-provided): labels `1g45…`, img `13T3…`, txt `1zGP…`, meta `1vyQ…`,
zip bundle `1Nfj…` (X_{train,test}_{img,txt,sim}.npy), unused all-identical `1LEm…`.
