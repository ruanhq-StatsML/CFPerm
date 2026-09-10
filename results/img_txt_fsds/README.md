# Image / Text FSDS · two datasets

FSDS attribution (not anomaly detection): **RF Domain Classifier** (covariate shift) + **PO-risk RF** (concept drift).  
Ranking = `np.argsort(-VIMP)[:20]` within each modality block after `X = [img | txt]`.

## Datasets

| name | n | d_img / d_txt | batch W |
|------|---|---------------|---------|
| `indiana_cxr` | 7426 | 512 / 512 | `projection`: Frontal vs Lateral |
| `mm_train_test` | 4034 | 512 / 512 | `split`: train vs test |

Image vs text assignment for Indiana: projection AUC ≈ 1 on one 512-block and ≈ 0.5 on the other; within-`uid` L2 is 0 on the text block (shared report across Frontal/Lateral).

## Board (ordered local indices)

### indiana_cxr

```python
feature_indices_covariate_shift = {
  "image": [330, 105, 235, 145, 278, 335, 494, 117, 293, 82, 111, 410, 154, 21, 115, 68, 0, 77, 431, 171],
  "text": [138, 187, 447, 38, 246, 423, 436, 476, 34, 426, 288, 115, 171, 444, 55, 310, 427, 364, 375, 210],
}
feature_indices_concept_drift = {
  "image": [133, 249, 146, 338, 57, 301, 170, 4, 370, 30, 511, 122, 248, 99, 28, 197, 199, 437, 268, 165],
  "text": [302, 241, 461, 489, 506, 196, 3, 257, 304, 256, 60, 72, 423, 428, 420, 268, 127, 361, 249, 0],
}
# CS modality mass ≈ image 0.996 / text 0.004  (RF domain AUC = 1.0)
# CD modality mass ≈ image 0.591 / text 0.409
```

### mm_train_test

```python
feature_indices_covariate_shift = {
  "image": [265, 377, 147, 124, 317, 108, 84, 316, 496, 459, 218, 348, 47, 18, 227, 38, 11, 114, 379, 479],
  "text": [261, 116, 414, 63, 244, 146, 504, 67, 505, 96, 469, 186, 164, 379, 412, 0, 259, 189, 211, 73],
}
feature_indices_concept_drift = {
  "image": [370, 39, 216, 175, 321, 304, 315, 489, 195, 360, 95, 28, 357, 336, 102, 66, 173, 209, 324, 501],
  "text": [115, 300, 501, 442, 113, 100, 23, 179, 96, 57, 348, 137, 453, 25, 1, 111, 419, 206, 67, 335],
}
# CS modality mass ≈ image 0.606 / text 0.394  (RF domain AUC ≈ 0.727)
# CD modality mass ≈ image 0.255 / text 0.745
```

## Artifacts

- `results/img_txt_fsds/*_fsds_feature_indices.json`
- `results/img_txt_fsds/img_txt_fsds_both_datasets_board.py`
- Runner: `python3 scripts/run_img_txt_fsds.py`
