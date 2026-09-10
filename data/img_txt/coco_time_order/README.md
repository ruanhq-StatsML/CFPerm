# COCO time-order + interpretable bbox features

Batch W: early vs late `image_id` (order proxy).

BBox board uses **handcrafted named features only**
(`bbox_named_feats.npy` + `bbox_named_feature_names.json`):
counts, geometry, supercategory area/count, key-class presence.
Do not rank CLIP coordinate indices for business claims.

```bash
python3 scripts/run_bbox_named_attribution.py
```
