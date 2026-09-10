# COCO time-order + bounding-box block

Same early/late `image_id` split as `coco_time_order`, with an added
`bbox_feats.npy` block joined via caption → COCO `image_id` → instances.

Feature layout (`bbox_feature_schema.json`):
- geo summaries (8)
- category histogram (80)
- top-5 boxes × (cx, cy, w, h, area_frac, cat_norm)

```bash
python3 scripts/run_bbox_attribution_board.py
```
