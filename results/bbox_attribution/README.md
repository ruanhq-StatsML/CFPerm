# Bounding-box hierarchical attribution board

Evidence board only. Business decides what the rankings mean.

## L1 · modality mass (COCO early/late + bbox block)

| Method | image | text | bbox | RF AUC |
|--------|-------|------|------|--------|
| RF | 0.525 | 0.442 | 0.033 | 0.570 |
| MMD | 0.472 | 0.478 | 0.050 | — |

## L2b · bbox subgroup mass (RF)

| geo | category_hist | top_boxes |
|-----|---------------|----------|
| 0.286 | 0.066 | 0.648 |

## BBox inject recovery (GT=bbox)

- RF: mass_on_bbox=0.9984, selAUC=0.8664, P@20=1.0
- MMD: mass_on_bbox=0.9583, selAUC=0.9335, P@20=1.0
