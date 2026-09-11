# Cross-modality + named bbox feature selection

Canonical LaTeX (paper-ready):

- `docs/method/CrossModality_BBox_SelectedFeatures_tables_only.tex`
- mirror: `results/bbox_attribution/CrossModality_BBox_SelectedFeatures_tables_only.tex`

Board JSON: `bbox_business_schema_board.json` (COCO early/late concat).

## Cross-modality mass (COCO early/late)

| Block | Mass share |
|---|---:|
| Image CLIP | 0.525 |
| Text CLIP | 0.449 |
| BBox named | 0.026 |
| RF Domain AUC | 0.600 |

CLIP blocks are modality-level only. BBox block is handcrafted named features (business-readable).

## Selected bbox features (top-20 RF Domain VIMP)

| Rank | Family | Feature | VIMP |
|---:|---|---|---:|
| 1 | position | `cx_std` | 0.00501 |
| 2 | aspect | `aspect_std` | 0.00269 |
| 3 | size | `h_mean` | 0.00240 |
| 4 | aspect | `aspect_mean` | 0.00219 |
| 5 | size | `w_mean` | 0.00163 |
| 6 | size | `h_std` | 0.00152 |
| 7 | area | `area_mean` | 0.00136 |
| 8 | position | `cx_mean` | 0.00127 |
| 9 | spatial rel. | `pairwise_iou_mean` | 0.00119 |
| 10 | category | `cat_hist_1 (person)` | 0.00108 |
| 11 | size | `w_std` | 0.00101 |
| 12 | spatial dist. | `nn_center_dist_mean` | 0.00096 |
| 13 | area | `area_std` | 0.00093 |
| 14 | position | `cy_mean` | 0.00088 |
| 15 | position | `cy_std` | 0.00076 |
| 16 | category | `cat_hist_42 (surfboard)` | 0.00019 |
| 17 | occlusion | `iscrowd_ratio` | 0.00012 |
| 18 | category | `cat_hist_16 (bird)` | 0.00010 |
| 19 | category | `cat_hist_15 (bench)` | 0.00008 |
| 20 | category | `cat_hist_32 (tie)` | 0.00007 |

## Inject consensus ($S^{\star}$, recovered on all 3 boards)

`n_objects`, `n_people`, `area_frac_person`, `area_frac_vehicle`, `region_center_area_prop`, `region_periphery_area_prop`, `quad_TL/TR/BL/BR_area_prop` — **10/10** on time-order / outdoor-indoor / center-split.

Details: `bbox_named_consensus_board.json`, `BBox_Named_Consensus_tables_only.tex`.

## COCO multi-layer logic (IoU / crowding / area)

Beyond conventional mean/std — hierarchical board:

- L1 modality → L2 bbox family (geometry / iou / crowding / area / semantic) → L3 named ranks
- LaTeX: `docs/method/COCO_Layered_BBox_tables_only.tex`
- Board: `bbox_layered_coco_board.json`
- Builder: `scripts/build_bbox_layered_coco.py`

| L2 family | mass share |
|---|---:|
| geometry | 0.440 |
| crowding | 0.253 |
| area | 0.190 |
| iou | 0.084 |
| semantic | 0.033 |

Emphasis tops: `pairwise_iou_mean`, `coverage_area` / `obj_density`, `area_gini` / `person_area_frac`.

## Post-hoc LOBO bounding-box attribution

Token-attribution analogue on COCO early/late:

1. Form: `{caption, bbx_text: {bbx_k: {box, text_description, category}}}`
2. Select named `f*` via RF / PO-risk VIMP (Y = rank-quantile → [0,1])
3. Locate high-`τ̂²` / high-`f*` images
4. Leave-one-box-out: `remain = boxes \ {k}` → `recalculate_bbox_features` → Δ PO rel

- Script: `scripts/run_coco_lobo_bbox_attribution.py`
- LaTeX: `docs/method/COCO_LOBO_BBox_Attribution_tables_only.tex`
- Board: `coco_lobo_bbox_attribution_board.json`
- Form sample: `coco_bbx_text_form_sample.jsonl` (full jsonl under `data/.../coco_bbx_text_form/`, gitignored)
