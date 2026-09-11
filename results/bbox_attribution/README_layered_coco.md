# COCO multi-layer bbox features (IoU / crowding / area)

RF Domain AUC = 0.584

## L1 modality mass

| block | share |
|---|---:|
| image_clip | 0.524 |
| text_clip | 0.447 |
| bbox_layered | 0.028 |

## L2 bbox family mass

| family | share |
|---|---:|
| geometry | 0.440 |
| crowding | 0.253 |
| area | 0.190 |
| iou | 0.084 |
| semantic | 0.033 |

## L3 top-20 named

| rank | family | feature | vimp |
|---:|---|---|---:|
| 1 | geometry | `cx_std` | 0.00327 |
| 2 | geometry | `aspect_std` | 0.00275 |
| 3 | area | `area_gini` | 0.00184 |
| 4 | geometry | `w_mean` | 0.00151 |
| 5 | geometry | `aspect_mean` | 0.00145 |
| 6 | crowding | `coverage_area` | 0.00133 |
| 7 | crowding | `obj_density` | 0.00106 |
| 8 | geometry | `cx_mean` | 0.00103 |
| 9 | crowding | `densest_quad_area_prop` | 0.00093 |
| 10 | crowding | `nn_center_dist_mean` | 0.00092 |
| 11 | geometry | `cy_mean` | 0.00089 |
| 12 | iou | `pairwise_iou_mean` | 0.00088 |
| 13 | crowding | `nn_center_dist_min` | 0.00088 |
| 14 | geometry | `h_mean` | 0.00085 |
| 15 | area | `person_area_frac` | 0.00078 |
| 16 | area | `area_std` | 0.00071 |
| 17 | crowding | `periphery_vs_center_area` | 0.00068 |
| 18 | geometry | `cy_std` | 0.00068 |
| 19 | crowding | `congest_score` | 0.00067 |
| 20 | area | `vehicle_area_frac` | 0.00065 |
