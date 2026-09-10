# Interpretable bbox attribution (named features)

Do **not** claim CLIP index importance. BBox evidence is handcrafted names.

## L1 modality mass (RF AUC=0.568)

| image_clip | text_clip | bbox_named |
|---:|---:|---:|
| 0.527 | 0.451 | 0.023 |

## BBox group mass

| group | mass |
|---|---:|
| geometry_counts | 0.759 |
| supercategory | 0.126 |
| key_class | 0.109 |
| largest_flags | 0.006 |

## Named feature ranking (top 15)

| rank | feature | vimp |
|---:|---|---:|
| 1 | `std_cx` | 0.00486 |
| 2 | `std_cy` | 0.00194 |
| 3 | `area_frac_vehicle` | 0.00153 |
| 4 | `mean_area_frac` | 0.00146 |
| 5 | `n_objects` | 0.00123 |
| 6 | `count_person` | 0.00094 |
| 7 | `largest_cy` | 0.00091 |
| 8 | `mean_cy` | 0.00090 |
| 9 | `largest_w_frac` | 0.00086 |
| 10 | `largest_cx` | 0.00083 |
| 11 | `n_people` | 0.00077 |
| 12 | `largest_h_frac` | 0.00076 |
| 13 | `coverage_union_proxy` | 0.00072 |
| 14 | `mean_cx` | 0.00068 |
| 15 | `area_frac_person` | 0.00062 |

## Inject recovery (GT=bbox_named)

- mass_on_bbox_named=0.997, AUC=1.000, P@20=1.000
