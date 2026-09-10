# Named bbox concat VIMP · cross-dataset consensus

S* = `['n_objects', 'n_people', 'area_frac_person', 'area_frac_vehicle', 'region_center_area_prop', 'region_periphery_area_prop', 'quad_TL_area_prop', 'quad_TR_area_prop', 'quad_BL_area_prop', 'quad_BR_area_prop']`

## Inject recovery by dataset

| dataset | AUC | mass_bbox | recover | recovered features |
|---|---:|---:|---:|---|
| coco_time_order | 1.000 | 0.912 | 1.00 | n_objects, n_people, area_frac_person, area_frac_vehicle, region_center_area_prop, region_periphery_area_prop, quad_TL_area_prop, quad_TR_area_prop, quad_BL_area_prop, quad_BR_area_prop |
| coco_outdoor_indoor | 1.000 | 0.616 | 1.00 | n_objects, n_people, area_frac_person, area_frac_vehicle, region_center_area_prop, region_periphery_area_prop, quad_TL_area_prop, quad_TR_area_prop, quad_BL_area_prop, quad_BR_area_prop |
| coco_center_split | 1.000 | 0.927 | 1.00 | n_objects, n_people, area_frac_person, area_frac_vehicle, region_center_area_prop, region_periphery_area_prop, quad_TL_area_prop, quad_TR_area_prop, quad_BL_area_prop, quad_BR_area_prop |

## Consensus (recovered on ≥2 datasets under inject)

**10/10**: n_objects, n_people, area_frac_person, area_frac_vehicle, region_center_area_prop, region_periphery_area_prop, quad_TL_area_prop, quad_TR_area_prop, quad_BL_area_prop, quad_BR_area_prop

## Baseline top-20 named ranks (for reference)

- `coco_time_order` AUC=0.584 mass_bbox=0.030: std_cx, area_frac_vehicle, mean_cy, count_person, mean_area_frac, n_people, person_to_total_area_prop, max_area_frac
- `coco_outdoor_indoor` AUC=1.000 mass_bbox=0.219: furniture_to_total_area_prop, count_furniture, area_frac_furniture, count_kitchen, largest_is_furniture, area_frac_kitchen, count_vehicle, area_frac_sports
- `coco_center_split` AUC=1.000 mass_bbox=0.762: region_center_area_prop, region_periphery_area_prop, region_center_count_prop, region_periphery_count_prop, largest_to_total_area_prop, std_cx, mean_area_frac, max_area_frac
