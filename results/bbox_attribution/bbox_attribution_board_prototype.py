# Hierarchical attribution board · image / text / bbox
# L1 = modality mass share; L2 = within-block ranking; L2b = bbox subgroups
modality_mass_share_rf = {'image': 0.5252620343942326, 'text': 0.44222916773549575, 'bbox': 0.03250879786927152}
feature_ranking_rf = {'image': [120, 148, 69, 30, 418, 169, 54, 332, 301, 213, 294, 286, 215, 363, 238, 21, 16, 346, 303, 272], 'text': [330, 376, 187, 365, 194, 59, 505, 61, 360, 195, 383, 501, 126, 118, 30, 202, 186, 45, 54, 456], 'bbox': [6, 93, 107, 95, 8, 90, 97, 2, 106, 98, 108, 111, 1, 0, 100, 103, 110, 102, 94, 5]}
bbox_subgroup_mass_share_rf = {'geo': 0.2861797737714111, 'category_hist': 0.06601810191179595, 'top_boxes': 0.6478021242860321}
bbox_geo_ranking_rf = ['std_cx', 'mean_area_frac', 'total_area_frac', 'n_boxes_log1p', 'mean_cy', 'std_cy', 'mean_cx', 'max_area_frac']
modality_mass_share_mmd = {'image': 0.4723165941507731, 'text': 0.4775969298173213, 'bbox': 0.05008647603154711}
bbox_inject_recovery_rf = {'selection_AUC': 0.8664, 'P@20': 1.0, 'mass_on_bbox': np.float64(0.9984)}
bbox_inject_recovery_mmd = {'selection_AUC': 0.9335, 'P@20': 1.0, 'mass_on_bbox': np.float64(0.9583)}
