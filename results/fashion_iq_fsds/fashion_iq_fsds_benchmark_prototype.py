# Fashion-IQ CLIP · FSDS benchmark_feature_selection (RF / MMD / PO-risk)
# Ranking = np.argsort(-VIMP)[:20] within modality; X=[img|txt], Y=last col
feature_indices_rf_domain = {'image': [381, 316, 496, 218, 124, 265, 226, 18, 379, 151, 393, 253, 492, 426, 459, 47, 108, 46, 355, 118], 'text': [326, 177, 212, 198, 209, 294, 53, 64, 35, 369, 415, 271, 119, 108, 0, 216, 295, 179, 55, 106]}
modality_vimp_share_rf_domain = {'image': 0.6969528319554447, 'text': 0.30304716804355525}
feature_indices_loco_mmd = {'image': [190, 496, 426, 84, 493, 331, 338, 133, 95, 490, 313, 47, 264, 362, 224, 507, 363, 176, 151, 20], 'text': [233, 186, 104, 211, 175, 378, 359, 201, 50, 10, 368, 271, 364, 90, 219, 375, 406, 60, 75, 471]}
modality_vimp_share_loco_mmd = {'image': 0.550107321396805, 'text': 0.44989267860290677}
feature_indices_po_risk = {'image': [39, 471, 489, 200, 499, 102, 358, 504, 195, 501, 370, 146, 71, 181, 317, 27, 345, 19, 239, 160], 'text': [315, 442, 443, 67, 117, 179, 369, 351, 111, 448, 15, 64, 115, 198, 23, 181, 57, 100, 159, 82]}
modality_vimp_share_po_risk = {'image': 0.2797336576809514, 'text': 0.7202663423180485}
