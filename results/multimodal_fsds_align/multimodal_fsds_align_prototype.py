# Multimodal FSDS alignment · RF / coord-MMD / PO-risk

# --- fashion_iq · train vs test ---
feature_indices_fashion_iq_rf_domain = {'image': [381, 316, 496, 218, 124, 265, 226, 18, 379, 151, 393, 253, 492, 426, 459, 47, 108, 46, 355, 118], 'text': [326, 177, 212, 198, 209, 294, 53, 64, 35, 369, 415, 271, 119, 108, 0, 216, 295, 179, 55, 106]}
modality_vimp_share_fashion_iq_rf_domain = {'image': 0.6969528319554447, 'text': 0.30304716804355525}
feature_indices_fashion_iq_coord_mmd = {'image': [190, 496, 426, 84, 493, 331, 338, 133, 95, 490, 313, 47, 264, 362, 224, 507, 363, 176, 151, 20], 'text': [233, 186, 104, 211, 175, 378, 359, 201, 50, 10, 368, 271, 364, 90, 219, 375, 406, 60, 75, 471]}
modality_vimp_share_fashion_iq_coord_mmd = {'image': 0.550107321396805, 'text': 0.44989267860290677}
feature_indices_fashion_iq_po_risk = {'image': [39, 471, 489, 200, 499, 102, 501, 504, 195, 239, 370, 27, 146, 181, 358, 345, 71, 317, 19, 465], 'text': [315, 442, 443, 67, 117, 179, 369, 351, 111, 448, 15, 115, 64, 198, 23, 181, 57, 100, 159, 82]}
modality_vimp_share_fashion_iq_po_risk = {'image': 0.2783185893121341, 'text': 0.7216814106868659}

# --- indiana_cxr · Frontal vs Lateral ---
feature_indices_indiana_cxr_rf_domain = {'image': [330, 235, 111, 335, 145, 117, 82, 293, 494, 154, 105, 278, 115, 281, 410, 298, 68, 77, 233, 0], 'text': [70, 49, 182, 285, 468, 77, 1, 395, 325, 250, 188, 322, 295, 90, 336, 138, 230, 449, 317, 301]}
modality_vimp_share_indiana_cxr_rf_domain = {'image': 0.9998551756367472, 'text': 0.00014482436225266647}
feature_indices_indiana_cxr_coord_mmd = {'image': [105, 330, 117, 278, 145, 410, 335, 235, 82, 293, 154, 77, 494, 327, 171, 68, 21, 0, 100, 111], 'text': [425, 424, 180, 492, 131, 171, 440, 99, 316, 346, 412, 191, 427, 168, 399, 485, 108, 446, 201, 31]}
modality_vimp_share_indiana_cxr_coord_mmd = {'image': 0.9848043343542222, 'text': 0.015195665645767492}
feature_indices_indiana_cxr_po_risk = {'image': [133, 4, 38, 28, 146, 446, 30, 301, 87, 151, 461, 112, 170, 77, 248, 93, 368, 291, 145, 437], 'text': [238, 424, 429, 69, 224, 291, 13, 485, 316, 483, 434, 222, 433, 102, 457, 237, 193, 128, 380, 376]}
modality_vimp_share_indiana_cxr_po_risk = {'image': 0.7210521340913945, 'text': 0.2789478659076055}

# --- coco_outdoor_indoor · outdoor vs indoor caption keywords ---
feature_indices_coco_outdoor_indoor_rf_domain = {'image': [149, 325, 484, 272, 447, 264, 24, 8, 90, 185, 401, 82, 326, 130, 53, 101, 75, 207, 510, 108], 'text': [149, 305, 325, 79, 24, 272, 462, 438, 484, 90, 99, 158, 491, 458, 413, 429, 53, 233, 218, 320]}
modality_vimp_share_coco_outdoor_indoor_rf_domain = {'image': 0.3235247474906695, 'text': 0.6764752525083304}
feature_indices_coco_outdoor_indoor_coord_mmd = {'image': [149, 325, 484, 272, 447, 24, 53, 264, 75, 130, 8, 427, 99, 101, 90, 59, 256, 210, 37, 123], 'text': [149, 305, 24, 325, 484, 79, 270, 99, 272, 130, 365, 158, 413, 53, 320, 491, 384, 429, 203, 462]}
modality_vimp_share_coco_outdoor_indoor_coord_mmd = {'image': 0.4209598107952401, 'text': 0.5790401892047393}
feature_indices_coco_outdoor_indoor_po_risk = {'image': [212, 179, 356, 20, 60, 66, 446, 248, 318, 141, 74, 467, 417, 445, 312, 437, 25, 335, 92, 497], 'text': [439, 8, 83, 128, 28, 260, 238, 293, 248, 504, 388, 243, 168, 124, 254, 24, 492, 441, 281, 236]}
modality_vimp_share_coco_outdoor_indoor_po_risk = {'image': 0.3802457061188839, 'text': 0.619754293880116}

# --- coco_time_order · early vs late image_id (order) ---
feature_indices_coco_time_order_rf_domain = {'image': [30, 303, 346, 332, 338, 238, 215, 169, 148, 480, 69, 132, 112, 401, 277, 164, 16, 416, 301, 101], 'text': [376, 330, 304, 338, 213, 117, 383, 79, 482, 156, 357, 480, 30, 106, 273, 133, 219, 194, 469, 303]}
modality_vimp_share_coco_time_order_rf_domain = {'image': 0.5337179169259694, 'text': 0.4662820830730305}
feature_indices_coco_time_order_coord_mmd = {'image': [442, 83, 165, 122, 329, 181, 123, 247, 222, 55, 22, 483, 469, 143, 400, 437, 493, 245, 81, 5], 'text': [1, 346, 338, 331, 39, 116, 442, 502, 299, 25, 505, 406, 182, 356, 22, 474, 408, 506, 441, 251]}
modality_vimp_share_coco_time_order_coord_mmd = {'image': 0.4972206229653164, 'text': 0.5027793770343063}
feature_indices_coco_time_order_po_risk = {'image': [355, 189, 191, 417, 373, 61, 77, 262, 256, 138, 392, 53, 190, 348, 52, 511, 203, 469, 108, 3], 'text': [436, 42, 264, 448, 353, 17, 58, 116, 153, 210, 360, 194, 11, 325, 478, 506, 450, 287, 388, 226]}
modality_vimp_share_coco_time_order_po_risk = {'image': 0.4210315412691836, 'text': 0.5789684587298164}

# --- microscopy_clip · short vs long caption ---
feature_indices_microscopy_clip_rf_domain = {'image': [454, 287, 298, 314, 284, 117, 426, 306, 28, 265, 249, 328, 141, 68, 76, 285, 329, 414, 442, 157], 'text': [287, 420, 292, 142, 342, 27, 349, 417, 52, 247, 190, 123, 106, 507, 14, 20, 492, 150, 235, 306]}
modality_vimp_share_microscopy_clip_rf_domain = {'image': 0.3248889000106188, 'text': 0.6751110999883813}
feature_indices_microscopy_clip_coord_mmd = {'image': [447, 431, 299, 102, 248, 212, 77, 508, 267, 216, 124, 90, 335, 107, 438, 69, 493, 423, 144, 476], 'text': [20, 27, 506, 54, 142, 188, 269, 243, 260, 247, 287, 325, 248, 493, 303, 399, 108, 189, 93, 374]}
modality_vimp_share_microscopy_clip_coord_mmd = {'image': 0.5522888894606603, 'text': 0.4477111105390633}
feature_indices_microscopy_clip_po_risk = {'image': [459, 404, 440, 225, 33, 32, 266, 15, 375, 402, 18, 154, 377, 386, 231, 406, 309, 112, 265, 65], 'text': [129, 154, 134, 479, 435, 3, 19, 91, 356, 145, 390, 457, 86, 473, 241, 1, 505, 149, 72, 322]}
modality_vimp_share_microscopy_clip_po_risk = {'image': 0.4374922525558728, 'text': 0.562507747443127}
