# AFFEC FSDS board · RF Domain Classifier rank

feature-indices (RF Domain Classifier rank): eye_tracking:[6, 15, 12, 11, 4, 5, 0, 1], pupil:[5, 4, 13, 11, 12, 15, 17, 10], cursor:[1, 0, 2, 3], gsr_eda:[15, 32, 31, 4, 3, 35, 30, 14], eeg:[50, 54, 23, 2, 56, 3, 47, 15]

- board order: `np.argsort(-RF_VIMP)[:8]` within each modality
- RF domain AUC: 0.999
- n=10000, p=144
- batch: W=0: run∈{0,1}; W=1: run∈{2,3}

| modality | ranked indices (high → low VIMP) | names |
|---|---|---|
| eye_tracking | [6, 15, 12, 11, 4, 5, 0, 1] | fixation_duration, validity, pupil_y_left, pupil_x_left, gaze_x_right, gaze_y_right, fixation_x, fixation_y |
| pupil | [5, 4, 13, 11, 12, 15, 17, 10] | pupil_latency, pupil_peak, eye_pos_y_raw, eye_pos_y, eye_pos_x_raw, eye_pos_velocity_y, eye_pos_range_x, eye_pos_x |
| cursor | [1, 0, 2, 3] | cursor_y, cursor_x, cursor_velocity, cursor_state |
| gsr_eda | [15, 32, 31, 4, 3, 35, 30, 14] | acc_y, gsr_feat_32, gsr_feat_31, gsr_scr_count, gsr_tonic, gsr_feat_35, gsr_feat_30, acc_x |
| eeg | [50, 54, 23, 2, 56, 3, 47, 15] | EEG_50, EEG_54, EEG_23, F3, EEG_56, F4, EEG_47, EEG_15 |
