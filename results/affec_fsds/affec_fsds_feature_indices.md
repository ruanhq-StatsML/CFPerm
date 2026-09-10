# AFFEC FSDS · selected feature indices

feature-indices: eye_tracking:[6, 9, 10, 15], pupil:[0, 2, 4, 5, 10, 11, 12, 13, 16, 17, 20], cursor:[1, 3], gsr_eda:[0, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 28, 31, 32, 33, 34, 35], eeg:[11, 19, 37, 46, 47, 52]

- n=10000, p=144, n_perm=40
- batch: W=0: run∈{0,1}; W=1: run∈{2,3}

| modality | indices | names |
|---|---|---|
| eye_tracking | [6, 9, 10, 15] | fixation_duration, eye_open_left, eye_open_right, validity |
| pupil | [0, 2, 4, 5, 10, 11, 12, 13, 16, 17, 20] | pupil_diameter, pupil_diameter_filt, pupil_peak, pupil_latency, eye_pos_x, eye_pos_y, eye_pos_x_raw, eye_pos_y_raw, eye_pos_dispersion, eye_pos_range_x, pupil_validity |
| cursor | [1, 3] | cursor_y, cursor_state |
| gsr_eda | [0, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 28, 31, 32, 33, 34, 35] | gsr_raw, gsr_tonic, gsr_scr_count, gsr_peaks, gsr_slope, gsr_std, gsr_range, body_temp, temp_mean, temp_std, temp_slope, acc_x, acc_y, acc_z, acc_magnitude, acc_slope, gsr_feat_20, gsr_feat_21, gsr_feat_22, gsr_feat_23, gsr_feat_24, gsr_feat_25, gsr_feat_28, gsr_feat_31, gsr_feat_32, gsr_feat_33, gsr_feat_34, gsr_feat_35 |
| eeg | [11, 19, 37, 46, 47, 52] | EEG_11, EEG_19, EEG_37, EEG_46, EEG_47, EEG_52 |
