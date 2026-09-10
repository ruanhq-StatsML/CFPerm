# AFFEC FSDS · covariate-shift vs concept-drift rankings

## Covariate shift (RF Domain Classifier · `argsort(-VIMP)[:k]`)

```python
feature_indices_covariate_shift = {'eye_tracking': [6, 15, 12, 11, 4, 5, 0, 1], 'pupil': [5, 4, 13, 11, 12, 15, 17, 10], 'cursor': [1, 0, 2, 3], 'gsr_eda': [15, 32, 31, 4, 3, 35, 30, 14], 'eeg': [50, 54, 23, 2, 56, 3, 47, 15]}
```

- domain AUC: 0.999

| modality | ranked indices | names |
|---|---|---|
| eye_tracking | [6, 15, 12, 11, 4, 5, 0, 1] | fixation_duration, validity, pupil_y_left, pupil_x_left, gaze_x_right, gaze_y_right, fixation_x, fixation_y |
| pupil | [5, 4, 13, 11, 12, 15, 17, 10] | pupil_latency, pupil_peak, eye_pos_y_raw, eye_pos_y, eye_pos_x_raw, eye_pos_velocity_y, eye_pos_range_x, eye_pos_x |
| cursor | [1, 0, 2, 3] | cursor_y, cursor_x, cursor_velocity, cursor_state |
| gsr_eda | [15, 32, 31, 4, 3, 35, 30, 14] | acc_y, gsr_feat_32, gsr_feat_31, gsr_scr_count, gsr_tonic, gsr_feat_35, gsr_feat_30, acc_x |
| eeg | [50, 54, 23, 2, 56, 3, 47, 15] | EEG_50, EEG_54, EEG_23, F3, EEG_56, F4, EEG_47, EEG_15 |

## Concept drift (PO-risk RF · `argsort(-VIMP)[:k]`)

```python
feature_indices_concept_drift = {'eye_tracking': [15, 9, 6, 0, 2, 1, 4, 3], 'pupil': [2, 0, 5, 4, 13, 10, 11, 15], 'cursor': [0, 1, 2, 3], 'gsr_eda': [0, 29, 2, 31, 1, 34, 14, 33], 'eeg': [62, 13, 46, 12, 52, 44, 27, 59]}
```

- PO-risk: 0.0000

| modality | ranked indices | names |
|---|---|---|
| eye_tracking | [15, 9, 6, 0, 2, 1, 4, 3] | validity, eye_open_left, fixation_duration, fixation_x, gaze_x_left, fixation_y, gaze_x_right, gaze_y_left |
| pupil | [2, 0, 5, 4, 13, 10, 11, 15] | pupil_diameter_filt, pupil_diameter, pupil_latency, pupil_peak, eye_pos_y_raw, eye_pos_x, eye_pos_y, eye_pos_velocity_y |
| cursor | [0, 1, 2, 3] | cursor_x, cursor_y, cursor_velocity, cursor_state |
| gsr_eda | [0, 29, 2, 31, 1, 34, 14, 33] | gsr_raw, gsr_feat_29, gsr_phasic, gsr_feat_31, gsr_filtered, gsr_feat_34, acc_x, gsr_feat_33 |
| eeg | [62, 13, 46, 12, 52, 44, 27, 59] | EEG_62, EEG_13, EEG_46, EEG_12, EEG_52, EEG_44, EEG_27, EEG_59 |
