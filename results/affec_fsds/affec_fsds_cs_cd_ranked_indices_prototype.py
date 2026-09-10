# Covariate-shift ranking (RF Domain Classifier, high→low)
feature_indices_covariate_shift = {'eye_tracking': [6, 15, 12, 11, 4, 5, 0, 1], 'pupil': [5, 4, 13, 11, 12, 15, 17, 10], 'cursor': [1, 0, 2, 3], 'gsr_eda': [15, 32, 31, 4, 3, 35, 30, 14], 'eeg': [50, 54, 23, 2, 56, 3, 47, 15]}

# Concept-drift ranking (PO-risk path, high→low)
feature_indices_concept_drift = {'eye_tracking': [15, 9, 6, 0, 2, 1, 4, 3], 'pupil': [2, 0, 5, 4, 13, 10, 11, 15], 'cursor': [0, 1, 2, 3], 'gsr_eda': [0, 29, 2, 31, 1, 34, 14, 33], 'eeg': [62, 13, 46, 12, 52, 44, 27, 59]}
