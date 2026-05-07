#PermuCATE vimp with the R-risk:
import os
import json
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.base import clone
from dataclasses import dataclass
from econml.dml import CausalForestDML
from abc import ABC, abstractmethod
from sklearn.pipeline import Pipeline
from scipy.stats import t as student_t
from scipy.sparse.linalg import spsolve
from scipy.sparse import diags, eye, csr_matrix
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, KFold, train_test_split
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier, GradientBoostingRegressor, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
'''
Utils functions:
- _as_1d(X):           return the 1d array - from shape (p, ) to (p)
- _safe_normalize(v):  return the normalized version of vector v: v/np.sum(v)
- make_folds:          return the n_folds indexes to store.
'''
def _as_1d(X):
    return np.asarray(X).reshape(-1)

def _safe_normalize(v):
    return v/np.sum(v)

def make_folds(n, n_folds = 5, seed = 1):
    rng = np.random.default_rng(seed)
    idx = np.arange(n)
    rng.shuffle(idx)
    return np.array_split(idx, n_folds)


model_factory = ModelRegistry(
  ntree = 150,
  ridge_alpha = 0.25,
  nthread = 1, maxit = 200, max_depth = 5,
  gamma = 0.25, eta = 0.15, mlp_hidden_size = 4,
  mlp_decay = 1e-5, mlp_max_iter = 500, mlp_trace = False,
  mlp_max_coef_reg = 10000, mlp_max_coef_clf = 10000,
  warn_xgb_labels = True, positive_class = 1
)
model_registry = model_factory.as_r_style_dict()

'''
Calculate the PO-risk:
E[(φ(Z) - τ̂(X))^2]
φ(Z) = (y - mu_{T})(T-pi(X))/(pi(X)(1-pi(X))) + (mu_1(X) - mu_0(X))
'''
def po_risk(tau, Y, T, mu0, mu1, pi, clip = 1e-3):
    pi_clip = np.clip(pi, clip, 1 - clip)
    pseudo = (Y - np.where(T == 1, mu1, mu0)) * (T - pi_clip)/(pi_clip * (1 - pi_clip)) + (mu1 - mu0)
    return np.mean((pseudo - tau) ** 2)


'''
Infer the type of each feature in X:
-----------
X:               Array of Shape (n_exist + n_new, p)

'''
def infer_feature_type(X):
    p = X.shape[1]
    feature_types = ["continuous"] * p
    for j in range(p):
        n_level = np.unique(X[:, j]).size
        if n_level == 2:
            feature_types[j] = 'binary'
        else if n_level < 5:
            feature_types[j] = 'categorical'
    return feature_types


'''
Fit the nuisance models: outcome model, propensity score model and the pseudo-outcome model
Parameters:
-----------
X:                    Array of Shape (n_exist + n_new, p)
Y:                    Array of Shape (n_exist + n_new, ) Response
W:                    Array of Shape (n_exist + n_new, ) Treatment, Currently only Support Binary
model_m(str):         Outcome Model
model_e(str):         Propensity Score Model
model_tau(str):       Pseudo Outcome Model(pseudo_outcome ~ X)
model_nu(str):        Conditional fitting model: one for each of the feature
model_nu_bi(str):     Conditional fitting model: binary classification model
model_nu_multi(str):  Conditional fitting model: mutlinomial classifaiction model
seed(int):            Random Seed
'''
def permuCATE_fit_nuisance(X, Y, W, 
    model_m = 'rf_regressor', 
    model_e = 'logistic_classifier',
    model_tau = 'rf_regressor', 
    model_nu = 'rf_regressor',
    model_nu_bi = 'logistic_classifier', 
    model_nu_multi = 'multinomial_classifier', 
    clip_e = 1e-2,
    seed = 2026):
    X = np.asarray(X)
    Y = _as_1d(Y)
    W = _as_1d(W).astype(int)
    n, p = X.shape
    #Specify the outcome model, propensity score model and the pseudo-outcome model
    model_outcome0 = model_registry[model_m]
    model_outcome1 = model_registry[model_m]
    model_propensity_score = model_registry[model_e]
    model_tau_spec = model_registry[model_tau]
    X_ctr = X[W == 0, ]
    X_trt = X[W == 1, ]
    Y_ctr = Y[W == 0]
    Y_trt = Y[W == 1]
    mu0_fit = model_outcome0['fit'](
        X_ctr, Y_ctr, seed = seed
    )
    mu1_fit = model_outcome1['fit'](
        X_trt, Y_trt, seed = seed
    )
    #mu0_est and mu1_est:
    mu0_est = model_outcome0['predict'](
        mu0_fit, X
    )
    mu1_est = model_outcome1['predict'](
        mu1_fit, X
    )
    #Propensity Score Models:
    model_e_fit = model_propensity_score['fit'](
        X, W, seed = seed + 2
    )
    pi_est = model_propensity_score['predict'](
        model_e_fit, X
    )
    pi_est = np.clip(pi_est, clip_e, 1 - clip_e)
    #DR pseudo outcome:
    pseudo_outcome = mu1_est - mu0_est + (Y - np.where(W == 1, mu1_est, mu0_est)) * (W - pi_est)/(pi_est * (1.0 - pi_est))
    model_tau_fit = model_tau_spec['fit'](
        X, pseudo_outcome, seed = seed + 3
    )
    tau_est = model_tau_spec['predict'](
        model_tau_fit, X
    )
    feature_types = infer_feature_type(X_train)#continuous, binary, categorical.
    nu_model = {} 
    for j in range(p):
        if feature_types[j] == 'continuous':
            model_nu_spec = model_registry[model_nu]
        else if feature_types[j] == 'binary':
            model_nu_spec = model_registry[model_nu_bi]
        else:
            model_nu_spec = model_registry[model_nu_multi]
        X_minus_j = np.delete(X, j, axis = 1)
        X_j = X[:, j]
        model_nu_fit = model_nu_spec['fit'](
            X_minus_j, X_j, seed = seed + j + 100
        )
        nu_model[j] = {
            'model_nu': model_nu_spec,
            'model_nu_fit': model_nu_fit
        }
    output = {
      #Model Specifications:
      'model_outcome0': model_outcome0,
      'model_outcome1': model_outcome1,
      'model_propensity_score': model_propensity_score,
      'model_tau': model_tau_spec,
      'model_nu': model_nu_spec,
      #Fitted Objects:
      'mu0_fit': mu0_fit,
      'mu1_fit': mu1_fit,
      'model_e_fit': model_e_fit,
      'model_tau_fit': model_tau_fit,
      #Estimation on the Training Set D_{train}:
      'mu0_est': mu0_est,
      'mu1_est': mu1_est,
      'pi_est':  pi_est,
      'tau_est': tau_est,
      'pseudo_outcome': pseudo_outcome,
      #Feature Level Nuisance Models
      'nu_model': nu_model
    }
    return output


'''
Fit the nuisance models: outcome model, propensity score model and the pseudo-outcome model
Parameters:
-----------
X:                    Array of Shape (n_exist + n_new, p)
Y:                    Array of Shape (n_exist + n_new, ) Response
W:                    Array of Shape (n_exist + n_new, ) Treatment, Currently only Support Binary
model_m(str):         Outcome Model
model_e(str):         Propensity Score Model
model_tau(str):       Pseudo Outcome Model(pseudo_outcome ~ X)
model_nu(str):        Conditional fitting model: one for each of the feature
model_nu_bi(str):     Conditional fitting model: binary classification model
model_nu_multi(str):  Conditional fitting model: mutlinomial classifaiction model
n_perm:               Number of permutation for the residual of the conditional fitted model
test_size:            Proportion of the number of observations in the test set
seed(int):            Random Seed(by default 2026)
clip_e:               Clipping value for the propensity score.
normalize:            Whether to normalize the variable importance or not.
'''
def permuCATE_vimp(
    X: np.ndarray, Y: np.ndarray, W: np.ndarray,
    model_m = 'rf_regressor', model_e = 'logistic_classifier', model_tau = 'rf_regressor', 
    model_nu = 'rf_regressor', model_nu_bi = 'logistic_classifier', model_nu_multi = 'multinomial_classifier',
    n_perm: int = 100, test_size = 0.5, *,
    seed: int = 0,  clip_e: float = 0.01, normalize: bool = False
):
    X = np.asarray(X)
    Y = _as_1d(Y)
    W = _as_1d(W).astype(int)
    #split into training set and test set:
    X_train, X_test, W_train, W_test, Y_train, Y_test = train_test_split(
        X, W, Y, test_size = 0.5, stratify = W
    )
    n, p = X.shape
    #specify the feature types:
    rng = np.random.default_rng(seed)
    outputs = permuCATE_fit_nuisance(X_train, Y_train, W_train, model_m, model_e, model_tau, model_nu, seed = seed)
    vimp_permucate = np.zeros(p, dtype = float)
    all_scores = []
    #The estimation of the parameters in the test set:
    tau_est = outputs['model_tau']['predict'](
        outputs['model_tau_fit'], X_new = X_test
    )
    mu0_test = outputs['model_outcome0']['predict'](
        outputs['mu0_fit'], X_new = X_test
    )
    mu1_test = outputs['model_outcome1']['predict'](
       outputs['mu1_fit'], X_new = X_test
    )
    pi_test = outputs['model_propensity_score']['predict'](
        outputs['model_e_fit'], X_new = X_test
    )
    pi_test = np.clip(pi_test, clip_e, 1 - clip_e)
    po_risk_original = po_risk(tau_est, Y_test, W_test, mu0_test, mu1_test, pi_test)
    for j in range(p):
        X_minus_j_test = np.delete(X_test, j, axis = 1)
        nu_spec = outputs['nu_model'][j]['model_nu']
        nu_fit = outputs['nu_model'][j]['model_nu_fit']
        nu_hat = nu_spec['predict'](
            nu_fit, X_minus_j_test
        )
        r_j = (X_test[:, j] - nu_hat).ravel()
        psi_k = np.zeros(n_perm, dtype = float)
        for k in range(n_perm):
            r_shuffle = rng.permutation(r_j)
            X_perm = X_test.copy()
            X_perm[:, j] = nu_hat + r_shuffle
            tau_perm = outputs['model_tau']['predict'](
                outputs['model_tau_fit'], X_perm
            )
            po_perm = po_risk(tau_perm, Y_test, W_test, mu0_test, mu1_test, pi_test)
            psi_k[k] = po_perm - po_risk_original
            print(k)
        vimp_permucate[j] = np.mean(psi_k)
        all_scores.append(psi_k)
    if normalize:
        return vimp_permucate/np.sum(vimp_permucate)
    else:
        return vimp_permucate



#test case:
model_factory = ModelRegistry(
ntree = 150, ridge_alpha = 0.25,
nthread = 1, maxit = 500,
max_depth = 5, gamma = 0.25,
eta = 0.15, mlp_hidden_size = 4,
mlp_decay = 1e-4, mlp_max_iter = 500)
seed = 2026
model_registry = model_factory.as_r_style_dict()
X = np.random.random((100, 20))
Y = np.random.random((100, ))
W = np.random.choice((0,1), 100)
model_m = 'rf_regressor'
model_e = 'rf_classifier'
model_tau = 'rf_regressor'
model_nu = 'rf_regressor'

permuCATE_vimp(X, Y, W, model_m = 'rf_regressor', model_e = 'rf_classifier',
    model_tau = 'rf_regressor', model_nu = 'rf_regressor', n_perm = 25)
























