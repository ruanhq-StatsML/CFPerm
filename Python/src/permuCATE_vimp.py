#PermuCATE vimp with the R-risk:
import os
import json
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.base import clone
from dataclasses import dataclass
from model_registry_class import *
from econml.dml import CausalForestDML
from abc import ABC, abstractmethod
from sklearn.pipeline import Pipeline
from scipy.stats import t as student_t
from scipy.sparse.linalg import spsolve
from scipy.sparse import diags, eye, csr_matrix
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, KFold
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


model_registry = default_model_registry(
  ntree = 150,
  ridge_alpha = 0.25,
  nthread = 1, maxit = 200, max_depth = 5,
  gamma = 0.25, eta = 0.15, mlp_hidden_size = 4,
  mlp_decay = 1e-5, mlp_max_iter = 500, mlp_trace = False,
  mlp_max_coef_reg = 10000, mlp_max_coef_clf = 10000,
  warn_xgb_labels = True, positive_class = 1
)


'''
Conditional Permutation Variable Importance for the PO-risk:
Fitting the conditional dependence models here:

'''
def permuCATE_fit_nuisance(X, Y, W, model_m, model_e, model_tau):
    X = np.asarray(X)
    Y = _as_1d(Y)
    W = _as_1d(W).astype(int)
    p = X.shape[0]
    model_outcome = model_registry[model_m]
    model_propensity_score = model_registry[model_e]
    model_tau = model_registry[model_tau]
    
    

def permuCATE(
    X: np.ndarray, Y: np.ndarray, W: np.ndarray,
    model_e, model_m,
    *, seed: int = 0, n_splits = 5, clip_e: float = 0.01,
    clip_wtilde: float = 1e-3, nomralize: bool = False
):
    X = np.asarray(X)
    Y = _as_1d(Y)
    W = _as_1d(W).astype(int)
    p = X.shape[0]




























