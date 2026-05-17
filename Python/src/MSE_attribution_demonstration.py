#Causal Attribution Demonstration:
import numpy as np
import pandas as pd
from dataclasses import dataclass
from scipy.stats import spearmanr, kendalltau
from sklearn.model_selection import train_test_split
from typing import Callable, Dict, List, Optional, Tuple
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, GradientBoostingClassifier, HistGradientBoostingClassifier
import xgboost
import lightgbm
from model_registry_class import *

model_registry = default_model_registry(
  ntree = 150,
  ridge_alpha = 0.25,
  nthread = 1, maxit = 200, max_depth = 5,
  gamma = 0.25, eta = 0.15, mlp_hidden_size = 4,
  mlp_decay = 1e-5, mlp_max_iter = 500, mlp_trace = False,
  mlp_max_coef_reg = 10000, mlp_max_coef_clf = 10000,
  warn_xgb_labels = True, positive_class = 1
)



def block_corr(p, p_signal, rho):
    Cov_Sigma = np.eye(p)
    if p_signal > 1:
        SigS = np.full((p_signal, p_signal), rho_signal, dtype = float)
        np.fill_diagonal(SigS, 1.0)
        Cov_Sigma[:p_signal, :p_signal] = SigS
        Cov_Sigma = np.tril(Cov_sigma, k = 3)
        Cov_Sigma = np.triu(Cov_sigma, k = -3)
    return Cov_Sigma


def conditional_gaussian_sampling(mu0, Sigma0, XS, S_idx, clip_val = 1e-5, rng):
    p = Sigma0.shape[0]
    S = np.array(S_idx, dtype = int)
    C = np.setdiff1d(np.arange(p), S)
    muS = mu0[S]
    muC = mu0[C]
    '''
    Fetch the off-diagonal and diagonal component of the Covariance Matrix, conduct the conditional sampling on the multivariate normal distribution.
    '''
    SigSS = Sigma0[muS, muS]
    SigCC = Sigma0[muC, muC]
    SigSC = Sigma0[muC, muS]
    SigCS = Sigma0[muS, muC]
    #Conditional Sampling on the Multivariate Conditional Distribution:
    invSigSS = np.linalg.inv(SigSS + clip_val * np.eye(SigSS.shape[0]))
    A = SigCS @ invSigSS
    #Conditional Mean and Conditional Variance for the Multivariate Normal Distribution:
    cond_mean = (XS - muS.reshape(1, -1)) @ A.T + muC.reshape(1, -1)
    cond_cov = SigCC - SigCS @ invSigSS @ SigSC
    cond_cov = cond_cov + clip_val * np.eye(cond_cov.shape[0])
    L = np.linalg.cholesky(cond_cov)
    Z = rng.standard_normal((XS.shape[0], len(C)))
    return cond_mean + Z @ L.T, C




'''
The orders of the rankings of when each feature is included?
Inputs include:
- order: order of when each variable is included in the model:
- 

Return:
- delta_MSE_path
- average ranking:
'''

def sigma_banded(p, rho, n_banded):
    Sigma = np.zeros((p,p))
    Sigma_rho = np.array([[rho**(i-j) for i in range(p)] for j in range(p)])
    Sigma_rho = np.tril(Sigma_rho, k =  n_banded)
    Sigma_rho = np.triu(Sigma_rho, k = -n_banded)
    return Sigma_rho

def sigma2(p, rho, n_banded):
    

#MSE-order-reank
def mse_order_rerank(model, df, B, order, n_banded = 3, seed = 2026):
    p = df.shape[1] - 1
    n = df.shape[0]
    #by default the mu0 and mu1 are 0 vectors:
    rng = np.random.default_rng(seed = seed)
    delta_MSE_path = np.zeros((p, B))
    order_path = np.zeros((p, B))
    Y1 = np.array(df[:, p])
    for i in range(B):
        random_shuffle_path = np.random.permutation(np.arange(p))
        cur_path = [random_shuffle_path[0]]
        #####
        #Conditional Sampling on top of the existing data:
        delta_MSE_path[0, i]
        for k in range(1, p):
            df1_X = df1[:, np.array(cur_path)]
            mu0 = np.zeros()
            Sigma0 = sigma_banded(p, rho, n_banded = 3)
            X_new = conditional_gaussian_sampling(mu0)
            #each of the time, randomly selecting another indice here:
            X_train, X_test, Y_train, Y_test = train_test_split(X_new, Y, test_size = 0.4)
            model0 = RandomForestRegressor().fit(X_train, df1)
            MSE = np.mean((model0.predict(X_new) - Y_test)**2)
            delta_MSE_path[k, i] = 


















'''
Benchmark Methods: SGShift, 
'''


#Implements a penalized regression here:




















































