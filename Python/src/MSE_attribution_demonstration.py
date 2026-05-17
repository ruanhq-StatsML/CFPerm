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
        #why not 

def block_corr(p, p_signal, rho):
    Cov_Sigma = np.eye(p)
    if p_signal > 1:
        SigS = np.full((p_signal, p_signal), rho_signal, dtype = float)
        np.fill_diagonal(SigS, 1.0)
        Cov_Sigma[:p_signal, :p_signal] = SigS
        Cov_Sigma = np.tril(Cov_sigma, k = 3)
        Cov_Sigma = np.triu(Cov_sigma, k = -3)
    return Cov_Sigma


def mvn_sample(n, mu, Sigma, rng):
    L = np.linalg.cholesky(Sigma)
    Z = rng.standard_normal((n, len(mu)))
    return mu.reshape(1, -1) + Z @ L.T

rng = np.random


def conditional_gaussian_sampling(mu0, Sigma0, XS, S_idx, clip_val = 1e-5, rng):

#Conditional gaussian sampling here:
p = Sigma0.shape[0]
S = np.array(S_idx, dtype = int)
C = np.setdiff1d(np.arange(p), S)
muS = mu0[S]
muC = mu0[C]
'''
Fetch the off-diagonal and diagonal component of the Covariance Matrix, conduct the conditional sampling on the multivariate normal distribution.
'''
SigSS = Sigma0[S, S]
SigCC = Sigma0[C, C]
SigSC = Sigma0[S, C]#
SigCS = Sigma0[C, S]
#Conditional Sampling on the Multivariate Conditional Distribution:
SigSS_inv = np.linalg.inv(SigSS + clip_val * np.eye(SigSS.shape[0])) #(n_s, n_s)
A = SigmaCS @ SigSS_inv #(n_c, n_s)
#Conditional Mean and Conditional Variance for the Multivariate Normal Distribution:
cond_mean = (XS - muS.reshape(1, -1)) @ A.T + muC.reshape(1, -1) #(n_c, 1)
cond_cov = SigCC - SigCS @ SigSS_inv @ SigSC
cond_cov = cond_cov + clip_val * np.eye(cond_cov.shape[0]) #(n_c, n_c)
L = np.linalg.cholesky(cond_cov)
Z = rng.standard_normal((XS.shape[0], len(C)))
return cond_mean + Z @ L.T, C

import numpy as np
rng = np.random.default_rng(seed = 2026)
Sigma = np.random.random((25, 25)) + np.eye(25) * 12
mu = np.zeros(25)
L = np.linalg.cholesky(Sigma)
Z = rng.standard_normal((200, 25))
X = mu.reshape(1, -1) + Z @ L.T
idx_S = np.random.choice(np.arange(25), 10)
idx_nS = np.setdiff1d(np.arange(25), idx_S)
next_ind = np.random.choice(idx_nS, 1).item()
SigSS = Sigma[np.ix_(idx_S, idx_S)]
SigSC = Sigma[np.ix_(idx_S, idx_nS)]
SigCS = Sigma[np.ix_(idx_nS, idx_S)]
muS = mu[idx_S]
XS = np.array(np.arange(len(idx_S)))
A = SigCS @ np.linalg.inv(SigSS + np.eye(SigSS.shape[0]) * clip_val) @ SigSC
cond_mean = (XS.reshape(1,-1) - muS.reshape(1, -1))@A.T + muC.reshape(1, -1)#(n_c, 1)
L = np.linalg.cholesky(cond_cov)
Z = rng.standard_normal((XS.shape[0], len(idx_nS)))
#Remember this! cond_mean + Z@L.T


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




#Unsupervised Contrastive Learning here:
#model0 = RandomForestRegressor().fit(X_train, df1)

def load_llm2vec_model():
    model_name = "McGill-NLP/LLM2Vec-Meta-Llama-3-8B-Instruct-mntp-unsup-simcse"
    # specify the original foundational model, it's a LoRA adapter with 4-bit model.
    base_model_name = "McGill-NLP/LLM2Vec-Meta-Llama-3-8B-Instruct-mntp"
    l2v_model = LLM2Vec.from_pretrained(
        base_model_name,
        peft_model_name_or_path = model_name,
        device_map = 'cuda' if torch.cuda.is_available() else 'cpu',
        torch_dtype = torch.bfloat16
    )
    instruction = 'Represent the following Sentence for Similarity Checking:'
    return l2v_model, instruction













'''
Benchmark Methods: SGShift, 
'''


#Implements a penalized regression here:














def train_model(n_train):
    Xtr, ytr





#sample conditional on S and then generate the next variable:
rng = np.random.default_rng(seed = 2026)
order = rng.permutation(p)
for m in range(n_orderings):
    order = rng.permutation(p)
    v = np.zeros(p + 1, dtype = float)
    for k in range(p + 1):
        S = order[:k]
        Xh = conditional_gaussian_sampling(S, rng)




def equicorr_block_cov(p_total = 25, p_signal = 10, rho_signal = 15):
    Sigma = np.eye(p_total)
    if p_signal > 1:
        SigS = np.full((p_signal, p_signal), rho_signal, dtype = float)
        np.fill_diagonal((SigS, 1.0))
        Sigma[:p_signal, :p_signal] = SigS
    return Sigma





#RBFSampler procedures and cond_cov = SigCC - A @ SigSC
import numba

@njit(parallel = True, nopython = True)
def gaussian_conditional_sampling(mu0, Sigma0, XS, S_idx, rng, jitter = 1e-5):
    A = np.linalg.solve(SigSS + jitter * np.eye(len), SigCS.T).T
    cond_cov = SigCC - A @ SigSC
    cond_cov += np.eye(cond_cov.shape[0]) * jitter
    n = Xs.shape[0]
    rng_vals = rng.standard_normal((n, len(C)))
    cond_mean = muC.reshape(1, -1) + (xS - muS.reshape(1,-1)) #(n, |C|)
    noise = rng_vals @ L.T#(lower = True)
    XC = cond_mean + noise
    return XC, C



invSigSS = np.linalg.inv(SigSS + jitter * np.eye(SigSS.shape[0]))
A = SigCS @ invSigSS

#develop the habit of taking initiative to own his own part.
#A = SigCS @ invSigSS
A = SigSC








from dataclass import dataclass

def topk_jaccard(a: np.ndarray, b: np.ndarray, k):
    






















