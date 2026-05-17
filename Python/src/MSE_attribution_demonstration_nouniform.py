#MSE_Attribution_Demonstration.py

import numpy as np 
import pandas as pd
from dataclasses import dataclass
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from scipy.linalg import solve, cholesky, LinAlgError
from scipy.stats import t

rng = np.random.default_rng(seed = seed)

def mvn_sample(n, mu, Sigma, rng):
    L = np.linalg.cholesky(Sigma)
    Z = rng.standard_normal((n, len(mu)))
    X = mu.reshape(1, -1) + Z @ L.T
    return X

def rank_data(values):
    order = np.argsort(-np.abs(values))
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, len(values) + 1)
    return ranks
   


#Instead of waiting passively for the solutions
#Importance of taking initiatives, encountered a block, 
#encouraged him not only to report this issue, 
#Helped him become more capable of -> instead of waiting passively for the instructions.

def conditional_sampling_gaussian(mu0, Sigma0, XS, S_idx, jitter = 1e-8):
    p = Sigma0.shape[0]
    S = np.asarray(S_idx, dtype = int)
    C = np.array([j for j in range(p) if j not in S], dtype = int)
    n = XS.shape[0]
    muS = mu0[S]
    muC = mu0[C]
    SigSS = Sigma0[np.ix_(S, S)]
    SigCC = Sigma0[np.ix_(C, C)]
    SigCS = Sigma0[np.ix_(C, S)]
    SigSC = Sigma0[np.ix_(S, C)]
    try: 
        A= solve(SigSS + jitter * np.eye(len(S)), SigCS.T).T
    except LinAlgError:
        A = SigCS @ np.linalg.pinv(SigSS + jitter * np.eye(len(S)))#(C, S)
    cond_mean = muC.reshape(1, -1) + (XS - muS.reshape(1, -1)) @ A.T
    cond_cov = SigCC - A @ SigSC
    cond_cov = cond_cov + jitter * np.eye(cond_cov.shape[0])
    #
    try:
    	L = cholesky(cond_cov)
    except LinAlgError:
    	L = cholesky(cond_cov + jitter * np.eye(cond_cov.shape[0]))
    Z = rng.standard_normal((n, len(C)))
    XC = cond_mean + Z @ L.T
    return XC, C

def conditional_sampling_

model_registry = default_model_registry(
  ntree = 150,
  ridge_alpha = 0.25,
  nthread = 1, maxit = 200, max_depth = 5,
  gamma = 0.25, eta = 0.15, mlp_hidden_size = 4,
  mlp_decay = 1e-5, mlp_max_iter = 500, mlp_trace = False,
  mlp_max_coef_reg = 10000, mlp_max_coef_clf = 10000,
  warn_xgb_labels = True, positive_class = 1
)

@dataclass
class ctype:
    kind: str


@dataclass:
class Batch:
    mu: np.ndarray
    Sigma: np.ndarray
    Y_generation: Callable[[np.ndarray, np.random.Generator], np.ndarray]#Function within this class.


#Very concise, very helpful:
def sample_batch(Batch: Batch, n, seed):
    rng = np.random.default_rng(seed)
    X = mvn_sample(n, env.mu, env.Sigma, rng)
    Y = Batch.Y_generation(X, rng)
    return X, Y

def make_model(beta1, beta2, noise_sd, p_signal, DGP = 'linear'):
    beta1 = np.asarray(beta, dtype = float).reshape(-1)
    beta2 = np.asarray(beta, dtype = float).reshape(-1)
    if DGP == 'linear':
        eps = rng.normal(0.0, noise_sd, size = X.shape[0])
        return X[:, :p_signal] @ beta + eps
    else:
        eps = t.rvs(df = 5, loc = 0, scale = 1, size = X.shape[0], rng = rng)
        return (X[:, :p_signal]**2) @ np.array(beta) + eps 


def LMGeneration(n = 1000, p = 20, covariance_matrix, beta = [1,1,1,1,1,1,1], variances = 1, cor = 0.3, eps_noi = 1, mean = 0):
    n_signal = len(beta)
    mean_vector = [mean] * n_signal
    cov_mat = np.zeros((n_signal, n_signal))
    for i in range(n_signal):
        for j in range(n_signal):
            cov_mat[i, j] = cor ** abs(i - j) * variances
    cov = np.dot(cov_mat, cov_mat.T)
    X_design = np.random.multivariate_normal(mean_vector, cov, size=n)
    #generate the response:
    Y = np.dot(X_design, beta) + np.random.normal(loc = 0, scale = eps_noi, size = n)
    X_noise = np.zeros((n, p - n_signal))
    for k in range(p - n_signal):
        X_noise[:, k] = np.random.normal(loc = 0, scale = 1, size = n)
    X_design = np.concatenate((X_design, X_noise), axis = 1)
    df_design = pd.DataFrame(np.concatenate((X_design, Y.reshape(-1, 1)), axis = 1))
    df_design.columns = ["X" + str(i) for i in range(p)] + ["Y"]
    return df_design





    




def row_replacement(X0_pool, X1_pool, switched_idx, rng):
    n, p = X1_pool.shape
    mask = np.zeros(p, dtype = bool)
    mask[switched_idx] = True
    j0 = rng.integers(0, X0_pool.shape[0], size = n)
    Xh = X0_pool[j0, :].copy()
    Xh[:, mask] = X1_pool[:, mask]
    return Xh


def hybrid_gaussian_conditional_replacement(
    env0, X1_pool, switched_idx,
    rng: np.random.Generator, jiter
):
    n, p = X1_pool.shape
    S = np.asarray(switched_idx, dtype = int)
    Xh = X1_pool.copy()
    if S.size < p:
        xS = X1_pool[:, S] if S.size > 0 else np.zeros((n, 0))
        XC, C = conditional_sampling_gaussian(
            env0.mu, env0.Sigma, xS, S, rng, jitter
        )
        if C.size > 0:
            Xh[:, C] = XC
    return Xh





#Abstract the data-generating process here:

def conditional_sampling_sequential(
	)














