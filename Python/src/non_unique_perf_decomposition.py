#benchmark model performance comparison:
import shap
import numpy as np
import pandas as pd
#import os
#os.chdir('Users/heqiaoruaniCloud Drive/Documents/Github 2/Causal_Objective_Permutation_Test/Python')
#from synthetic_DGP import *
from joblib import Parallel, delayed, parallel_backend
from sklearn.model_selection import train_test_split, KFold
from sklearn.linear_model import LogisticRegression, LassoCV
from sklearn.datasets import make_classification, make_regression
from sklearn.metrics import log_loss, roc_auc_score, mean_squared_error, mean_absolute_error
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, GradientBoostingClassifier, HistGradientBoostingClassifier
from sklearn.ensemble import RandomForestRegressor, AdaBoostRegressor, GradientBoostingRegressor, HistGradientBoostingRegressor

'''
We illustrate that there's not only one unique aspect for the attribution of the drop of performance --
much more beyond than just concept drift and covariate shift -- thus the observed drop of the 
performance can not be attributed to any of the individual source of the distributional changes.
'''


N = 1000
SEED = 1250
BETA_SHIFT_GRID = np.linspace(0.05, 0.5, 10)
GAMMA_SHIFT_GRID = np.linspace(0, 0.1, 21)
EPS_SHIFT_GRID = np.linspace(0, 0.1, 21)
MSE_drop_grid = np.zeros_like(BETA_SHIFT_GRID)
MSE_drop_search_mean = np.zeros((len(GAMMA_SHIFT_GRID), len(EPS_SHIFT_GRID)))
MSE_drop_search_sd = np.zeros((len(GAMMA_SHIFT_GRID), len(EPS_SHIFT_GRID)))
RHO = 0.3
B = 20

def banded_toeplitz_corr(rho, p=20, bandwidth=6):
    from scipy.linalg import toeplitz
    corr = toeplitz(rho ** np.arange(p))
    return np.triu(np.tril(corr, bandwidth), -bandwidth)

def generate_distributionshift_dgp(n = 1000, p = 20, delta_beta = 0.0,
    delta_gamma = 0.0, rho = 0.0, epsilon_shift = 0.0, seed = 2026):
    rng = np.random.default_rng(seed)
    corr = banded_toeplitz_corr(float(rho), p = p)
    X1 = rng.multivariate_normal(mean = np.zeros(p), cov = corr, size = n)
    mean_X2 = np.zeros(p)
    mean_X2[:p] = delta_gamma
    X2 = rng.multivariate_normal(mean = mean_X2, cov = corr, size = n)
    beta1 = np.zeros(p)
    beta1[:10] = 1
    beta2 = np.zeros(p)
    beta2[:10] = 1 + delta_beta
    #Then putting tat:
    Y1 = (X1 @ beta1 + rng.normal(0, 1, n)).reshape(-1, 1)
    Y2 = (X2 @ beta2 + rng.normal(0, 1 + epsilon_shift, n)).reshape(-1, 1)
    df1 = np.hstack([X1, Y1])
    df2 = np.hstack([X2, Y2])
    return df1, df2

def parallel_MSE_drop(delta_beta = 0, 
	                gamma_ind = j,
                    rho = 0.3, epsilon_ind = k,
                    seed1 = 2021, 
                    seed2 = 2022):
    MSE_list = np.zeros(B)
    delta_gamma = GAMMA_SHIFT_GRID[j]
    delta_eps = EPS_SHIFT_GRID[k]
    for b in range(B):
        df1, df2 = generate_distributionshift_dgp(
            n = N, p = P, delta_beta = delta_beta,
            delta_gamma = delta_gamma, rho = 0.3,
            epsilon_shift = delta_eps,
            seed = seed1
        )
        df1_X = df1[:, :(df1.shape[1] - 1)]
        df1_Y = df1[:, (df1.shape[1] - 1)].ravel()
        df2_X = df2[:, :(df2.shape[1] - 1)]
        df2_Y = df2[:, (df2.shape[1] - 1)].ravel()
        df1_X_train, df1_X_val, df1_Y_train, df1_Y_val = train_test_split(
        df1_X, df1_Y, test_size = 0.3
            )
        #Then evaluate the difference of the test MSE and the validation MSE 
        #before trying to replicate the effects:
        rf_model = RandomForestRegressor(
            max_depth = 5,
            n_estimators = 100,
            min_samples_leaf = max(1, round(np.sqrt(2*N)//2)),
            n_jobs = 1,
            random_state = seed2
	    )
        rf_model.fit(df1_X, df1_Y)
	    MSE_list[b] = np.mean((rf_model.predict(df2_X) - df2_Y)**2) - 
	                      np.mean((rf_model.predict(df1_X) - df1_Y)**2)
	return j, k, np.mean(MSE_list), np.std(MSE_list)

    
for i, delta_beta in enumerate(BETA_SHIFT_GRID):
    df1, df2 = generate_distributionshift_dgp(
        n = 1000, p = 20, delta_beta = delta_beta)
    df1_X = df1[:, :(df1.shape[1] - 1)]
    df2_X = df2[:, :(df2.shape[1] - 1)]
    df1_Y = df1[:, (df1.shape[1] - 1)].ravel()
    df2_Y = df2[:, (df2.shape[1] - 1)].ravel()
    #Then combine them together:
    df1_X_train, df1_X_val, df1_Y_train, df1_Y_val = train_test_split(
        df1_X, df1_Y, test_size = 0.3
    )
    #Then evaluate the difference of the test MSE and the validation MSE 
    #before trying to replicate the effects:
    rf_model = RandomForestRegressor(
        max_depth = 5,
        n_estimators = 100,
        min_samples_leaf = max(1, round(np.sqrt(2*N)//2)),
        n_jobs = 1,
        random_state = SEED - i * i
    )
    rf_model.fit(df1_X, df1_Y)
    MSE_drop_grid[i] = np.mean((rf_model.predict(df2_X) - df2_Y)**2) - 
                      np.mean((rf_model.predict(df1_X) - df1_Y)**2)
    #Then parallel for each of the group of observation here:
    result = Parallel(n_jobs = N_JOBS, verbose = 5){
        delayed(parallel_MSE_drop)(delta_beta = delta_beta, 
            gamma_ind = j,
            rho = 0.3, epsilon_ind = k,
            seed1 = SEED + j * k,
            seed2 = SEED + 2 * j * k
        )
        for j, gamma in enumerate(GAMMA_SHIFT_GRID)
        for k, epsilon in enumerate(EPS_SHIFT_GRID)
    }
    for j, k, mean_v, sd_v in result:
        MSE_drop_search_mean[j, k] = mean_v
        MSE_drop_search_sd[j, k] = sd_v
    np.savetxt("MSE_drop_search_mean_gammaShift" + str(delta_gamma) +
        "epsShift" + str(delta_eps) + ".txt", MSE_drop_search_mean, delimiter = ",")
    np.savetxt("MSE_drop_search_std_gammaShift" + str(delta_gamma) +
        "epsShift" + str(delta_eps) + ".txt", MSE_drop_search_sd, delimiter = ",")
    np.savetxt('MSE_drop_search_cd' + str(delta_beta) + ".txt",
        MSE_drop_grid, delimiter = ",")

#Illustration & visualization for the coverage of the reproducing the other type of distribution shift through this type of shifts.



#Looking for the specifi range of the degree of mean shift + epsilon shift that correpsonding to the degree
#of the concept drift
import os
import shap
import numpy as np
import pandas as pd
from joblib import Parallel, delayed, parallel_backend
os.chdir("/Users/heqiaoruan/Documents/GitHub 2/CFPerm/CFPerm_Python/results/mse_drop_mimic/")

grid_mean = np.loadtxt("MSE_drop_search_mean_gamma0.3-0.5_eps0.3-0.5_n11_B50.txt", delimiter = ",")
grid_mean = pd.DataFrame(grid_mean, columns = np.linspace(0.3, 0.5, 11))
grid_mean.index = np.linspace(0.3, 0.5, 11)
grid_std = np.loadtxt('MSE_drop_search_std_gamma0.3-0.5_eps0.3-0.5_n11_B50.txt', delimiter = ',')
grid_std = pd.DataFrame(grid_std, columns = np.linspace(0.3, 0.5, 11))
grid_std.index = np.linspace(0.3, 0.5, 11)
MSE_orig_sd = np.loadtxt('MSE_drop_seach_cd_beta0.500.txt')
MSE_orig_mean = np.loadtxt('MSE_drop_search_std_beta0.500.txt')
########################
MSE_orig_mean = np.loadtxt('MSE_drop_search_mean_gamma_eps_only.txt', delimiter = ",")
MSE_orig_sd = np.loadtxt('MSE_drop_search_std_gamma_eps_only.txt', delimiter = ",")


#CS orig:
cs_orig = 
result = {}
B = 50
#check whether they are in the 95% confidence interval: mean +/- 1.95 * \Delta -> would be fine here already.
for gamma in GRID_GAMMA:
    for eps in EPS_GRID:
        tuple0 = (gamma, eps, beta)
        result.append({
            'gamma': gamma,
            'eps': eps,
            'beta': beta,
            'MSE_orig': MSE_orig
            'lower_CI': MSE - MSE_sd * 1.96/np.sqrt(B),
            'high_CI':  MSE + MSE_sd * 1.96/np.sqrt(B),
            'coverage': coverage_result
        })
pd.DataFrame(result).to_csv('MSE_perf_whole.csv')














