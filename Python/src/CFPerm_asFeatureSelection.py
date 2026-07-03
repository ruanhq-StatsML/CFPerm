#Feature Selection CFPerm:

'''
We justify the procedure of leveraging causal testing procedure 
as a feature selection problem:
(1) Feature Selection for the Distribution Shift here:
(2) Feature Selection for the Distribution Shift here:

WeightedSim_{ij} = w
'''
from dataclasses import dataclass
import numpy as np
import pandas as pd
import shap
from sklearn.datasets import make_classification, make_regression
from sklearn.preprocessing import StandardScaler, FunctionTransformer, PolynomialFeatures
from sklearn.utils import check_random_state
from sklearn.model_selection import KFold
from sklearn.metrics import log_loss, roc_auc_score, mean_squared_error, mean_absolute_error
from sklearn.linear_model import LassoCV
from scipy.special import expit
from sklearn.utils import resample
from copy import deepcopy
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, GradientBoostingClassifier, HistGradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier,MLPRegressor
from sklearn.kernel_approximation import RBFSampler
from sklearn.neighbors import KNeighborsRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor, AdaBoostRegressor, GradientBoostingRegressor, HistGradientBoostingRegressor
from sklearn.svm import SVR
from joblib import Parallel, delayed
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet
from benchmark_config import MODEL_REGISTRY
from LOCO_vimp_r_risk import * 
from permuCATE_vimp import * 
from grf_vimp_causalForest import *
from utils import _as_1d


BETA_ORIGINAL = np.array([1,1,1,1,1,1,1,1])
BETA_NEW = np.array([1,1,1,1,0,0,0,0])
DIF_SEQ = np.linspace(0, 1, 21)
SIG_IND = np.array([4,5,6,7])#The degree of recovery here for these 4 features
#Then visualize the trend of the p-value procedures here:

def lm_generation(n=1000, p=20, beta=None, cor=0.3, eps_noi=1.0, mean_vector = None,
    variance_vector = None):
    if beta is None:
        beta = [1, 1, 1, 1, 1, 1, 1, 1]
    n_signal = len(beta)
    if mean_vector is None:
        mean_vector = np.zeros(n_signal, dtype=float)
    if variance_vector is None:
        variance_vector = np.ones(n_signal, dtype = float)
    base = np.fromfunction(lambda i, j: cor ** np.abs(i - j) * np.sqrt(variance_vector[i]) * np.sqrt(variance_vector[j]), (n_signal, n_signal))
    cov = base @ base.T
    X_signal = np.random.multivariate_normal(mean_vector, cov, size=n)
    Y = X_signal @ np.asarray(beta, dtype=float) + np.random.normal(0, eps_noi, size=n)
    X_noise = np.random.normal(0, 1, size=(n, p - n_signal))
    X = np.concatenate([X_signal, X_noise], axis=1)
    df = pd.DataFrame(np.concatenate([X, Y.reshape(-1, 1)], axis=1))
    df.columns = [f"X{i+1}" for i in range(p)] + ["Y"]
    return df




#Average Ranking here for the 
rejected_dict = {}
for i, dif in enumerate(DIF_SEQ):
    df1 = lm_generation(n = 1000, p = 20, beta = BETA_ORIGINAL)
    BETA_NEW = BETA_ORIGINAL
    for j in range(4, 7):
        BETA_NEW[j] = BETA_NEW[j] - dif
    df2 = lm_generation(n = 1000, p = 20, beta = BETA_NEW)
    X = np.vstack(df1[:, :])
    W = np.concatenate([np.zeros(n), np.ones(n)])
    Y = np.concatenate()
    feature_specific_pval_permu_res = cfperm_feature_pvals(
        X, Y, W,
        model_registry = MODEL_REGISTRY,
        vimp = 'permucate',
        n_perm = 100,
        n_perm_cate = 25,
        level_across_feature = 0.95,
        model_m = 'rf_regressor',
        model_e = 'logistic_classifier',
        parallel = True
    )
    feature_specific_pval_permu = feature_specific_pval_permu_res['feature_specific_pval']
    perm_rej = feature_specific_pval_permu['rejected_feature']
    feature_specific_pval_loco_res = cfperm_feature_pvals(
        X, Y, W,
        vimp = 'loco',
        model_registry = MODEL_REGISTRY,
        model_m = 'rf_regressor',
        model_e = 'logistic_classifier',
        parallel = True
    )
    feature_specific_pval_loco = feature_specific_pval_loco_res['feature_specific_pval']
    loco_rej = feature_specific_pval_loco['rejected_feature']
    feature_specific_pval_grf_res = cfperm_feature_pvals(
        X, Y, W,
        vimp = 'grf',
        model_registry = MODEL_REGISTRY,
        model_m = 'rf_regressor',
        model_e = 'logistic_classifier',
        parallel = True        
    )
    feature_specific_pval_grf = feature_specific_pval_grf_res['feature_specific_pval']
    grf_rej = feature_specific_pval_grf_res['rejected_feature']
    #Then evaluate this: 
    if feature_specific_pval_permu_res['rejected']:
        rejected_dict[dif].append(['permu', feature_specific_pval_loco_res])
    if feature_specific_pval_loco_res['rejected']:
        rejected_dict[dif].append(['loco', feature_specific_pval_loco_res])
    if feature_specific_pval_grf_res['rejected']:
        rejected_dict[dif].append(['grf', feature_specific_pval_grf_res])        
    for k in range(len(SIG_IND)):
        top_ind_rank_loco[k, i] = np.argsort(-feature_specific_pval_loco)[4+k]
        top_ind_rank_permu[k, i] = np.argsort(-feature_specific_pval_permu)[4+k]
        feature_specific_pval_grf[k, i] = np.argsort(-feature_specific_pval_grf)[4+k]
    #Then return all of them: Can you select them efficient and look at their trend here:




#The sample-specific p-value with the bootstrap is helpful here.
#rejected_dict[dif].append(['permu', feature_specific_pval_loco_res])

#feature_specific_pval_grf[k, i] = np.argsort(-feature_specific_pval_grf)[4+k]
#top_ind_rank_loco[k, i] = np.argsort(-feature_specific_pval_loco)[4+k]




'''
Youtube launch a new feature, how to design the experiment to measure the 
user's watching length?
Multicollinearity, ridge regression and the regression basics here are essential.
1000 students running, calculate the grade, the first 100 and the last 100

SELECT user_id
FROM signup table AS A
JOIN session table AS b ON A.user_id = B.user_id
GROUP BY user_id
HAVING(duration < 2) AND (b.date - a.signup_date) >= 7

hypothesis tesitng ->> power -> 
favorite google product, how to improve on top of this(A/B)

How to explain algorithm or statistics to Non-Tech?
How to resolve conflict?
n = m fair, n > m -> depend on the number here.

Calculate the median from a sample(no distributional assumption)
to use the bootstrap.
n = m fair, n > m

The variable selection path for lasso here may not be consistent here.

Regression Regularization, explain the lasso in plain English here!




y = ax1+bx2 
y = a(x1 + x2) + b(x1 - x2)

n = 1000, p = 900 ~ if works well in the training ->
median is more important than the mean

Calculate the median/mean from the bootstrap here:
for i in range(100):
    


median for the bootstrap confidence interval here -> 
youtube launch a new feature how to design the experiemtn to measure the
watching length:

regression regularization and explain lasso in plain english here
launch a new feature, how to design the experiment to measure the
watching legnth?

The groupby, window function, join and other case problems;;

-> 

What's MLE, what's logisitc regression
U(x,y)= 0/1 given x, y 
U(X,Y) -> estimate the area under the curve here, whats' the variance
MC -> median, pooled sample-median


U(x,y) = 0/1 given x and y point
How to use U(x, y) estimate the area under the curve?
The sample median and the median for the pooled sample:


X1_List <- rnorm(1000)
n <- length(X1_List)
B <- 200
median_val <- seq_len(B)
for(i in 1:B){
  X_resample = sample(X1_List, n, replace = TRUE)
}

Talking about the data intuition here with the groupby, window function
as well as the join & case in problems for others:

U(X, Y) -> 0/1 given X and Y.

Binomial Distribution, Margin of Error:
How to design survey, higher probability to get users' thinking

Search Analytics Brand _-. Correlation -> multicolineary ->
Which features to select? Stepwise Regression here?
Which features to remove? Malicious or 

Regression Regularization: explaining the LASSO in plain English
Launch a new feature, how to design the experiment to measure the 
user watching length? 
median from a sample 
How to design the experiment to measure the user watching length?
Regression Regularization: explaining the LASSO in plain English?


delta = E[Y|T=1] - phi -> E[Y|T=1] can be easily estimated
with the sample mean of the outcome among the exposed units.

delta = E[Y|T=1] - phi -> E[Y|T=1] -> easily estimated 
with the  sample-mean of the outcome among the exposed units.

Go through every detail here:

delta = E[Y|T=1] - phi -> E[Y|T=1] - phi -> easily estimated
can be

10 million frames from waymo fleet, only 500 frames is there a
child running across the street, how evaluate whether your detection
model is good at detecting this specific scenario?

Yes how would you quantify this efficiecny here?
Suppose 500 frames are leveraged, is there a child running across the street,
how to evaluate whether your detection model is good enough 
at detecting this specific scenario?
Delta = E[Y|T=1] - phi here?

same rate -> multiple items -> repeated measures here: Autocorrelation here!
Introduce the effects that the standard error is highly underestimated here!
'''















