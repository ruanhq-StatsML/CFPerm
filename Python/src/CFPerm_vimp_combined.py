#CFPerm VIMP:

@dataclass
class CFPermResult:
    feature_names: list,
    imp_observed: pd.Series
    p_values: pd.Series
    q1_upper: pd.Series
    threshold: float
    rejected_feature_indice: list
    perm_imp: np.ndarray
    method: str

model_factory = ModelRegistry(
ntree = 150, ridge_alpha = 0.25,
nthread = 1, maxit = 500,
max_depth = 5, gamma = 0.25,
eta = 0.15, mlp_hidden_size = 4,
mlp_decay = 1e-4, mlp_max_iter = 500
)
model_registry = model_factory.as_r_style_dict()



#Return another data class that satisfies this condition here:
def cfperm_feature_pvals(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray, 
    vimp: str = 'loco',
    seed: int = 0,
    n_perm: int = 200,
    n_perm_cate: int = 50,
    top_k: int = 1,
    level_feature: float = 0.99,
    level_across_feature: float = 0.99,
    clip_e: float = 0.01,
    clip_wtilde: float = 1e-3,
    model_m: str = 'rf_regressor',
    model_e: str = 'rf_classifier',
    model_tau: str = 'rf_regressor',
    model_nu: str = 'rf_regressor',
    n_estimators: int = 120
):
    rng = np.random.default_rng(seed)
    X = np.asarray(X)
    Y = _as_1d(Y)
    W = _as_1d(W).astype(int)
    n, p = X.shape        
    if vimp == 'permucate':
        imp = permuCATE_vimp(X, Y, W, model_m, model_e, model_tau, model_nu, n_perm = n_perm_cate)
    elif vimp == 'loco':
        imp = vimp_loco_r_risk(X, Y, W, model_m, model_e)
    elif vimp == 'grf':
        imp = cf_variable_importance(X, Y, W, model_psm = 'logistic', model_outcome = 'rf',
            n_estimators = n_estimators)
    else:
        raise ValueError("Wrong Value, please return the appropriate version here!")
    #Permute the refit the variable importance:
    perm_imp = np.zeros((p, n_perm), dtype = float)
    for b in range(n_perm):
        W_b = rng.permutation(W)
        if vimp == 'permucate':
            perm_imp[:, b] = permuCATE_vimp(X, Y, W_b, model_m, model_e, model_tau, model_nu, n_perm = n_perm_cate)
        elif vimp == 'loco':
            perm_imp[:, b] = vimp_loco_r_risk(X, Y, W_b, model_m, model_e)
        else:
            perm_imp[:, b] = cf_variable_importance(X, Y, W_b, model_psm = 'logistic', model_outcome = 'rf',
                n_estimators = n_estimators)
    feature_specific_pval = (1 + np.sum(perm_imp >= imp[:, None], axis = 1))/(1 + n_perm)
    q1_upper = np.quantile(perm_imp, q = level_feature, axis = 1)
    threshold = float(np.quantile(q1_upper, q = level_across_feature))
    rejected_feature = [j for j in range(p) if imp[j] > threshold]
    rejected = int(len(rejected_feature) >= top_k)
    return {
      'rejected': rejected,
      'rejected_feature': rejected_feature,
      'feature_specific_pval': feature_specific_pval,
      'threshold': threshold,
      'imp': imp
    }


#Testing 
rng = np.random.default_rng(2023)
n = 400
p = 16
X = rng.normal(size = (n, p))
Y = 10.0 * X[:, 0] + 2.0 * X[:, 1] + 1.0 * X[:, 2] + rng.normal(scale = 1.0, size = n)
W = rng.binomial(1, 0.5, size = n)
cfperm_feature_pvals(X, Y, W, n_perm = 30, 
    level_feature = 1.0, level_across_feature = 1.0, vimp = 'grf')

cfperm_feature_pvals(X, Y, W, n_perm = 30, vimp = 'grf')
cfperm_feature_pvals(X, Y, W, level_feature = 0.99, level_across_feature = 0.95, vimp = 'loco')
cfperm_feature_pvals(X, Y, W, level_feature = 0.99, level_across_feature = 0.95, vimp = 'permucate', n_perm_cate = 25)









