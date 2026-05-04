#GRF Variable Importance:
'''
'''
from typing import Optional, Tuple, List
from econml.dml import CausalForestDML
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.ensemble import RandomForestRegressor

def grf_vimp(
    X_exist, Y_exist,
    X_new,   Y_new,
    n_estimators = 100, random_state = None, resampling = 20):
    '''
    Impurity Based variable importance for the Causal Forest
    - resampling:   20(number of resamples)
    - n_estimators: number of estimators for the Causal Forests

    Returns: 
      imp: array for the variable importance across each of the feature
    '''
    n_exist, n_new = X_exist.shape[0], X_new.shape[0]
    p = X_exist.shape[1]
    T_exist = np.zeros(n_exist, dtype = np.int64)
    T_new = np.ones(n_new, dtype = np.int64)
    X_both = np.vstack([X_exist, X_new])
    Y_both = np.concatenate([Y_exist, Y_new])
    T_both = np.concatenate([T_exist, T_new])
    rng = np.random.default_rng(random_state)
    cf = CausalForestDML(
        RandomForestRegressor(n_estimators = 100, max_depth = 5),
        model_t = 'auto',
        n_estimators = n_estimators,
        discrete_treatment = True,
        random_state = int(rng.integers(0, 2 ** 31)) if random_state is not None else None
    )
    cf.fit(Y_both, T_both, X = X_both)
    tau_te = cf.const_marginal_effect(X_new)
    return cf.feature_importance_
    #Return the array of feature importance for each of the individual feature.
    #forest = getattr(cf, '_model_final', None) or getattr(cf, 'model_final_', None)
    #if forest is not None:
    #    inner = getattr(forest, '_model', None) or getattr(forest, 'model_', None)
    #    if inner is not None and hasattr(inner, 'feature_importances_'):
    #        imp = np.asarray(inner.feature_importances_).ravel()
    #        if len(imp) == X_exist.shape[1]:
    #            return imp, None  
    #except Exception:
    #    pass
    #fallback to the permutation variable importance:













        
