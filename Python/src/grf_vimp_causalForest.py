#Causal Forest based Variable Importance:
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from econml.dml import CausalForestDML
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LogisticRegression


'''
Fit the nuisance models: outcome model, propensity score model and the pseudo-outcome model
Parameters:
-----------
X:                Array of Shape (n_exist + n_new, p)
Y:                Array of Shape (n_exist + n_new, ) Response
W:                Array of Shape (n_exist + n_new, ) Treatment, Currently only Support Binary
model_psm:        Propensity Score Model
model_outcome:    Outcome Model
seed(int):        Random Seed
n_estimators:     Number of Trees
min_samples_leaf: Minimal number of samples in one leaf 
max_depth:        Maximal depth of causal tree
normalize:        Whether to normalize the variable importance to make them sum up to 1
'''
def cf_variable_importance(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    model_psm, model_outcome, 
    *, seed: int = 0,
    n_estimators: int = 120,
    min_samples_leaf: int = 5,
    max_depth = 4,
    normalize = True
):
    X = np.asarray(X)
    Y = _as_1d(Y)
    W = _as_1d(W).astype(object)
    p = X.shape[1]
    if model_psm == 'logistic':
        model_e = LogisticRegression(solver = 'lbfgs', max_iter = 5000, random_state = seed)
    else:
        model_e = RandomForestClassifier(n_estimators = n_estimators,
            max_depth = max_depth, min_samples_leaf = min_samples_leaf,
            max_features = 'sqrt', n_jobs = 1, random_state = seed)
    if model_outcome == 'lm':
        model_t = LinearRegression()
    elif model_outcome == 'rf':
        model_t = RandomForestRegressor(n_estimators = n_estimators,
            max_depth = max_depth, min_samples_leaf = min_samples_leaf,
            max_features = 0.33, n_jobs = 1, random_state = seed)
    else:
        model_t = Ridge(alpha = 1.0)
    cf_model = CausalForestDML(
        model_y = model_t,
        model_t = model_e,
        n_estimators = n_estimators,
        min_samples_leaf = min_samples_leaf,
        max_depth = max_depth,
        discrete_treatment = True,
        cv = 4,
        random_state = seed
    )
    cf_model.fit(Y = Y, T = W, X = X)
    imp = np.asarray(cf_model.feature_importances_, dtype = float)
    if normalize:
        return imp/np.sum(imp)
    else:
        return imp


#model_factory = ModelRegistry(
#ntree = 150, ridge_alpha = 0.25,
#nthread = 1, maxit = 500,
#max_depth = 5, gamma = 0.25,
#eta = 0.15, mlp_hidden_size = 4,
#mlp_decay = 1e-4, mlp_max_iter = 500)
#model_registry = model_factory.as_r_style_dict()
#X = np.random.random((100, 20))
#Y = np.random.random((100, ))
#W = np.random.choice((0,1), 100)
#cf_variable_importance(X, Y, W, model_psm = 'rf', model_outcome = 'rf')















