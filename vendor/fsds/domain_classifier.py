#VIMP configuration:
import math
import torch
import numpy as np
import pandas as pd
import torch.nn as nn
import xgboost as xgb
import torch.optim as optim
import torch.nn.functional as F
from dataclasses import dataclass
from sklearn.pipeline import Pipeline
from sklearn.linear_model import Ridge
from sklearn.impute import SimpleImputer
from sklearn.utils import check_random_state
from torch.utils.data import Dataset, DataLoader
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier, MLPRegressor
from typing import Any, Callable, Dict, Optional, Tuple, Union
from torch.nn import TransformerEncoder, TransformerEncoderLayer
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.linear_model import LinearRegression, LogisticRegression, Ridge

#Model Registry Class:
FitFn = Callable[..., Any]
PredictFn = Callable[..., np.ndarray]

#The model adapter class:
@dataclass(frozen = True)
class ModelAdapter:
    """ 
    The model adapter contract aligned with: name, fit, predict,
    For the entire model factory, return {'name', 'fit', 'predict'} enhanced together here.
    """
    name: str
    fit: FitFn
    predict: PredictFn


#Building the dataset here:
class DC_Dataset(Dataset):
    def __init__(self, X, Y):
        self.X = torch.tensor(X, dtype = torch.float32)
        self.Y = torch.tensor(Y, dtype = torch.long)
    def __len__(self):
        return len(self.Y)
    def __getitem__(self, i):
        X = self.X[i]
        Y = self.Y[i]
        return X, Y


#Building the domain regressor here - 
class domain_classifier(nn.Module):
    def __init__(self, input_dim, n_classes = 2, n_layers = 2, dropout = 0.15):
        super().__init__()
        dims = np.round(
            np.exp(np.linspace(np.log(max(input_dim, 2)), 0, n_layers))
        ).astype(int)
        dims = np.maximum(dims, 1)
        layers_dict = []
        for i in range(n_layers - 1):
            in_dim = dims[i]
            out_dim = dims[i + 1]
            layers_dict.append(nn.Linear(in_dim, out_dim))
            layers_dict.append(nn.BatchNorm1d(out_dim))
            layers_dict.append(nn.ReLU())
            if dropout > 0:
                layers_dict.append(nn.Dropout(dropout))
        layers_dict.append(nn.Linear(dims[-1], n_classes))
        self.net = nn.Sequential(*layers_dict)
    def forward(self, X):
        return torch.softmax(self.net(torch.tensor(X, dtype = torch.float32)), dim = -1)


'''
X = torch.rand((1000, 100))
input_dim = 100
model0 = domain_classifier(input_dim = 100, n_layers = 2)
model0.forward(X)
'''


'''
Incorporate the model agnostic nature for the variable importance:
'''
def _rng_seed(seed):
    if seed is None:
        return np.random.randint(0, 2 ** 30)
    return int(seed) % (2 ** 30)

def _as_2d_float(X):
    X = np.asarray(X, dtype = np.float64)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    return X


class ModelVIMPRegistry:
    """
    A factory class that registers and instantiates various classification models.
    Currently supports:
        - Random Forest (binary / multi-class)
        - XGBoost (binary / multi-class)
        - Logistic Regression (binary)
        - MLP Classifier (sklearn)
        - Custom domain classifier (PyTorch-based simple NN)
    All models are wrapped with a ModelAdapter to unify the fit/predict interface.
    """
    def __init__(
        self,
        ntree=150,                     # Number of trees for RF / XGBoost
        ridge_alpha=0.25,              # (Reserved, currently unused)
        nthread=1,                     # Number of threads for parallel computing
        maxit=200,                     # Max iterations for logistic regression
        max_depth=5,                   # Max depth for XGBoost trees
        gamma=0.05,                    # Minimum loss reduction for XGBoost leaf split
        eta=0.1,                       # XGBoost learning rate
        mlp_hidden_size=7,             # Number of neurons in the MLP hidden layer
        mlp_decay=1e-5,                # MLP regularization parameter (alpha)
        mlp_max_iter=500,              # MLP maximum number of iterations
        mlp_trace=False,               # (Reserved, currently unused)
        mlp_max_coef_reg=10000,        # (Reserved, currently unused)
        mlp_max_coef_clf=10000,        # (Reserved, currently unused)
        warn_xgb_labels=True,          # (Reserved, currently unused)
        positive_class=1,              # The positive class label for binary classification
        activation='relu'):            # Activation function for MLP):
        # Store parameters as instance attributes with proper type casting
        self.ntree = int(ntree)
        self.positive_class = positive_class
        self.ridge_alpha = float(ridge_alpha)
        self.nthread = int(nthread)
        self.maxit = int(maxit)
        self.max_depth = int(max_depth)
        self.gamma = float(gamma)
        self.eta = float(eta)
        self.mlp_hidden_size = int(mlp_hidden_size)
        self.mlp_decay = float(mlp_decay)
        self.mlp_max_iter = int(mlp_max_iter)
        self.warn_xgb_labels = bool(warn_xgb_labels)
        self.ridge_alpha = float(ridge_alpha) # Redundant assignment, kept for compatibility
        self.activation = activation
        self.mlp_max_coef_reg = float(mlp_max_coef_reg)
        self.mlp_max_coef_clf = float(mlp_max_coef_clf)
    #Binary:
    @staticmethod
    def _brier_score(y_true, y_probs):
        n = y_true.shape[0]
        n_classes = y_probs.shape[1]
        one_hot = np.zeros_like(y_probs, dtype = np.float64)
        one_hot[np.arange(n), y_true] = 1.0
        return np.mean((y_probs - one_hot) ** 2)
    # ---------- XGBoost Multi-class Classifier ----------
    def make_xgb_multiclassifier(self):
        #Incorporate the model dictionary here.
        def fit(X, y, seed):
            y_arr = np.asarray(y).ravel()
            _, y_enc = np.unique(y_arr, return_inverse=True)
            n_class = len(np.unique(y_arr))
            dtrain = xgb.DMatrix(_as_2d_float(X), label=y_enc.astype(int))
            unique_class = np.unique(y_arr)
            params = dict(
                max_depth = int(self.max_depth),
                gamma = float(self.gamma),
                eta = float(self.eta),
                objective = 'multi:softprob',
                eval_metric = 'mlogloss',
                num_class = len(unique_class)
            )
            booster = xgb.train(params, dtrain, num_boost_round = self.ntree)
            return booster, y_enc.astype(int)
        def predict(fit_obj, X_new):
            booster, label = fit_obj
            d = xgb.DMatrix(_as_2d_float(X_new))
            #output all of the probability here:
            out = np.asarray(booster.predict(d), dtype = float)
            return out#You need the two-dimensional array: 
        return ModelAdapter(name = 'xgb_multiclassifier', 
            fit = fit, predict = predict)
    # ---------- XGBoost Binary Classifier ----------
    def make_xgb_classifier(self):
        def fit(X, y, seed=None):
            y_arr = np.asarray(y).ravel()
            labels = np.unique(y_arr)
            y_binary = (y_arr == self.positive_class).astype(float)
            dtrain = xgb.DMatrix(_as_2d_float(X), label=y_binary)
            params = {
                'max_depth': self.max_depth,
                'gamma': 0.075,
                'eta': 0.075,
                'objective': 'binary:logistic',   
                'nthread': self.nthread,
                'eval_metric': 'logloss'
            }
            if seed is not None:
                params['seed'] = int(seed) % (2**10)
            booster = xgb.train(params, dtrain, num_boost_round=self.ntree)
            return booster, labels.astype(int)
        def predict(fit_obj, X_new):
            booster, labels = fit_obj
            d = xgb.DMatrix(_as_2d_float(X_new))
            p_pos = np.asarray(booster.predict(d), dtype = float)
            probs = np.column_stack([1-p_pos, p_pos])
            return probs
        return ModelAdapter(name = 'xgb_classifier',
            fit = fit, predict = predict)
    # ---------- Random Forest Binary Classifier ----------
    def make_rf_classifier(self):
        def fit(X, Y, seed):
            Xa = np.asarray(X, dtype = float)
            Ya = np.asarray(Y, dtype = int).ravel()
            n, p = Xa.shape
            return RandomForestClassifier(
                n_estimators = int(self.ntree),
                max_features = 'sqrt',
                min_samples_leaf = round(np.sqrt(n)/2),
                random_state = seed,
                n_jobs = (self.nthread if self.nthread > 1 else None)
            ).fit(Xa, Ya)
        def predict(fit_obj, X_new):
            probs = np.asarray(
                fit_obj.predict_proba(_as_2d_float(X_new)),
                dtype = float
            )
            return probs
        return ModelAdapter(name = 'rf_classifier',
            fit = fit, predict = predict)
    # ---------- Random Forest Multi-class Classifier ----------
    def make_rf_multinomial(self):
        def fit(X, Y, seed):
            n_rounds = int(self.ntree)
            s = _rng_seed(seed)
            _, y_enc = np.unique(np.asarray(Y).ravel(), return_inverse = True)
            Xa = np.asarray(X, dtype = float)
            Ya = np.asarray(Y, dtype = float).ravel()
            n, p = Xa.shape
            return RandomForestClassifier(
                n_estimators = n_rounds,
                max_features = 'sqrt',
                min_samples_leaf = round(np.sqrt(n)/2),
                random_state = seed,
                n_jobs = (self.nthread if self.nthread > 1 else None)
            ).fit(Xa, Ya)
        def predict(fit_obj, X_new):
            probs = np.asarray(
                fit_obj.predict_proba(_as_2d_float(X_new)),
                dtype = float
            )
            return probs
        return ModelAdapter(name = 'rf_multinomial', fit = fit, predict = predict)  
    # ---------- Logistic Regression Classifier ----------
    def make_logistic_classification(self):
        maxit = self.maxit
        def fit(X, y, seed = None):
            s = _rng_seed(seed)
            clf = Pipeline([
                ('impute', SimpleImputer(strategy='median')),
                ('scale', StandardScaler()),
                ('logistic', LogisticRegression(
                    solver='lbfgs',
                    max_iter=maxit,
                    random_state=s,
                )),
            ])
            return clf.fit(_as_2d_float(X), np.asarray(y).ravel())
        def predict(fit_obj, X_new):
            probs = np.asarray(
                fit_obj.predict_proba(_as_2d_float(X_new)),
                dtype = float
            )
            return probs
        return ModelAdapter(name = 'logistic_classifier', fit = fit, predict = predict)
    # ---------- Custom Domain Classifier (PyTorch NN) ----------
    def make_domain_classifier(self):
        def fit(X, Y, batch_size = 16, n_epochs = 5, seed = 2026):
            s = _rng_seed(seed)
            n_classes = len(set(np.asarray(Y)))
            model = domain_classifier(input_dim = X.shape[1], n_classes = n_classes, n_layers = 2, dropout = 0.15)
            optimizer = optim.AdamW(model.parameters(), lr = 0.001, betas = (0.9, 0.999), eps = 1e-7)            
            dataset_dc = DC_Dataset(X, Y)
            data_loader = DataLoader(dataset_dc, batch_size = batch_size, shuffle = True)
            loss_fn = nn.CrossEntropyLoss()
            total_loss = []
            train_loss = 0
            optimizer = torch.optim.AdamW(model.parameters(), lr = 0.0001, betas = (0.9, 0.999), eps = 1e-7)
            for epoch in range(n_epochs):
                model.train()
                for x, y in data_loader:
                    loss = loss_fn(model.forward(x), y)
                    optimizer.zero_grad()
                    loss.backward()
                    optimizer.step()
                print(f"Epoch {epoch + 1}/{n_epochs}, The training loss is: {train_loss:.6f}")
                with torch.no_grad():
                    train_loss = loss_fn(model.forward(x), y).item()
                    total_loss.append(train_loss)
                #return the predictive model:
            return model
        def predict(fit_obj, X_new):
            y_output = fit_obj.forward(X_new)#(n, k)
            return np.asarray(y_output.detach().cpu().numpy())
            #output the 2d probability matrix, the performance here is not that important.
        return ModelAdapter(name = 'nn_clf', fit = fit, predict = predict)
    def make_mlp_classifier(self):
        """
        Build an MLP classifier adapter using sklearn's MLPClassifier.
        Includes imputation and scaling in a preprocessing pipeline.
        Hidden layer size, regularization alpha, and max iterations are controlled
        by constructor parameters.
        predict outputs a probability matrix of shape (n_samples, n_classes).
        """
        hidden = self.mlp_hidden_size
        max_iter = self.mlp_max_iter
        alpha = self.mlp_decay
        activation = self.activation
        def fit(X, Y, seed=None):
            s = _rng_seed(seed)
            clf = Pipeline([
                ('impute', SimpleImputer(strategy='median')),
                ('scale', StandardScaler()),
                ('mlp', MLPClassifier(
                    hidden_layer_sizes = (hidden, ),
                    max_iter = max_iter,
                    alpha = alpha,
                    random_state = s,
                    early_stopping = True,
                    )),
                ]
            )
            model = clf.fit(_as_2d_float(X), np.asarray(Y).ravel().astype(int))
            return model
        def predict(fit_obj, X_new):
            probs = np.asarray(
                fit_obj.predict_proba(_as_2d_float(X_new)),
                dtype=float,
            )
            return probs
        return ModelAdapter(name="mlp_classifier", fit=fit, predict=predict)
    def adapters(self):
        """
        Return a dictionary mapping model names to their corresponding ModelAdapter instances.
        This provides a central registry for all available models.
        """
        return {
        'rf_classifier': self.make_rf_classifier(),
        'rf_multinomial': self.make_rf_multinomial(),
        'neural_network_clf': self.make_domain_classifier(),
        'mlp_classifier': self.make_mlp_classifier(),
        'xgb_classifier': self.make_xgb_classifier(),
        'xgb_multiclassifier': self.make_xgb_multiclassifier(),
        'logistic_clf': self.make_logistic_classification()
        }
    def as_r_style_dict(self):
        """
        Convert the adapter registry into an R-style dictionary format.
        Each entry contains 'name', 'fit', and 'predict' keys.
        Useful for interoperability with R-like data structures or external frameworks.
        """
        mapping = self.adapters()
        registry_map = {k: {'name': v.name, 
        'fit': v.fit, 'predict': v.predict}
        for k, v in mapping.items()}
        return registry_map

#MODEL_REGISTRY = ModelVIMPRegistry()
#model_registry_map = MODEL_REGISTRY.as_r_style_dict()


#Calculate the Brier(L2) Risk from the predicted probabilities and the ground-truth value:
#The gt is (1,0,0,0,...,0) with the 1 as the ground-truth category.
def _brier_risk(Y_true, probs, n_classes):
    Y_true = np.asarray(Y_true, dtype=int).ravel()
    probs = np.asarray(probs, dtype=float)
    if probs.ndim == 1:
        probs = probs.reshape(-1, 1)
    n = len(Y_true)
    onehot = np.zeros((n, n_classes), dtype=np.float64)
    onehot[np.arange(n), Y_true] = 1.0
    if probs.shape[1] == 1 and n_classes == 2:
        probs = np.column_stack([1.0 - probs[:, 0], probs[:, 0]])
    return float(np.mean((probs - onehot) ** 2))



'''
Calculate the variable importance from an arbitrary model adapters:
The notions of the variable importance include LOCO(leave-one-covariate-out) and 
CondPerm(Conditional Permutation Variable Importance).
'''
def compute_vimp_from_adapter(
    adapter, fitted_obj,
    X_train, Y_train,
    X_test, Y_test, vimp_type = 'LOCO',
    n_classes = 2, random_state = 42,
    n_perm_reps = 10,
    oob_config: Optional[OOBVIMPConfig] = None,
):
    n, p = X_test.shape
    predictions_test = adapter['predict'](fitted_obj, X_test)
    if predictions_test.ndim == 1:
        predictions_test = predictions_test.reshape(-1, 1)
    #Getting the One-Hot results and return to the Overall L2-loss:
    n_classes = len(set(Y_train).union(set(Y_test)))
    Y_onehot = np.zeros((n, n_classes), dtype = np.float64)
    Y_onehot[np.arange(n), Y_test] = 1.0
    baseline_risk = np.mean((predictions_test - Y_onehot) ** 2)
    vimp_scores = np.zeros(p)
    if vimp_type == 'LOCO':
        for j in range(p):
            keep_cols = [c for c in range(p) if c != j]
            adapter_new = adapter.copy()
            X_train_j = X_train[:, keep_cols]
            X_test_j  = X_test[:, keep_cols]
            fitted_j = adapter_new['fit'](X_train_j, Y_train, seed = random_state + j)
            prob_j = adapter_new['predict'](fitted_j, X_test_j)
            #It's not necessarily -> to make them aligned..
            if prob_j.ndim == 1:
                prob_j = prob_j.reshape(-1, 1)
            risk_j = np.mean((prob_j - Y_onehot) ** 2)
            vimp_scores[j] = risk_j - baseline_risk
    elif vimp_type == 'condperm':
        # Conditional permutation: average risk increase over n_perm_reps draws.
        rng = np.random.default_rng(seed = random_state)
        n_reps = max(1, int(n_perm_reps))
        for j in range(p):
            rep_scores = np.zeros(n_reps, dtype = np.float64)
            for rep in range(n_reps):
                X_test_perm = X_test.copy()
                X_test_perm[:, j] = rng.permutation(X_test[:, j])
                prob_perm = adapter['predict'](fitted_obj, X_test_perm)
                if prob_perm.ndim == 1:
                    prob_perm = prob_perm.reshape(-1, 1)
                risk_perm = np.mean((prob_perm - Y_onehot) ** 2)
                rep_scores[rep] = risk_perm - baseline_risk
            vimp_scores[j] = float(np.mean(rep_scores))
    elif vimp_type == 'oob':
        cfg = oob_config or OOBVIMPConfig()
        oob_mask = _domain_oob_mask(
            X_train, Y_train, seed=random_state, oob_config=cfg
        )
        X_oob = np.asarray(X_train, dtype=float)[oob_mask]
        Y_oob = np.asarray(Y_train, dtype=int).ravel()[oob_mask]
        p = X_oob.shape[1]
        n_classes_oob = len(set(Y_train).union(set(Y_oob)))
        probs_base = adapter['predict'](fitted_obj, X_oob)
        baseline_risk = _brier_risk(Y_oob, probs_base, n_classes_oob)
        rng = np.random.default_rng(seed=random_state)
        n_reps = max(1, int(cfg.n_perm_reps if oob_config is not None else n_perm_reps))
        if oob_config is None and n_perm_reps != cfg.n_perm_reps:
            n_reps = max(1, int(n_perm_reps))
        vimp_scores = np.zeros(p, dtype=np.float64)
        for j in range(p):
            rep_scores = np.zeros(n_reps, dtype=np.float64)
            for rep in range(n_reps):
                X_perm = X_oob.copy()
                X_perm[:, j] = rng.permutation(X_oob[:, j])
                probs_perm = adapter['predict'](fitted_obj, X_perm)
                rep_scores[rep] = _brier_risk(Y_oob, probs_perm, n_classes_oob) - baseline_risk
            vimp_scores[j] = float(np.mean(rep_scores))
    else:
        raise ValueError(
            f"Unknown vimp_type: {vimp_type}. Use 'LOCO', 'condperm', or 'oob'."
        )
    # Sort descending: higher VIMP = more important
    rank_desc = np.argsort(-vimp_scores)
    return vimp_scores, rank_desc



'''
Optional:
Benchmark Boruta/ShadowVIMP/VitaPIMP + XGB/MLP/Logistic Regression
'''
DC_VIMP_MODEL_KEYS = {
    "rf": "rf_classifier",
    "xgb": "xgb_classifier",
    "mlp": "mlp_classifier",
    "mlp_classifier": "mlp_classifier",
    "logistic": "logistic_clf",
    "logistic_clf": "logistic_clf",
}


def compute_dc_vimp_batches(
    df1_X,
    df2_X,
    model="rf",
    vimp_type="LOCO",
    seed=0,
    test_size=0.3,
    registry=None,
    n_perm_reps=10,
):
    """
    Fit domain classifier on train split; LOCO/condperm VIMP on eval split.

    This function computes the variable importance 
    Parameters
    ----------
    df1_X : array-like, shape (n_samples1, n_features)
        Feature matrix for domain 1.
    df2_X : array-like, shape (n_samples2, n_features)
        Feature matrix for domain 2.
    model : str, default="rf"
        Model key to use. Must correspond to a key in `DC_VIMP_MODEL_KEYS` mapping
        or directly match a key in the registry dictionary.
    vimp_type : str, default="LOCO"
        Type of variable importance to compute. Options depend on the underlying
        `compute_vimp_from_adapter` function (e.g., "LOCO", "condperm").
    seed : int, default=0
        Random seed for reproducibility of train/test split and model fitting.
    test_size : float, default=0.3
        Proportion of the combined data to use as the evaluation set.
    registry : dict, optional
        A dictionary mapping model names to adapters (with 'fit' and 'predict' keys).
        If None, a default `ModelVIMPRegistry` is created with `ntree=100, nthread=1`.
    n_perm_reps : int, default=10
        Number of permutation repetitions for VIMP computation (for the conditional permutation variable importances).
    **registry_kwargs : additional keyword arguments
        Passed to the registry constructor when `registry` is None.

    Returns
    -------
    vimp_result : Any
        The result from `compute_vimp_from_adapter` (structure depends on the VIMP
        implementation, typically a DataFrame or array of importance scores).
    """
    X = np.vstack([np.asarray(df1_X, dtype=float), np.asarray(df2_X, dtype=float)])
    W = np.concatenate(
        [np.zeros(len(df1_X), dtype=int), np.ones(len(df2_X), dtype=int)]
    )
    X_train, X_test, W_train, W_test = train_test_split(
        X,
        W,
        test_size=test_size,
        random_state=seed,
        stratify=W,
    )
    if registry is None:
        registry = ModelVIMPRegistry(ntree=100, nthread=1).as_r_style_dict()
    adapter_key = DC_VIMP_MODEL_KEYS.get(model, model)
    adapter = registry[adapter_key]
    fitted = adapter["fit"](X_train, W_train, seed=seed)
    return compute_vimp_from_adapter(
        adapter,
        fitted,
        X_train,
        W_train,
        X_test,
        W_test,
        vimp_type=vimp_type,
        random_state=seed,
        n_perm_reps=n_perm_reps,
    )


'''
MODEL_REGISTRY = ModelVIMPRegistry()
model_registry_map = MODEL_REGISTRY.as_r_style_dict()

model1 = model_registry_map['logistic_clf']
X_train = np.random.random((1000, 20))
Y_train = np.random.choice([0,1], 1000, replace = True)
fitted_obj = model1['fit'](X_train, Y_train, seed = 42)
model1['predict'](fitted_obj, X_train)

#fitted_obj.predict_proba(X_train)

X_train = np.random.random((1000, 20))
Y_train = np.random.choice([0,1,2,3], 1000, replace = True)
X_test = np.random.random((1000, 20))
Y_test = np.random.choice([0,1,2,3], 1000, replace = True)
#Would need StandardScaler() here.
adapter = model_registry_map['rf_multinomial']
fitted_obj = adapter['fit'](X_train, Y_train, seed = 42)#You need to fit first on the evaluation set!
vimp_scores, top_features = compute_vimp_from_adapter(
    adapter=adapter,
    fitted_obj=fitted_obj,
    X_train=X_train,
    Y_train=Y_train,
    X_test=X_test,
    Y_test=Y_test,
    vimp_type='condperm',
    n_classes=2
)
'''



























































