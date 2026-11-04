
# Meta-Learner Based Hypothesis Testing Procedure for Distribution Shift via Variable Importance

Implementation for the "Leveraging Causal Inference for Detecting Distribution Shift": It provides implementation for leveraging causal inference methods to detect distribution shift via the varying types of the causal models including Causal Forests, DR-learner and Dragon-Net with permutation of the Variable Importance. The variable importance notions include PermuCATE(conditional permutation variable importance), LOCO(leave-one-covariate-out) and the Impurity Gain/Split Frequency Based Variable Importance for the Causal Forest.

The hypothesis testing procedure is based on permuting the batch assignment among the existing group and new group with the high empirical quantile on the variable as the decision threshold. The hypothesis testing procedure is conducted via looking at whether there exist features can be extracted, the pseudo-code for the algorithm is presented as below:
<img width="1326" height="1392" alt="CFPerm_MetaLearner" src="https://github.com/user-attachments/assets/07c9f1d7-e01d-4f33-b922-8e6f4b033424" />


Thus the post-hoc subgroup analysis can be developed for the distribution shift localization - to calculate the mean difference between the two groups for discretization - below is one of the illustration examples. 
<img width="1630" height="1296" alt="SubsetAnalysis" src="https://github.com/user-attachments/assets/8fdd99f6-b3ef-4c06-af73-f31f17c2ff01" />


### Motivation: Meta-Learner and Distribution Shift
We leverage the meta-learner(causal forest, doubly-robust pseudo-outcome learner and the R-learner) to quantify the distance between the two batches of datasets, then the variable importance for this version of meta-learner is a proxy for the distribution shift - from another aspect, so called OOD Variable Importance. The illustration of the feature selection as the validity of this proposed framework in both Covariate Shift(P(X) shift) and Concept Drift(P(Y|X) shift).

- For Covariate Shift(the feature selection consistency in terms of kendall's tau correlation is leveraged for efficiency of the ranking of the feature importance in the covariate shift).
<img width="900" height="970" alt="vecshift_lambda_cor06_polished" src="https://github.com/user-attachments/assets/0275705c-7a57-4edb-9832-6356cc2c541b" />

- For Concept Drift(the feature selection consistency in terms of kendall's tau correlation is leveraged for efficiency of the ranking of the feature importance in the concept drift).
<img width="900" height="970" alt="uq_vimpood_cd_on_cs_allmethods_polished" src="https://github.com/user-attachments/assets/65925cae-d308-4651-8ba1-979d122e5dae" />

This framework possess high flexibility both in terms of the models as well as the notions of variable importances. 
<img width="1700" height="900" alt="metalearner_flexibility" src="https://github.com/user-attachments/assets/dd3f66ed-4d54-4f57-b5ba-58760b803d3a" />

The variable importance here can be denoted as a proxy of the Out-of-Distribution(OOD) Variable Importance - characterizing each feature's contribution to the distribution shift. Thus the procedure is highly flexible with different notions of the variable importances and the meta-learners(as long as a valid objective function exist for the meta-learner to estimate the discrepancy between the two batches of data). ***Emphasize: The meta-learner procedure here won't possess causal interpretation as under the alternative hypothesis, the strong ignorability condition won't be satisfied - the pseudo-outcome learner itself is primarily served as the distance estimator instead of the treatment effect estimator, the Causal Forest(CF-Split) procedure can serve as a parsimonious empirical tool for the hypothesis testing procedure as the validity of its variable importance requires causal interpretation.***

## Python version:
```python
!pip install cfperm
import cfperm
@dataclass
class CFPermResult:
    feature_names: list
    imp_observed: pd.Series
    p_values: pd.Series
    q1_upper: pd.Series
    threshold: float
    rejected_feature_indice: list
    perm_imp: np.ndarray
    method: str

'''
model_factory = ModelRegistry(
ntree = 150, ridge_alpha = 0.25,
nthread = 1, maxit = 500,
max_depth = 5, gamma = 0.25,
eta = 0.15, mlp_hidden_size = 4,
mlp_decay = 1e-4, mlp_max_iter = 500
)
model_registry = model_factory.as_r_style_dict()
'''

'''
Hypothesis Testing for Distribution Shift via Variable Importance:
- X:               Covariates array - (n_exist + n_new, p)
- Y:               Responses - (n_exist + n_new, )
- W:               Treatment(Batch) Assignment - (n_exist + n_new, )
- seed:            Random State Seeds
- n_perm:          Number of permutations in the variable importance
- n_perm_cate:     Number of permutations for the conditional permutation variable importance
- n_split:         Number of splits for the cross-fitting procedure
- clip_e:          Clipping value for the propensity score
- alpha:           Significance Level for the test
- model_m(str):    Outcome Model
- model_e(str):    Propensity Score Model
- model_tau(str):  Model for treatment effect
- model_nu(str):   Model for the fitting covariates in X
Returns the cross-fitted value (m_hat, e_hat) evaluated on the whole data
as well as the fitted outcome model and propensity score models.
'''
def cfperm_feature_pvals(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray, 
    vimp: str = 'loco',
    seed: int = 0,
    n_perm: int = 200,
    n_perm_cate: int = 50,
    n_splits: int = 5,
    top_k: int = 1,
    level_feature: float = 0.99,
    level_across_feature: float = 0.99,
    clip_e: float = 0.01,
    clip_wtilde: float = 1e-3,
    model_m: str = 'rf_regressor',
    model_e: str = 'rf_classifier',
    model_tau: str = 'rf_regressor',
    model_nu: str = 'rf_regressor',
    model_nu_binaryclassifier: str = 'logistic_classifier',
    model_nu_multiclassifier: str = 'multinomial_classifier',
    n_estimators: int = 120
):
    rng = np.random.default_rng(seed)
    X = np.asarray(X)
    Y = _as_1d(Y)
    W = _as_1d(W).astype(int)
    n, p = X.shape        
    if vimp == 'permucate':
        imp = permuCATE_vimp(X, Y, W, MODEL_REGISTRY, model_m, model_e, model_tau, model_nu, n_perm = n_perm_cate)
    elif vimp == 'loco':
        imp = vimp_loco_r_risk(X, Y, W, MODEL_REGISTRY, model_m, model_e)
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

'''
#Testing 
rng = np.random.default_rng(2023)
n = 400
p = 16
X = rng.normal(size = (n, p))
Y = 10.0 * X[:, 0] + 2.0 * X[:, 1] + 1.0 * X[:, 2] + rng.normal(scale = 1.0, size = n)
W = rng.binomial(1, 0.5, size = n)
cfperm_feature_pvals(X, Y, W, n_perm = 30, 
    level_feature = 1.0, level_across_feature = 1.0, vimp = 'grf')
cfperm_feature_pvals(X, Y, W, n_perm = 30, 
    level_feature = 1.0, level_across_feature = 1.0, vimp = 'loco')
cfperm_feature_pvals(X, Y, W, n_perm = 2, 
    level_feature = 1.0, level_across_feature = 1.0, vimp = 'permucate')

#cfperm_feature_pvals(X, Y, W, n_perm = 30, vimp = 'grf')
#cfperm_feature_pvals(X, Y, W, level_feature = 0.99, level_across_feature = 0.95, vimp = 'loco')
#cfperm_feature_pvals(X, Y, W, level_feature = 0.99, level_across_feature = 0.95, vimp = 'permucate', n_perm_cate = 25)
#cfperm_feature_pvals(X, Y, W, level_feature = 1.0, level_across_feature = 0.95, vimp = 'permucate', n_perm_cate = 25)
'''
```

## R version:
#### Installation (local)

```r
install.packages(c("devtools", "roxygen2", "testthat", "grf", "MASS"))
devtools::install_local("path/to/CFPerm")
```

## Quick example

```r
library(CFPerm)
set.seed(1)
#Testing Case I, simulate with no distribution shift:
set.seed(2026)
n0 <- 200
n1 <- 200
n <- n0 + n1
p <- 10
X <- matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)
colnames(X) <- paste0("X", seq_len(p))
T <- c(rep(0, n0), rep(1, n1))
Y <- 1 + 0.8 * X[, 1] - 0.5 * X[, 2] + rnorm(n, sd = 1)
#The test result via the permuCATE variable importance
es_null_permucate <- cfperm(
  X = X,
  Y = Y,
  T = T,
  n_perm = 30,
  vimp = "permuCATE",
  seed = 2026,
  level_feature = 0.05,
  level_across_feature = 0.05,
  top_k = 1
)
#The test result via the LOCO variable importance
es_null_loco <- cfperm(
  X = X,
  Y = Y,
  T = T,
  n_perm = 30,
  vimp = "loco",
  seed = 2026,
  level_feature = 0.05,
  level_across_feature = 0.05,
  top_k = 1
)
#The test result via the GRF variable importance
es_null_grf <- cfperm(
  X = X,
  Y = Y,
  T = T,
  n_perm = 200,
  vimp = 'GRF',
  seed = 2020,
  level_feature = 0.005,
  level_across_feature = 0.005
)
#A reasonable amount of the permutation is required.
```

## Development

```r
devtools::document()
devtools::test()
devtools::check()
```

## Demonstration for the Attribution Analysis: No Uniform Attribution for the MSE and the benchmark Comparisons:

https://colab.research.google.com/drive/1t12mtdzDb9pouSae2bvrSjFcm19miFK2

### Impossibility Issue issue in the Decomposition of Distribution Shift into Individual Feature or Separate Sources

- **Upon detecting the change of the data, one crucial yet important aspect would be - is that possible to attribute the root-cause of the change of data into individual sources? Our answer is that, it's impossible! Thus The impossibility issue is proposed with the following three levels.**
-  Level1: The drop of the model performance(increase of MSE), there are infinite number of the distribution shifts to characterize or achieve this degree of performance drop.
-  Level2: Upon observing the distribution shift, there are infinite number of the Covariate Shift and Concept Drift to characterize or achieve this degree of distribution shift.
-  Level3: Upon observing the change of the variable importance, there are infinite number of change of the feature and the dependence structure to characterize or achieve this degree of variable importance change. 

- **For any observed distribution shift or performance drop, disentangle it into separate sources would be impossible, there are infinite number of ways to replicate this distribution shift or this degree of the performance degradation** The mean shift + epsilon shift(data quality degradation) can mimic(with 50 bootstrap resamples and the confidence interval is established) the degree of the performance degradation like that from the concept drift.
<img width="1796" height="552" alt="Screenshot 2026-08-03 at 09 57 54" src="https://github.com/user-attachments/assets/bd74444f-992d-4928-ac54-df080e397cf7" />
- **Consequently, upon observing a notable drop in model performance, we prioritize post-hoc feature selection or localization of distribution-shift drivers over disentangling the shift into concept drift versus covariate shift, as such decomposition is not identifiable, Subset Localization is all you need!**
- It gives people concise proxy for efficiently dealing with the model performance degradation in the deployed ML model - distribution shift driver localization is what you will need.

***The Impossibility Issue for the Disentanglement has deep implications in areas including Econometrics and Quantitative Finance***
- ****Econometrics**** ~~ For econometrics, the impossibility issue means that observed changes in policy outcomes cannot be uniquely decomposed into covariate shift versus concept drift, nor into individual feature contributions—infinitely many mixtures reproduce the same shift. Rather than pursuing an unidentifiable decomposition, the practical priority is post-hoc subset localization: identifying the subgroups or covariates driving the change. This yields robust, assumption-light diagnostics for policy adjustment.
- ****Quantitative Finance**** ~~ In quantitative finance, a decline in factor performance or IC likewise cannot be uniquely attributed to exposure drift versus alpha decay, since infinite combinations of universe changes, style rotations, and regime shifts produce identical degradation. Factor diagnostics should therefore focus on subset localization—pinpointing the specific factors, sectors, or regimes responsible—rather than attempting a non-identifiable decomposition. Subset localization is all you need.
- ****Credit Risk Modeling**** ~~ In the regime of credit risk modeling, an observed rise in delinquency cannot be uniquely split into covariate shift versus concept drift, or into individual features. Income shifts, borrower-mix changes, macro effects, and approval-policy changes can all raise delinquency, and infinite mixtures can produce the same rate. Instead of an unidentifiable decomposition, the practical goal is subset localization: find the income bands, borrower segments, channels, or risk buckets driving the change. This gives robust, assumption-light diagnostics for credit-policy adjustment.
- ****Manufacturing Maintenance**** ~~ For manufacturing predictive maintenance, a rise in false alarms or missed failures cannot be uniquely split into covariate shift versus concept drift, or into sensor/feature contributions. Operating conditions, sensor drift, maintenance policy, and failure physics can all degrade performance, and infinite mixtures can produce the same alarm rate. Instead of an unidentifiable decomposition, the practical goal is subset localization: pinpoint the equipment units, operating regimes, sensor channels, or production lines responsible. This gives robust, assumption-light diagnostics for maintenance and calibration.





