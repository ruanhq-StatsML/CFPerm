
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

Thus the procedure is highly flexible with different notions of the variable importances and the meta-learners(as long as a valid objective function exist for the meta-learner to estimate the discrepancy between the two batches of data). ***Emphasize: The meta-learner procedure here won't possess causal interpretation as under the alternative hypothesis, the strong ignorability condition won't be satisfied - the pseudo-outcome learner itself is primarily served as the distance estimator instead of the treatment effect estimator, the Causal Forest(CF-Split) procedure can serve as a parsimonious empirical tool for the hypothesis testing procedure.***

## Python version:
```python
!pip install cfperm
import cfperm

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

### Identifiability issue in the Decomposition of Distribution Shift into Individual Feature or Separate Sources

- **For any observed distribution shift or performance drop, disentangle it into separate sources would be impossible, there are infinite number of ways to replicate this distribution shift or this degree of the performance degradation** The mean shift + epsilon shift(data quality degradation) can mimic(with 50 bootstrap resamples and the confidence interval is established) the degree of the performance degradation like that from the concept drift.
<img width="1796" height="552" alt="Screenshot 2026-08-03 at 09 57 54" src="https://github.com/user-attachments/assets/bd74444f-992d-4928-ac54-df080e397cf7" />
- **Consequently, upon observing a notable drop in model performance, we prioritize post-hoc feature selection or localization of distribution-shift drivers over disentangling the shift into concept drift versus covariate shift, as such decomposition is not identifiable**
- It gives people concise proxy for efficiently dealing with the model performance degradation in the deployed ML model - distribution shift driver localization is what you will need.
