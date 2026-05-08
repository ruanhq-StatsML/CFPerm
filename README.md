# CFPerm

Implementation for the "Leveraging Causal Inference for Detecting Distribution Shift": It provides implementation for leveraging causal inference methods to detect distribution shift via the varying types of the causal models including Causal Forests, DR-learner and Dragon-Net with permutation of the Variable Importance. The variable importance notions include PermuCATE(conditional permutation variable importance), LOCO(leave-one-covariate-out) and the Impurity Gain Variable Importance for the Causal Forest.


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
)#A reasonable amount of the permutation is required.

```

## Development

```r
devtools::document()
devtools::test()
devtools::check()
```
