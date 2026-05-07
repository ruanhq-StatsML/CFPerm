# CFPerm

Implementation for the "Leveraging Causal Inference for Detecting Distribution Shift": It provides implementation for causal machine learning model to detect distribution shfit via permutation test on the variable importances.


## Python version:
```python
!pip install cfperm_vimp
from cfperm_vimp import RRPerm, DRPerm

```

## R version:
#### Installation (local)

```r
install.packages(c("devtools", "roxygen2", "testthat", "grf", "MASS", ""))
devtools::install_local("path/to/CFPerm")
```

## Quick example

```r
library(CFPerm)

set.seed(1)
sim <- LM_generation(
  n = 200,
  beta_hat = c(1, -1, 0.5),
  mean_shift = 0,
  var_shift = 1,
  cor = 0.3,
  n_nuisance = 5,
  eps = 1
)

df <- sim$df_return
df_train <- df[1:100, c(paste0("X", 1:3), paste0("X_nuis", 1:5), "Y")]
df_test  <- df[101:200, c(paste0("X", 1:3), paste0("X_nuis", 1:5), "Y")]

res <- cfperm(df_train, df_test, n_perm = 50, num.trees = 150, seed = 123)
res
```

## Development

```r
devtools::document()
devtools::test()
devtools::check()
```
