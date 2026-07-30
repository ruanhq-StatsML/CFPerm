# Meta-Learner Based Hypothesis Testing Procedure for Distribution Shift via Variable Importance

Implementation for the "Leveraging Causal Inference for Detecting Distribution Shift": It provides implementation for leveraging causal inference methods to detect distribution shift via the varying types of the causal models including Causal Forests, DR-learner and Dragon-Net with permutation of the Variable Importance. The variable importance notions include PermuCATE(conditional permutation variable importance), LOCO(leave-one-covariate-out) and the Impurity Gain/Split Frequency Based Variable Importance for the Causal Forest.

The hypothesis testing procedure is based on permuting the batch assignment among the existing group and new group, where the 
$$
\begin{algorithm}[tb]
\caption{Testing Distribution Shift via Meta-Learner and variable importance - \textbf{CFPerm}}.\label{alg:alg_cfperm}
\\
\textbf{BEGIN}
Start with existing batch of data $\mathscr{D}_{exist} =  \{ (\mathbf{X}_{1},Y_{1}),...,(\mathbf{X}_{n_{exist}},Y_{n_{exist}})\}$ where $(\mathbf{X}_{i},Y_{i})$ and new batch of data $\mathscr{D}_{new} = \{ (\mathbf{X}_{n_{exist}+1},Y_{n_{exist}+1}),...,(\mathbf{X}_{n_{exist}+n_{new}},Y_{n_{exist}+n_{new}})\}$ where $(\mathbf{X}_{i}, Y_{i}) $, with each $\mathbf{X}_{j} = (x_{1},...,x_{p})$. Assign (control) treatment $T=0$ to existing batch of data and $T=1$ (treatment) to new batch of data. Define hypotheses $H_{0}:P(\mathbf{X}_{new}) = P(\mathbf{X}_{exist}) \ and \ P(Y_{new}|\mathbf{X}_{new}) = P(Y_{exist}|\mathbf{X}_{exist}) $ v.s. $ H_{a}: P(\mathbf{X}_{new}) \neq P(\mathbf{X}_{exist}) \ or \ P(Y_{new}|\mathbf{X}_{new}) \neq P(Y_{exist}|\mathbf{X}_{exist})$.
\begin{enumerate}
\item[\textbf{S1}] Fit the causal forest on $(\mathbf{X},Y,T)$ and record the variable importance (VIMP) for each feature $v_{1},...,v_{p}$ The notion of variable importance include \textbf{PermuCATE}, \textbf{LOCO} and \textbf{CF-split}.
\item[\textbf{S2}] Permute the treatment assignment T randomly across the data a total of $B$ times (by default $B$ = 500).  Each time, fit the causal forest on the resulting permuted batches and calculate the variable importance, recorded as $v_{b,1},...,v_{b,p}$ for $b=1,..,B$.
\item[\textbf{S3}] Construct an interval for the variable importance of each of the feature with the lower bound as the smallest value and the higher bound as the largest value across the permutations $I_{i} = (v_{i,min},v_{i,max}), i = 1,2,...,p$
\item[\textbf{S4-1}] \textbf{CFPerm0}: The maximal value of the permuted confidence interval of variable importance for each feature is recorded as $I_{max} = \{v_{max,1},...,v_{max,p} \}$.  Reject $H_{0}$ if there exists a feature whose original variable importance is larger than the maximal of upper bound among all of the features $max(I_{max})$ served as the decision threshold; otherwise, do not reject.
\item[\textbf{S4-2}] \textbf{CFPerm1}: The 99\% upper quantile for the permuted confidence interval of variable importance in each of the feature is recorded: $I_{99\%} = \{v_{1,99\%},...,v_{p,99\%}\}$ Reject $H_{0}$ based on whether there exists a feature whose original variable importance is larger than the 95\% upper quantile $q_{0.95}(I_{99\%})$ of the recorded feature-wise 99\% quantile across the features served as the decision threshold. Otherwise the null hypothesis would not be rejected. 
\end{enumerate}
\textbf{END}
\end{algorithm}\label{alg: algorithm1}
$$

### Motivation: Meta-Learner and Distribution Shift
We leverage the meta-learner(causal forest, doubly-robust pseudo-outcome learner and the R-learner) to quantify the distance between the two batches of datasets, then the variable importance for this version of meta-learner is a proxy for the distribution shift - from another aspect, so called OOD Variable Importance. The illustration of the feature selection as the validity of this proposed framework in both Covariate Shift(P(X) shift) and Concept Drift(P(Y|X) shift).

- For Covariate Shift(the feature selection consistency in terms of kendall's tau correlation is leveraged for efficiency of the ranking of the feature importance in the covariate shift).
<img width="900" height="970" alt="vecshift_lambda_cor06_polished" src="https://github.com/user-attachments/assets/0275705c-7a57-4edb-9832-6356cc2c541b" />

- For Concept Drift(the feature selection consistency in terms of kendall's tau correlation is leveraged for efficiency of the ranking of the feature importance in the concept drift).
<img width="900" height="970" alt="uq_vimpood_cd_on_cs_allmethods_polished" src="https://github.com/user-attachments/assets/65925cae-d308-4651-8ba1-979d122e5dae" />

Then we adapt the second-stage feature selection followed by the distribution shift localization procedure - where the test itself is built upon whether significant features can be selected from the hypothesis testing procedure.


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
