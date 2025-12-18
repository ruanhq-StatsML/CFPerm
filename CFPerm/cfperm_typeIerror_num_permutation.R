#Type-I error evaluation for the CFPerm procedure with respect to varying number of permutations:
setwd("/ihome/lmentch/her76/causalforestProject")
source("permutation_covariateShift_utils.R")
source("vimp_cf_variance_utils.R")
source("CFPerm_func.R")
library(grf)
set.seed(42)
corr_seq <- seq(0, 0.7, 0.1)
beta_train <- c(1, 1, 1, 1, 1, 1, 1, 1)
beta_test <- c(1, 1, 1, 1, 0, 0, 0, 0)
eps_seq <- c(0.1, 0.2, 0.5, 1, 2, 3, 5, 7.5, 10, 15, 20, 30, 50, 100)
mean_shift_seq <- seq(0.2, 2, 0.2)
B_seq <- c(300, 500, 1000, 2000, 3000, 5000, 10000, 20000)
power_df <- data.frame(
  n_nuisance = B_seq,
  power0 = numeric(6),
  power1 = numeric(6)
  )
B <- 500
power_list0 <- numeric(B)
power_list1 <- numeric(B)
for(j in 1:6){
  n_B <- B_seq[j]
  for(i in 1:B){
    df_train1 <- LM_generation(500, beta_hat = beta_test, cor = 0.3, n_nuisance = 12, eps = 5)$df_return
    df_test1 <- LM_generation(500, beta_hat = beta_test, cor = 0.3, n_nuisance = 12, eps = 5)$df_return
    #Conduct the causal forest for testing:
    power_list0 <- cfperm(df_train1, df_test1, n_perm = n_B, q1 = 0.99, level = 0.05, num.trees = 150)$result
    power_list1 <- cfperm(df_train1, df_test1, n_perm = n_B, q1 = 1, level = 0, num.trees = 150)$result
  }
  power_df$power0[j] <- mean(power_list0)
  power_df$power1[j] <- mean(power_list1)
  write.csv(power_df, "cfperm_typeIerror_scale_n_permutation.csv")
}


