source("CFPerm_func.R")
source("data_generating_processes.R")
library(grf)
set.seed(42)

corr_seq <- seq(0, 0.7, 0.1)
beta_train <- c(1, 1, 1, 1, 1, 1, 1, 1)
beta_test <- c(1, 1, 0, 0, 0, 0, 0, 0)
scaler <- c(0.1, 0.2, 0.5, 0.75, 0.95, 0.98, 1.05, 1.1, 1.25, 1.5, 2, 3, 5)
eps_seq <- c(0.1, 0.2, 0.5, 1, 2, 3, 5, 7.5, 10, 15, 20, 30, 50, 100)
B <- 500
df_vimp_cover <- matrix(0, 20, B + 1)
#repeating the power for 15 times:
power_df <- data.frame(
  eps_seq = eps_seq,
  power_cfperm1 = numeric(14),
  power_cfperm2 = numeric(14)
  )
for(j in 1:8){
  for(k in 1:14){
    pval_list1 <- numeric(B)
    pval_list2 <- numeric(B)
    for(i in 1:B){
      df_train1 <- LM_generation(1000, beta_hat = beta_train, cor = corr_seq[j], n_nuisance = 12, eps = eps_seq[k])$df_return
      df_test1 <- LM_generation(1000, beta_hat = beta_test, cor = corr_seq[j], n_nuisance = 12, eps = eps_seq[k])$df_return
      #Conduct the causal forest for testing:
      pval_list1[i] <- cfperm(df_train1, df_test1, q1 = 1, level = 0, n_perm = 500)$result
      pval_list2[i] <- cfperm(df_train1, df_test1, q1 = 0.99, level = 0.05, n_perm = 500)$result
    }
    power_df$power_cfperm1[k] <- mean(pval_list1)
    power_df$power_cfperm2[k] <- mean(pval_list2)
  }
  write.csv(power_df, paste("power_conceptdrift_cfperm_corr", corr_seq[j] , ".csv",
    sep= "_"))
}
