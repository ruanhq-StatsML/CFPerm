#This script include the power evaluation for robustness when both covariate shift and concept drift exists:
source("CFPerm_func.R")
source("data_generating_processes.R")
library(grf)
set.seed(42)

beta_train <- c(1, 1, 1, 1, 1, 1, 1, 1)
beta_test <- c(1, 1, 0, 0, 0, 0, 0, 0)
scaler <- c(0.1, 0.2, 0.5, 0.75, 0.95, 0.98, 1.05, 1.1, 1.25, 1.5, 2, 3, 5)
eps_seq <- c(0.1, 0.2, 0.5, 1, 2, 3, 5, 7.5, 10, 15, 20, 30, 50, 100)
B <- 200
df_vimp_cover <- matrix(0, 20, B + 1)
mean_shift_seq <- seq(0, 2, 0.2)
power_df <- data.frame(
  mean_shift_seq = mean_shift_seq,
  power0 = numeric(11),
  power1 = numeric(11)
  )
n_perm <- 500
for(k in c(4, 7, 9)){
  for(j in 1:11){
  	power_list0 <- numeric(B)
	power_list1 <- numeric(B)
	for(i in 1:B){
	  df_train1 <- LM_generation(500, beta_hat = beta_train, mean_shift = 0, cor = 0.3, n_nuisance = 12, eps = eps_seq[k])$df_return
	  df_test1 <- LM_generation(500, beta_hat = beta_test, mean_shift = mean_shift_seq[j], cor = 0.3, n_nuisance = 12, eps = eps_seq[k])$df_return
      power_list0 <- cfperm(df_train1, df_test1, n_perm = 1000, q1 = 0.99, level = 0.05, num.trees = 150)$result
      power_list1 <- cfperm(df_train1, df_test1, n_perm = 1000, q1 = 1, level = 0, num.trees = 150)$result
	}
	power_df$power0[j] <- mean(power_list0)
	power_df$power1[j] <- mean(power_list1)
  }
  write.csv(power_df, paste("CFPerm_Power_Sigma", k, "SensitivityMeanShift_ConceptDrift1.csv", sep = "_"))  
}

