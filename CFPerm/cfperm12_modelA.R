#paper1_ModelA:
setwd("/ihome/lmentch/her76/causalforestProject")
source("permutation_covariateShift_utils.R")
source("vimp_cf_variance_utils.R")
source("CFPerm_func.R")
library(grf)
set.seed(42)
library(grf)
sample_size_seq <- c(200, 500, 1000, 2000)

B <- 500
df_vimp_cover <- matrix(0, 5, B + 1)
power_list <- numeric(4)
power_df <- data.frame(
  sample_size = c(200, 500, 1000, 2000),
  power_null1 = numeric(4),
  power_alt1 = numeric(4),
  power_null2 = numeric(4),
  power_alt2 = numeric(4))
for(j in 1:4){
  power_cover_null1 <- numeric(B)
  power_cover_alt1 <- numeric(B)
  power_cover_null2 <- numeric(B)
  power_cover_alt2 <- numeric(B)
  for(i in 1:B){
    n <- sample_size_seq[j]
    random <- sample(c(1, -1), 5, replace = TRUE, prob = c(0.5, 0.5))
    beta1 <- random * as.matrix(c(1,1,1,1,1))
    df_train1 <- model_A_DGP(n, alpha = 0, train = TRUE, beta = beta1)
    df_train2 <- model_A_DGP(n, alpha = 0, train = TRUE, beta = beta1)
    df_test1 <- model_A_DGP(n, alpha = 0, train = FALSE, beta = beta1)
    df_test2 <- model_A_DGP(n, alpha = 0.5, train = FALSE, beta = beta1)
    ####
    #Conduct the causal forest for testing:
    ncol_df <- ncol(df_train1)
    df_train_1 <- df_train1[,2:ncol_df]
    df_test_1 <- df_test1[,2:ncol_df]
    power_cover_null1[i] <- cfperm1(df_train_1, df_test_1, n_perm = B)
    power_cover_null1[i] <- cfperm2(df_train_1, df_test_1, n_perm = B)     
    ##########
    ncol_df <- ncol(df_train2)
    df_train_2 <- df_train2[,2:ncol_df]
    df_test_2 <- df_test2[,2:ncol_df]
    n_train <- nrow(df_train_2)
    n_test <- nrow(df_test_2)
    n_feature <- ncol(df_train_2) - 1
    power_cover_alt1[i] <- cfperm1(df_train_2, df_test_2, n_perm = B)
    power_cover_alt2[i] <- cfperm2(df_train_2, df_test_2, n_perm = B)    
  }
  power_df$power_null1[j] <- mean(power_cover_null1)
  power_df$power_null2[j] <- mean(power_cover_null2)
  power_df$power_alt1[j] <- mean(power_cover_alt1)
  power_df$power_alt2[j] <- mean(power_cover_alt2)
  write.csv(power_df, "power_df_modelA_CFPerm12.csv")
}

