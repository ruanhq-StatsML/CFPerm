#This script evaluate the power for covariate shift only:
source("CFPerm_func.R")
source("data_generating_processes.R")
library(grf)
set.seed(42)

B <- 500
#Variance Shift on top of No Mean Shift:
dif_seq <- seq(0, 1.5, 0.05)
power_df <- data.frame(
  dif_seq = dif_seq,
  power1 = numeric(31),
  power2 = numeric(31)
  )
for(j in 1:31){
  power_list1 = numeric(B)
  power_list2 = numeric(B)
  for(i in 1:B){
    dif <- dif_seq[j]
    var_seq_train <- c(1 + 5 * dif, 1 + 4 * dif, 1 + 3 * dif, 1 + 2 * dif, 1 + dif)
    var_seq_test <- rev(var_seq_train)
    df1 <- mars_model_df_normal(n = 1000, p_nuisance = 15, sigma = 1, mean_seq = c(1,1,1,1,1),
      var_seq = var_seq_train)
    df2 <- mars_model_df_normal(n = 1000, p_nuisance = 15, sigma = 1, mean_seq = c(1,1,1,1,1),
      var_seq = var_seq_test)
    power_list1[i] <- cfperm(df1, df2, q1 = 1, level = 0, n_perm = 500)$result
    power_list2[i] <- cfperm(df1, df2, q1 = 0.99, level = 0.95, n_perm = 500)$result
  }
  power_df$power1[j] <- mean(power_list1)
  power_df$power2[j] <- mean(power_list2) 
  write.csv(power_df, "cfperm_covariateshiftonly1.csv")
}

#Variance Shift on top of Small Mean Shift:
dif_seq <- seq(0, 1.5, 0.05)
power_df <- data.frame(
  dif_seq = dif_seq,
  power1 = numeric(31),
  power2 = numeric(31)
  )
for(j in 1:31){
  power_list1 = numeric(B)
  power_list2 = numeric(B)
  for(i in 1:B){
    dif <- dif_seq[j]
    var_seq_train <- c(1 + 5 * dif, 1 + 4 * dif, 1 + 3 * dif, 1 + 2 * dif, 1 + dif)
    var_seq_test <- rev(var_seq_train)
    df1 <- mars_model_df_normal(n = 1000, p_nuisance = 15, sigma = 1, mean_seq = c(1,1,1,1,1),
      var_seq = var_seq_train)
    df2 <- mars_model_df_normal(n = 1000, p_nuisance = 15, sigma = 1, mean_seq = c(1.05, 1.05, 1.05, 1.05, 1.05),
      var_seq = var_seq_test)
    power_list1[i] <- cfperm(df1, df2, q1 = 1, level = 0, n_perm = 500)$result
    power_list2[i] <- cfperm(df1, df2, q1 = 0.99, level = 0.95, n_perm = 500)$result
  }
  power_df$power1[j] <- mean(power_list1)
  power_df$power2[j] <- mean(power_list2) 
  write.csv(power_df, "cfperm_covariateshiftonly2.csv")
}

#Variance Shift on top of Moderate Mean Shift:
dif_seq <- seq(0, 1.5, 0.05)
power_df <- data.frame(
  dif_seq = dif_seq,
  power1 = numeric(31),
  power2 = numeric(31)
  )
for(j in 1:31){
  power_list1 = numeric(B)
  power_list2 = numeric(B)
  for(i in 1:B){
    dif <- dif_seq[j]
    var_seq_train <- c(1 + 5 * dif, 1 + 4 * dif, 1 + 3 * dif, 1 + 2 * dif, 1 + dif)
    var_seq_test <- rev(var_seq_train)
    df1 <- mars_model_df_normal(n = 1000, p_nuisance = 15, sigma = 1, mean_seq = c(1,1,1,1,1),
      var_seq = var_seq_train)
    df2 <- mars_model_df_normal(n = 1000, p_nuisance = 15, sigma = 1, mean_seq = c(1.2,1.2,1.2,1.2,1.2),
      var_seq = var_seq_test)
    power_list1[i] <- cfperm(df1, df2, q1 = 1, level = 0, n_perm = 500)$result
    power_list2[i] <- cfperm(df1, df2, q1 = 0.99, level = 0.95, n_perm = 500)$result
  }
  power_df$power1[j] <- mean(power_list1)
  power_df$power2[j] <- mean(power_list2) 
  write.csv(power_df, "cfperm_covariateshiftonly3.csv")
}

#Variance Shift on top of Large Mean Shift:
dif_seq <- seq(0, 1.5, 0.05)
power_df <- data.frame(
  dif_seq = dif_seq,
  power1 = numeric(31),
  power2 = numeric(31)
  )
for(j in 1:31){
  power_list1 = numeric(B)
  power_list2 = numeric(B)
  for(i in 1:B){
    dif <- dif_seq[j]
    var_seq_train <- c(1 + 5 * dif, 1 + 4 * dif, 1 + 3 * dif, 1 + 2 * dif, 1 + dif)
    var_seq_test <- rev(var_seq_train)
    df1 <- mars_model_df_normal(n = 1000, p_nuisance = 15, sigma = 1, mean_seq = c(1,1,1,1,1),
      var_seq = var_seq_train)
    df2 <- mars_model_df_normal(n = 1000, p_nuisance = 15, sigma = 1, mean_seq = c(2,2,2,2,2),
      var_seq = var_seq_test)
    power_list1[i] <- cfperm(df1, df2, q1 = 1, level = 0, n_perm = 500)$result
    power_list2[i] <- cfperm(df1, df2, q1 = 0.99, level = 0.95, n_perm = 500)$result
  }
  power_df$power1[j] <- mean(power_list1)
  power_df$power2[j] <- mean(power_list2) 
  write.csv(power_df, "cfperm_covariateshiftonly4.csv")
}


#Mean Shift on top of no Variance Shift
dif_seq <- seq(0, 1.5, 0.05)
power_df <- data.frame(
  dif_seq = dif_seq,
  power1 = numeric(31),
  power2 = numeric(31)
  )
for(j in 1:31){
  power_list1 = numeric(B)
  power_list2 = numeric(B)
  for(i in 1:B){
    dif <- dif_seq[j]
    mean_seq_train <- c(1, 1, 1, 1, 1)
    mean_seq_test <- c(1 + 5 * dif, 1 + 4 * dif, 1 + 3 * dif, 1 + 2 * dif, 1 + dif)
    df1 <- mars_model_df_normal(n = 1000, p_nuisance = 15, sigma = 1, mean_seq = mean_seq_train,
      var_seq = c(1, 1, 1, 1, 1))
    #df2_1 <- correlation_model_df(n = 2000, p = 20, snr = 0.8, rho = 0.3)
    df2 <- mars_model_df_normal(n = 500, p_nuisance = 15, sigma = 1, mean_seq = mean_seq_test,
      var_seq = c(1, 1, 1, 1, 1))
    power_list1[i] <- cfperm(df1, df2, q1 = 1, level = 0, n_perm = 500)$result
    power_list2[i] <- cfperm(df1, df2, q1 = 0.99, level = 0.95, n_perm = 500)$result
  }
  power_df$power1[j] <- mean(power_list1)
  power_df$power2[j] <- mean(power_list2) 
  write.csv(power_df, "cfperm_covariateshiftonly5.csv")
}


#Mean Shift on top of Moderate Variance Shift
dif_seq <- seq(0, 1.5, 0.05)
power_df <- data.frame(
  dif_seq = dif_seq,
  power1 = numeric(31),
  power2 = numeric(31)
  )
for(j in 1:31){
  power_list1 = numeric(B)
  power_list2 = numeric(B)
  for(i in 1:B){
    dif <- dif_seq[j]
    mean_seq_train <- c(1, 1, 1, 1, 1)
    mean_seq_test <- c(1 + 5 * dif, 1 + 4 * dif, 1 + 3 * dif, 1 + 2 * dif, 1 + dif)
    df1 <- mars_model_df_normal(n = 1000, p_nuisance = 15, sigma = 1, mean_seq = mean_seq_train,
      var_seq = c(6,5,4,3,2))
    #df2_1 <- correlation_model_df(n = 2000, p = 20, snr = 0.8, rho = 0.3)
    df2 <- mars_model_df_normal(n = 500, p_nuisance = 15, sigma = 1, mean_seq = mean_seq_test,
      var_seq = c(2,3,4,5,6))
    power_list1[i] <- cfperm(df1, df2, q1 = 1, level = 0, n_perm = 500)$result
    power_list2[i] <- cfperm(df1, df2, q1 = 0.99, level = 0.95, n_perm = 500)$result
  }
  power_df$power1[j] <- mean(power_list1)
  power_df$power2[j] <- mean(power_list2) 
  write.csv(power_df, "cfperm_covariateshiftonly6.csv")
}
