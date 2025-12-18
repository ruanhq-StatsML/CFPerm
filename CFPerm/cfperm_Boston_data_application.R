source("CFPerm_func.R")
source("data_generating_processes.R")
library(grf)
set.seed(42)
library(mltools)
library(caret)

df <- read.csv("real_datasets/Boston.csv")[, 2:15]
colnames(df)[14] <- "Y"
power_list <- list()
#split with respect to rm, lstat, crim and nox to perform the hypothesis test:
########
mc <- 500
result_matrix <- matrix(0, mc, 2)
for(i in 1:mc){
  med <- median(df$rm)
  df_train <- df[which(df$rm <= med), ]
  df_test <- df[which(df$rm > med), ]
  df_train1 <- df_train[sample(c(1:nrow(df_train)), nrow(df_train), replace = TRUE), ]
  df_test1 <- df_test[sample(c(1:nrow(df_test)), nrow(df_test), replace = TRUE), ]
  result_matrix[i, 1] <- cfperm(df_train1, df_test1, q1 = 1, level = 0, n_perm = 500)$result
  result_matrix[i, 2] <- cfperm(df_train1, df_test1, q1 = 0.99, level = 0.05, n_perm = 500)$result
}
power_list[["rm"]] <- result_matrix

result_matrix <- matrix(0, mc, 2)
for(i in 1:mc){
  med <- median(df$lstat)
  df_train <- df[which(df$lstat <= med), ]
  df_test <- df[which(df$lstat > med), ]
  df_train1 <- df_train[sample(c(1:nrow(df_train)), nrow(df_train), replace = TRUE), ]
  df_test1 <- df_test[sample(c(1:nrow(df_test)), nrow(df_test), replace = TRUE), ]
  result_matrix[i, 1] <- cfperm(df_train1, df_test1, q1 = 1, level = 0, n_perm = 500)$result
  result_matrix[i, 2] <- cfperm(df_train1, df_test1, q1 = 0.99, level = 0.05, n_perm = 500)$result
}
power_list[["lstat"]] <- result_matrix

result_matrix <- matrix(0, mc, 2)
for(i in 1:mc){
  med <- median(df$crim)
  df_train <- df[which(df$crim <= med), ]
  df_test <- df[which(df$crim > med), ]
  df_train1 <- df_train[sample(c(1:nrow(df_train)), nrow(df_train), replace = TRUE), ]
  df_test1 <- df_test[sample(c(1:nrow(df_test)), nrow(df_test), replace = TRUE), ]
  result_matrix[i, 1] <- cfperm(df_train1, df_test1, q1 = 1, level = 0, n_perm = 500)$result
  result_matrix[i, 2] <- cfperm(df_train1, df_test1, q1 = 0.99, level = 0.05, n_perm = 500)$result
}
power_list[["crim"]] <- result_matrix

result_matrix <- matrix(0, mc, 2)
for(i in 1:mc){
  med <- median(df$nox)
  df_train <- df[which(df$nox <= med), ]
  df_test <- df[which(df$nox > med), ]
  df_train1 <- df_train[sample(c(1:nrow(df_train)), nrow(df_train), replace = TRUE), ]
  df_test1 <- df_test[sample(c(1:nrow(df_test)), nrow(df_test), replace = TRUE), ]
  result_matrix[i, 1] <- cfperm(df_train1, df_test1, q1 = 1, level = 0, n_perm = 500)$result
  result_matrix[i, 2] <- cfperm(df_train1, df_test1, q1 = 0.99, level = 0.05, n_perm = 500)$result
}
power_list[["nox"]] <- result_matrix

saveRDS(power_list, "Boston_cfperm_procedure_comparison.rds")
