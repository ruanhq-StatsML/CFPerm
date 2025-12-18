#cfperm power evaluation for the auto-price dataset:

source("CFPerm_func.R")
source("data_generating_processes.R")
library(grf)
set.seed(42)

df <- read.csv("real_datasets/auto_price.csv", sep = ",")[,2:17]
colnames(df)[16] <- "Y"


power_list <- list()
#split with respect to rm, lstat, crim and nox
B <- 500
power_df <- data.frame(
  power_cfperm0 = numeric(4),
  power_cfperm1 = numeric(4)
)
power_list1 <- numeric(B)
power_list2 <- numeric(B)
for(i in 1:B){
  med <- median(df[,5])
  df_train <- df[which(df[,5] <= med), ]
  df_test <- df[which(df[,5] > med), ]
  power_list1[i] <- cfperm(df_train, df_test, q1 = 1, level = 0, n_perm = 500)$result
  power_list2[i] <- cfperm(df_train, df_test, q1 = 0.99, level = 0.05, n_perm = 500)$result   
  print(i)
}
power_df$power_cfperm0[1] <- mean(power_list1)
power_df$power_cfperm1[1] <- mean(power_list2)

for(i in 1:B){
  med <- median(df[,7])
  df_train <- df[which(df[,7] <= med), ]
  df_test <- df[which(df[,7] > med), ]
  power_list1[i] <- cfperm(df_train, df_test, q1 = 1, level = 0, n_perm = 500)$result
  power_list2[i] <- cfperm(df_train, df_test, q1 = 0.99, level = 0.05, n_perm = 500)$result   
}
power_df$power_cfperm0[2] <- mean(power_list1)
power_df$power_cfperm1[2] <- mean(power_list2)

for(i in 1:B){
  med <- median(df[,8])
  df_train <- df[which(df[,8] <= med), ]
  df_test <- df[which(df[,8] > med), ]
  power_list1[i] <- cfperm(df_train, df_test, q1 = 1, level = 0, n_perm = 500)$result
  power_list2[i] <- cfperm(df_train, df_test, q1 = 0.99, level = 0.05, n_perm = 500)$result   
}
power_df$power_cfperm0[3] <- mean(power_list1)
power_df$power_cfperm1[3] <- mean(power_list2)

for(i in 1:B){
  med <- median(df[,10])
  df_train <- df[which(df[,10] <= med), ]
  df_test <- df[which(df[,10] > med), ]
  power_list1[i] <- cfperm(df_train, df_test, q1 = 1, level = 0, n_perm = 500)$result
  power_list2[i] <- cfperm(df_train, df_test, q1 = 0.99, level = 0.05, n_perm = 500)$result   
}
power_df$power_cfperm0[4] <- mean(power_list1)
power_df$power_cfperm1[4] <- mean(power_list2)
