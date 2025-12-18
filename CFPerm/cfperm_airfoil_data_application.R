source("CFPerm_func.R")
source("data_generating_processes.R")
library(grf)
set.seed(42)

air_foil_dt <- read.table("airfoil_self_noise.dat", header = FALSE)
colnames(air_foil_dt) <- c("X1", "X2", "X3", "X4", "X5", "Y")

set.seed(42)
power_df <- data.frame(
  power_cfperm0 = numeric(5),
  power_cfperm1 = numeric(5)
)
power_list1 <- numeric(B)
power_list2 <- numeric(B)

for(i in 1:200){
  n <- nrow(air_foil_dt)
  n_col <- ncol(air_foil_dt)
  sample_idx <- sample(c(1:n), 751, replace = FALSE)
  df_train <- air_foil_dt[sample_idx, ]
  rownames(df_train) <- c(1:751)
  df_test <- air_foil_dt[-sample_idx, ]
  rownames(df_test) <- c(1:752)
  power_list1[i] <- cfperm(df_train, df_test, q1 = 1, level = 0, n_perm = 500)$result
  power_list2[i] <- cfperm(df_train, df_test, q1 = 0.99, level = 0.05, n_perm = 500)$result   
}
power_df$power_cfperm0[1] <- mean(power_list1)
power_df$power_cfperm1[1] <- mean(power_list2)


power_list1 <- numeric(B)
power_list2 <- numeric(B)
for(i in 1:B){
  n <- nrow(air_foil_dt)
  n_col <- ncol(air_foil_dt)
  #########
  med <- quantile(air_foil_dt[,3], 0.5)
  df_train <- air_foil_dt[air_foil_dt[,3] <= med, ]
  df_test <- air_foil_dt[air_foil_dt[,3] > med,]  
  power_list1[i] <- cfperm(df_train, df_test, q1 = 1, level = 0, n_perm = 500)$result
  power_list2[i] <- cfperm(df_train, df_test, q1 = 0.99, level = 0.05, n_perm = 500)$result  
}
power_df$power_cfperm0[2] <- mean(power_list1)
power_df$power_cfperm1[2] <- mean(power_list2)


power_list1 <- numeric(B)
power_list2 <- numeric(B)
for(i in 1:B){
  n <- nrow(air_foil_dt)
  n_col <- ncol(air_foil_dt)
  #########
  med <- quantile(air_foil_dt[,4], 0.5)
  df_train <- air_foil_dt[air_foil_dt[,4] <= med, ]
  df_test <- air_foil_dt[air_foil_dt[,4] > med,]
  power_list1[i] <- cfperm(df_train, df_test, q1 = 1, level = 0, n_perm = 500)$result
  power_list2[i] <- cfperm(df_train, df_test, q1 = 0.99, level = 0.05, n_perm = 500)$result  
}
power_df$power_cfperm0[3] <- mean(power_list1)
power_df$power_cfperm1[3] <- mean(power_list2)


power_list1 <- numeric(B)
power_list2 <- numeric(B)
for(i in 1:B){
  #########
  med <- quantile(air_foil_dt[,6], 0.5)
  df_train <- air_foil_dt[air_foil_dt[,6] <= med, ]
  df_test <- air_foil_dt[air_foil_dt[,6] > med,] 
  power_list1[i] <- cfperm(df_train, df_test, q1 = 1, level = 0, n_perm = 500)$result
  power_list2[i] <- cfperm(df_train, df_test, q1 = 0.99, level = 0.05, n_perm = 500)$result
}
power_df$power_cfperm0[4] <- mean(power_list1)
power_df$power_cfperm1[4] <- mean(power_list2)


power_list1 <- numeric(B)
power_list2 <- numeric(B)
for(i in 1:B){
  n <- nrow(air_foil_dt)
  n_col <- ncol(air_foil_dt)
  sample_idx <- sample(c(1:n), 301, replace = FALSE)
  df_train <- air_foil_dt[sample_idx, ]
  n_train <- nrow(df_train)
  rownames(df_train) <- c(1:301)
  df_test1 <- air_foil_dt[-sample_idx, ]
  rownames(df_test1) <- c(1:nrow(df_test1))
  sample_prob <- apply(df_test1, 1, function(x){
    exp(-x[1] + x[5])
    })
  sample_prob <- sample_prob/sum(sample_prob)
  sample_idx_test <- sample(c(1:nrow(df_test1)), 301, replace = FALSE)
  df_test <- df_test1[sample_idx_test,]
  power_list1[i] <- cfperm(df_train, df_test, q1 = 1, level = 0, n_perm = 500)$result
  power_list2[i] <- cfperm(df_train, df_test, q1 = 0.99, level = 0.05, n_perm = 500)$result  
}
power_df$power_cfperm0[5] <- mean(power_list1)
power_df$power_cfperm1[5] <- mean(power_list2)

write.csv(power_df, "power_cfperm_airfoil_application.csv")

