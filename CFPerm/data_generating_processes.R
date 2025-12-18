#This script include all of the data-generating processes:
library(glmnet)
library(ranger)
library(densratio)
library(knockoff)
library(MASS)
library(ranger)
library(randomForest)
library(pROC)
library(radiant.model)
library(tableone)
library(LaplacesDemon)
mars_model_df <- function(n, p_nuisance, sigma){
  X1 <- runif(n, 0, 1)
  X2 <- runif(n, 0, 1)
  X3 <- runif(n, 0, 1)
  X4 <- runif(n, 0, 1)
  X5 <- runif(n, 0, 1)
  Y1 <- 0.1 * exp(4 * X1) + 
    4/(1+exp(-20 * (X2 - 0.5))) + 3 * X3 +
    2 * X4 + X5
  random_noise <- rnorm(n, 0, sigma)
  Y <- Y1 + random_noise
  if(p_nuisance > 0){
    #nuisance parameters with standard normal random noise:
    X_nuiss <- matrix(NA, n, p_nuisance)
    for(i in 1:p_nuisance){
      X_nuiss[,i] <- rnorm(n, 0, 1)
    }
    data_mat <- cbind(X1, X2, X3, X4, X5, X_nuiss, Y)
    colnames(data_mat) <-
      c(paste("X", c(1:(5+as.numeric(p_nuisance))), sep = ""), "Y")
  }
  else{
    data_mat <- cbind(X1, X2, X3, X4, X5, Y)
    colnames(data_mat) <- 
      c(paste("X", c(1:5), sep = ""), "Y")
  }
  return(data.frame(data_mat))
}
mars_model_df_normal <- function(n, p_nuisance, sigma, mean_seq = c(0,0,0,0,0), var_seq = c(1,1,1,1,1)){
  X1 <- rnorm(n, mean_seq[1], var_seq[1])
  X2 <- rnorm(n, mean_seq[2], var_seq[2])
  X3 <- rnorm(n, mean_seq[3], var_seq[3])
  X4 <- rnorm(n, mean_seq[4], var_seq[4])
  X5 <- rnorm(n, mean_seq[5], var_seq[5])
  #Y1 <- 0.1 * exp(4 * X1) +
  #  4/(1+exp(-20 * (X2 - 0.5))) + 3 * X3 +
  #  2 * X4 + X5
  Y1 <- 50 * sin(pi * X1 * X2) + 5 * (X3 - 0.05)^2 + 5 * X4 + 5 * X5
  random_noise <- rnorm(n, 0, sigma)
  Y <- Y1 + random_noise
  if(p_nuisance > 0){
    #nuisance parameters with standard normal random noise:
    X_nuiss <- matrix(NA, n, p_nuisance)
    for(i in 1:p_nuisance){
      X_nuiss[,i] <- rnorm(n, 0, 1)
    }
    data_mat <- cbind(X1, X2, X3, X4, X5, X_nuiss, Y)
    colnames(data_mat) <-
      c(paste("X", c(1:(5+as.numeric(p_nuisance))), sep = ""), "Y")
  }
  else{
    data_mat <- cbind(X1, X2, X3, X4, X5, Y)
    colnames(data_mat) <-
      c(paste("X", c(1:5), sep = ""), "Y")
  }
  return(data.frame(data_mat))
}
correlation_model_df <- function(n, p, snr, rho = 0.35){
  corr_matrix <- matrix(NA, p, p)
  for(i in 1:p){
    for(j in 1:p){
      corr_matrix[i, j] <- rho ^ (abs(i - j))
    }
  }
  X_design <- mvrnorm(n, mu = rep(0, p), corr_matrix)
  beta <- as.matrix(rep(1, p))
  Y1 <- X_design %*% beta
  sigma_noise_square <- t(beta) %*% corr_matrix %*% beta / snr
  Y <- Y1 + rnorm(n, 0, sqrt(sigma_noise_square))
  data_mat <- cbind(X_design, Y)
  colnames(data_mat) <- c(paste("X", c(1:as.numeric(p)), sep = ""), "Y")
  return(data.frame(data_mat))
}
#Linear model data-generating process:
LM_generation <- function(n, beta_hat, cor, n_nuisance, eps){
  p <- length(beta_hat)
  corr_matrix <- matrix(NA, p, p)
  for(i in 1:p){
    for(j in 1:p){
      corr_matrix[i, j] <- cor ^ (abs(i - j))
    }
  }
  X_design <- mvrnorm(n, mu = rep(0, p), corr_matrix)
  Y1 <- as.matrix(X_design) %*% as.matrix(beta_hat, nrow = p)
  random_error <- rnorm(n, 0, eps)
  #adding the nuisance random noise features:
  X_nuiss <- matrix(0, n, n_nuisance)
  for(i in 1:n_nuisance){
    X_nuiss[,i] <- rnorm(n, 0, 1)
  }
  Y <- Y1 + random_error
  df_return <- data.frame(cbind(X_design, X_nuiss, Y))
  ncol_df <- ncol(df_return)
  colnames(df_return) <- c(paste("X", c(1:p), sep = ""), paste("X_nuis", c(1:n_nuisance), sep = ""), "Y")
  X_return <- df_return[,1:(ncol_df - 1)]
  return(list(df_return = df_return, X_return = X_return))
}
#Model A/B/C from Hu_JASA_2023:
model_A_DGP <- function(n, alpha, train = TRUE, beta){
  if(train == TRUE){
    X1 <- mvrnorm(n, mu = c(0,0,0,0,0), Sigma = diag(c(1,1,1,1,1)))
  }
  else{
    X1 <- mvrnorm(n, mu = c(1,1,-1,-1,0), Sigma = diag(c(1,1,1,1,1)))
  }
  Y <- X1 %*% beta + rnorm(n, 0, 1) + alpha
  df_mat <- data.frame(cbind(X1, Y))
  colnames(df_mat) <- c("X1", "X2", "X3", "X4", "X5", "Y")
  return(df_mat)
}
#Model B from Hu_JASA_2023
model_B_DGP <- function(n, alpha, train = TRUE, beta){
  if(train == TRUE){
    X <- 0.5 *  mvrnorm(n, mu = c(0,0,0,0,0), Sigma = diag(c(1,1,1,1,1))) + 
      0.5 * mvrnorm(n, mu = c(0.5, 0.5, -0.5, -0.5, 0), Sigma = diag(c(1,1,1,1,1)))
  }
  else{
    X <- 0.5 *  mvrnorm(n, mu = c(0,0,0,0,0), Sigma = diag(c(1,1,1,1,1))) + 
      0.5 * mvrnorm(n, mu = c(0,0,0,0,0), Sigma = 1.5 * diag(c(1,1,1,1,1)))
  }
  Y <- alpha + X[,1] + X[,2] + X[,3]^2 + X[,4]^2 + X[,5]^3 + rt(n, 5)
  df_mat <- data.frame(cbind(X, Y))
  colnames(df_mat) <- c("X1", "X2", "X3", "X4", "X5", "Y")
  return(df_mat)
}