#CFPerm_Unittest1

library(nnet)
library(glmnet)
library(xgboost)
library(randomForest)
library(grf)
library(testthat)

test_that('cfperm returns TRUE/FALSE and validates inputs', 
  {
  skip_if_not_installed(c('grf', 'nnet', 'randomForest', 'glmnet', 'xgboost'))
  set.seed(123)
  n0 <- 200
  n1 <- 200
  p <- 10
  X <- matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)
  colnames(X) <- paste0("X", seq_len(p))
  T <- c(rep(0, n0), rep(1, n1))
  Y <- 1 + 0.8 * X[, 1] - 0.5 * X[, 2] + rnorm(n, sd = 1)
  res_permu <- cfperm(
  X = X,
  Y = Y,
  T = T,
  n_perm = 30,
  vimp = "permuCATE",
  seed = 2026,
  level_feature = 0.05,
  level_across_feature = 0.05,
  top_k = 1)
  res_loco <- cfperm(
  X = X,
  Y = Y,
  T = T,
  n_perm = 30,
  vimp = "loco",
  seed = 2026,
  level_feature = 0.05,
  level_across_feature = 0.05,
  top_k = 1)
  res_grf <- cfperm(
  X = X,
  Y = Y,
  T = T,
  n_perm = 200,
  vimp = 'GRF',
  seed = 2020,
  level_feature = 0.005,
  level_across_feature = 0.005)
  expect_true(res_permu %in% c(0L, 1L))
  expect_true(res_loco %in% c(0L, 1L))
  expect_true(res_grf %in% c(0L, 1L))
})











