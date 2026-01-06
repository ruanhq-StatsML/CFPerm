test_that("cfperm returns 0/1 and validates inputs", {
  skip_if_not_installed("grf")

  set.seed(123)

  # small synthetic dataset where Y is present
  n <- 40
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1] + 0.1 * rnorm(n)

  df_all <- data.frame(X, Y = Y, check.names = FALSE)
  colnames(df_all) <- c(paste0("X", 1:p), "Y")

  df_train <- df_all[1:20, , drop = FALSE]
  df_test  <- df_all[21:40, , drop = FALSE]

  res <- cfperm(
    df_train = df_train,
    df_test = df_test,
    n_perm = 3,
    num.trees = 20,
    seed = 999
  )

  expect_true(res %in% c(0L, 1L))

  # input validation
  expect_error(cfperm(df_train[, 1:p, drop = FALSE], df_test)) # missing Y
})
