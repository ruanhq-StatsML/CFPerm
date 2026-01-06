test_that("LM_generation returns correct shapes and names", {
  skip_if_not_installed("MASS")

  set.seed(1)
  out <- LM_generation(
    n = 30,
    beta_hat = c(1, -1, 0.5),
    mean_shift = 0,
    var_shift = 1,
    cor = 0.3,
    n_nuisance = 2,
    eps = 1
  )

  expect_true(is.list(out))
  expect_true(all(c("df_return", "X_return") %in% names(out)))

  df <- out$df_return
  X <- out$X_return

  expect_s3_class(df, "data.frame")
  expect_s3_class(X, "data.frame")

  # df columns: Y1 + p signal + n_nuis + Y
  expect_equal(nrow(df), 30)
  expect_equal(ncol(df), 1 + 3 + 2 + 1)

  expect_equal(colnames(df), c("Y1", "X1", "X2", "X3", "X_nuis1", "X_nuis2", "Y"))
  expect_equal(ncol(X), 3 + 2)
  expect_equal(colnames(X), c("X1", "X2", "X3", "X_nuis1", "X_nuis2"))

  expect_true(is.numeric(df$Y))
  expect_true(is.numeric(df$Y1))
})

test_that("Mars_generation returns correct shapes and names", {
  set.seed(1)
  df <- Mars_generation(
    n = 25,
    beta_hat = rep(1, 7),
    mean_shift_vec = rep(1, 8),
    var_shift_vec = rep(1, 8),
    n_nuisance = 3,
    eps = 0.5
  )

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 25)
  expect_equal(ncol(df), 1 + 8 + 3 + 1)
  expect_equal(
    colnames(df),
    c("Y1", paste0("X", 1:8), paste0("X_nuis", 1:3), "Y")
  )
})
