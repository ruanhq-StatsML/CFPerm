test_that("aggregate_power returns expected indicator on toy example", {
  # 3 features, 4 perms + 1 original => 5 cols
  # Construct so that feature 2 original is clearly above perm quantiles
  df <- rbind(
    c(0.1, 0.2, 0.15, 0.18, 0.05),   # feature1 original small
    c(0.05, 0.08, 0.07, 0.06, 1.00), # feature2 original large
    c(0.2, 0.25, 0.22, 0.24, 0.10)   # feature3 original small
  )

  res <- aggregate_power(
    df = df,
    top_k = 1,
    level_feature = 0.25,          # upper quantile among perms
    level_across_feature = 0.25,   # global threshold among features
    normalization = FALSE
  )

  expect_true(res %in% c(0L, 1L))
  expect_equal(res, 1L)
})

test_that("aggregate_power errors on invalid input", {
  expect_error(aggregate_power(matrix(1, nrow = 2, ncol = 1)))
})
