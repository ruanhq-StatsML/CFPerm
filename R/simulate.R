#' Linear-model data generator with correlated design
#'
#' Generates a response from a linear model with AR(1)-style correlation among
#' covariates, plus optional nuisance (noise) features.
#'
#' @param n Integer; number of samples.
#' @param beta_hat Numeric vector; coefficients for signal features.
#' @param mean_shift Numeric; mean shift applied to all signal covariates.
#' @param var_shift Numeric; variance scale applied to the covariance matrix.
#' @param cor Numeric in (-1, 1); AR(1) correlation parameter.
#' @param n_nuisance Integer; number of nuisance features.
#' @param eps Numeric > 0; noise standard deviation for the outcome.
#'
#' @return A list with:
#' \describe{
#'   \item{df_return}{Data frame containing Y1, signal covariates, nuisance covariates, and Y.}
#'   \item{X_return}{Data frame of all covariates (signal + nuisance), excluding Y1 and Y.}
#' }
#'
#' @importFrom MASS mvrnorm
#' @export
LM_generation <- function(n,
                          beta_hat,
                          mean_shift = 0,
                          var_shift = 1,
                          cor,
                          n_nuisance,
                          eps) {
  # Basic validation
  if (!is.numeric(n) || length(n) != 1 || n <= 0) stop("`n` must be a positive integer-like value.")
  n <- as.integer(n)

  if (!is.numeric(beta_hat) || length(beta_hat) < 1) stop("`beta_hat` must be a non-empty numeric vector.")
  p <- length(beta_hat)

  if (!is.numeric(mean_shift) || length(mean_shift) != 1) stop("`mean_shift` must be a numeric scalar.")
  if (!is.numeric(var_shift) || length(var_shift) != 1 || var_shift <= 0) stop("`var_shift` must be a positive numeric scalar.")
  if (!is.numeric(cor) || length(cor) != 1 || abs(cor) >= 1) stop("`cor` must be a numeric scalar with |cor| < 1.")
  if (!is.numeric(n_nuisance) || length(n_nuisance) != 1 || n_nuisance < 0) stop("`n_nuisance` must be a non-negative integer-like value.")
  n_nuisance <- as.integer(n_nuisance)
  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0) stop("`eps` must be a positive numeric scalar.")

  # Correlation / covariance
  corr_matrix <- matrix(0, p, p)
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      corr_matrix[i, j] <- cor^(abs(i - j))
    }
  }

  # If var_shift is intended as a variance multiplier, covariance = corr * var_shift
  Sigma <- corr_matrix * var_shift
  X_design <- MASS::mvrnorm(n, mu = rep(mean_shift, p), Sigma = Sigma)

  Y1 <- as.numeric(as.matrix(X_design) %*% matrix(beta_hat, nrow = p))
  random_error <- stats::rnorm(n, mean = 0, sd = eps)

  X_nuiss <- if (n_nuisance > 0) {
    matrix(stats::rnorm(n * n_nuisance, mean = 0, sd = 1), nrow = n, ncol = n_nuisance)
  } else {
    matrix(numeric(0), nrow = n, ncol = 0)
  }

  Y <- Y1 + random_error

  df_return <- data.frame(
    Y1 = Y1,
    X_design,
    X_nuiss,
    Y = Y,
    check.names = FALSE
  )

  colnames(df_return) <- c(
    "Y1",
    paste0("X", seq_len(p)),
    if (n_nuisance > 0) paste0("X_nuis", seq_len(n_nuisance)) else character(0),
    "Y"
  )

  ncol_df <- ncol(df_return)
  X_return <- df_return[, 2:(ncol_df - 1), drop = FALSE]

  return(list(df_return = df_return, X_return = X_return))
}

#' MARS-style nonlinear data generator (8 signal covariates)
#'
#' Generates covariates X1..X8 and a nonlinear response, plus optional nuisance features.
#'
#' @param n Integer; number of samples.
#' @param beta_hat Numeric vector of length 7; coefficients used in the response model.
#' @param mean_shift_vec Numeric vector of length 8; per-feature multiplicative mean shifts.
#' @param var_shift_vec Numeric vector of length 8; per-feature variance scales (via upper bound sqrt(var)).
#' @param n_nuisance Integer; number of nuisance features.
#' @param eps Numeric > 0; noise standard deviation for the outcome.
#'
#' @return Data frame with columns Y1, X1..X8, X_nuis1.., Y.
#' @export
Mars_generation <- function(n,
                            beta_hat,
                            mean_shift_vec = rep(1, 8),
                            var_shift_vec = rep(1, 8),
                            n_nuisance,
                            eps) {
  if (!is.numeric(n) || length(n) != 1 || n <= 0) stop("`n` must be a positive integer-like value.")
  n <- as.integer(n)

  if (!is.numeric(beta_hat) || length(beta_hat) < 7) stop("`beta_hat` must be numeric with length >= 7.")
  if (!is.numeric(mean_shift_vec) || length(mean_shift_vec) != 8) stop("`mean_shift_vec` must have length 8.")
  if (!is.numeric(var_shift_vec) || length(var_shift_vec) != 8 || any(var_shift_vec <= 0)) stop("`var_shift_vec` must have length 8 and be positive.")
  if (!is.numeric(n_nuisance) || length(n_nuisance) != 1 || n_nuisance < 0) stop("`n_nuisance` must be a non-negative integer-like value.")
  n_nuisance <- as.integer(n_nuisance)
  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0) stop("`eps` must be a positive numeric scalar.")

  X <- vector("list", 8)
  for (j in 1:8) {
    X[[j]] <- stats::runif(n, min = 0, max = sqrt(var_shift_vec[j])) * mean_shift_vec[j]
  }
  X1 <- X[[1]]; X2 <- X[[2]]; X3 <- X[[3]]; X4 <- X[[4]]
  X5 <- X[[5]]; X6 <- X[[6]]; X7 <- X[[7]]; X8 <- X[[8]]

  Y1 <- beta_hat[1] * exp(0.2 * X1) +
    beta_hat[2] * (1 / (1 + exp(-20 * (X2 - 0.5)))) +
    beta_hat[3] * X3 +
    beta_hat[4] * X4 +
    beta_hat[5] * X5 +
    beta_hat[6] * (X6)^2 +
    beta_hat[7] * sin(pi * X7 * X8)

  random_error <- stats::rnorm(n, mean = 0, sd = eps)
  Y <- Y1 + random_error

  X_nuiss <- if (n_nuisance > 0) {
    matrix(stats::rnorm(n * n_nuisance, mean = 0, sd = 1), nrow = n, ncol = n_nuisance)
  } else {
    matrix(numeric(0), nrow = n, ncol = 0)
  }

  df_return <- data.frame(
    Y1 = Y1,
    X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5, X6 = X6, X7 = X7, X8 = X8,
    X_nuiss,
    Y = Y,
    check.names = FALSE
  )

  colnames(df_return) <- c(
    "Y1",
    paste0("X", 1:8),
    if (n_nuisance > 0) paste0("X_nuis", seq_len(n_nuisance)) else character(0),
    "Y"
  )

  df_return
}
