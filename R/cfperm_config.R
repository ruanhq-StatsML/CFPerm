#Causal Inference Permutation Test here:


#' Configuration for PermuCATE-based CFPerm
#'
#' @param clip_e Propensity-score clipping value.
#' @param n_permutation Number of residual permutations used inside PermuCATE VIMP.
#' @param lambda Regularization or mixing parameter used by the PermuCATE routine.
#' @param model_mu outcome model.
#' @param model_nu conditional covariate model for residual permutation.
#' @param model_tau CATE/pseudo-outcome regression model.
#'
#' @return A PermuCATE configuration list
permucate_config <- function(
  clip_e = 0.05,
  n_permutation = 30,
  lambda = 0.25,
  model_mu = 'xgb_regression',
  model_nu = 'xgb_regression',
  model_tau = 'xgb_regression'
){
  cfg <- list(
    method = 'permuCATE',
    clip_e = clip_e,
    n_permutation = n_permutation,
    lambda = lambda,
    model_mu = model_mu,
    model_nu = model_nu,
    model_tau = model_tau
  )
  class(cfg) <- c('permucate_config', 'cfperm_vimp_config')
  cfg
}

#' Configuration for LOCO-based CFPerm
#'
#' @param clip_e Propensity-score clipping value.
#' @param n_permutation Number of residual permutations used inside PermuCATE VIMP.
#' @param lambda Regularization or mixing parameter used by the PermuCATE routine.
#' @param model_mu outcome model.
#' @param model_nu conditional covariate model for residual permutation.
#' @param model_tau CATE/pseudo-outcome regression model.
#'
#' @return A PermuCATE configuration list
loco_config <- function(
  clip_e = 0.05, clip_wtilde = 1e-3,
  lambda = 0.25, n_splits = 5,
  m_model = 'xgb_regression',
  e_model = 'rf_classification'
){
  cfg <- list(
    method = 'LOCO',
    clip_e = clip_e,
    clip_wtilde = clip_wtilde,
    lambda = lambda,
    n_splits = n_splits,
    m_model = m_model,
    e_model = e_model
  )
  class(cfg) <- c('loco_config', 'cfperm_vimp_config')
  cfg
}


#' Configuration for grf based CFPerm
#'
#' @param n_estimators Number of trees used in the causal forest.
#' @param n_bootstrap  Number of bootstrap repetitions for uncertainty estimation(setting as 0 means no repeatance)
#' @param ci_level     Confidence level used when bootstrap inference is requested.
#'
#' @return A GRF configuration list
grf_config <- function(
  n_estimators = 200, n_bootstrap = 0,
  ci_level = 0.95){
  cfg <- list(method = 'GRF',
    n_estimators = n_estimators,
    n_bootstrap = n_bootstrap,
    ci_level = ci_level)
  class(cfg) <- c('grf_config', 'cfperm_vimp_config')
  cfg
}

























