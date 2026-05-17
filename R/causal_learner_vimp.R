#LOCO variable importance:
source('model_registry.R')
########################################################################################
#' Leave-one-Covariate-Out Variable Importance for the R-risk:
#'@param X,Y,W:               Numeric matrix/data.frame of dimension (n_exist + n_new) * (p + 2), p is the number of features, with the response column(the second last column) and the batch assignment column(the last column)
#'                            For the batch assignment column T, the first n_exist value represent batch 0 and the next n_new value represent batch 1        
#'@param outcome_model:       Type of the outcome model
#'@param n_folds:             Number of folds of cross-fitting
#'@param propensity_model:    Type of the propensity score model
#'@param binary_outcome:      Whether the response(Y) is binary or not
#'@param clip_e:              The clip for the propensity score here.
#'@param clip_wtilde:         The clip for the weights in the R-learner
#'@return vimp:               Numeric vector-list of (p) for each of the feature
########################################################################################
importance_LOCO_r_risk <- function(
X, Y, W, seed = 0, n_folds = 5, binary_outcome = FALSE, lambda = 0.1,
clip_e = 0.01, clip_wtilde = 1e-3){
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  W <- as.numeric(W)
  nuisance_fit <- cross_nuisance_fit(X, Y, W,
    seed = seed, n_fold = n_folds, binary_outcome = binary_outcome)
  Y_tilde = Y - nuisance_fit$m_hat
  W_tilde = W - nuisance_fit$e_hat
  tau_full_model <- fit_tau_rlearner_weighted(X, Y_tilde, W_tilde,
    ridge_lambda = lambda, 
    seed = seed, clip_wtilde = clip_wtilde
  )
  #outcome model with the regression:
  tau_full <- tau_full_model$predict(X)
  risk_full <- r_risk(Y_tilde, W_tilde, tau_full)
  p = ncol(X)
  imps = rep(0, p)
  for(j in seq_len(p)){
    X_minus = X[, -j, drop = FALSE]
    tau_j_model <- fit_tau_rlearner_weighted(X_minus, Y_tilde, W_tilde,
      ridge_lambda = lambda,
      seed = 2026 + seed + j, clip_wtilde = clip_wtilde
    )
    tau_j <- tau_j_model$predict(X_minus)
    risk_j <- r_risk(Y_tilde, W_tilde, tau_j)
    imps[j] <- risk_j - risk_full
  }
  return(abs(imps)/sum(abs(imps)))
}

#r_risk:
#'@param Y_tilde       :Y-m(X)
#'@param W_tilde       :(T-e(X))
#'@param tau           :tau(X)
r_risk <- function(Y_tilde, W_tilde, tau){
  mean((Y_tilde - W_tilde * tau)^2)
}


########################################################################################
#' R-learner fitting via ridge regression
#'
########################################################################################
fit_tau_rlearner_weighted <- function(X, Y_tilde, W_tilde, ridge_lambda,
  clip_wtilde = 1e-4, seed = 2026){
  mask <- abs(W_tilde) > clip_wtilde
  if(!is.null(seed)) set.seed(as.integer(seed))
  X0 = X[mask, , drop = FALSE]
  z <- Y_tilde[mask]/W_tilde[mask]
  w <- (W_tilde[mask]^2)
  #standardize the
  mu <- colMeans(X0)
  sd <- apply(X0, 2, sd)
  sd[sd == 0 | !is.finite(sd)] <- 1
  #Standardization here:
  Xs <- sweep(sweep(X0, 2, mu, '-'), 2, sd, '/')
  X_design <- cbind(rep(1, nrow(Xs)), Xs)
  #sqrt(w)
  sw <- sqrt(w)#no major difference here.
  Xw <- X_design * sw
  yw <- z * sw#Weighted Ridge Regression here
  #Xw %*% t(Xw)
  XtX <- t(Xw) %*% Xw
  #Xw %*% yw
  Xty <- t(Xw) %*% yw
  #Ridge Regression:
  diag_mat <- diag(c(0, rep(1, ncol(X_design) - 1L)))#diag(c(0, rep(1, ncol(X_design) - 1L)))
  beta <- solve(XtX + ridge_lambda * diag_mat, Xty)
  pred_fun <- function(Xnew){
    Xnew <- as.matrix(Xnew)
    #do the standardization here, it's more efficient than the scale here.
    Xs_new <- sweep(sweep(Xnew, 2, mu, '-'), 2, sd, '/')
    as.numeric(cbind(1, Xs_new) %*% beta)
  }
  list(coef = beta, mu = mu, sd = sd, predict = pred_fun)
}

#PO risk:
po_risk <- function(tau_pred, Y, T, mu0, mu1, e, clip_eps = 1e-3){
  pi_clip <- pmax(pmin(e, 1 - clip_eps), clip_eps)
  tau_po <- (Y - ifelse(T == 1, mu1, mu0)) * (T - pi_clip)/(pi_clip * (1 - pi_clip))
  return(mean((tau_po - tau_pred)^2))
}

########################################################################################################
#' PermuCATE Variable Importance for the PO-risk, shuffle the residual for the covariate fitting model
#' Then return the variable importances and the combined permutation score:
#' @param X:             Covariate Space
#' @param Y:             Response
#' @param T:             Treatment Batch Label
#' @param lambda:        L2 Penalty for the Regularized Logistic Regression
#' @param model_mu:      Outcome Model
#' @param model_nu:      Conditional Model(X_j ~ X_{-j}) for continuous covariates X_j
#' @param model_nu_bina: Conditional Model(X_j ~ X_{-j}) for binary covariates X_j
#' @param model_nu_disc: Conditional Model(X_j ~ X_{-j}) for discrete covariates X_j
#' @param model_tau:     Pseudo-Outcome Model
#' @param clip_e:        The clipping value for the propensity score
########################################################################################################
permucate_fit <- function(X, Y, T, lambda = 0.01,
  model_mu, model_nu, 
  model_nu_bina, model_nu_disc,
  model_tau, 
  seed = 2026, clip_e = 1e-3) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  idx0 <- which(T == 0)
  idx1 <- which(T == 1)
  if (length(idx0) == 0 || length(idx1) == 0) {
    stop("PermuCATE requires both A=0 and A=1 in combined data")
  }
  infer_type <- function(X){
    p <- ncol(X)
    feature_type = rep("", p)
    for(j in 1:p){
      if(len(unique(X[,j])) <= 2){
        feature_type[j] = "binary"
      }
      else if(len(unique(X[,j])) <= 5){
        feature_type[j] = 'discrete'
      }
      else{
        feature_type[j] = 'continuous'
      }
    }
  }
  model_regis <- default_model_registry()
  model_mu0 <- model_regis[[model_mu]]
  model_mu1 <- model_regis[[model_mu]]
  ########
  model_mu0_fit <- model_mu0$fit(X[idx0, , drop = FALSE], Y[idx0])
  model_mu1_fit <- model_mu1$fit(X[idx1, , drop = FALSE], Y[idx1])
  # Propensity: regularized logistic (glmnet) to handle possible separation
  pi_fit <- glmnet::glmnet(X, T, family = "binomial", alpha = 0, lambda = 0.1)
  pi_pred <- function(newx, lambda = 0.1) {
    p <- predict(pi_fit, as.matrix(newx), s = lambda, type = "response")
    pmin(pmax(as.numeric(p), clip_e), 1 - clip_e)
  }
  pi_vals <- pi_pred(X, lambda = 0.1)
  mu0_est <- model_mu0$predict(fit = model_mu0_fit, X_new = X)
  mu1_est <- model_mu1$predict(fit = model_mu1_fit, X_new = X)
  pseudo <- (Y - ifelse(T == 1, mu1_est, mu0_est)) * (T - pi_vals) /
    (pi_vals * (1 - pi_vals)) + (mu1_est - mu0_est)
  tau_model <- model_regis[[model_tau]]
  tau_fit <- tau_model$fit(X, pseudo)
  tau_est <- tau_model$predict(fit = tau_fit, X_new = X)
  nu_estimators <- vector("list", d)
  #Specify the list of the feature types:
  feature_type_list <- infer_type(X)
  for (j in seq_len(d)) {
    cols <- setdiff(seq_len(d), j)
    X_minus_j <- X[, cols, drop = FALSE]
    y_j <- X[, j]
    if(feature_type_list[j] == 'continuous'){
      model_nu_estimator <- model_regis[[model_nu]]
    } 
    else if(feature_type_list[j] == 'binary'){
      model_nu_estimator <- model_regis[[model_nu_bina]]
    }
    else{
      model_nu_estimator <- model_regis[[model_nu_disc]]
    }
    nu_j <- model_nu_estimator$fit(X_minus_j, y_j, seed = 2000 + j)
    nu_estimators[[j]] <- list(cols = cols, nu = nu_j, model_nu = model_nu_estimator)
  }
  tau_estimator = list(tau_est = as.numeric(tau_est), tau_model = tau_model, tau_fit = tau_fit)
  pi_estimator = list(pi_fit = pi_fit, pi_pred = pi_pred, pi_est = as.numeric(pi_vals))
  list(
    tau_estimator = tau_estimator,
    pi_estimator = pi_estimator,
    mu0_est = mu0_est,
    mu1_est = mu1_est,
    nu_estimators = nu_estimators,
    d = d
  )
}

########################################################################################################
#' PermuCATE Variable Importance for the PO-risk, shuffle the residual for the covariate fitting model
#' Then return the variable importances and the combined permutation score:
#' @param pc:            List of the cross-fitted nuisance components.
#' @param X:             Covariate Space
#' @param Y:             Response
#' @param T:             Treatment Batch Label
#' @param clip_e:        The clipping value for the propensity score
#' @param n_permutation: Number of permutations
########################################################################################################
permucate_vimp <- function(
    pc,
    X, Y, T, clip_e,
    n_permutation, lambda = 0.5){
  X <- as.matrix(X)
  T <- as.numeric(T)
  Y <- as.numeric(Y)
  pi_vals <- pmin(1 - clip_e, pmax(clip_e, pc$pi_estimator$pi_est))
  mu0 <- pc$mu0_est
  mu1 <- pc$mu1_est
  tauX = pc$tau_estimator$tau_est
  Risk_baseline = po_risk(tauX, Y, T, mu0, mu1, pi_vals)
  d <- ncol(X)
  importance = numeric(d)
  all_scores = vector('list', d)
  for(j in 1:d){
    #Fitting X_{j}~X_{-j} on the predictive model:
    nu_info <- pc$nu_estimators[[j]]
    cols <- nu_info$cols
    X_minus_j <- X[, cols, drop = FALSE]
    nu_hat <- nu_info$model_nu$predict(fit = nu_info$nu, X_new = X_minus_j)
    r_j <- X[, j] - nu_hat
    psi_k <- numeric(n_permutation)
    for(k in 1:n_permutation){
      r_shuffle <- sample(r_j)#sampling the residuals!
      X_perm <- X
      X_perm[, j] <- nu_hat + r_shuffle
      tau_perm <- pc$tau_estimator$tau_model$predict(fit = pc$tau_estimator$tau_fit, X_new = X_perm)
      Risk_perm <- po_risk(tau_perm, Y, T, mu0, mu1, pi_vals)#permute the residuals and fit the PO-risk
      psi_k[k] <- (Risk_perm - Risk_baseline)
    }
    importance[j] <- mean(psi_k)
    all_scores[[j]] <- psi_k
  }
  score_per_perm <- do.call(rbind, all_scores)#return (d * n_permutation)
  return(list(importance = importance, score_per_perm = score_per_perm))
}

########################################################################################################
#' GRF variable importance with the Bootstrap CI(Impurity Gain)
#' Then return the variable importances and the combined permutation score:
#' @param X:             Covariate Space
#' @param Y:             Response
#' @param W:             Treatment Batch Label
#' @param n_bootstrap:   The number of bootstraps for estimate the confidence interval.
#' @param ci_level:      Confidence Interva Levels
#' @param n_estimators:  Number of Trees in Causal Forest
#' @param n_permutation: Number of permutations
########################################################################################################
grf_vimp_with_ci <- function(X, Y, W,
    n_bootstrap = 100, ci_level = 0.95, n_estimators = 100, 
    random_state = NULL){
  d <- ncol(X)
  X_tr <- X[W==1, , drop = FALSE]
  Y_tr <- Y[W==1]
  X_te <- X[W==0, , drop = FALSE]
  Y_te <- Y[W==0]
  W_tr <- W[W==1]
  W_te <- W[W==0]
  n_exist <- nrow(X_tr)
  n_new <- nrow(X_te)
  if(!is.null(random_state)) set.seed(random_state)
  alpha <- 1 - ci_level
  low_pct <- 100 * alpha / 2
  high_pct <- 100 * (1 - alpha / 2)
  if(n_bootstrap > 0){
    boot_imp <- matrix(0, n_bootstrap, d)
    for(b in seq_len(n_bootstrap)){
      idx_e <- sample(n_exist, n_exist, replace = TRUE)
      X_tr_b <- X_tr[idx_e, ]
      W_tr_b <- W_tr[idx_e]
      Y_tr_b <- Y_tr[idx_e]
      idx_n <- sample(n_new, n_new, replace = TRUE)
      X_te_b <- X_te[idx_n, ]
      W_te_b <- W_te[idx_n]
      Y_te_b <- Y_te[idx_n]
      cf <- grf::causal_forest(
        X = rbind(X_tr_b, X_te_b),
        Y = c(Y_tr_b, Y_te_b),
        W = c(W_tr_b, W_te_b),
        num.trees = n_estimators,
        sample.fraction = 0.5,
        honesty = TRUE,
        seed = random_state,
        mtry = round(ncol(X_tr)/2),
        min.node.size = round(nrow(X_tr)^(1/2)/2)
      )
      imp <- grf::variable_importance(cf)
      boot_imp[b, ] <- abs(imp)/sum(abs(imp))
    }
    importance <- colMeans(boot_imp)
    variance <- apply(boot_imp, 2, var)
    std_error <- sqrt(variance/ n_bootstrap)
    ci_low <- apply(boot_imp, 2, quantile, probs = alpha / 2)
    ci_high <- apply(boot_imp, 2, quantile, probs = 1 - alpha / 2)
    ci_width <- ci_high - ci_low
    return(
      list(
        importance = importance,
        std_error = std_error,
        variance = variance,
        ci_low = ci_low,
        ci_high = ci_high,
        ci_width = ci_width
      )
    )
  }
  else{
    cf <- grf::causal_forest(
        X = X, Y = Y, W = W,
        num.trees = n_estimators, 
        sample.fraction = 0.5,
        honesty = TRUE,
        seed = random_state,
        mtry = round(ncol(X_tr)/2),
        min.node.size = round(nrow(X_tr)^(1/2)/2)
    )
    imp <- grf::variable_importance(cf)
    return(imp)
  }
}
####################################################################################################
#' Estimate the uncertainty & confidence for the VIMP scores:
#'@param scores_per_perm:    Variable Importance for each of the variable and d times of permutations
#'@param lower_q:            Lower quantile for each of the feature
#'@param higher_q:           Higher quantile for each of the feature
#'@return: list of importance, std_error, variance, ci_low and ci_high
####################################################################################################
permutation_importance_uncertainty <- function(scores_per_perm,
  lower_q = 0.025, higher_q = 0.975){
  importances <- rowMeans(scores_per_perm)
  variance <- apply(scores_per_perm, 1, var)
  n_perm <- ncol(scores_per_perm)
  std_error <- sqrt(variance/n_perm)
  ci_low <- apply(scores_per_perm, 1, quantile, probs = lower_q, names = FALSE)
  ci_high <- apply(scores_per_perm, 1, quantile, probs = higher_q, names = FALSE)
  ci_width <- ci_high - ci_low
  return(list(
    importance = importances,
    std_error = std_error,
    variance = variance,
    ci_low = ci_low,
    ci_high = ci_high,
    ci_width = ci_width
  ))
}



































