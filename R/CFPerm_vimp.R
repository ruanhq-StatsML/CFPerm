#Causal Inference Permutation Test here:

source('cfperm_config.R')

#Call the config here:
default_vimp_config <- function(vimp){
  if(vimp == 'permuCATE'){
    return(permucate_config())
  }
  if(vimp == 'LOCO'){
    return(loco_config())
  }
  if(vimp == 'GRF'){
    return(grf_config())
  }
}
########################################################################################
#' CFPerm main function to conduct the permutation testing for causal inference,
#' The propensity score model is the regularized logistic regression by default.
#' @param X                         :Covariates, (n_exist + n_new, p)
#' @param Y                         :Responses vector, (n_exist + n_new, )
#' @param T                         :Treatment(Batch) Assignment Labels, (n_exist + n_new, ) 
#' @param n_perm                    :Number of Permutations for the batch assignment in the test
#' @param n_folds                   :The number of folds for cross-fitting
#' @param model_mu                  :The types of the outcome model
#' @param model_nu                  :The types of the conditional fitting model, by default random forest.
#' @param model_tau                 :The pseudo-outcome model
#' @param vimp                      :The notion of the variable importance
#' @param level_feature             :The level for the featurewise confidence interval
#' @param level_across_feature      :The level across the quantile for each of the feature
#' @param top_k                     :top_k: the decision rule regarding how many features are significantly higher?
#' @param n_permutation             :Number of Permutations for the conditional permutation variable importance(PermuCATE)
########################################################################################
cfperm <- function(X, Y, T, 
    n_perm = 200,
    vimp = 'permuCATE',
    level_feature = 0.01,
    level_across_feature = 0.05,
    top_k = 1,
    clip_e = 0.01,
    clip_wtilde = 1e-3,
    n_permutation = 30,
    lambda = 0.25,
    n_folds = 5,
    model_mu = 'xgb_regression',
    model_nu = 'xgb_regression',
    model_tau = 'rf_regression',
    seed = 2026
    ){
    n <- nrow(X)
    p <- ncol(X)
    if(vimp == 'permuCATE'){
      pc <- permucate_fit(X, Y, T, lambda = 0.01,
        model_mu = model_mu, model_nu = model_nu,
        model_tau = model_tau, clip_e = clip_e)
      imp <- permucate_vimp(
        pc, X, Y, T, clip_e = clip_e,
        n_permutation = n_permutation,
        lambda = lambda
      )$importance
    }
    else if(vimp == 'loco'){
      imp <- importance_LOCO_r_risk(X, Y, T, seed = seed, n_folds = n_folds,
        binary_outcome = TRUE, lambda = 0.1, clip_e = 0.05, clip_wtilde = 1e-3)
    }
    else if(vimp == 'GRF'){
      imp <- grf_vimp_with_ci(X, Y, T, n_bootstrap = 2, ci_level = 0.95,
        n_estimators = 100, random_state = 2026)$importance
    }
    else{
      stop('Only support permuCATE, loco and GRF variable importance notions.')
    }
    imp_perm <- matrix(NA, nrow = p, ncol = n_perm)
    for(b in seq_len(n_perm)){
      T_resampled <- sample(T, size = n, replace = FALSE)
      if(vimp == 'permuCATE'){
        pc <- permucate_fit(X, Y, T_resampled, lambda = 0.01,
          model_mu = model_mu, model_nu = model_nu,
          model_tau = model_tau, clip_e = clip_e)
        imp_perm[,b] <- permucate_vimp(
          pc, X, Y, T_resampled, clip_e = clip_e,
          n_permutation = n_permutation,
          lambda = lambda)$importance
        }
      else if(vimp == 'loco'){
        imp_perm[,b] <- importance_LOCO_r_risk(X, Y, T_resampled, seed = seed, n_folds = n_folds,
        binary_outcome = TRUE, lambda = 0.1, clip_e = 0.05, clip_wtilde = 1e-3)
      }
      else if(vimp == 'GRF'){
        imp_perm[,b] <- grf_vimp_with_ci(X, Y, T_resampled, n_bootstrap = 2, ci_level = 0.95,
          n_estimators = 100, random_state = 2026)$importance
      }
      else{
        stop('Only support permuCATE, loco and GRF variable importance notions.')
      }
    }
    pvals <- rowMeans(imp_perm >= imp)#feature-specific p-value
    q1_upper <- apply(imp_perm, 1, quantile, probs = 1 - level_feature, names = FALSE)
    threshold <- as.numeric(quantile(q1_upper, probs = 1 - level_across_feature, names = FALSE))
    rejected <- (sum(imp > threshold) >= top_k)
    return(rejected)
}


