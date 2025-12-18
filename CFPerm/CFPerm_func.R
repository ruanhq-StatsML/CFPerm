#The proposed CFPerm Procedure:
library(grf)
source("vimp_causal_forest_refit.R")
########
#Decision Rule Evaluation:
aggregate_power_ <- function(df, top_k = 1, q1 = 0.95, level = 0.05, normalization = TRUE){
  B <- ncol(df) - 1
  #getting the variable importance confidence intervals:
  vimp_ci <- t(apply(df[,1:B], 1,
               function(x){
                 c(min(x),
                   quantile(x, q1))
               }))
  vimp_ci <- as.data.frame(vimp_ci)
  if(normalization == TRUE){
    vimp_ci$original <- (df[, B+1])/sum(df[, B+1])
  }
  else{
    vimp_ci$original <- df[, B+1]
  }
  #Then we evaluate the result for the individual test:
  max_perm <- quantile(vimp_ci[,2], 1 - level)
  result1 <- apply(vimp_ci, 1, function(x){
    ifelse((x[3]>max_perm), 1, 0)
    })
  return(as.numeric(sum(result1) >= top_k))
}

#The conservative procedure with the most conservative decision rule -
#the decision threshold is the maximal VIMP across features for the maximal value in 
#the feature-specific interval.
#The method included include frequency based and the refitting based procedure:
cfperm <- function(df_train, df_test, n_perm = 200, q1 = 0.99, level = 0.05, num.trees = 150,
  vimp = "freq"){
  ncol_df <- ncol(df_train)
  n_train <- nrow(df_train)
  n_test <- nrow(df_test)
  n_feature <- ncol(df_train) - 1
  df_vimp_cover <- matrix(0, n_feature, n_perm + 1)
  df_train$trt <- 0
  df_test$trt <- 1
  combined_df <- rbind(df_train, df_test)
  ####
  n <- nrow(combined_df)
  mtry_now <- round(ncol_df/2)
  min_node_size_now <- round(n^(1/2)/2)
  cf_original <- causal_forest(X = combined_df[1:n, 1:n_feature], Y = combined_df$Y[1:n],
    W = combined_df$trt, num.trees = 150, sample.fraction = 0.5, honesty = TRUE,
    mtry = mtry_now, min.node.size = min_node_size_now
    )
  if(vimp == "freq"){
    df_vimp_cover[, n_perm + 1] <- variable_importance(cf_original)
  }
  else{
    df_vimp_cover[, n_perm + 1] <- vimp_causal_forests(cf_original)
  }
  for(b in 1:n_perm){
    combined_df1 <- combined_df
    combined_df1$trt <- sample(combined_df$trt, n, replace = FALSE)
    #rerun the two versions of variable importance for this causal forest:
    cf_new <- causal_forest(X = combined_df1[1:n, 1:n_feature], Y = combined_df1$Y[1:n],
      W = combined_df1$trt, num.trees = 150, sample.fraction = 0.5, honesty = TRUE,
      mtry = mtry_now, min.node.size = min_node_size_now)
    if(vimp == "freq"){
      df_vimp_cover[, b] <- variable_importance(cf_new)
    }
    else{
      df_vimp_cover[, b] <- vimp_causal_forests(cf_new)
    }
  }
  pvalue_list <- apply(df_vimp_cover, 1, function(x){1 - ecdf(x[1:n_perm])(x[n_perm+1])})
  result <- aggregate_power_(df_vimp_cover, q1 = q1, level = level)
  return(list(result = result, pvalue_list = pvalue_list))
}




