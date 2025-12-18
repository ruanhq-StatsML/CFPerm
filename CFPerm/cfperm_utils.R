#cfperm utils:
library(grf)

#This script summarizes the 
#aggregate the power, operation on the matrix with B(number of permutations) * p(number of features)
aggregate_power1 <- function(df, level = 0.5, top_k = 1, normalization = TRUE){
  B <- ncol(df) - 1
  vimp_ci <- t(apply(df[,1:B], 1,
    function(x){
      c(min(x), max(x))
      }))
  vimp_ci <- as.data.frame(vimp_ci)
  if(normalization == TRUE){
    vimp_ci$original <- (df[, B+1])/sum(df[,B+1])
  }
  else{
    vimp_ci$original <- df[,B+1]
  }
  max_perm = quantile(vimp_ci[,2], 1 - level)
  result1 <- apply(vimp_ci, 1, function(x){
    ifelse((x[3] > max_perm), 1, 0)
    })
  return(as.numeric(sum(result1) >= top_k))
}

#The generalized procedure where q_1 represent the feature-level quantile and q_2 represent the 
#across feature-level quantile for aggregation:
cfperm <- function()