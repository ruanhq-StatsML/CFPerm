#' CFPerm: causal-forest permutation test using variable importance
#'
#' Trains a causal forest on pooled train/test with a treatment indicator (trt),
#' computes the original variable importance, and compares it to permutation
#' distribution obtained by permuting the treatment indicator.
#'
#' @param df_train Data frame containing features and a response column.
#' @param df_test Data frame with the same columns as df_train.
#' @param n_perm Integer; number of permutations.
#' @param level_feature Numeric in (0,1); feature-wise tail level for aggregation.
#' @param level_across_feature Numeric in (0,1); across-feature tail level for aggregation.
#' @param num.trees Integer; number of trees for the causal forest.
#' @param y_col Character; name of the response column. Default "Y".
#' @param seed Optional integer seed for reproducibility.
#'
#' @return Integer 0 or 1 (reject indicator).
#' @importFrom grf causal_forest variable_importance
#' @export
cfperm <- function(df_train,
                   df_test,
                   n_perm = 200,
                   level_feature = 0.01,
                   level_across_feature = 0.05,
                   num.trees = 150,
                   y_col = "Y",
                   seed = NULL) {
  if (!is.data.frame(df_train) || !is.data.frame(df_test)) stop("`df_train` and `df_test` must be data.frames.")
  if (!identical(colnames(df_train), colnames(df_test))) stop("`df_train` and `df_test` must have identical columns (same order).")
  if (!y_col %in% colnames(df_train)) stop("Response column `y_col` not found in df_train/df_test.")
  if (!is.numeric(n_perm) || length(n_perm) != 1 || n_perm < 1) stop("`n_perm` must be a positive integer-like value.")
  n_perm <- as.integer(n_perm)
  if (!is.numeric(num.trees) || length(num.trees) != 1 || num.trees < 1) stop("`num.trees` must be a positive integer-like value.")
  num.trees <- as.integer(num.trees)

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1) stop("`seed` must be a scalar integer-like value.")
    set.seed(as.integer(seed))
  }

  n_train <- nrow(df_train)
  n_test <- nrow(df_test)
  if (n_train < 2 || n_test < 2) stop("Need at least 2 rows in both df_train and df_test.")

  df_train$trt <- 0L
  df_test$trt <- 1L
  combined_df <- rbind(df_train, df_test)

  feature_cols <- setdiff(colnames(df_train), c(y_col))
  X <- as.matrix(combined_df[, feature_cols, drop = FALSE])
  Y <- combined_df[[y_col]]
  W <- combined_df$trt
  n <- nrow(combined_df)
  p <- ncol(X)

  # heuristic settings (as in the original script, but guarded)
  mtry_now <- max(1L, round((p + 1) / 2))
  min_node_size_now <- max(1L, round(sqrt(n) / 2))

  df_vimp_cover <- matrix(0, nrow = p, ncol = n_perm + 1)

  cf_original <- grf::causal_forest(
    X = X, Y = Y, W = W,
    num.trees = num.trees,
    sample.fraction = 0.5,
    honesty = TRUE,
    mtry = mtry_now,
    min.node.size = min_node_size_now
  )

  df_vimp_cover[, n_perm + 1] <- grf::variable_importance(cf_original)

  for (b in seq_len(n_perm)) {
    W_perm <- sample(W, size = n, replace = FALSE)
    cf_new <- grf::causal_forest(
      X = X, Y = Y, W = W_perm,
      num.trees = num.trees,
      sample.fraction = 0.5,
      honesty = TRUE,
      mtry = mtry_now,
      min.node.size = min_node_size_now
    )
    df_vimp_cover[, b] <- grf::variable_importance(cf_new)
  }

  aggregate_power(
    df = df_vimp_cover,
    level_feature = level_feature,
    level_across_feature = level_across_feature
  )
}
