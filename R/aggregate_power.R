#' Aggregate variable-importance evidence across features
#'
#' Given a matrix/data.frame of permutation-based variable importance (columns 1..B)
#' plus the original (non-permuted) importance in the last column, this function
#' performs a two-level aggregation:
#' (i) feature-wise upper quantile thresholding; (ii) across-feature thresholding.
#'
#' @param df Numeric matrix/data.frame of dimension (n_features) x (B + 1).
#' @param top_k Integer; return 1 if at least top_k features exceed the global threshold.
#' @param level_feature Numeric in (0,1); feature-wise tail level.
#' @param level_across_feature Numeric in (0,1); across-feature tail level.
#' @param normalization Logical; if TRUE, normalize original importances to sum to 1.
#' @param typeIerror Unused placeholder for future extensions (kept for compatibility).
#'
#' @return Integer 0 or 1.
#' @export
aggregate_power <- function(df,
                            top_k = 1,
                            level_feature = 0.01,
                            level_across_feature = 0.05,
                            normalization = TRUE,
                            typeIerror = TRUE) {
  if (is.data.frame(df)) df <- as.matrix(df)
  if (!is.matrix(df) || !is.numeric(df)) stop("`df` must be a numeric matrix or data.frame.")
  if (ncol(df) < 2) stop("`df` must have at least 2 columns (permutations + original).")

  n_features <- nrow(df)
  B <- ncol(df) - 1L

  if (!is.numeric(top_k) || length(top_k) != 1 || top_k < 1) stop("`top_k` must be a positive integer-like value.")
  top_k <- as.integer(top_k)
  if (top_k > n_features) stop("`top_k` cannot exceed the number of features (nrow(df)).")

  if (!is.numeric(level_feature) || length(level_feature) != 1 || level_feature <= 0 || level_feature >= 1)
    stop("`level_feature` must be in (0,1).")
  if (!is.numeric(level_across_feature) || length(level_across_feature) != 1 || level_across_feature <= 0 || level_across_feature >= 1)
    stop("`level_across_feature` must be in (0,1).")

  perm_mat <- df[, 1:B, drop = FALSE]
  original <- df[, B + 1]

  vimp_min <- apply(perm_mat, 1, min)
  vimp_q <- apply(perm_mat, 1, stats::quantile, probs = 1 - level_feature, names = FALSE)

  if (normalization) {
    s <- sum(original)
    if (s == 0) stop("Cannot normalize because sum of original importances is 0.")
    original <- original / s
  }

  max_perm <- stats::quantile(vimp_q, probs = 1 - level_across_feature, names = FALSE)

  hits <- original > max_perm
  as.integer(sum(hits) >= top_k)
}
