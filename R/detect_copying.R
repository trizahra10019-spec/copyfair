#' Compute item x ability-bin option-choice probability table
#'
#' Internal helper: for each item and each ability bin, estimates the
#' population probability of choosing each answer option, using
#' Laplace (add-one) smoothing so probabilities are never exactly 0 or
#' 1 -- avoiding degenerate variances downstream when an ability bin
#' has few examinees for a given item.
#'
#' @param R Examinee x item response matrix (NA = not administered).
#' @param ability_bin Integer vector of ability-bin membership, one
#'   per row of \code{R}.
#' @param answer_levels Numeric vector of possible option codes.
#' @param min_bin_n Minimum examinees required in a bin before its
#'   empirical probabilities are trusted (otherwise the uniform prior
#'   \eqn{1/K} is kept).
#'
#' @return A 3-D array [item, option, bin] of probabilities.
#' @keywords internal
#' @noRd
build_option_prob_table <- function(R, ability_bin, answer_levels = 1:4,
                                      min_bin_n = 5) {
  n_items <- ncol(R)
  n_bins <- max(ability_bin, na.rm = TRUE)
  n_opt <- length(answer_levels)

  prob_tab <- array(1 / n_opt, dim = c(n_items, n_opt, n_bins),
                     dimnames = list(colnames(R), as.character(answer_levels), NULL))

  for (b in seq_len(n_bins)) {
    rows_b <- which(ability_bin == b)
    if (length(rows_b) < min_bin_n) next
    Rb <- R[rows_b, , drop = FALSE]
    for (k in seq_len(n_items)) {
      col <- Rb[, k]
      valid_n <- sum(!is.na(col))
      if (valid_n == 0) next
      for (vi in seq_along(answer_levels)) {
        v <- answer_levels[vi]
        prob_tab[k, vi, b] <- (sum(col == v, na.rm = TRUE) + 1) / (valid_n + n_opt)
      }
    }
  }
  prob_tab
}

#' Detect ability- and item-conditioned answer-copying signal
#'
#' The core \code{copyfair} routine. For every pair of examinees who
#' share at least \code{min_common} administered items, computes the
#' observed number of matching WRONG answers, and compares it against
#' an EXPECTED count derived from the population probability of
#' choosing that specific option, conditional on both the item and the
#' copier's ability bin. This jointly accounts for (a) ability-driven
#' baseline similarity and (b) item-specific "attractive distractor"
#' effects that a single ability-only or fully-global baseline would
#' conflate with genuine collusion.
#'
#' @param responses Examinee x item response matrix of chosen option
#'   codes (e.g. 1:4 for A-D). Use \code{NA} for items not administered
#'   to a given examinee (e.g. under a rotated-booklet design).
#' @param correct_key Named numeric vector or single-row-matching
#'   vector of correct-option codes, one per item, in the same column
#'   order as \code{responses}.
#' @param min_common Minimum number of jointly-administered items
#'   required to test a pair.
#' @param n_ability_bins Number of ability bins used to estimate
#'   item-and-ability-conditioned option probabilities.
#' @param answer_levels Numeric vector of possible option codes.
#'
#' @return A data frame with one row per directed pair (source =
#'   examinee whose wrong answers are being matched, copier = examinee
#'   being compared against), and columns \code{common_n}, \code{denom}
#'   (source's wrong-and-jointly-administered item count), \code{matches}
#'   (observed matching-wrong count), \code{expected}, \code{variance},
#'   \code{z} (standardized residual), and \code{mean_ability}.
#'
#' @examples
#' exam <- build_copyfair_demo(n_examinees = 150, n_items = 25,
#'                              n_copy_pairs = 5, seed = 1)
#' key <- setNames(rep(1, ncol(exam$responses)), colnames(exam$responses))
#' pairs <- copyfair_pairs(exam$responses, key, min_common = 5)
#' head(pairs[order(-pairs$z), ])
#'
#' @export
copyfair_pairs <- function(responses, correct_key, min_common = 5,
                             n_ability_bins = 5, answer_levels = 1:4) {
  R <- as.matrix(responses)
  if (is.null(colnames(R))) colnames(R) <- paste0("item", seq_len(ncol(R)))
  if (is.null(rownames(R))) rownames(R) <- paste0("s", seq_len(nrow(R)))
  users <- rownames(R); n_u <- nrow(R); n_i <- ncol(R)

  key <- correct_key[colnames(R)]
  if (any(is.na(key))) {
    stop("`correct_key` must have names matching every column of `responses` ",
         "(after auto-naming, if `responses` had no colnames). Missing key for: ",
         paste(colnames(R)[is.na(key)], collapse = ", "))
  }

  valid <- !is.na(R)
  wrong <- valid & sweep(R, 2, key, `!=`); wrong[is.na(wrong)] <- FALSE

  ability <- rowSums(valid & sweep(R, 2, key, `==`), na.rm = TRUE) / pmax(rowSums(valid), 1)
  ability_bin <- as.integer(cut(ability, breaks = n_ability_bins,
                                  include.lowest = TRUE, labels = FALSE))
  ability_bin[is.na(ability_bin)] <- 1L

  prob_tab <- build_option_prob_table(R, ability_bin, answer_levels)

  valid_num <- valid * 1
  common_n <- valid_num %*% t(valid_num)
  wrong_num <- wrong * 1
  denom_mat <- wrong_num %*% t(valid_num)

  matches_mat <- matrix(0, n_u, n_u)
  expected_mat <- matrix(0, n_u, n_u)
  var_mat <- matrix(0, n_u, n_u)

  for (v_i in seq_along(answer_levels)) {
    v <- answer_levels[v_i]
    Xv <- (wrong & (R == v)) * 1; Xv[is.na(Xv)] <- 0
    Yv <- (valid & (R == v)) * 1; Yv[is.na(Yv)] <- 0
    matches_mat <- matches_mat + Xv %*% t(Yv)

    Pv_by_bin <- t(prob_tab[, v_i, ])
    Pv_j <- Pv_by_bin[ability_bin, , drop = FALSE]
    Pv_j[is.na(Pv_j)] <- 0
    Pv_j_valid <- Pv_j * valid_num

    expected_mat <- expected_mat + Xv %*% t(Pv_j_valid)
    var_mat <- var_mat + Xv %*% t(Pv_j_valid * (1 - Pv_j_valid))
  }

  diag(common_n) <- 0; diag(denom_mat) <- 0
  diag(matches_mat) <- 0; diag(expected_mat) <- 0; diag(var_mat) <- 0

  eligible <- which(common_n >= min_common & denom_mat > 0, arr.ind = TRUE)
  eligible <- eligible[eligible[, 1] != eligible[, 2], , drop = FALSE]
  if (nrow(eligible) == 0) {
    return(data.frame(source = character(0), copier = character(0),
                       common_n = numeric(0), denom = numeric(0),
                       matches = numeric(0), expected = numeric(0),
                       variance = numeric(0), z = numeric(0),
                       mean_ability = numeric(0)))
  }

  z_raw <- (matches_mat[eligible] - expected_mat[eligible]) /
    sqrt(pmax(var_mat[eligible], 1e-6))

  data.frame(
    source = users[eligible[, 1]], copier = users[eligible[, 2]],
    common_n = common_n[eligible], denom = denom_mat[eligible],
    matches = matches_mat[eligible], expected = expected_mat[eligible],
    variance = var_mat[eligible], z = z_raw,
    mean_ability = (ability[eligible[, 1]] + ability[eligible[, 2]]) / 2
  )
}

#' Flag statistically significant answer-copying pairs
#'
#' Converts the per-pair standardized residual \code{z} from
#' \code{\link{copyfair_pairs}} into a one-sided p-value (excess
#' similarity only), applies Benjamini-Hochberg FDR correction across
#' all tested pairs, and flags pairs below \code{q_threshold}.
#'
#' @param pair_df Output of \code{\link{copyfair_pairs}}.
#' @param q_threshold FDR-adjusted significance threshold.
#' @param min_denom Minimum source-wrong-and-jointly-administered item
#'   count required to test a pair (small denom -> unstable z; such
#'   pairs are dropped before FDR correction so they do not dilute the
#'   correction for well-supported pairs).
#'
#' @return \code{pair_df} restricted to tested pairs, with added
#'   columns \code{p_value}, \code{q_value}, and \code{flagged}.
#'
#' @examples
#' exam <- build_copyfair_demo(n_examinees = 150, n_items = 25,
#'                              n_copy_pairs = 5, seed = 1)
#' key <- setNames(rep(1, ncol(exam$responses)), colnames(exam$responses))
#' pairs <- copyfair_pairs(exam$responses, key, min_common = 5)
#' flagged <- copyfair_flag(pairs)
#' sum(flagged$flagged)
#'
#' @export
copyfair_flag <- function(pair_df, q_threshold = 0.05, min_denom = 5) {
  pair_df <- pair_df[pair_df$denom >= min_denom, ]
  pair_df$p_value <- 1 - stats::pnorm(pair_df$z)
  pair_df$q_value <- stats::p.adjust(pair_df$p_value, method = "BH")
  pair_df$flagged <- pair_df$q_value < q_threshold
  pair_df
}

#' Wild-bootstrap omnibus test for localized copying disparity
#'
#' Tests whether the standardized residuals \code{z} show a
#' significant excess concentrated in any one ability bin (H1),
#' against the null of no localized disparity (H0), using a
#' Rademacher wild-bootstrap on the pair-level \code{z} values,
#' following the resampling logic of Wu (1986).
#'
#' @param pair_df Output of \code{\link{copyfair_pairs}}.
#' @param n_bins Number of ability bins for the omnibus statistic.
#' @param B Number of bootstrap replications.
#' @param seed RNG seed.
#'
#' @return A list with \code{observed_stat}, \code{p_value}, and
#'   \code{boot_dist}.
#'
#' @examples
#' exam <- build_copyfair_demo(n_examinees = 300, n_items = 30,
#'                              n_copy_pairs = 10, seed = 1)
#' key <- setNames(rep(1, ncol(exam$responses)), colnames(exam$responses))
#' pairs <- copyfair_pairs(exam$responses, key, min_common = 5)
#' test <- copyfair_omnibus_test(pairs, B = 200)
#' test$p_value
#'
#' @export
copyfair_omnibus_test <- function(pair_df, n_bins = 5, B = 1000, seed = 1) {
  if (nrow(pair_df) == 0) {
    stop("`pair_df` has no rows to test -- check that `min_common`/`min_denom` ",
         "are not too strict for your data.")
  }
  set.seed(seed)
  z <- pair_df$z
  z[is.na(z)] <- 0
  bins <- cut(pair_df$mean_ability, breaks = n_bins)
  bin_int <- as.integer(bins)
  n <- length(z)
  n_bin_levels <- max(bin_int, na.rm = TRUE)
  bin_count <- as.numeric(table(factor(bin_int, levels = seq_len(n_bin_levels))))

  stat_from_binsum <- function(binsum_vec) {
    m <- binsum_vec / bin_count
    max(abs(m), na.rm = TRUE)
  }

  obs_binsum <- rowsum(z, group = bin_int, reorder = TRUE)[, 1]
  T_obs <- stat_from_binsum(obs_binsum)

  boot_stats <- numeric(B)
  for (b in seq_len(B)) {
    w <- sample(c(-1, 1), size = n, replace = TRUE)
    binsum_b <- rowsum(z * w, group = bin_int, reorder = TRUE)[, 1]
    boot_stats[b] <- stat_from_binsum(binsum_b)
  }

  p_value <- mean(boot_stats >= T_obs)
  list(observed_stat = T_obs, p_value = p_value, boot_dist = boot_stats)
}

#' Run the full copyfair pipeline
#'
#' Convenience wrapper chaining \code{\link{copyfair_pairs}},
#' \code{\link{copyfair_flag}}, and \code{\link{copyfair_omnibus_test}}.
#'
#' @inheritParams copyfair_pairs
#' @param q_threshold FDR threshold passed to \code{\link{copyfair_flag}}.
#' @param min_denom Passed to \code{\link{copyfair_flag}}.
#' @param B Bootstrap replications passed to \code{\link{copyfair_omnibus_test}}.
#'
#' @return A list with \code{pairs} (flagged pair-level results) and
#'   \code{omnibus_test}.
#'
#' @examples
#' exam <- build_copyfair_demo(n_examinees = 200, n_items = 25,
#'                              n_copy_pairs = 8, seed = 1)
#' key <- setNames(rep(1, ncol(exam$responses)), colnames(exam$responses))
#' result <- copyfair_analyze(exam$responses, key)
#' result$omnibus_test$p_value
#'
#' @export
copyfair_analyze <- function(responses, correct_key, min_common = 5,
                               n_ability_bins = 5, answer_levels = 1:4,
                               q_threshold = 0.05, min_denom = 5, B = 1000) {
  pairs <- copyfair_pairs(responses, correct_key, min_common,
                            n_ability_bins, answer_levels)
  flagged <- copyfair_flag(pairs, q_threshold, min_denom)
  test <- copyfair_omnibus_test(pairs, n_bins = n_ability_bins, B = B)
  list(pairs = flagged, omnibus_test = test)
}
