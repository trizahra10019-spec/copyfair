#' Estimate examinee ability via joint MLE under a known 2PL model
#'
#' A lightweight Newton-Raphson maximum-likelihood ability estimator,
#' used when item parameters (\code{a}, \code{b}) are already known or
#' calibrated externally (e.g. via \pkg{mirt} or \pkg{TAM}). Included
#' so \code{copyfair} has no hard IRT-estimation dependency for its
#' simulation examples and vignette.
#'
#' @param scored An n_examinees x n_items 0/1 matrix.
#' @param item_params A data frame with columns \code{a} (discrimination)
#'   and \code{b} (difficulty), one row per item, in the same column
#'   order as \code{scored}.
#' @param max_iter Maximum Newton-Raphson iterations per examinee.
#' @param tol Convergence tolerance on the ability update step.
#'
#' @return A numeric vector of ability estimates, length
#'   \code{nrow(scored)}.
#'
#' @examples
#' exam <- simulate_exam(n_examinees = 100, n_items = 20, seed = 1)
#' theta_hat <- estimate_ability(exam$scored, exam$item_params)
#' cor(theta_hat, exam$theta)
#'
#' @export
estimate_ability <- function(scored, item_params, max_iter = 50, tol = 1e-5) {
  a <- item_params$a; b <- item_params$b
  n <- nrow(scored)
  theta <- rep(0, n)

  for (i in seq_len(n)) {
    x <- scored[i, ]
    th <- 0
    for (it in seq_len(max_iter)) {
      p <- 1 / (1 + exp(-a * (th - b)))
      score <- sum(a * (x - p))
      info <- sum(a^2 * p * (1 - p))
      if (info < 1e-8) break
      step <- score / info
      th <- th + step
      th <- max(min(th, 4), -4)
      if (abs(step) < tol) break
    }
    theta[i] <- th
  }
  theta
}
