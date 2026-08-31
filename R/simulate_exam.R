#' Simulate a multiple-choice exam response matrix with distractors
#'
#' Simulates examinee ability, item parameters (2PL), and a
#' distractor-selection process in which lower-ability examinees are
#' more likely to be pulled toward a specific "common misconception"
#' distractor -- so that error-similarity between examinees is a
#' meaningful, ability-dependent signal even absent any collusion.
#'
#' @param n_examinees Integer, number of examinees.
#' @param n_items Integer, number of multiple-choice items.
#' @param n_options Integer, number of options per item (incl. the key).
#' @param theta_sd Standard deviation of the ability distribution,
#'   \eqn{\theta \sim N(0, \code{theta_sd})}.
#' @param seed Integer RNG seed for reproducibility.
#'
#' @return A list with components \code{theta} (true generating
#'   ability), \code{item_params} (data frame of discrimination \code{a}
#'   and difficulty \code{b}), \code{responses} (examinee x item matrix
#'   of chosen option labels; 1 = key), and \code{scored} (0/1 matrix).
#'
#' @examples
#' exam <- simulate_exam(n_examinees = 100, n_items = 20, seed = 1)
#' dim(exam$responses)
#'
#' @export
simulate_exam <- function(n_examinees = 500,
                           n_items = 40,
                           n_options = 4,
                           theta_sd = 1,
                           seed = 2026) {
  set.seed(seed)

  theta <- stats::rnorm(n_examinees, mean = 0, sd = theta_sd)
  a <- stats::runif(n_items, min = 0.7, max = 1.8)
  b <- stats::runif(n_items, min = -2.2, max = 2.2)

  distractor_bias <- matrix(stats::rnorm(n_items * (n_options - 1), sd = 0.6),
                             nrow = n_items)
  misconception_slope <- matrix(stats::runif(n_items * (n_options - 1), 0.1, 0.5),
                                 nrow = n_items)

  responses <- matrix(NA_integer_, nrow = n_examinees, ncol = n_items,
                       dimnames = list(paste0("s", seq_len(n_examinees)),
                                        paste0("item", seq_len(n_items))))
  P_correct <- matrix(NA_real_, nrow = n_examinees, ncol = n_items)

  for (j in seq_len(n_items)) {
    p_correct <- 1 / (1 + exp(-a[j] * (theta - b[j])))
    P_correct[, j] <- p_correct
    correct_draw <- stats::runif(n_examinees) < p_correct

    logits <- outer(-theta, misconception_slope[j, ]) +
      matrix(distractor_bias[j, ], nrow = n_examinees,
             ncol = n_options - 1, byrow = TRUE)
    expl <- exp(logits)
    probs <- expl / rowSums(expl)

    distractor_choice <- apply(probs, 1, function(p) {
      sample(2:n_options, size = 1, prob = p)
    })
    responses[, j] <- ifelse(correct_draw, 1L, distractor_choice)
  }

  scored <- (responses == 1L) * 1L

  list(
    theta = theta,
    item_params = data.frame(item = seq_len(n_items), a = a, b = b),
    distractor_bias = distractor_bias,
    misconception_slope = misconception_slope,
    responses = responses,
    scored = scored,
    P_correct = P_correct
  )
}

#' Inject an ability-concentrated answer-copying signal
#'
#' For each designated (source, copier) pair, on items the source got
#' WRONG, the copier's distractor choice is overwritten to match the
#' source's with elevated probability -- but only when the source's
#' ability falls inside \code{ability_band}. This creates a disparity
#' concentrated at one part of the ability continuum, the scenario
#' \code{copyfair}'s conditional test is designed to catch.
#'
#' @param exam A list produced by \code{\link{simulate_exam}}.
#' @param n_copy_pairs Number of colluding (source, copier) pairs to inject.
#' @param copy_rate Probability that a copier's wrong-answer match is
#'   overwritten to equal the source's, when inside the ability band.
#' @param ability_band Numeric length-2 vector, the theta range within
#'   which copying is active for the source examinee.
#' @param seed Integer RNG seed.
#'
#' @return \code{exam} with an added \code{copy_pairs} data frame
#'   (ground truth) and modified \code{responses}/\code{scored}.
#'
#' @examples
#' exam <- simulate_exam(n_examinees = 200, n_items = 30, seed = 1)
#' exam <- inject_copying(exam, n_copy_pairs = 5, ability_band = c(-0.5, 0.5))
#' nrow(exam$copy_pairs)
#'
#' @export
inject_copying <- function(exam,
                            n_copy_pairs = 15,
                            copy_rate = 0.85,
                            ability_band = c(-0.5, 0.5),
                            seed = 2027) {
  set.seed(seed)
  n <- nrow(exam$responses)

  eligible <- which(exam$theta >= ability_band[1] & exam$theta <= ability_band[2])
  if (length(eligible) < n_copy_pairs * 2) {
    stop("Not enough examinees in the ability band to draw non-overlapping pairs.")
  }

  chosen <- sample(eligible, size = n_copy_pairs * 2, replace = FALSE)
  sources <- chosen[seq_len(n_copy_pairs)]
  copiers <- chosen[(n_copy_pairs + 1):(2 * n_copy_pairs)]

  responses <- exam$responses
  for (k in seq_len(n_copy_pairs)) {
    s <- sources[k]; cpr <- copiers[k]
    wrong_items <- which(responses[s, ] != 1L)
    overwrite <- wrong_items[stats::runif(length(wrong_items)) < copy_rate]
    responses[cpr, overwrite] <- responses[s, overwrite]
  }

  exam$responses <- responses
  exam$scored <- (responses == 1L) * 1L
  ids <- rownames(exam$responses)
  exam$copy_pairs <- data.frame(
    source = ids[sources], copier = ids[copiers],
    source_theta = exam$theta[sources], copier_theta = exam$theta[copiers]
  )
  exam
}

#' Simulate a demonstration exam with an injected copying signal
#'
#' Convenience wrapper combining \code{\link{simulate_exam}} and
#' \code{\link{inject_copying}} in one call.
#'
#' @inheritParams simulate_exam
#' @inheritParams inject_copying
#'
#' @return A list as returned by \code{\link{inject_copying}}.
#'
#' @examples
#' exam <- build_copyfair_demo(n_examinees = 200, n_items = 30, n_copy_pairs = 5)
#' nrow(exam$copy_pairs)
#'
#' @export
build_copyfair_demo <- function(n_examinees = 500,
                                 n_items = 40,
                                 n_options = 4,
                                 n_copy_pairs = 15,
                                 copy_rate = 0.85,
                                 ability_band = c(-0.5, 0.5),
                                 seed = 2026) {
  exam <- simulate_exam(n_examinees, n_items, n_options, seed = seed)
  inject_copying(exam, n_copy_pairs, copy_rate, ability_band, seed = seed + 1)
}
