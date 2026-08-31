## =====================================================================
## analysis/simulation_study.R
##
## Reproduces the simulation study reported in the copyfair manuscript:
## injects an ability-concentrated answer-copying signal into simulated
## exam data and confirms that copyfair's wild-bootstrap omnibus test
## detects it, while contrasting with what a global (unconditioned)
## index would show.
##
## Self-contained -- generates its own data, no external files needed.
## Usage: Rscript analysis/simulation_study.R
## =====================================================================

library(copyfair)

cat("=== copyfair simulation study ===\n\n")

exam <- build_copyfair_demo(
  n_examinees  = 500,
  n_items      = 40,
  n_options    = 4,
  n_copy_pairs = 15,
  copy_rate    = 0.85,
  ability_band = c(-0.5, 0.5),
  seed         = 2026
)

cat("Ground truth:", nrow(exam$copy_pairs),
    "colluding pairs injected, source ability in [-0.5, 0.5]\n\n")

key <- setNames(rep(1, ncol(exam$responses)), colnames(exam$responses))
result <- copyfair_analyze(exam$responses, key, min_common = 5, B = 2000)

cat("Wild-bootstrap omnibus test:\n")
cat("  observed statistic:", round(result$omnibus_test$observed_stat, 4), "\n")
cat("  p-value            :", round(result$omnibus_test$p_value, 4), "\n\n")

truth_keys <- paste(exam$copy_pairs$source, exam$copy_pairs$copier)
result$pairs$is_truth <- paste(result$pairs$source, result$pairs$copier) %in% truth_keys

cat("Mean standardized residual (z), true copy pairs   :",
    round(mean(result$pairs$z[result$pairs$is_truth], na.rm = TRUE), 3), "\n")
cat("Mean standardized residual (z), innocent pairs     :",
    round(mean(result$pairs$z[!result$pairs$is_truth], na.rm = TRUE), 3), "\n\n")

cat("Pairs flagged (q < .05):", sum(result$pairs$flagged), "\n")
cat("  of which true copy pairs:",
    sum(result$pairs$flagged & result$pairs$is_truth), "/", nrow(exam$copy_pairs), "\n")
