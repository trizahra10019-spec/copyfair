test_that("simulate_exam produces correctly-shaped output", {
  exam <- simulate_exam(n_examinees = 50, n_items = 10, seed = 1)
  expect_equal(dim(exam$responses), c(50, 10))
  expect_equal(dim(exam$scored), c(50, 10))
  expect_equal(length(exam$theta), 50)
})

test_that("inject_copying adds detectable ground-truth pairs", {
  exam <- simulate_exam(n_examinees = 100, n_items = 20, seed = 1)
  exam <- inject_copying(exam, n_copy_pairs = 5, ability_band = c(-1, 1), seed = 2)
  expect_equal(nrow(exam$copy_pairs), 5)
  expect_true(all(exam$copy_pairs$source_theta >= -1 & exam$copy_pairs$source_theta <= 1))
})

test_that("copyfair_pairs recovers an injected ability-concentrated signal", {
  exam <- build_copyfair_demo(n_examinees = 300, n_items = 30,
                               n_copy_pairs = 10, copy_rate = 0.9,
                               ability_band = c(-0.5, 0.5), seed = 1)
  key <- setNames(rep(1, ncol(exam$responses)), colnames(exam$responses))
  pairs <- copyfair_pairs(exam$responses, key, min_common = 5)

  expect_true(nrow(pairs) > 0)
  expect_true(all(c("source", "copier", "z", "expected", "variance") %in% names(pairs)))

  truth_keys <- paste(exam$copy_pairs$source, exam$copy_pairs$copier)
  pairs$is_truth <- paste(pairs$source, pairs$copier) %in% truth_keys

  # pasangan menyontek sungguhan harus punya z rata-rata jauh lebih
  # tinggi daripada pasangan tak bersalah
  expect_true(mean(pairs$z[pairs$is_truth], na.rm = TRUE) >
                mean(pairs$z[!pairs$is_truth], na.rm = TRUE))
})

test_that("copyfair_flag applies FDR correction and respects min_denom", {
  exam <- build_copyfair_demo(n_examinees = 150, n_items = 25,
                               n_copy_pairs = 5, seed = 1)
  key <- setNames(rep(1, ncol(exam$responses)), colnames(exam$responses))
  pairs <- copyfair_pairs(exam$responses, key, min_common = 5)
  flagged <- copyfair_flag(pairs, q_threshold = 0.05, min_denom = 5)

  expect_true(all(flagged$denom >= 5))
  expect_true(all(c("p_value", "q_value", "flagged") %in% names(flagged)))
  expect_true(all(flagged$q_value[flagged$flagged] < 0.05))
})

test_that("copyfair_omnibus_test returns a valid p-value", {
  exam <- build_copyfair_demo(n_examinees = 200, n_items = 25,
                               n_copy_pairs = 8, seed = 1)
  key <- setNames(rep(1, ncol(exam$responses)), colnames(exam$responses))
  pairs <- copyfair_pairs(exam$responses, key, min_common = 5)
  test <- copyfair_omnibus_test(pairs, B = 200)

  expect_true(test$p_value >= 0 && test$p_value <= 1)
  expect_true(test$observed_stat >= 0)
})

test_that("copyfair_analyze wrapper returns pairs and omnibus_test", {
  exam <- build_copyfair_demo(n_examinees = 150, n_items = 20,
                               n_copy_pairs = 5, seed = 1)
  key <- setNames(rep(1, ncol(exam$responses)), colnames(exam$responses))
  result <- copyfair_analyze(exam$responses, key, B = 200)

  expect_true(is.data.frame(result$pairs))
  expect_true(is.list(result$omnibus_test))
  expect_true("p_value" %in% names(result$omnibus_test))
})

test_that("estimate_ability correlates with true theta", {
  exam <- simulate_exam(n_examinees = 200, n_items = 30, seed = 1)
  theta_hat <- estimate_ability(exam$scored, exam$item_params)
  expect_true(cor(theta_hat, exam$theta) > 0.7)
})
