# =============================================================================
# Selectivity_block / Q_block time-block columns in index_data / catch_data.
#
# - Selectivity_block is read ONLY for Time_varying_{sel,q} == "Block"; every
#   other configuration ignores it. It is now OPTIONAL: clean_data() default-fills
#   a missing column with 1 (a single block = no time-blocking), so a user need
#   only supply it for Block-mode fleets.
# - Q_block is VESTIGIAL: never read (q time-blocking reuses Selectivity_block).
#   data_check() warns that a supplied Q_block is deprecated/ignored.
# =============================================================================

testthat::test_that("clean_data() default-fills a missing Selectivity_block with 1", {
  testthat::skip_if_not_installed("Rceattle")
  d <- make_test_data()
  d$index_data$Selectivity_block <- NULL
  d$catch_data$Selectivity_block <- NULL
  out <- Rceattle:::clean_data(d)
  testthat::expect_true(all(out$index_data$Selectivity_block == 1L))
  testthat::expect_true(all(out$catch_data$Selectivity_block == 1L))
})

testthat::test_that("clean_data() leaves an existing Selectivity_block untouched", {
  testthat::skip_if_not_installed("Rceattle")
  d <- make_test_data()
  d$index_data$Selectivity_block <- 2L   # a non-default value must survive
  out <- Rceattle:::clean_data(d)
  testthat::expect_true(all(out$index_data$Selectivity_block == 2L))
})

testthat::test_that("a model with no Selectivity_block column fits", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")
  d <- make_test_data()
  d$index_data$Selectivity_block <- NULL
  d$catch_data$Selectivity_block <- NULL
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0))))
  testthat::expect_true(is.finite(fit$obj$fn()))
})

testthat::test_that("a supplied Q_block warns (deprecated / ignored)", {
  testthat::skip_if_not_installed("Rceattle")
  d <- make_test_data()
  d$index_data$Q_block <- 1L
  dc <- get("data_check", asNamespace("Rceattle"))
  w <- character(0)
  withCallingHandlers(
    tryCatch(suppressMessages(dc(suppressMessages(Rceattle::switch_check(d)))),
             error = function(e) NULL),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  testthat::expect_true(any(grepl("Q_block", w)))
})
