# The estimability table shown when sdreport fails is gated by `verbose`.
# TMBhelper print()s it to stdout, which suppressMessages() cannot catch.

# Same return shape as TMBhelper::check_estimability(), and it print()s.
fake_estimability <- function(obj) {
  bp <- data.frame(Param       = c("rec_pars", "beta_linkage_re_pen"),
                   MLE         = c(1.1, -2.2),
                   Param_check = c(1L, 2L))
  print(bp)
  list(BadParams = bp, WhichBad = 2L)
}

# iter.max = 1 stops far from the optimum, so sdreport fails and the branch runs.
run_and_capture <- function(verbose) {
  so <- textConnection("so_out", "w", local = TRUE)
  se <- textConnection("se_out", "w", local = TRUE)
  sink(so); sink(se, type = "message")
  on.exit({
    suppressWarnings(try(sink(type = "message"), silent = TRUE))
    suppressWarnings(try(sink(), silent = TRUE))
    close(so); close(se)
  }, add = TRUE)

  testthat::local_mocked_bindings(.check_estimability = fake_estimability,
                                  .package = "Rceattle")
  # The non-positive-definite Hessian is the condition being induced.
  fit <- suppressWarnings(try(Rceattle::fit_mod(
    Rceattle::BS2017SS, file = NULL, estimateMode = 1, msmMode = 0,
    random_rec = FALSE,
    fit_control = Rceattle::fit_control(
      getsd = TRUE, verbose = verbose, phase = FALSE, newtonsteps = 0,
      nlminb_control = list(eval.max = 2, iter.max = 1, trace = 0))),
    silent = TRUE))

  suppressWarnings(try(sink(type = "message"), silent = TRUE))
  suppressWarnings(try(sink(), silent = TRUE))
  list(fit = fit, stdout = so_out, stderr = se_out)
}

testthat::test_that("the estimability table is silent at verbose = 0", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  r <- run_and_capture(0)

  testthat::expect_true(is.null(r$fit$opt$SD))   # branch actually reached
  testthat::expect_false(any(grepl("Param_check", r$stdout)))
  testthat::expect_false(any(grepl("Param_check", r$stderr)))
  testthat::expect_false(any(grepl("did not converge", r$stderr)))
})

testthat::test_that("the estimability table is silent at verbose = 1", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # From 5.26.0 the table is one row per parameter, each named after its block,
  # and `fit$convergence` reports the flagged ones by fleet, bin and year
  # instead. The rest of the verbose = 1 reporting is unchanged.
  r <- run_and_capture(1)

  testthat::expect_true(is.null(r$fit$opt$SD))
  testthat::expect_false(any(grepl("Param_check", r$stderr)))
  testthat::expect_false(any(grepl("Param_check", r$stdout)))
  testthat::expect_true(any(grepl("did not converge", r$stderr)))
})

testthat::test_that("the estimability table is shown at verbose > 1", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  r <- run_and_capture(2)

  testthat::expect_true(is.null(r$fit$opt$SD))
  # Re-emitted as a message, so still suppressible and never on stdout.
  testthat::expect_true(any(grepl("Param_check", r$stderr)))
  testthat::expect_false(any(grepl("Param_check", r$stdout)))
})

testthat::test_that("the estimability verdict is kept whatever verbose says", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  r <- run_and_capture(0)
  testthat::expect_false(is.null(r$fit$identified))
  testthat::expect_equal(r$fit$identified$WhichBad, 2L)
  testthat::expect_true("beta_linkage_re_pen" %in%
                          as.character(r$fit$identified$BadParams$Param))
})
