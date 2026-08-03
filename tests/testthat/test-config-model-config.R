# model_config() slot and fit_mod()'s missing()-based override. Constructor
# tests are fast/unguarded; anything invoking fit_mod() (even estimateMode = 3,
# which builds MakeADFun) is skipped on CRAN.

testthat::test_that("model_config() returns a classed object with fit_mod defaults", {
  cfg <- model_config()
  testthat::expect_s3_class(cfg, "Rceattle_model_config")
  testthat::expect_identical(cfg$msmMode, 0)
  testthat::expect_identical(cfg$initMode, "NonEquilibrium")
  testthat::expect_identical(cfg$niter, 3)
  # Carries the build_*() process objects.
  for (f in c("HCR", "recFun", "M1Fun", "growthFun", "qFun", "selFun", "compFun"))
    testthat::expect_false(is.null(cfg[[f]]), info = f)
})

testthat::test_that("model_config() validates its scalar switches", {
  testthat::expect_error(model_config(msmMode = c(0, 1)), "single value")
  testthat::expect_error(model_config(niter = 0), "positive integer")
})

testthat::test_that("print.Rceattle_model_config is a compact summary", {
  out <- utils::capture.output(print(model_config(msmMode = 1)))
  testthat::expect_true(any(grepl("Rceattle model_config", out)))
  testthat::expect_true(any(grepl("msmMode", out)))
})

testthat::test_that("build_data() accepts a model_config slot", {
  testthat::skip_if_not_installed("Rceattle")
  d <- suppressWarnings(suppressMessages(
    build_data(base = BS2017SS, model_config = model_config(msmMode = 0),
               .check = FALSE)))
  testthat::expect_s3_class(d$model_config, "Rceattle_model_config")
})

testthat::test_that("fit_mod honours the slot only for omitted arguments", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("Rceattle")

  # Config value is used when the caller omits the argument...
  dm <- BS2017MS
  dm$model_config <- model_config(msmMode = 1)
  b_omit <- suppressWarnings(suppressMessages(
    Rceattle::fit_mod(dm, file = NULL, estimateMode = 3)))
  testthat::expect_equal(b_omit$data_list$msmMode, 1)

  # ...and an explicitly-passed argument (even a default value) overrides it.
  b_arg <- suppressWarnings(suppressMessages(
    Rceattle::fit_mod(dm, file = NULL, estimateMode = 3, msmMode = 0)))
  testthat::expect_equal(b_arg$data_list$msmMode, 0)

  # initMode likewise flows from the slot.
  ds <- BS2017SS
  ds$model_config <- model_config(msmMode = 0, initMode = "Equilibrium")
  b_init <- suppressWarnings(suppressMessages(
    Rceattle::fit_mod(ds, file = NULL, estimateMode = 3)))
  testthat::expect_equal(b_init$data_list$initMode, "Equilibrium")
})

testthat::test_that("write_data warns that model_config will not persist", {
  testthat::skip_if_not_installed("Rceattle")
  d <- BS2017SS
  d$model_config <- model_config(msmMode = 0)
  tmp <- tempfile(fileext = ".xlsx")
  testthat::expect_warning(
    suppressMessages(Rceattle::write_data(d, tmp)),
    "not written to the xlsx")
})

testthat::test_that("model_config is scrubbed from the TMB data but kept on the fit", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("Rceattle")
  d <- BS2017SS
  d$model_config <- model_config(msmMode = 0)
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    d, file = NULL, estimateMode = 3)))
  # Not in the taped TMB data (a nested spec list must not ride into MakeADFun)...
  testthat::expect_false("model_config" %in% names(fit$obj$env$data))
  # ...but retained on the returned object so print()/summary() can show it.
  testthat::expect_s3_class(fit$data_list$model_config, "Rceattle_model_config")
})

testthat::test_that("a model_config slot resolving to defaults does not move the fit", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("Rceattle")
  fc <- function() fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  ss0 <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    BS2017SS, file = NULL, inits = NULL, estimateMode = 0, msmMode = 0,
    fit_control = fc())))
  d <- BS2017SS
  d$model_config <- model_config(msmMode = 0)
  ssc <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    d, file = NULL, inits = NULL, estimateMode = 0, fit_control = fc())))
  # The slot must not leak into the TMB data: bit-identical objective and pars.
  testthat::expect_equal(ssc$opt$objective, ss0$opt$objective, tolerance = 1e-10)
  testthat::expect_equal(ssc$opt$par, ss0$opt$par, tolerance = 1e-10)
})
