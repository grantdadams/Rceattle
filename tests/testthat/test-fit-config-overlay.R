# fit_mod(config=) overlay. A run configuration from load_config() / run_config()
# fills only the arguments the caller did NOT pass (an explicit argument always
# wins), and every fit records the run config it resolved so save_config(fit)
# round-trips. config = NULL (default) is a complete no-op.

.build3 <- function(d, ...) {
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, file = NULL, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0), ...)))
}

testthat::test_that("fit records its resolved run_config; config fills omitted args, explicit args win", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  d <- make_test_data()

  # Every fit records the run configuration it used.
  f <- .build3(d)
  testthat::expect_s3_class(f$run_config, "Rceattle_run_config")
  testthat::expect_equal(f$run_config$estimateMode, 3)
  testthat::expect_equal(f$run_config$model_config$initMode, "NonEquilibrium")

  # A config carrying non-default fields fills them when the caller omits them
  # (a model_config structure field and a top-level estimation control). Call
  # fit_mod() directly, omitting initMode and random_rec, so the config -- not a
  # helper default -- supplies them.
  cfg <- Rceattle:::.rce_run_config(
    mc = Rceattle::model_config(initMode = "OffsetEquilibrium"),
    estimateMode = 3, random_rec = TRUE)
  f_cfg <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, file = NULL, estimateMode = 3, msmMode = 0, config = cfg,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))
  testthat::expect_equal(f_cfg$run_config$model_config$initMode, "OffsetEquilibrium")
  testthat::expect_true(f_cfg$run_config$random_rec)

  # An explicitly-passed argument beats the config (explicit initMode + the
  # helper's explicit random_rec = FALSE both win over the config's values).
  f_exp <- .build3(d, config = cfg, initMode = "NonEquilibrium")
  testthat::expect_equal(f_exp$run_config$model_config$initMode, "NonEquilibrium")
  testthat::expect_false(f_exp$run_config$random_rec)

  # A non-run-config object is rejected.
  testthat::expect_error(.build3(d, config = list(estimateMode = 1)),
                         "Rceattle_run_config")
})

testthat::test_that("load_config(save_config(fit)) reproduces a non-default fit bit-identical", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not(exists("BS2017SS"))

  # The end-to-end proof: fit a NON-default configuration, round-trip its config
  # through YAML, and drive a fresh fit purely from the loaded config. The two
  # fits must be bit-identical (objective and SSB), the standard from the plan.
  fc <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  fit1 <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    BS2017SS, file = NULL, estimateMode = 1, initMode = "OffsetEquilibrium",
    random_rec = TRUE, fit_control = fc)))

  tf <- tempfile(fileext = ".yaml")
  suppressMessages(Rceattle::save_config(fit1, tf))
  cfg <- Rceattle::load_config(tf)

  fit2 <- suppressMessages(suppressWarnings(Rceattle::fit_mod(BS2017SS, file = NULL, config = cfg)))

  testthat::expect_identical(fit2$opt$objective, fit1$opt$objective)
  testthat::expect_identical(fit2$quantities$ssb, fit1$quantities$ssb)
})


# An absent DSEM must round-trip through YAML as ABSENT, not as a default one.
#
# Regression, and a silently-wrong-MODEL one. save_config() writes `dsem: []`
# for a fit with no DSEM, and load_config() materializes every builder field
# present in the YAML through its constructor. That is a no-op for the other
# seven -- build_srr(), build_M1() and the rest return the ordinary
# parameterization when called bare -- but build_DSEM()'s default makes the
# recruitment deviations the latent states of a GMRF. So the empty block
# reloaded as a REAL DSEM spec, fit_mod()'s `if (is.null(dsem))` guard never
# fired, and the refit declared dsem_x_tj random instead of rec_dev: the
# deviations stopped being integrated out. The objective moved only 5.4e-08 on
# this fixture, which is why it reads as rounding rather than as a different
# model.
#
# Assert the random block directly, not just the objective. A model that lands
# in the same place by luck would pass a tolerance check.
testthat::test_that("a fit with no DSEM does not reload as one", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  fc <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    Rceattle::BS2017SS, file = NULL, estimateMode = 1,
    initMode = "OffsetEquilibrium", random_rec = TRUE, fit_control = fc)))
  testthat::expect_null(fit[["dsem"]])
  testthat::expect_identical(fit[["random_vars"]], "rec_dev")

  tf <- tempfile(fileext = ".yaml")
  suppressMessages(Rceattle::save_config(fit, tf))
  cfg <- Rceattle::load_config(tf)

  # The config itself must not carry a materialized DSEM...
  testthat::expect_true(is.null(cfg$model_config$dsem) ||
                        length(cfg$model_config$dsem) == 0L)

  # ...and neither must the fit driven from it.
  again <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(Rceattle::BS2017SS, file = NULL, config = cfg)))
  testthat::expect_null(again[["dsem"]])
  testthat::expect_null(again[["data_list"]][["dsem_settings"]])
  testthat::expect_identical(again[["random_vars"]], "rec_dev")
})
