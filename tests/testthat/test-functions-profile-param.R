# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

testthat::test_that("profile: 1-D sigmaR profile fixes R_log_sd at the grid value", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  library(Rceattle)

  data(BS2017SS)

  ss_run <- fit_mod(
    data_list   = BS2017SS,
    inits       = NULL,
    file        = NULL,
    estimateMode = 1,        # Hindcast only -- keeps the smoke test fast
    random_rec  = FALSE,
    msmMode     = 0,
    fit_control = fit_control(phase = TRUE, getsd = FALSE, verbose = 0)
  )

  grid_vals <- c(0.2, 0.6, 1.2)
  prof <- profile(
    fitted   = ss_run,
    param    = "R_log_sd",
    slots    = list(1),
    values   = list(grid_vals),
    cores    = 1
  )

  # Structural expectations
  testthat::expect_named(prof, c("Rceattle_list", "grid", "nll", "param", "slots"))
  testthat::expect_equal(prof$param, "R_log_sd")
  testthat::expect_equal(nrow(prof$grid), length(grid_vals))
  testthat::expect_equal(length(prof$Rceattle_list), length(grid_vals))
  testthat::expect_equal(length(prof$nll), length(grid_vals))
  testthat::expect_equal(prof$grid$slot_1, grid_vals)

  # The profiled cell should be log(grid_val) in each fit's estimated params
  for (i in seq_along(grid_vals)) {
    fit <- prof$Rceattle_list[[i]]
    if (is.null(fit)) next
    testthat::expect_equal(
      as.numeric(fit$estimated_params$R_log_sd[1]),
      log(grid_vals[i]),
      tolerance = 1e-6
    )
  }
})


testthat::test_that("profile: sigmaR alias matches raw R_log_sd call", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  library(Rceattle)
  data(BS2017SS)

  ss_run <- fit_mod(
    data_list   = BS2017SS,
    inits       = NULL,
    file        = NULL,
    estimateMode = 1,
    random_rec  = FALSE,
    msmMode     = 0,
    fit_control = fit_control(phase = TRUE, getsd = FALSE, verbose = 0)
  )

  grid_vals <- c(0.3, 0.9)

  # Natural-scale alias
  p_alias <- profile(ss_run, param = "sigmaR",
                           slots = list(1), values = list(grid_vals),
                           cores = 1)

  # Equivalent raw form: pass log-scale values with identity transform
  p_raw <- profile(ss_run, param = "R_log_sd",
                         slots = list(1), values = list(log(grid_vals)),
                         transform = "identity", cores = 1)

  for (i in seq_along(grid_vals)) {
    expected <- log(grid_vals[i])
    if (!is.null(p_alias$Rceattle_list[[i]])) {
      testthat::expect_equal(
        as.numeric(p_alias$Rceattle_list[[i]]$estimated_params$R_log_sd[1]),
        expected, tolerance = 1e-6
      )
    }
    if (!is.null(p_raw$Rceattle_list[[i]])) {
      testthat::expect_equal(
        as.numeric(p_raw$Rceattle_list[[i]]$estimated_params$R_log_sd[1]),
        expected, tolerance = 1e-6
      )
    }
  }

  # Same NLL up to optimiser noise
  testthat::expect_equal(p_alias$nll, p_raw$nll, tolerance = 1e-3)
})


testthat::test_that("profile: alpha alias fills in rec_pars column", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  library(Rceattle)
  data(BS2017SS)

  ss_run <- fit_mod(
    data_list   = BS2017SS,
    inits       = NULL,
    file        = NULL,
    estimateMode = 1,
    random_rec  = FALSE,
    msmMode     = 0,
    fit_control = fit_control(phase = TRUE, getsd = FALSE, verbose = 0)
  )

  alpha_vals <- c(5, 15)

  prof <- profile(ss_run, param = "alpha",
                        slots = list(1), values = list(alpha_vals),
                        cores = 1)

  for (i in seq_along(alpha_vals)) {
    fit <- prof$Rceattle_list[[i]]
    if (is.null(fit)) next
    # Column 2 of rec_pars is alpha; alias should have set [1, 2]
    testthat::expect_equal(
      as.numeric(fit$estimated_params$rec_pars[1, 2]),
      log(alpha_vals[i]), tolerance = 1e-6
    )
  }
})


testthat::test_that("profile: alias warns when transform is overridden", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  library(Rceattle)
  data(BS2017SS)

  ss_run <- fit_mod(
    data_list   = BS2017SS,
    inits       = NULL,
    file        = NULL,
    estimateMode = 1,
    random_rec  = FALSE,
    msmMode     = 0,
    fit_control = fit_control(phase = TRUE, getsd = FALSE, verbose = 0)
  )

  testthat::expect_warning(
    profile(ss_run, param = "sigmaR",
                  slots = list(1), values = list(0.5),
                  transform = "identity", cores = 1),
    "ignoring the supplied `transform`"
  )
})


testthat::test_that("profile: defaults slots to species 1 with a warning", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  library(Rceattle)
  data(BS2017SS)

  ss_run <- fit_mod(
    data_list   = BS2017SS,
    inits       = NULL,
    file        = NULL,
    estimateMode = 1,
    random_rec  = FALSE,
    msmMode     = 0,
    fit_control = fit_control(phase = TRUE, getsd = FALSE, verbose = 0)
  )

  # sigmaR alias: default slot should be list(1) (species 1)
  prof <- testthat::expect_warning(
    profile(ss_run, param = "sigmaR",
                  values = list(c(0.4, 0.8)), cores = 1),
    "defaulting to species 1"
  )
  testthat::expect_equal(prof$slots, list(1L))

  # alpha alias: user slot dim is 1 (column appended); default still list(1)
  prof_a <- testthat::expect_warning(
    profile(ss_run, param = "alpha",
                  values = list(c(5, 10)), cores = 1),
    "defaulting to species 1"
  )
  testthat::expect_equal(prof_a$slots, list(c(1L, 2L)))

  # Default requires a single grid; reject if values has length > 1
  testthat::expect_error(
    profile(ss_run, param = "sigmaR",
                  values = list(c(0.4, 0.8), c(0.4, 0.8))),
    "species-1 default"
  )
})


testthat::test_that("profile: rec_pars alias rejects multi-element slots", {
  testthat::skip_if_not_installed("Rceattle")
  library(Rceattle)
  data(BS2017SS)

  ss_run <- fit_mod(
    data_list   = BS2017SS,
    inits       = NULL,
    file        = NULL,
    estimateMode = 1,
    random_rec  = FALSE,
    msmMode     = 0,
    fit_control = fit_control(phase = TRUE, getsd = FALSE, verbose = 0)
  )

  testthat::expect_error(
    profile(ss_run, param = "alpha",
                  slots = list(c(1, 2)), values = list(0.5),
                  cores = 1),
    "single species index"
  )
})


testthat::test_that("profile: input validation", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  library(Rceattle)

  # Not an Rceattle object
  testthat::expect_error(
    profile(fitted = list(), param = "R_log_sd",
            slots = list(1), values = list(1)),
    "no applicable method"
  )

  data(BS2017SS)
  ss_run <- fit_mod(
    data_list   = BS2017SS,
    inits       = NULL,
    file        = NULL,
    estimateMode = 1,
    random_rec  = FALSE,
    msmMode     = 0,
    fit_control = fit_control(phase = TRUE, getsd = FALSE, verbose = 0)
  )

  # Unknown parameter
  testthat::expect_error(
    profile(ss_run, param = "no_such_param",
                  slots = list(1), values = list(0.5)),
    "not found"
  )

  # Slot length mismatch: R_log_sd is 1-D, c(1, 2) is 2-D
  testthat::expect_error(
    profile(ss_run, param = "R_log_sd",
                  slots = list(c(1, 2)), values = list(0.5)),
    "dimension"
  )

  # values length must match slots length
  testthat::expect_error(
    profile(ss_run, param = "R_log_sd",
                  slots = list(1), values = list(0.5, 0.6)),
    "same length"
  )

  # Bad transform
  testthat::expect_error(
    profile(ss_run, param = "R_log_sd",
                  slots = list(1), values = list(0.5),
                  transform = "exp"),
    "transform"
  )
})
