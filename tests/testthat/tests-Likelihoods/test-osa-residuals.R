# One-step-ahead (OSA) residuals -- Phase 1 (aggregate catch + index).
#
# A gold-standard "joint nll is numerically unchanged by the obsvec/keep
# refactor" check (comparing against a DLL built from the pre-change source)
# lives in R/dev/osa_phase1_check.R -- it is too slow for routine testing
# because it recompiles a second DLL. Here we test the parts that guarantee
# correctness and are fast/CI-appropriate: the obsvec/obs_ctl construction and
# the end-to-end osa_residuals() machinery.

testthat::test_that("build_osa_data lays observations into obsvec correctly", {
  testthat::skip_if_not_installed("Rceattle")

  dat <- make_test_data(nyrs = 8, nages = 5, seed = 123)
  dl  <- Rceattle::rearrange_data(dat)

  # obs_ctl must be a DATA FRAME (regression guard: rearrange_data's
  # data.frame->matrix coercion must not sweep it up).
  testthat::expect_s3_class(dl$obs_ctl, "data.frame")
  testthat::expect_true(all(c("obs_pos", "type", "data_row", "fleet_code",
                              "year", "is_last_bin") %in% names(dl$obs_ctl)))

  # Defaults and types.
  testthat::expect_identical(dl$osa_mode, 0L)
  testthat::expect_type(dl$obsvec, "double")
  testthat::expect_type(dl$index_obsvec_idx, "integer")
  testthat::expect_type(dl$catch_obsvec_idx, "integer")

  endyr <- dl$endyr
  ft    <- dl$flt_type

  # Index: included rows are exactly those passing the TMB guard, and obsvec
  # stores log(observation) at the mapped position.
  inc_i <- which(dl$index_ctl[, 3] > 0 & dl$index_ctl[, 3] <= endyr &
                   ft[dl$index_ctl[, 1]] > 0 & dl$index_obs[, 1] > 0)
  testthat::expect_true(all(dl$index_obsvec_idx[inc_i] >= 0))
  testthat::expect_equal(dl$obsvec[dl$index_obsvec_idx[inc_i] + 1L],
                         log(dl$index_obs[inc_i, 1]))

  # Catch: same, with the fishery-only guard (flt_type == 1).
  inc_c <- which(dl$catch_ctl[, 3] > 0 & dl$catch_ctl[, 3] <= endyr &
                   ft[dl$catch_ctl[, 1]] == 1 & dl$catch_obs[, 1] > 0)
  testthat::expect_equal(dl$obsvec[dl$catch_obsvec_idx[inc_c] + 1L],
                         log(dl$catch_obs[inc_c, 1]))

  # One obs_ctl row per included observation.
  testthat::expect_equal(nrow(dl$obs_ctl), length(inc_i) + length(inc_c))
})


testthat::test_that("obsvec/keep refactor leaves the fitted objective finite", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  dat <- make_test_data(nyrs = 8, nages = 5, seed = 123)
  res <- Rceattle::fit_mod(data_list = dat, inits = NULL, estimateMode = 3,
                           fit_control = fit_control(phase = FALSE, verbose = 0))

  testthat::expect_false(is.null(res$obj))
  testthat::expect_true(is.finite(res$obj$fn()))
  jc <- res$obj$report(res$obj$env$last.par.best)$jnll_comp
  testthat::expect_true(is.finite(sum(jc[1, ])))   # index slot
  testthat::expect_true(is.finite(sum(jc[2, ])))   # catch slot

  # obs_ctl is carried onto the fitted object as a data frame.
  testthat::expect_s3_class(res$obs_ctl, "data.frame")
})


testthat::test_that("osa_residuals() runs end-to-end on a converging model", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOApollock", package = "Rceattle")
  fit <- Rceattle::fit_mod(
    data_list = GOApollock, inits = NULL, file = NULL,
    estimateMode = 1, random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(verbose = 0, phase = TRUE))

  osa <- osa_residuals(fit, types = c("index", "catch"))

  testthat::expect_s3_class(osa, "rceattle_osa")
  testthat::expect_true(all(c("type", "fleet", "year", "residual") %in% names(osa)))
  testthat::expect_true(all(is.finite(osa$residual)))
  testthat::expect_setequal(unique(osa$type), c("index", "catch"))

  # One residual per included aggregate observation.
  testthat::expect_equal(nrow(osa), nrow(fit$obs_ctl))

  # Deterministic with a fixed seed.
  osa2 <- osa_residuals(fit, types = c("index", "catch"))
  testthat::expect_equal(osa$residual, osa2$residual)

  # Diagnostics return one row per source plus an overall row, with the null
  # intervals populated.
  diag <- osa_diagnostics(osa)
  testthat::expect_true(nrow(diag) >= 2)
  testthat::expect_true(all(c("sdnr", "sdnr_lo", "sdnr_hi") %in% names(diag)))
  testthat::expect_true(is.finite(diag$sdnr[diag$source == "all"]))

  # residuals(type = "osa") reshapes into the common residual schema.
  r <- residuals(fit, type = "osa")
  testthat::expect_true(all(c("Source", "Fleet_code", "Year", "Residual") %in% names(r)))
  testthat::expect_equal(nrow(r), nrow(osa))

  # OSA must refuse a debug (estimateMode >= 3) fit.
  fit_dbg <- Rceattle::fit_mod(data_list = GOApollock, inits = NULL,
                               estimateMode = 3,
                               fit_control = fit_control(phase = FALSE, verbose = 0))
  testthat::expect_error(osa_residuals(fit_dbg), "estimateMode")
})
