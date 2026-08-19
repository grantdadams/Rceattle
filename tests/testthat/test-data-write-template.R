# write_template() emits a minimal, structurally complete single-species
# workbook that round-trips through read_data()/data_check() and builds.

testthat::test_that("write_template writes a workbook that round-trips to canonical names", {
  testthat::skip_if_not_installed("writexl")
  testthat::skip_if_not_installed("readxl")
  testthat::skip_if_not_installed("Rceattle")

  f <- withr::local_tempfile(fileext = ".xlsx")
  d <- suppressMessages(suppressWarnings(Rceattle::write_template(f)))
  testthat::expect_true(file.exists(f))

  back <- suppressMessages(suppressWarnings(Rceattle::read_data(f)))
  # Canonical column names throughout, and no deprecated spelling survives.
  fc <- back$fleet_control
  testthat::expect_true(all(c("Catchability_index", "Catchability_init",
                              "Comp_distribution", "Sel_norm_bin",
                              "Observation_units", "Proj_F_proportion") %in% names(fc)))
  testthat::expect_false(any(c("Q_index", "Comp_loglike", "Sel_norm_bin1",
                               "Weight1_Numbers2", "proj_F_prop") %in% names(fc)))
  testthat::expect_equal(back$nspp, 1)
  testthat::expect_equal(nrow(fc), 2)
})

testthat::test_that("the write_template workbook builds under fit_mod(estimateMode = 3)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  f <- withr::local_tempfile(fileext = ".xlsx")
  suppressMessages(suppressWarnings(Rceattle::write_template(f)))
  back <- suppressMessages(suppressWarnings(Rceattle::read_data(f)))

  # End-to-end: fit_mod() runs clean_data -> switch_check -> data_check ->
  # build, so a finite objective proves the template passes validation and
  # builds under TMB.
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = back, file = NULL, estimateMode = 3, random_rec = FALSE,
    msmMode = 0, fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))
  obj <- fit$opt$objective %||% sum(fit$quantities$jnll_comp)
  testthat::expect_true(is.finite(obj))
})
