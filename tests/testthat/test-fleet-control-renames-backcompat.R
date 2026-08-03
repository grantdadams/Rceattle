# =============================================================================
# Back-compat for the PR-5 intuitive-naming column renames.
#
# Every fleet_control column (and the data_list-level sigma_rec / Diet
# elements) renamed in the intuitive-naming pass keeps a deprecation alias that
# is upgraded in place by the schema-driven `.rce_upgrade_*_aliases()` helpers.
# These tests pin the released-package contract: an existing data list, .rda, or
# workbook that still uses the OLD names keeps reading and fitting IDENTICALLY.
# =============================================================================

# New canonical name -> a representative old name accepted on read.
.pr5_fc_new_to_old <- c(
  Catchability_index    = "Q_index",
  Catchability_init     = "Q_prior",
  Catchability_prior_sd = "Q_sd_prior",
  Comp_distribution     = "Comp_loglike",
  CAAL_distribution     = "CAAL_loglike",
  Index_distribution    = "Index_loglike",
  Sel_norm_bin          = "Sel_norm_bin1",
  Sel_norm_bin_upper    = "Sel_norm_bin2",
  Observation_units     = "Weight1_Numbers2",
  Proj_F_proportion     = "proj_F_prop"
)

# Rename canonical fleet_control columns (and sigma_rec) back to their old names.
.to_old_names <- function(d) {
  fc <- d$fleet_control
  for (new in names(.pr5_fc_new_to_old)) {
    j <- which(names(fc) == new)
    if (length(j)) names(fc)[j] <- unname(.pr5_fc_new_to_old[new])
  }
  d$fleet_control <- fc
  j <- which(names(d) == "sigma_rec")
  if (length(j)) names(d)[j] <- "sigma_rec_prior"
  d
}

testthat::test_that("switch_check upgrades every deprecated PR-5 column name", {
  testthat::skip_if_not_installed("Rceattle")
  # Only the canonical columns actually present on the fixture can be renamed
  # back (e.g. Index_distribution is defaulted at runtime, not stored).
  present_new <- intersect(names(Rceattle::BS2017SS$fleet_control),
                           names(.pr5_fc_new_to_old))
  old_present  <- unname(.pr5_fc_new_to_old[present_new])

  d_old <- .to_old_names(Rceattle::BS2017SS)
  testthat::expect_true(all(old_present %in% names(d_old$fleet_control)))

  out <- suppressMessages(suppressWarnings(Rceattle::switch_check(d_old)))
  fc <- out$fleet_control
  # Every canonical name present (switch_check also defaults the absent ones),
  # and NO deprecated name survives.
  testthat::expect_true(all(names(.pr5_fc_new_to_old) %in% names(fc)))
  testthat::expect_false(any(unname(.pr5_fc_new_to_old) %in% names(fc)))
  testthat::expect_true("sigma_rec" %in% names(out))
  testthat::expect_false("sigma_rec_prior" %in% names(out))
})

testthat::test_that("an old-name data list fits identically to the canonical one", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  ctl <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  fit_new <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, file = NULL, estimateMode = 3,
    random_rec = FALSE, msmMode = 0, fit_control = ctl)))
  fit_old <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = .to_old_names(Rceattle::BS2017SS), file = NULL, estimateMode = 3,
    random_rec = FALSE, msmMode = 0, fit_control = ctl)))

  obj_new <- sum(fit_new$quantities$jnll_comp)
  obj_old <- sum(fit_old$quantities$jnll_comp)
  testthat::expect_true(is.finite(obj_new))
  # Bit-identical: the alias upgrade is numerically inert.
  testthat::expect_equal(obj_old, obj_new, tolerance = 0)
})

testthat::test_that("a workbook with legacy fleet_control column names reads to canonical", {
  testthat::skip_if_not_installed("readxl")
  testthat::skip_if_not_installed("writexl")
  testthat::skip_if_not_installed("Rceattle")

  d <- make_test_data()
  f_new <- withr::local_tempfile(fileext = ".xlsx")
  suppressMessages(suppressWarnings(Rceattle::write_data(d, f_new)))

  # Rebuild the workbook verbatim but rename the fleet_control columns to their
  # deprecated spellings, mimicking a workbook saved before the rename.
  sheets <- readxl::excel_sheets(f_new)
  wb <- stats::setNames(
    lapply(sheets, function(s) as.data.frame(readxl::read_xlsx(f_new, sheet = s))),
    sheets)
  fc <- wb[["fleet_control"]]
  for (new in names(.pr5_fc_new_to_old)) {
    j <- which(names(fc) == new)
    if (length(j)) names(fc)[j] <- unname(.pr5_fc_new_to_old[new])
  }
  wb[["fleet_control"]] <- fc
  f_old <- withr::local_tempfile(fileext = ".xlsx")
  writexl::write_xlsx(wb, f_old)

  back <- suppressMessages(suppressWarnings(Rceattle::read_data(f_old)))
  fcb <- back$fleet_control
  # read_data upgraded every legacy name to canonical.
  testthat::expect_true(all(names(.pr5_fc_new_to_old) %in% names(fcb)))
  testthat::expect_false(any(unname(.pr5_fc_new_to_old) %in% names(fcb)))
})
