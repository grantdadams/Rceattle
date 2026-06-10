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

  # One obs_ctl row per included aggregate observation.
  testthat::expect_equal(sum(dl$obs_ctl$type %in% c("index", "catch")),
                         length(inc_i) + length(inc_c))

  # Composition: bin counts = (proportion + 1e-5) * Sample_size, with one
  # decomposition group per row and the final bin flagged is_last_bin.
  if (nrow(dl$comp_obs) > 0) {
    r1 <- which(dl$comp_obsvec_idx >= 0)[1]
    sp <- dl$comp_ctl[r1, 2]
    n_comp <- dl$nages[sp]   # make_test_data: single-sex age comps
    start  <- dl$comp_obsvec_idx[r1]
    expected <- (as.numeric(dl$comp_obs[r1, seq_len(n_comp)]) + 0.00001) *
      dl$comp_n[r1, 2]
    testthat::expect_equal(dl$obsvec[start + seq_len(n_comp)], expected)

    comp_ctl_rows <- dl$obs_ctl[dl$obs_ctl$type == "comp", ]
    testthat::expect_true("length" %in% names(dl$obs_ctl))
    # exactly one last-bin (dropped) per composition group
    testthat::expect_equal(sum(comp_ctl_rows$is_last_bin),
                           length(unique(comp_ctl_rows$group_id)))
  }
})


testthat::test_that("build_osa_data lays CAAL observations into obsvec correctly", {
  testthat::skip_if_not_installed("Rceattle")

  # Parametric growth populates conditional age-at-length (CAAL) data.
  dat <- make_test_data(nyrs = 8, nages = 5, seed = 123,
                        growth = "vonBertalanffy")
  dl  <- Rceattle::rearrange_data(dat)
  testthat::skip_if(nrow(dl$caal_obs) == 0, "no CAAL data in fixture")

  testthat::expect_true("length" %in% names(dl$obs_ctl))

  r1     <- which(dl$caal_obsvec_idx >= 0)[1]
  sp     <- dl$caal_ctl[r1, 2]
  n_caal <- dl$nages[sp]                 # CAAL bins span ages (no joint-sex)
  start  <- dl$caal_obsvec_idx[r1]
  expected <- (as.numeric(dl$caal_obs[r1, seq_len(n_caal)]) + 0.00001) *
    dl$caal_n[r1, 1]
  testthat::expect_equal(dl$obsvec[start + seq_len(n_caal)], expected)

  caal_rows <- dl$obs_ctl[dl$obs_ctl$type == "caal", ]
  testthat::expect_equal(sum(caal_rows$is_last_bin),
                         length(unique(caal_rows$group_id)))
  testthat::expect_true(all(!is.na(caal_rows$length)))           # conditioning length
  testthat::expect_true(all(!is.na(caal_rows$age_or_length)))    # age bin
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

  osa <- osa_residuals(fit, types = c("index", "catch", "comp"))

  testthat::expect_s3_class(osa, "rceattle_osa")
  testthat::expect_true(all(c("type", "fleet", "year", "residual") %in% names(osa)))
  testthat::expect_true(all(is.finite(osa$residual)))
  testthat::expect_setequal(unique(osa$type), c("index", "catch", "comp"))

  # Aggregate: one residual per included catch/index observation.
  agg <- osa[osa$type %in% c("index", "catch"), ]
  testthat::expect_equal(nrow(agg),
                         sum(fit$obs_ctl$type %in% c("index", "catch")))

  # Composition: each composition drops its sum-to-N last bin, so there is one
  # fewer residual than bins per composition (here per fleet x year group).
  n_comp_bins   <- sum(fit$obs_ctl$type == "comp")
  n_comp_groups <- length(unique(fit$obs_ctl$group_id[fit$obs_ctl$type == "comp"]))
  testthat::expect_equal(sum(osa$type == "comp"), n_comp_bins - n_comp_groups)
  testthat::expect_true(all(is.finite(osa$residual[osa$type == "comp"])))

  # Deterministic with a fixed seed (cheap aggregate path).
  a1 <- osa_residuals(fit, types = "catch")
  a2 <- osa_residuals(fit, types = "catch")
  testthat::expect_equal(a1$residual, a2$residual)

  # Diagnostics return one row per source plus an overall row, with the null
  # intervals populated.
  diag <- osa_diagnostics(osa)
  testthat::expect_true(nrow(diag) >= 2)
  testthat::expect_true(all(c("sdnr", "sdnr_lo", "sdnr_hi") %in% names(diag)))
  testthat::expect_true(is.finite(diag$sdnr[diag$source == "all"]))

  # residuals(type = "osa") reshapes into the common residual schema (aggregate
  # subset to keep the test cheap).
  r <- residuals(fit, type = "osa", types = c("index", "catch"))
  testthat::expect_true(all(c("Source", "Fleet_code", "Year", "Residual") %in% names(r)))
  testthat::expect_equal(nrow(r), sum(fit$obs_ctl$type %in% c("index", "catch")))

  # OSA must refuse a debug (estimateMode >= 3) fit.
  fit_dbg <- Rceattle::fit_mod(data_list = GOApollock, inits = NULL,
                               estimateMode = 3,
                               fit_control = fit_control(phase = FALSE, verbose = 0))
  testthat::expect_error(osa_residuals(fit_dbg), "estimateMode")
})
