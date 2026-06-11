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
  dl  <- Rceattle::rearrange_data(dat, build_osa = TRUE)

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


testthat::test_that("default fit skips composition OSA data (fast path)", {
  testthat::skip_if_not_installed("Rceattle")

  dat  <- make_test_data(nyrs = 8, nages = 5, seed = 123)
  fast <- Rceattle::rearrange_data(dat)                    # build_osa = FALSE (default)
  full <- Rceattle::rearrange_data(dat, build_osa = TRUE)

  # Fast path: aggregate index/catch entries are still built (the template
  # always reads them), but no composition / caal / diet metadata.
  testthat::expect_true(all(fast$comp_obsvec_idx < 0))
  testthat::expect_true(all(fast$caal_obsvec_idx < 0))
  testthat::expect_false(any(fast$obs_ctl$type %in% c("comp", "caal", "diet")))

  # Aggregate obsvec entries are identical between the fast and full builds, so
  # the fitted objective is unaffected by the toggle.
  ii <- fast$index_obsvec_idx[fast$index_obsvec_idx >= 0] + 1L
  testthat::expect_equal(fast$obsvec[ii], full$obsvec[ii])

  # The full build does produce composition entries when comp data are present.
  if (nrow(full$comp_obs) > 0) {
    testthat::expect_true(any(full$comp_obsvec_idx >= 0))
    testthat::expect_true(any(full$obs_ctl$type == "comp"))
  }
})


testthat::test_that("build_osa_data lays CAAL observations into obsvec correctly", {
  testthat::skip_if_not_installed("Rceattle")

  # Parametric growth populates conditional age-at-length (CAAL) data.
  dat <- make_test_data(nyrs = 8, nages = 5, seed = 123,
                        growth = "vonBertalanffy")
  dl  <- Rceattle::rearrange_data(dat, build_osa = TRUE)
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


testthat::test_that("build_osa_data lays diet observations into obsvec correctly", {
  testthat::skip_if_not_installed("Rceattle")

  # Minimal data_list with two stomachs of predator species 1 (suitMode > 0):
  # stomach 0 has 2 prey items, stomach 1 has 1. Aggregate/comp/caal are empty.
  empty_i <- function(nc) matrix(integer(0), 0, nc)
  empty_n <- function(nc) matrix(numeric(0), 0, nc)
  dl <- list(
    endyr = 5, flt_type = c(1L, 2L), nages = c(3L, 3L), nlengths = c(3L, 3L),
    index_ctl = empty_i(3), index_obs = empty_n(2),
    catch_ctl = empty_i(3), catch_obs = empty_n(2),
    comp_ctl = empty_i(5), comp_obs = empty_n(3), comp_n = empty_n(2),
    caal_ctl = empty_i(5), caal_obs = empty_n(3), caal_n = empty_n(1),
    n_stomach_obs = 2L,
    stomach_id = c(0L, 0L, 1L),
    diet_ctl = matrix(c(1, 1, 0, 0, 1, 1, 3,
                        1, 2, 0, 0, 1, 1, 3,
                        1, 1, 0, 0, 2, 1, 4), nrow = 3, byrow = TRUE),
    diet_obs = matrix(c(50, 0.6, 50, 0.3, 40, 0.7), nrow = 3, byrow = TRUE),
    suitMode = c(2L, 0L))

  dl <- Rceattle:::build_osa_data(dl, build_osa = TRUE)

  testthat::expect_length(dl$diet_obsvec_idx, 2L)
  testthat::expect_true(all(dl$diet_obsvec_idx >= 0))

  diet_rows <- dl$obs_ctl[dl$obs_ctl$type == "diet", ]
  # stomach 0: 2 prey + other = 3 bins; stomach 1: 1 prey + other = 2 bins.
  testthat::expect_equal(nrow(diet_rows), 5L)
  testthat::expect_equal(sum(diet_rows$is_last_bin), 2L)   # one "other prey" per stomach
  testthat::expect_setequal(unique(diet_rows$stomach_id), c(0L, 1L))

  # Counts = (proportions + other + 1e-5) normalized, scaled to the sample size.
  v <- c(0.6, 0.3, 1 - 0.9) + 0.00001
  expected0 <- v / sum(v) * 50
  start0 <- dl$diet_obsvec_idx[1]
  testthat::expect_equal(dl$obsvec[start0 + seq_len(3L)], expected0)
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


testthat::test_that("process_residuals() runs on a converging model", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOApollock", package = "Rceattle")
  fit <- Rceattle::fit_mod(
    data_list = GOApollock, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(verbose = 0, phase = TRUE, getsd = TRUE))

  pr <- process_residuals(fit, process = "recruitment")
  testthat::expect_s3_class(pr, "rceattle_osa")
  testthat::expect_equal(unique(pr$type), "recruitment")
  testthat::expect_true(all(is.finite(pr$residual)))
  # One recruitment residual per hindcast year.
  testthat::expect_equal(nrow(pr), length(GOApollock$styr:GOApollock$endyr))
  testthat::expect_setequal(range(pr$year), c(GOApollock$styr, GOApollock$endyr))

  # Deterministic with a fixed seed.
  pr2 <- process_residuals(fit, process = "recruitment")
  testthat::expect_equal(pr$residual, pr2$residual)

  # process = "all" returns every supported process present, all finite. On
  # GOApollock the catchability deviates use a random-walk prior, so the
  # marginal-SD standardization is approximate and emits a warning.
  testthat::expect_warning(pr_all <- process_residuals(fit, process = "all"),
                           "approximate")
  testthat::expect_true("recruitment" %in% pr_all$type)
  testthat::expect_true(all(is.finite(pr_all$residual)))

  # Recruitment-only residuals use an exact iid prior -> no approximation warning.
  testthat::expect_no_warning(process_residuals(fit, process = "recruitment"))

  # residuals(type = "process") reshapes into the common residual schema.
  r <- residuals(fit, type = "process")
  testthat::expect_true(all(c("Source", "Species", "Year", "Residual") %in% names(r)))

  # This fit used the default osa = FALSE, so composition OSA data was not built;
  # osa_residuals() for composition must fail with an actionable message.
  testthat::expect_error(osa_residuals(fit, types = "comp"), "osa = TRUE")
})


testthat::test_that("osa_residuals() runs end-to-end on a converging model", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOApollock", package = "Rceattle")
  fit <- Rceattle::fit_mod(
    data_list = GOApollock, inits = NULL, file = NULL,
    estimateMode = 1, random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(verbose = 0, phase = TRUE, osa = TRUE))

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

  # residuals(type = "osa") reshapes into the common residual schema; `source`
  # selects the data source(s) (aggregate subset to keep the test cheap).
  r <- residuals(fit, type = "osa", source = c("index", "catch"))
  testthat::expect_true(all(c("Source", "Fleet_code", "Year", "Residual") %in% names(r)))
  testthat::expect_equal(nrow(r), sum(fit$obs_ctl$type %in% c("index", "catch")))

  # glm-style residual kinds: response / pearson cover all sources by default.
  rr <- residuals(fit, type = "response")
  testthat::expect_true(all(c("index", "catch", "comp") %in% rr$Source))
  testthat::expect_true(all(is.finite(rr$Residual)))
  rp <- residuals(fit, type = "pearson", source = "index")
  testthat::expect_setequal(unique(rp$Source), "index")
  testthat::expect_true(all(is.finite(rp$Residual)))

  # Legacy type = <source> still works with a deprecation warning.
  testthat::expect_warning(rl <- residuals(fit, type = "index"), "source")
  testthat::expect_setequal(unique(rl$Source), "index")

  # OSA must refuse a debug (estimateMode >= 3) fit.
  fit_dbg <- Rceattle::fit_mod(data_list = GOApollock, inits = NULL,
                               estimateMode = 3,
                               fit_control = fit_control(phase = FALSE, verbose = 0))
  testthat::expect_error(osa_residuals(fit_dbg), "estimateMode")
})
