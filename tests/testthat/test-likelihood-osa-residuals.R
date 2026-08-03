# One-step-ahead (OSA) residuals -- Phase 1 (aggregate catch + index).
#
# A gold-standard "joint nll is numerically unchanged by the obsvec/keep
# refactor" check (comparing against a DLL built from the pre-change source)
# lives in dev/osa_phase1_check.R -- it is too slow for routine testing
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
  testthat::expect_true(all(c("obs_pos", "source", "data_row", "fleet_code",
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
  testthat::expect_equal(sum(dl$obs_ctl$source %in% c("index", "catch")),
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

    comp_ctl_rows <- dl$obs_ctl[dl$obs_ctl$source == "comp", ]
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
  testthat::expect_false(any(fast$obs_ctl$source %in% c("comp", "caal", "diet")))

  # Aggregate obsvec entries are identical between the fast and full builds, so
  # the fitted objective is unaffected by the toggle.
  ii <- fast$index_obsvec_idx[fast$index_obsvec_idx >= 0] + 1L
  testthat::expect_equal(fast$obsvec[ii], full$obsvec[ii])

  # The full build does produce composition entries when comp data are present.
  if (nrow(full$comp_obs) > 0) {
    testthat::expect_true(any(full$comp_obsvec_idx >= 0))
    testthat::expect_true(any(full$obs_ctl$source == "comp"))
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

  caal_rows <- dl$obs_ctl[dl$obs_ctl$source == "caal", ]
  testthat::expect_equal(sum(caal_rows$is_last_bin),
                         length(unique(caal_rows$group_id)))
  testthat::expect_true(all(!is.na(caal_rows$length)))           # conditioning length
  testthat::expect_true(all(!is.na(caal_rows$age_length_bin)))    # age bin
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

  diet_rows <- dl$obs_ctl[dl$obs_ctl$source == "diet", ]
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


testthat::test_that("fitted objective is unchanged (golden jnll on BS2017SS)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Total negative log-likelihood at the initial parameters of a fixed bundled
  # model. This guards that the fitted objective is unchanged -- in particular
  # that the OSA obsvec/keep machinery leaves normal fitting bit-for-bit
  # identical. estimateMode = 3 builds the model and reports jnll_comp without
  # optimizing, so the value is deterministic. Update the golden value ONLY when
  # the likelihood is intentionally changed. Uses the default comp_offset = 1e-5
  # (the historical composition proportion offset); set comp_offset = 0 for a
  # WHAM-style multinomial.
  data("BS2017SS", package = "Rceattle")
  fit <- Rceattle::fit_mod(BS2017SS, estimateMode = 3, msmMode = 0,
                           fit_control = fit_control(phase = FALSE, verbose = 0))
  jc <- fit$obj$report(fit$obj$par)$jnll_comp
  testthat::expect_equal(sum(jc[is.finite(jc)]), 1537036.2876293703,
                         tolerance = 1e-6)
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
  testthat::expect_equal(unique(pr$source), "recruitment")
  testthat::expect_true(all(is.finite(pr$residual)))
  # One recruitment residual per hindcast year.
  testthat::expect_equal(nrow(pr), length(GOApollock$styr:GOApollock$endyr))
  testthat::expect_setequal(range(pr$year), c(GOApollock$styr, GOApollock$endyr))

  # Deterministic with a fixed seed.
  pr2 <- process_residuals(fit, process = "recruitment")
  testthat::expect_equal(pr$residual, pr2$residual)

  # Does not clobber the caller's global RNG: the posterior draw is seeded on a
  # local stream and .Random.seed is restored on exit.
  set.seed(101L)
  before <- get(".Random.seed", envir = .GlobalEnv)
  invisible(process_residuals(fit, process = "recruitment", seed = 7))
  testthat::expect_identical(before, get(".Random.seed", envir = .GlobalEnv))

  # process = "all" returns every supported process present, all finite. On
  # GOApollock the catchability deviates use a random-walk prior, so the
  # marginal-SD standardization is approximate and emits a warning.
  testthat::expect_warning(pr_all <- process_residuals(fit, process = "all"),
                           "approximate")
  testthat::expect_true("recruitment" %in% pr_all$source)
  testthat::expect_true(all(is.finite(pr_all$residual)))

  # Recruitment-only residuals use an exact iid prior -> no approximation warning.
  testthat::expect_no_warning(process_residuals(fit, process = "recruitment"))

  # residuals(type = "process") reshapes into the common residual schema.
  r <- residuals(fit, type = "process")
  testthat::expect_true(all(c("Source", "Species", "Year", "Residual") %in% names(r)))

  # Composition OSA data is built on demand, so comp residuals work from any fit
  # -- no fit_control(osa = TRUE) needed.
  osa_comp <- osa_residuals(fit, source = "comp")
  testthat::expect_s3_class(osa_comp, "rceattle_osa")
  testthat::expect_true(all(is.finite(osa_comp$residual)))
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

  osa <- osa_residuals(fit, source = c("index", "catch", "comp"))

  testthat::expect_s3_class(osa, "rceattle_osa")
  testthat::expect_true(all(c("source", "fleet", "year", "residual") %in% names(osa)))
  testthat::expect_true(all(is.finite(osa$residual)))
  testthat::expect_setequal(unique(osa$source), c("index", "catch", "comp"))

  # Full observation map, regenerated on demand (as osa_residuals() does
  # internally); the stored fit$obs_ctl carries only the aggregate series.
  ctl <- Rceattle:::build_osa_data(fit$obj$env$data, build_osa = TRUE)$obs_ctl

  # Aggregate: one residual per included catch/index observation.
  agg <- osa[osa$source %in% c("index", "catch"), ]
  testthat::expect_equal(nrow(agg),
                         sum(ctl$source %in% c("index", "catch")))

  # Composition: each composition drops its sum-to-N last bin, so there is one
  # fewer residual than bins per composition (here per fleet x year group).
  n_comp_bins   <- sum(ctl$source == "comp")
  n_comp_groups <- length(unique(ctl$group_id[ctl$source == "comp"]))
  testthat::expect_equal(sum(osa$source == "comp"), n_comp_bins - n_comp_groups)
  testthat::expect_true(all(is.finite(osa$residual[osa$source == "comp"])))

  # Deterministic with a fixed seed (cheap aggregate path).
  a1 <- osa_residuals(fit, source = "catch")
  a2 <- osa_residuals(fit, source = "catch")
  testthat::expect_equal(a1$residual, a2$residual)

  # Diagnostics return one row per source plus an overall row, with the null
  # intervals populated.
  diag <- osa_diagnostics(osa)
  testthat::expect_true(nrow(diag) >= 2)
  testthat::expect_true(all(c("sdnr", "sdnr_lo", "sdnr_hi") %in% names(diag)))
  testthat::expect_true(is.finite(diag$sdnr[diag$group == "all"]))

  # residuals(type = "osa") reshapes into the common residual schema; `source`
  # selects the data source(s) (aggregate subset to keep the test cheap).
  r <- residuals(fit, type = "osa", source = c("index", "catch"))
  testthat::expect_true(all(c("Source", "Fleet_code", "Year", "Residual") %in% names(r)))
  testthat::expect_equal(nrow(r), sum(ctl$source %in% c("index", "catch")))

  # glm-style residual kinds: response / pearson cover all sources by default.
  rr <- residuals(fit, type = "response")
  testthat::expect_true(all(c("index", "catch", "comp") %in% rr$Source))
  testthat::expect_true(all(is.finite(rr$Residual)))
  rp <- residuals(fit, type = "pearson", source = "index")
  testthat::expect_setequal(unique(rp$Source), "index")
  testthat::expect_true(all(is.finite(rp$Residual)))

  # `type` selects the residual kind only; passing a data-source name (use
  # `source` for that) is an error.
  testthat::expect_error(residuals(fit, type = "index"))

  # species filter (GOApollock is a single-species model -> species 1).
  testthat::expect_setequal(unique(residuals(fit, species = 1)$Species), 1)
  testthat::expect_equal(nrow(residuals(fit, species = 999)), 0L)

  # plot() builds a separate aggregate (Q-Q only) and composition (Q-Q + OSA and
  # Pearson bubbles) figure. The osa object carries fleet_name and the Pearson
  # residuals needed for the composition bubbles.
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    testthat::expect_true("fleet_name" %in% names(osa))
    testthat::expect_false(is.null(attr(osa, "pearson")))
    pf <- tempfile(fileext = ".pdf"); grDevices::pdf(pf)
    pl    <- plot(osa)
    p_idx <- plot(osa, source = "index")                # source filter
    p_sep <- plot(osa, source = "comp", combine = FALSE) # split age/length
    p_sp  <- plot(osa, species = 1)                     # species filter
    pc    <- plot_comp(fit)                             # ggplot2 comp figures
    grDevices::dev.off()
    testthat::expect_true(all(c("aggregate", "composition") %in% names(pl)))
    testthat::expect_equal(names(p_idx), "aggregate")
    testthat::expect_true("composition_age" %in% names(p_sep))
    testthat::expect_false("composition" %in% names(p_sep))
    testthat::expect_true(length(p_sp) >= 1L)
    testthat::expect_true("pearson" %in% names(pc))
    testthat::expect_true(any(grepl("^annual", names(pc))))
    testthat::expect_true(any(grepl("^aggregated", names(pc))))
  }

  # OSA must refuse a debug (estimateMode >= 3) fit.
  fit_dbg <- Rceattle::fit_mod(data_list = GOApollock, inits = NULL,
                               estimateMode = 3,
                               fit_control = fit_control(phase = FALSE, verbose = 0))
  testthat::expect_error(osa_residuals(fit_dbg), "estimateMode")
})


testthat::test_that("residuals(source = 'diet') computes diet Pearson residuals", {
  testthat::skip_if_not_installed("Rceattle")

  # Minimal fitted-model stand-in carrying just the diet pieces residuals() uses.
  dd <- data.frame(Pred = c(1L, 1L, 2L), Pred_sex = 0L, Prey = c(1L, 2L, 1L),
                   Prey_sex = 0L, Pred_age = c(3L, 3L, 4L), Prey_age = c(1L, 2L, 1L),
                   Year = c(2000L, 2000L, 2001L),
                   Stomach_proportion_by_weight = c(0.6, 0.3, 0.5),
                   Sample_size = c(50, 50, 40))
  fit <- structure(list(
    data_list  = list(diet_data = dd, spnames = c("A", "B")),
    quantities = list(diet_hat = cbind(NA_real_, c(0.55, 0.35, 0.45)))),
    class = "Rceattle")

  r <- residuals(fit, type = "pearson", source = "diet")
  testthat::expect_setequal(unique(r$Source), "diet")
  testthat::expect_true(all(c("Pred_sex", "Prey", "Prey_sex", "Pred_age", "Prey_age",
                              "Observed", "Fitted", "Residual") %in% names(r)))
  hat <- fit$quantities$diet_hat[, 2]
  testthat::expect_equal(r$Residual, (dd$Stomach_proportion_by_weight - hat) /
                           sqrt(hat * (1 - hat) / dd$Sample_size))
  testthat::expect_equal(residuals(fit, type = "response", source = "diet")$Residual,
                         dd$Stomach_proportion_by_weight - hat)
  # diet uses a predator/prey schema, so it must be requested on its own.
  testthat::expect_error(residuals(fit, source = c("diet", "comp")), "on its own")
  # species filter acts on the predator species.
  testthat::expect_setequal(unique(residuals(fit, source = "diet", species = 1)$Species),
                            1L)
})


testthat::test_that("diet residuals and plot_diet_comp run on a fitted diet model", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Small multispecies model with diet (same fixture as test-diet-likelihood.R);
  # estimateMode = 3 builds the report (diet_hat) without optimizing.
  nyrs <- 10L; nspp <- 2L
  Fmort  <- c(seq(0.02, 0.3, length.out = nyrs / 2), seq(0.3, 0.05, length.out = nyrs / 2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  sim <- make_msm_test_data(
    years = seq_len(nyrs),
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),
    gam_a = c(1, 0.1), gam_b = rep(0.15, nspp),
    log_phi = matrix(c(-5, 0.5, -10, -2), nspp, nspp, byrow = TRUE))
  fit <- Rceattle::fit_mod(sim$data_list, estimateMode = 3, msmMode = 1,
                           suitMode = 4, niter = 5, initMode = "NonEquilibrium",
                           fit_control = fit_control(phase = FALSE, verbose = 0))

  testthat::expect_false(is.null(fit$quantities$diet_hat))
  r <- residuals(fit, type = "pearson", source = "diet")
  testthat::expect_setequal(unique(r$Source), "diet")
  testthat::expect_equal(nrow(r), nrow(fit$data_list$diet_data))
  # The diet Pearson matches the proportion formula (finite rows).
  dd  <- fit$data_list$diet_data
  hat <- fit$quantities$diet_hat[, 2]
  fin <- is.finite(r$Residual)
  testthat::expect_equal(r$Residual[fin],
    ((dd$Stomach_proportion_by_weight - hat) /
       sqrt(hat * (1 - hat) / dd$Sample_size))[fin])

  # plot_diet_comp() now sources its residuals from residuals(source = "diet").
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    pf <- tempfile(fileext = ".pdf"); grDevices::pdf(pf)
    testthat::expect_error(plot_diet_comp(fit), NA)
    grDevices::dev.off()
  }
})


# A small multispecies diet fixture reused by the two OSA-diet tests below.
# Diet is multinomial (Diet_distribution = 0) with unit weights, so the OSA
# (conditional-binomial) decomposition of the diet likelihood must reproduce the
# ordinary multinomial diet likelihood slot exactly (see test (a)).
.make_diet_osa_fixture <- function() {
  nyrs <- 10L; nspp <- 2L
  Fmort  <- c(seq(0.02, 0.3, length.out = nyrs / 2),
              seq(0.3, 0.05, length.out = nyrs / 2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  make_msm_test_data(
    years   = seq_len(nyrs),
    Fmort   = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),
    gam_a   = c(1, 0.1), gam_b = rep(0.15, nspp),
    log_phi = matrix(c(-5, 0.5, -10, -2), nspp, nspp, byrow = TRUE))
}


testthat::test_that("OSA mode reproduces the diet fitting likelihood (decomposition invariant)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Diet has no WHAM reference (WHAM has no diet data), so we validate the diet
  # OSA path against the model itself. With keep == 1 the conditional-binomial
  # decomposition (dmultinom_osa) sums back to the joint multinomial density, so
  # building the model in OSA mode (osa_mode = 1, reading the diet composition
  # from obsvec) must reproduce the diet likelihood slot of ordinary fitting
  # (osa_mode = 0, reading diet_obs). This is a deterministic check that the diet
  # obsvec plumbing -- diet_obsvec_idx, the n_prey + 1 "other prey" segment, and
  # the obs<->pred alignment -- is wired correctly. The fixture is multinomial
  # diet with unit weights, so the weighted slot equals the unweighted joint.
  sim <- .make_diet_osa_fixture()

  # estimateMode = 3 builds and reports jnll_comp without optimizing; the
  # decomposition identity holds at any parameter values, so no fit is needed.
  # The diet obsvec segment (read when osa_mode = 1) is rebuilt on demand.
  fit <- Rceattle::fit_mod(sim$data_list, estimateMode = 3, msmMode = 1,
                           suitMode = 4, niter = 5, initMode = "NonEquilibrium",
                           fit_control = fit_control(phase = FALSE, verbose = 0))
  ctl <- Rceattle:::build_osa_data(fit$obj$env$data, build_osa = TRUE)$obs_ctl
  testthat::skip_if(!any(ctl$source == "diet"), "no diet observations")

  diet_row <- 19L                       # jnll_comp C++ slot 18 (diet) -> R row 19
  jc0  <- fit$obj$report(fit$obj$par)$jnll_comp                 # osa_mode = 0
  obj1 <- Rceattle:::.osa_build_obj(fit)                        # rebuild osa_mode = 1
  jc1  <- obj1$report(obj1$par)$jnll_comp

  testthat::expect_gte(nrow(jc0), diet_row)
  # The diet likelihood is active (non-trivial) ...
  testthat::expect_true(sum(jc0[diet_row, ]) != 0)
  # ... and the OSA decomposition reproduces it to ~machine precision.
  testthat::expect_equal(sum(jc1[diet_row, ]), sum(jc0[diet_row, ]),
                         tolerance = 1e-8)
})


testthat::test_that("osa_residuals(source = 'diet') runs end-to-end on a fitted diet model", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # End-to-end exercise of oneStepPredict() over the diet obsvec positions (never
  # run elsewhere in the suite). osa_residuals requires an optimized
  # (estimateMode < 3) fit; the simulated data are generated from the model with
  # a large effective sample size, so estimateMode = 1 converges. The diet
  # composition obsvec is rebuilt on demand by osa_residuals(). getsd = FALSE
  # skips sdreport -- OSA residuals use the fitted object, not the standard
  # errors, and this tiny fixture has a non-positive-definite Hessian that would
  # otherwise warn.
  sim <- .make_diet_osa_fixture()
  fit <- Rceattle::fit_mod(sim$data_list, estimateMode = 1, msmMode = 1,
                           suitMode = 4, niter = 5, initMode = "NonEquilibrium",
                           fit_control = fit_control(phase = FALSE, verbose = 0,
                                                     getsd = FALSE))
  ctl <- Rceattle:::build_osa_data(fit$obj$env$data, build_osa = TRUE)$obs_ctl
  testthat::skip_if(!any(ctl$source == "diet"), "no diet observations")

  osa <- osa_residuals(fit, source = "diet")

  testthat::expect_s3_class(osa, "rceattle_osa")
  testthat::expect_setequal(unique(osa$source), "diet")

  # Each stomach drops its sum-to-one "other prey" last bin, so the residual
  # count equals the kept (non-last) diet bins in obs_ctl.
  n_keep <- sum(ctl$source == "diet" & !ctl$is_last_bin)
  testthat::expect_equal(nrow(osa), n_keep)
  testthat::expect_true(all(is.finite(osa$residual)))

  # Deterministic for a fixed seed (the default continuous method has no RNG).
  osa2 <- osa_residuals(fit, source = "diet")
  testthat::expect_equal(osa$residual, osa2$residual)

  # Diagnostics run on the diet source and return a finite SDNR. The fixture's
  # diet proportions are noiseless (expected values, not multinomial draws), so
  # an SDNR ~ 1 calibration is not expected here -- that needs sampled data and
  # belongs in an offline validation; here we only guard that the machinery runs
  # and produces a finite statistic.
  dg <- osa_diagnostics(osa)
  testthat::expect_true(is.finite(dg$sdnr[dg$group == "all"]))
})

testthat::test_that("OSA supports an MVN index fleet without warning (per-family coverage in test-likelihood-osa-index-families.R)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # OSA residuals now cover every index family: the MVN covariance block is
  # whitened by the Cholesky of its covariance so oneStepPredict() can residualize
  # it (the correctness oracle lives in test-likelihood-osa-index-families.R). Here
  # we only pin that build_osa_data() no longer excludes/ warns for an MVN fleet
  # and that its index observations are laid into the OSA obsvec.
  nyrs <- 8; nages <- 5
  dat  <- make_test_data(nyrs = nyrs, nages = nages, seed = 42)
  sds  <- rep(20, nyrs); Rho <- matrix(0.3, nyrs, nyrs); diag(Rho) <- 1
  Sigma <- diag(sds) %*% Rho %*% diag(sds)
  srv  <- dat$fleet_control$Fleet_name == "Survey"       # Survey == Fleet_code 1
  dat$fleet_control$Index_distribution[srv] <- "MVN"
  dat$fleet_control$Catchability[srv]       <- "AnalyticalArith"
  dat$index_cov <- list(Survey = Sigma)

  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0,
    fit_control = fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))

  # The OSA build no longer warns, and the MVN Survey's fitted observations get
  # whitened obsvec entries (one per fitted year).
  testthat::expect_no_warning(
    osa_dat <- Rceattle:::build_osa_data(fit$obj$env$data, build_osa = TRUE))
  testthat::expect_equal(sum(osa_dat$obs_ctl$source == "index"), nyrs)
  testthat::expect_true(all(osa_dat$index_obsvec_idx >= 0L))

  osa <- suppressWarnings(osa_residuals(fit, source = c("index", "comp")))
  testthat::expect_s3_class(osa, "rceattle_osa")
  testthat::expect_true(all(is.finite(osa$residual)))
  testthat::expect_true(any(osa$source == "index"))
})
