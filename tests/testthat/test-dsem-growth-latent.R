# A DSEM growth latent must land on the right growth_linkage_offset cell: right
# species, right sex, right growth parameter, right YEAR. Nothing else can check
# that -- the R-side index map (test-dsem-latent-registry.R) pins only that the
# values handed to the template are the intended ones, not where they arrive.
#
# Provenance of the per-parameter expectations, which are NOT the same for every
# growth parameter and were measured rather than assumed:
#
#   * L1    enters the LINEAR RAMP (growth.hpp, b_len = (l1 - Lmin)/age_L1), so
#           a spike at year t moves AGE 0 in year t itself, then one age per
#           year as that cohort ages.
#   * Linf  enters the PLUS GROUP (growth.hpp, diff = linf - current_size), so a
#           spike at year t moves the OLDEST age in year t, then the on-curve
#           ages from t+1 via laglinf.
#   * K     would move NO age in year t -- the Jan-1 recursion reads lagkappa --
#           which is why the plan's "first moves at t+1" is a K-specific
#           statement. K is not drivable, so it is not tested here.
#
# In every case successive year effects decay by exactly exp(-K), which is the
# acceptance criterion that catches a one-year index shift: the recursion already
# lags by one, so "the offset moved" alone cannot see it.

testthat::skip_on_cran()

.grow_fit <- function(sem, nyrs = 12L, nages = 8L) {
  e <- new.env(parent = asNamespace("Rceattle"))
  sys.source(testthat::test_path("helpers.R"), envir = e)
  d <- e$make_test_data(nyrs = nyrs, nprojyrs = 3L, nages = nages, minage = 1L,
                        growth = "vonBertalanffy", seed = 1)
  d$random_rec <- TRUE
  set.seed(4)
  d$env_data <- data.frame(Year = seq(d$styr, d$endyr), BT = stats::rnorm(nyrs))
  suppressWarnings(Rceattle::fit_mod(
    data_list = d, estimateMode = 3, msmMode = 0,
    growthFun = Rceattle::build_growth(fun = "vonBertalanffy"),
    dsem = Rceattle::build_DSEM(sem = sem),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0)))
}

# dsem_x_tj is column-major, so cell (t, j) is parameter (j-1)*nrow + t.
.spike <- function(fit, column, t, delta = 0.5) {
  obj <- fit$obj
  p <- obj$env$last.par.best
  if (is.null(p)) p <- obj$env$par
  xi <- which(names(p) == "dsem_x_tj")
  nr <- nrow(fit$dsem$tmb_inputs$data$y_tj)
  j  <- which(colnames(fit$dsem$tmb_inputs$data$y_tj) == column)
  testthat::expect_length(j, 1L)
  k <- xi[(j - 1L) * nr + t]
  p[k] <- p[k] + delta
  obj$report(p)
}

SEM3 <- paste0("recdevs1 <-> recdevs1, 0, sigmaR1, 1\n",
               "l1devs1 <-> l1devs1, 0, sigL1, 0.3\n",
               "linfdevs1 <-> linfdevs1, 0, sigLinf, 0.3\n")


testthat::test_that("a growth latent reaches only its own parameter slot", {
  testthat::skip_if_not_installed("dsem")
  fit <- .grow_fit(SEM3)
  r0  <- fit$obj$report(fit$obj$env$par)
  TS  <- 6L

  # The discriminating check: with L1 and Linf both driven, a wrong index shows
  # up as the WRONG PARAMETER moving, which a single-latent test cannot see.
  targets <- list(l1devs1 = 2L, linfdevs1 = 3L)         # L1 -> slot 1, Linf -> 2
  for (nm in names(targets)) {
    d_off <- .spike(fit, nm, TS)$growth_linkage_offset[1, 1, , ] -
      r0$growth_linkage_offset[1, 1, , ]
    moved <- targets[[nm]]
    # The spiked parameter moves by exactly the spike, in exactly that year.
    testthat::expect_equal(d_off[TS, moved], 0.5, tolerance = 1e-10)
    testthat::expect_equal(sum(abs(d_off[-TS, moved])), 0)
    # Every other parameter slot is untouched -- no leakage into K or m.
    testthat::expect_equal(sum(abs(d_off[, -moved])), 0)
  }
})


testthat::test_that("L1 moves the ramp age in its own year, then one age per year", {
  testthat::skip_if_not_installed("dsem")
  fit <- .grow_fit(SEM3)
  r0 <- fit$obj$report(fit$obj$env$par)
  TS <- 6L
  dL <- .spike(fit, "l1devs1", TS)$length_hat - r0$length_hat

  # Year t: the ramp age only (age index 1 = the minage cohort).
  testthat::expect_gt(abs(dL[1, 1, 1, TS]), 0)
  testthat::expect_equal(sum(abs(dL[1, 1, -1, TS])), 0)
  # Year t+1: that cohort, one age older, and nothing else.
  testthat::expect_gt(abs(dL[1, 1, 2, TS + 1L]), 0)
  testthat::expect_equal(sum(abs(dL[1, 1, -2, TS + 1L])), 0)
  # Nothing before the spiked year.
  testthat::expect_equal(sum(abs(dL[1, 1, , seq_len(TS - 1L)])), 0)
})


testthat::test_that("Linf moves the plus group in its own year, the curve after", {
  testthat::skip_if_not_installed("dsem")
  fit <- .grow_fit(SEM3)
  r0 <- fit$obj$report(fit$obj$env$par)
  TS <- 6L
  na <- fit$data_list$nages[1]
  dL <- .spike(fit, "linfdevs1", TS)$length_hat - r0$length_hat

  # Year t: the plus group (oldest age) only -- diff = linf - current_size.
  testthat::expect_gt(abs(dL[1, 1, na, TS]), 0)
  testthat::expect_equal(sum(abs(dL[1, 1, -na, TS])), 0)
  # Year t+1: the on-curve ages move through laglinf.
  testthat::expect_gt(sum(abs(dL[1, 1, 2:(na - 1L), TS + 1L])), 0)
  testthat::expect_equal(sum(abs(dL[1, 1, , seq_len(TS - 1L)])), 0)
})


testthat::test_that("successive year effects decay by exactly exp(-K)", {
  testthat::skip_if_not_installed("dsem")
  # This is what catches a one-year index shift. A shifted map still produces a
  # decaying series, but the FIRST responding year moves, so the ratio sequence
  # is checked against the parameter the recursion actually uses.
  fit <- .grow_fit(SEM3)
  r0 <- fit$obj$report(fit$obj$env$par)
  K  <- r0$growth_parameters[1, 1, 1, 1]
  TS <- 6L
  for (nm in c("l1devs1", "linfdevs1")) {
    dL <- apply(.spike(fit, nm, TS)$length_hat - r0$length_hat, 4,
                function(z) max(abs(z)))
    nz <- dL[dL > 0]
    testthat::expect_gt(length(nz), 5L)
    ratios <- unname(nz[-1] / nz[-length(nz)])
    # Steady state is the MIDDLE of the series, and both ends are excluded for
    # stated reasons rather than convenience. The first step superposes the
    # current-year channel (ramp or plus group) on the propagated one; the last
    # breaks when the perturbed cohort reaches the plus group, whose mean length
    # is an average over ages (measured 0.488 against exp(-K) = 0.741).
    testthat::expect_equal(ratios[2:4], rep(exp(-K), 3), tolerance = 1e-6)
  }
})


testthat::test_that("a recruitment-only DSEM leaves the growth offset at zero", {
  testthat::skip_if_not_installed("dsem")
  # The no-op guarantee for every existing DSEM fit: with no growth latent the
  # index map is empty, so section 5.5b-growth's loop never executes.
  fit <- .grow_fit("recdevs1 <-> recdevs1, 0, sigmaR1, 1\n")
  testthat::expect_length(fit$dsem$tmb_inputs$data$growth_dev_col, 0L)
  testthat::expect_equal(
    sum(abs(fit$obj$report(fit$obj$env$par)$growth_linkage_offset)), 0)
})


testthat::test_that("the growth latent is bias-corrected on the marginal variance", {
  testthat::skip_if_not_installed("dsem")

  # The spike tests above CANNOT see this: dsem_margvar_tj depends only on
  # beta_z, so it is unchanged by a spike to x_tj and cancels in the difference.
  # Review demonstrated that `expect_equal(d_off[TS, moved], 0.5)` passes
  # identically with the -margvar/2 term present, absent, or sign-flipped. So
  # this asserts the ABSOLUTE offset against the identity the template claims:
  #
  #   growth_linkage_offset = z_tj - bias_adjust_proc * margvar_tj / 2
  #
  # which is what makes E[exp(latent)] = 1, i.e. the PARAMETER's mean unbiased.
  # It does NOT make mean weight-at-age unbiased: weight is a convex power of
  # length, so Jensen leaves that biased upward.
  fit <- .grow_fit(SEM3)
  r   <- fit$obj$report(fit$obj$env$par)
  dt  <- fit$dsem$tmb_inputs$data
  nyd <- dt$nyrs_dsem
  ba  <- fit$data_list$bias_adjust_proc

  testthat::expect_gt(nyd, 0L)
  for (g in seq_along(dt$growth_dev_col)) {
    col <- dt$growth_dev_col[g] + 1L
    sp  <- dt$growth_dev_sp[g]  + 1L
    sx  <- dt$growth_dev_sex[g] + 1L
    par <- dt$growth_dev_par[g] + 1L
    got <- r$growth_linkage_offset[sp, sx, seq_len(nyd), par]
    want <- r$dsem_z_tj[seq_len(nyd), col] -
      ba * r$dsem_margvar_tj[seq_len(nyd), col] / 2
    testthat::expect_equal(as.numeric(got), as.numeric(want), tolerance = 1e-12)
  }

  # And the correction is actually DOING something here -- otherwise the
  # identity above would hold trivially and this test would prove nothing.
  mv <- r$dsem_margvar_tj[seq_len(nyd), dt$growth_dev_col + 1L, drop = FALSE]
  testthat::expect_gt(ba, 0)
  testthat::expect_gt(max(abs(mv)), 0)
})


# ---- process error: the DSEM field is drawn, selectively -------------------
# A DSEM's latent field is one joint GMRF, so sim_mod() draws it here rather
# than in the template. These pin that growth participates, that a caller can
# ask for one process without redrawing another (a cell left out is CONDITIONED
# on, not ignored), and that the drawn field comes back as the truth.

testthat::test_that("sim_mod() draws the growth field, and only what was asked for", {
  testthat::skip_if_not_installed("dsem")
  fit <- .grow_fit(SEM3)
  if (is.null(fit$obj$env$last.par.best))
    fit$obj$env$last.par.best <- fit$obj$env$par
  cols <- colnames(fit$dsem$tmb_inputs$data$y_tj)
  X0 <- as.matrix(fit$estimated_params$dsem_x_tj)
  moved <- function(st) {
    set.seed(3)
    d <- Rceattle:::.dsem_draw_process(fit, Rceattle:::.sim_state_codes(st))
    vapply(seq_along(cols), function(j) !isTRUE(all.equal(X0[, j], d[, j])),
           logical(1))
  }
  jR <- which(cols == "recdevs1"); jL <- which(cols == "l1devs1")

  all_ <- moved(TRUE)
  testthat::expect_true(all_[jR]); testthat::expect_true(all_[jL])
  # Asking for recruitment must NOT redraw growth -- before this it did, which
  # would corrupt a self-test's attribution without saying so.
  rec_ <- moved("recruitment")
  testthat::expect_true(rec_[jR]); testthat::expect_false(rec_[jL])
  # And asking for growth must not redraw recruitment.
  gro_ <- moved("growth")
  testthat::expect_false(gro_[jR]); testthat::expect_true(gro_[jL])
})


testthat::test_that("a growth-only sem can draw its process, and reports the truth", {
  testthat::skip_if_not_installed("dsem")
  # The gate used to require the SEM to drive RECRUITMENT, so a growth-only sem
  # could not draw its process at all; and the drawn field was attached to
  # attr(, "process_sim") only alongside rec_dev, so a growth-only draw lost its
  # truth silently and compare_sim() would have measured against the wrong one.
  fit <- .grow_fit("l1devs1 <-> l1devs1, 0, sigL1, 0.25\n")
  if (is.null(fit$obj$env$last.par.best))
    fit$obj$env$last.par.best <- fit$obj$env$par
  testthat::expect_false(Rceattle:::.dsem_drives(fit, "recruitment"))
  testthat::expect_true(Rceattle:::.dsem_drives_growth(fit))

  set.seed(8)
  ds <- suppressWarnings(Rceattle::sim_mod(fit, simulate = TRUE, process = TRUE))
  truth <- attr(ds, "process_sim")
  testthat::expect_false(is.null(truth$dsem_x_tj))
  # The field was actually drawn, not handed back at its fitted (zero) values.
  testthat::expect_gt(stats::sd(as.matrix(truth$dsem_x_tj)[, 1]), 0.05)
  testthat::expect_true(all(truth$dsem_x_tj_drawn))
})


testthat::test_that("a DSEM growth latent is not reported as having no distribution", {
  testthat::skip_if_not_installed("dsem")
  # .sim_warn_process_absent() tested only for a LINKAGE-GRAMMAR growth random
  # effect, so on a DSEM growth latent it said "no random linkage on a growth
  # parameter ... Nothing was redrawn and the fitted values stand" while the
  # field was in fact being drawn. That false alarm cost two wrong diagnoses.
  fit <- .grow_fit(SEM3)
  obj <- Rceattle:::.sim_obj(fit)
  w <- character(0)
  withCallingHandlers(
    Rceattle:::.sim_warn_process_absent(obj, Rceattle:::.sim_state_codes(TRUE)),
    warning = function(x) { w <<- c(w, conditionMessage(x))
                            invokeRestart("muffleWarning") })
  testthat::expect_false(any(grepl("growth", w)))
  # The M arm must still fire -- M1_re = 0 here, so M really cannot be drawn.
  testthat::expect_true(any(grepl("M1_re = 0", w)))
})


testthat::test_that("a growth latent is refused beside a RANDOM linkage on its parameter", {
  testthat::skip_if_not_installed("dsem")

  # Two random fields on one growth parameter, each scored by its own density.
  # A FIXED linkage composes, and a random one on another parameter is apart.
  e <- new.env(parent = asNamespace("Rceattle"))
  sys.source(testthat::test_path("helpers.R"), envir = e)
  d <- e$make_test_data(nyrs = 12L, nprojyrs = 3L, nages = 8L, minage = 1L,
                        growth = "vonBertalanffy", seed = 1)
  d$random_rec <- TRUE
  set.seed(4)
  d$env_data <- data.frame(Year = seq(d$styr, d$endyr), BT = stats::rnorm(12))
  fit_with <- function(linkages) suppressWarnings(Rceattle::fit_mod(
    data_list = d, estimateMode = 3, msmMode = 0,
    growthFun = Rceattle::build_growth(fun = "vonBertalanffy", linkages = linkages),
    dsem = Rceattle::build_DSEM(sem = paste0(
      "recdevs1 <-> recdevs1, 0, sigmaR1, 1\n",
      "linfdevs1 <-> linfdevs1, 0, sigLinf, 0.3\n")),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0)))

  testthat::expect_error(
    fit_with(list(Linf = Rceattle::linkage_spec(~ (1 | Year)))),
    "RANDOM Linf growth linkage")
  testthat::expect_error(
    fit_with(list(Linf = Rceattle::linkage_spec(~ (1 | Year),
                                                species = d$spnames[1]))),
    "RANDOM Linf growth linkage")
  for (lk in list(list(Linf = Rceattle::linkage_spec(~ BT)),
                  list(K = Rceattle::linkage_spec(~ (1 | Year))))) {
    testthat::expect_length(fit_with(lk)$dsem$tmb_inputs$data$growth_dev_col, 1L)
  }
})


testthat::test_that("pre-built objects are refused under a growth model that never reads them", {
  testthat::skip_if_not_installed("dsem")

  # Built on a data_list with no growth_model, then fitted under empirical
  # weight-at-age: the check has to run against the fitted settings. The msm
  # fixture fits no CAAL data (every Sample_size is 0), which empirical growth
  # would refuse first.
  e <- new.env(parent = asNamespace("Rceattle"))
  sys.source(testthat::test_path("helpers.R"), envir = e)
  sys.source(testthat::test_path("helpers-make-msm-data.R"), envir = e)
  set.seed(3)
  d <- e$make_msm_test_data()$data_list
  d$random_rec <- TRUE
  d$growth_model <- NULL
  d$env_data <- data.frame(Year = seq(d$styr, d$endyr),
                           BT = stats::rnorm(d$endyr - d$styr + 1L))
  obj <- build_dsem_objects(dsem_settings = Rceattle::build_DSEM(sem = paste0(
    "recdevs1 <-> recdevs1, 0, sR1, 1\nrecdevs2 <-> recdevs2, 0, sR2, 1\n",
    "linfdevs1 <-> linfdevs1, 0, sigLinf, 0.3\n")), data_list = d)
  testthat::expect_error(suppressWarnings(Rceattle::fit_mod(
    data_list = d, estimateMode = 3, msmMode = 0, dsem = obj,
    growthFun = Rceattle::build_growth(fun = "empirical"),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))),
    "empirical weight-at-age")
})
