# self_test() on a DSEM.
#
# WHAT A SELF TEST MEANS HERE, because it is easy to over-read: sim_mod() draws
# new OBSERVATIONS in R from the fitted quantities. It does not redraw the
# process random effects -- dsem.hpp's SIMULATE blocks were commented out when
# the header was vendored (they need `this`, and calculate_dsem() is a free
# function), and ceattle.cpp has none. So the fitted recruitment deviations are
# held fixed and the test asks whether the model recovers itself under fresh
# observation error. Under a DSEM those deviations come from the latent states,
# which are likewise held fixed, so the SEM's own process is never re-realized.
# That is the same limitation the standard recruitment path already has. A clean
# self test therefore says nothing about how well the SEM structure is
# identified across recruitment realizations.
#
# What it does establish, and what these tests pin: the simulated data really is
# different, the refit re-estimates the DSEM rather than echoing the parent back,
# and the replicates converge.

testthat::test_that("self_test() refits a DSEM on simulated data", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d <- Rceattle::BS2017SS
  yrs <- d$styr:d$endyr
  set.seed(1)
  d$env_data <- data.frame(
    Year = yrs, temp = as.numeric(scale(cumsum(stats::rnorm(length(yrs))))))
  sem <- "
    recdevs1 -> recdevs1, 1, rho_R,     0.3
    temp     -> recdevs1, 0, temp_to_R, 0.2
    recdevs1 <-> recdevs1, 0, sigmaR1,  0.6
    recdevs2 <-> recdevs2, 0, sigmaR2,  0.6
    recdevs3 <-> recdevs3, 0, sigmaR3,  0.6
  "
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = TRUE, msmMode = 0,
    dsem = Rceattle::build_DSEM(sem = sem, family = "fixed"),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))

  # The simulated data must actually differ, or the "self test" refits the very
  # data set it started from and passes for free.
  sim <- suppressWarnings(suppressMessages(
    Rceattle::sim_mod(fit, simulate = TRUE)))
  testthat::expect_false(isTRUE(all.equal(sim$index_data$Observation,
                                          d$index_data$Observation)))

  # self_test() returns the model list DIRECTLY, named Sim_1..Sim_k -- not
  # wrapped in $Rceattle_list the way retrospective() and jitter() are.
  st <- suppressWarnings(suppressMessages(
    Rceattle::self_test(fit, nsim = 2, seed = 11, cores = 1, getsd = FALSE)))
  testthat::expect_gt(length(st), 0L)

  bp <- fit$estimated_params$dsem_beta_z
  for (nm in names(st)) {
    b <- st[[nm]]$estimated_params$dsem_beta_z
    testthat::expect_equal(length(b), length(bp))
    # Re-estimated, not carried over: a refit that silently reused the parent's
    # DSEM would return these unchanged.
    testthat::expect_gt(max(abs(b - bp)), 1e-6)
    testthat::expect_true(is.finite(st[[nm]]$opt$objective))
  }
})

# `process` accepts FALSE/TRUE, "none"/"all"/"dynamics"/"observation", and any
# subset of the five process names. An earlier guard tested isTRUE(process),
# which is FALSE for every character spelling -- so process = "recruitment" and
# process = "all" walked past it and returned a data set whose recruitment was
# never redrawn, silently. Under "all" it was worse than silent: the warning
# enumerated the OTHER un-drawn processes and so read as confirmation that
# recruitment HAD been drawn.

testthat::test_that("redrawing recruitment on a DSEM is refused by every spelling", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d <- Rceattle::BS2017SS
  d$env_data <- data.frame(Year = d$styr:d$endyr, BT = 0)
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 3,
    random_rec = TRUE, msmMode = 0, dsem = Rceattle::build_DSEM(),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))

  for (p in list(TRUE, "all", "dynamics", "recruitment",
                 c("recruitment", "M"))) {
    testthat::expect_error(
      suppressWarnings(Rceattle::sim_mod(fit, simulate = TRUE, process = p)),
      "not supported on a DSEM fit",
      info = paste(deparse(p), collapse = " "))
  }

  # ...and the spellings that do NOT ask for recruitment must still work.
  for (p in list(FALSE, "none", "observation", "M")) {
    testthat::expect_error(
      suppressWarnings(suppressMessages(
        Rceattle::sim_mod(fit, simulate = TRUE, process = p))),
      NA, info = paste(deparse(p), collapse = " "))
  }
})

testthat::test_that("a DSEM fit warns that R0 is not comparable at the default bias correction", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  # Lognormal bias correction is not applied to DSEM deviations. SSB absorbs the
  # offset and looks fine; R0 does not -- 20-51% below the non-DSEM fit on
  # BS2017SS with an IID sem -- and dynamic B0 and the B40% proxy are keyed to
  # R0. The number must not go out unmarked.
  d <- Rceattle::BS2017SS
  d$env_data <- data.frame(Year = d$styr:d$endyr, BT = 0)
  mk <- function(bias) Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                             verbose = 1, bias_adjust_proc = bias)
  testthat::expect_warning(
    suppressMessages(Rceattle::fit_mod(
      data_list = d, inits = NULL, file = NULL, estimateMode = 3,
      random_rec = TRUE, msmMode = 0, dsem = Rceattle::build_DSEM(),
      fit_control = mk(TRUE))),
    "not comparable to a non-DSEM fit")

  # Silent where the two paths ARE comparable.
  testthat::expect_warning(
    suppressMessages(Rceattle::fit_mod(
      data_list = d, inits = NULL, file = NULL, estimateMode = 3,
      random_rec = TRUE, msmMode = 0, dsem = Rceattle::build_DSEM(),
      fit_control = mk(FALSE))),
    NA)
})

# The diagnostics that refit through `inits` were all silently broken by the
# same defect: fit_mod() dropped the dsem_* blocks out of a warm start, so each
# rebuilt the model with the recruitment SD at its start value (0.707107). They
# work now; these pin that the DSEM actually survives the round trip rather than
# merely that the call returns.

testthat::test_that("remove_F(), reweight_comps() and osa_residuals() keep the DSEM", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d <- Rceattle::BS2017SS
  d$env_data <- data.frame(Year = d$styr:d$endyr, BT = 0)
  fc <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0,
                              bias_adjust_proc = FALSE)
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = TRUE, msmMode = 0, dsem = Rceattle::build_DSEM(),
    fit_control = fc)))
  bp <- as.numeric(fit$estimated_params$dsem_beta_z)
  # The value the bug produced, so the assertions below name what they exclude.
  START <- 0.707107

  rf <- suppressWarnings(suppressMessages(Rceattle::remove_F(fit)))
  # remove_F() rebuilds at estimateMode = 3 and re-estimates nothing, so the
  # DSEM must come back EXACTLY, not merely close.
  testthat::expect_equal(as.numeric(rf$estimated_params$dsem_beta_z), bp,
                         tolerance = 1e-10)

  rw <- suppressWarnings(suppressMessages(
    Rceattle::reweight_comps(fit, n_iter = 2, verbose = FALSE)))
  br <- as.numeric(rw$estimated_params$dsem_beta_z)
  testthat::expect_equal(length(br), length(bp))
  # Re-estimated against the new weights, so not equal to the parent -- but it
  # must not have collapsed onto the start value for every species.
  testthat::expect_gt(max(abs(br - START)), 1e-3)
  testthat::expect_gt(stats::sd(br), 1e-6)

  # osa_residuals(): one-step-ahead residuals for the DATA. The random effects
  # are integrated out whatever their structure, so a DSEM must give the same
  # index residuals as the non-DSEM fit of the same model.
  fp <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = TRUE, msmMode = 0, fit_control = fc)))
  od <- suppressWarnings(suppressMessages(
    Rceattle::osa_residuals(fit, source = "index", parallel = FALSE)))
  op <- suppressWarnings(suppressMessages(
    Rceattle::osa_residuals(fp, source = "index", parallel = FALSE)))
  testthat::expect_equal(nrow(od), nrow(op))
  testthat::expect_true(all(is.finite(od$residual)))
  testthat::expect_equal(od$residual, op$residual, tolerance = 1e-3)
})
