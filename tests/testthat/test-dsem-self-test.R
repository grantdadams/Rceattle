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
