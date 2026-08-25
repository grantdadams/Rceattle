# Over-dispersion of the Dirichlet-multinomial diet draw.
#
# Split out of test-functions-sim-mod-diet.R: that file made ~1200 sim_mod()
# calls and was 66% of the coverage run's wall clock. Same assertions, on their
# own runner. Shared fixture: helpers-diet-sim.R.
testthat::skip_on_cran()

testthat::test_that("the Dirichlet-multinomial diet draw is over-dispersed", {
  testthat::skip_if_not_installed("TMB")

  # A Dirichlet-multinomial at theta = 0.1 has sd sqrt((N + N*theta)/(1 + N*theta))
  # times the multinomial's -- about 2.7 at this sample size. Comparing the two
  # draws against each other needs no theoretical constant, and it is the check
  # that has power: at the concentration this model fits to (theta ~ 1e7) the
  # Dirichlet-multinomial is numerically a multinomial, so a test run there would
  # pass whatever the simulator did.
  mn <- .diet_sim_fit("Multinomial", 1)
  dm <- .diet_sim_fit("DirichletMultinomial", log(0.1))

  tgt <- .diet_pred(mn)
  big <- which(tgt > 0.05)
  testthat::expect_gt(length(big), 20)

  sd_mn <- apply(.diet_reps(mn, nrep = 150)[big, , drop = FALSE], 1, stats::sd)
  sd_dm <- apply(.diet_reps(dm, nrep = 150)[big, , drop = FALSE], 1, stats::sd)

  N <- as.numeric(Rceattle::rearrange_data(mn$data_list)$diet_obs[1, 1])
  expected <- sqrt((N + N * 0.1) / (1 + N * 0.1))
  testthat::expect_equal(mean(sd_dm / sd_mn), expected, tolerance = 0.15)
})

