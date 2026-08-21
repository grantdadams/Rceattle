# Diet_comp_weights and the multinomial diet draw.
#
# Split out of test-functions-sim-mod-diet.R: that file made ~1200 sim_mod()
# calls and was 66% of the coverage run's wall clock. Same assertions, on their
# own runner. Shared fixture: helpers-diet-sim.R.
testthat::skip_on_cran()

testthat::test_that("the multinomial diet draw honours Diet_comp_weights", {
  testthat::skip_if_not_installed("TMB")

  # For diet_ll_type == 0 the weight multiplies the log-likelihood, which makes
  # it an effective sample size. Drawing at the nominal N regardless would give
  # a down-weighted predator stomachs the likelihood treats as w times less
  # informative than they are, and self_test() would report recovery against
  # data of the wrong precision. (Note the weight lives in the PARAMETER
  # diet_comp_weights, not the data column, which fit_mod() ignores when inits
  # are supplied.)
  w <- 4
  mn1 <- .diet_sim_fit("Multinomial", 1)
  mnw <- .diet_sim_fit("Multinomial", w)
  testthat::expect_equal(as.numeric(mnw$estimated_params$diet_comp_weights)[1], w)

  tgt <- .diet_pred(mn1)
  big <- which(tgt > 0.02)
  testthat::expect_gt(length(big), 20)

  sd1 <- apply(.diet_reps(mn1, nrep = 150)[big, , drop = FALSE], 1, stats::sd)
  sdw <- apply(.diet_reps(mnw, nrep = 150)[big, , drop = FALSE], 1, stats::sd)

  # Weighting by w shrinks the spread by sqrt(w).
  testthat::expect_equal(mean(sdw / sd1), 1 / sqrt(w), tolerance = 0.12)
})

