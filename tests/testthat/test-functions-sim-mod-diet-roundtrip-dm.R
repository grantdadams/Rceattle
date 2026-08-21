# Round trip for the Dirichlet-multinomial diet draw.
#
# Split out of test-functions-sim-mod-diet.R: that file made ~1200 sim_mod()
# calls and was 66% of the coverage run's wall clock. Same assertions, on their
# own runner. Shared fixture: helpers-diet-sim.R.
testthat::skip_on_cran()

testthat::test_that("what a Dirichlet-multinomial refit sees is centred on the density", {
  testthat::skip_if_not_installed("TMB")
  .diet_roundtrip_check("DirichletMultinomial")
})
