# Age/length compositions and conditional age-at-length are drawn by the TMB
# model, in RAW bin space (ceattle.cpp, slots 2 and 3). This file holds the two
# cheap contract checks; the fold, weighting and CAAL claims are in
# test-functions-sim-mod-comp-{fold,weights,caal}.R, split so the shard planner
# can put them on different runners. Shared fixtures: helpers-comp-sim.R.
testthat::skip_on_cran()

testthat::test_that("an empty composition row is reported, not passed through", {
  testthat::skip_if_not_installed("TMB")

  # A row sums to zero when the model predicts no composition for it -- a fleet
  # with no catch that year, or one switched off -- or when Sample_size times the
  # composition weight rounds away. run_mse() reads a zero row as "not sampled"
  # and drops the sample size with it, so the row must come back zero rather than
  # holding its original proportions (which sum to 1 and look entirely
  # plausible). Say so rather than let a self_test() score a replicate whose
  # data quietly went missing.
  testthat::expect_warning(
    .sim_warn_empty_comp(matrix(c(0, 0, 1, 2), nrow = 2, byrow = TRUE),
                         c(1, 2), "composition"),
    "came back empty")
  testthat::expect_no_warning(
    .sim_warn_empty_comp(matrix(c(3, 1, 1, 2), nrow = 2, byrow = TRUE),
                         c(1, 2), "composition"))
})


testthat::test_that("simulate = FALSE returns the predicted composition and draws nothing", {
  testthat::skip_if_not_installed("TMB")

  fit  <- .comp_fixture()
  cols <- grep("^Comp_", names(fit$data_list$comp_data))

  set.seed(2); before <- .Random.seed
  sim <- suppressWarnings(Rceattle::sim_mod(fit, simulate = FALSE))
  testthat::expect_identical(.Random.seed, before)
  testthat::expect_equal(unname(as.matrix(sim$comp_data[, cols])),
                         unname(fit$quantities$comp_hat[, seq_along(cols)]),
                         tolerance = 1e-10)
})
