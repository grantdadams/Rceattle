# Regression test for the `Time_varying_q` overload.
#
# `Time_varying_q` holds a tv_q_map mode code for most fleets, but when
# `Catchability` is "AR1" or "Environmental" it instead holds 1-based
# `env_data` column indices, and convert_switches() deliberately leaves
# those rows unconverted. "Environmental" (Catchability = 5) may carry a
# comma-separated list ("1,3"), which as.integer() turned into NA with a
# coercion warning. TMB reads the resulting `index_varying_q` vector at
# ceattle_v01_11.cpp:3042/3049/3058 and process_residuals() tests it
# against 4, so the NA silently poisoned a branch test.
#
# Environmental fleets must now emit an explicit 0 ("not applicable") --
# their env indices reach the model via the index_q_beta map instead.

.rearranged_index_varying_q <- function(catchability, time_varying_q) {
  dat <- make_test_data(seed = 1)
  dat$fleet_control$Catchability   <- catchability
  dat$fleet_control$Time_varying_q <- time_varying_q
  dat <- Rceattle::switch_check(dat)
  Rceattle::rearrange_data(dat)$index_varying_q
}


testthat::test_that("a comma-separated Environmental q list does not become NA", {
  n_flt <- nrow(make_test_data(seed = 1)$fleet_control)
  testthat::skip_if(n_flt < 1L)

  ivq <- testthat::expect_no_warning(
    .rearranged_index_varying_q(
      catchability   = rep("Environmental", n_flt),
      time_varying_q = rep("1,3", n_flt)
    )
  )

  testthat::expect_false(anyNA(ivq))
  testthat::expect_type(ivq, "integer")
  # 0 = "not applicable"; must not collide with the mode codes the cpp
  # tests for (1 IID, 2 AR1, 4 RandomWalk).
  testthat::expect_true(all(ivq == 0L))
})


testthat::test_that("non-environmental fleets keep their tv_q mode code", {
  n_flt <- nrow(make_test_data(seed = 1)$fleet_control)
  testthat::skip_if(n_flt < 1L)

  ivq <- .rearranged_index_varying_q(
    catchability   = rep("Estimated", n_flt),
    time_varying_q = rep("RandomWalk", n_flt)
  )

  testthat::expect_false(anyNA(ivq))
  # tv_q_map: RandomWalk = 4
  testthat::expect_true(all(ivq == 4L))
})
