# In single-species mode a species with input numbers-at-age (estDynamics > 0)
# has rec_pars fixed, so its equilibrium SB0 and B0 are built on placeholder
# recruitment (R0 = exp(9)), and with DynamicHCR = FALSE ssb_depletion and
# biomass_depletion divide by them. Those four are NA for such a species.
# Under DynamicHCR = TRUE the depletions are DynamicSB0-based, which for input
# numbers is the numbers relative to themselves, and stay reported. Under
# predation MSSB0 replaces SB0 (test-functions-fixed-n-hcr-outputs.R).
# A model whose only species is fixed has nothing to estimate, so the fixture is
# the two-species one run without predation.

testthat::skip_on_cran()

testthat::test_that("placeholder SB0, B0 and depletion are NA in single-species mode", {
  m <- fixed_natage_build(fixed_natage_data(c(1, 0)), msmMode = 0)
  q <- m$quantities
  for (nm in c("SB0", "B0", "ssb_depletion", "biomass_depletion", "R", "R0")) {
    testthat::expect_true(all(is.na(q[[nm]][1, ])), info = nm)
    testthat::expect_true(all(is.finite(q[[nm]][2, ])), info = nm)
  }
  for (nm in c("R_init", "steepness", "SPR0")) testthat::expect_true(is.na(q[[nm]][1]), info = nm)
  testthat::expect_true(all(is.finite(q$ssb[1, ])))
  testthat::expect_true(all(is.finite(q$DynamicSB0[1, ])))

  # sample_rec() reads R and R0 for every species it projects; the fixed one is skipped.
  s <- sample_rec(m, sample_rec = FALSE, update_model = FALSE)
  testthat::expect_true(all(is.finite(s$estimated_params$rec_dev)))
})

testthat::test_that("estDynamics = 2 in single-species mode is announced as fitting as 1", {
  d <- fixed_natage_data(c(2, 0))
  d$msmMode <- 0
  d <- suppressMessages(Rceattle::switch_check(d))
  testthat::expect_message(suppressWarnings(Rceattle:::data_check(d)),
                           "estimated only under predation")
})

testthat::test_that("under DynamicHCR = TRUE the depletions stay reported", {
  m <- fixed_natage_build(fixed_natage_data(c(1, 0)), msmMode = 0,
                          HCR = build_hcr(HCR = "ConstantF", DynamicHCR = TRUE, Ftarget = 0.1))
  q <- m$quantities
  testthat::expect_true(all(is.na(q$SB0[1, ])))
  testthat::expect_true(all(is.finite(q$ssb_depletion[1, ])))
})
