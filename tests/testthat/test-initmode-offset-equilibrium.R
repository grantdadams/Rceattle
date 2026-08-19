# =============================================================================
# initMode = "OffsetEquilibrium" (5): an F = 0 equilibrium age-structure seeded
# by the FIRST-YEAR recruitment exp(rec_pars + rec_dev[year 1]) decayed by M1,
# with initial deviates turned off and no init-dev penalty. This is the Cole
# Monnahan / AFSC GOA pollock convention -- the initial numbers track the
# year-1 recruitment deviation under free estimation, unlike initMode = 1
# (which seeds off the mean-recruitment equilibrium R0). The name refers to
# that offset: the equilibrium is unfished, but displaced from mean recruitment
# by the year-1 deviation.
#
# NOTE: initMode is a fit_mod() ARGUMENT; it overwrites data_list$initMode, so
# these tests pass it as the argument (setting data_list$initMode has no effect).
# =============================================================================

.oe_fit <- function(mode, inits = NULL) {
  data(BS2017SS, package = "Rceattle")
  suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = BS2017SS, inits = inits, initMode = mode, estimateMode = 3,
    msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, phase = FALSE, verbose = 0))))
}

testthat::test_that("OffsetEquilibrium is a recognised initMode string alias (= 5)", {
  testthat::expect_true("OffsetEquilibrium" %in% names(Rceattle:::initMode_map))
  testthat::expect_identical(unname(Rceattle:::initMode_map["OffsetEquilibrium"]), 5)
  # The pre-release spelling was renamed and is no longer accepted.
  testthat::expect_false("FishedEquilibrium" %in% names(Rceattle:::initMode_map))
})

testthat::test_that("OffsetEquilibrium initN equals R1 * exp(-cumsum M1) closed form", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  data(BS2017SS, package = "Rceattle")
  oe <- .oe_fit("OffsetEquilibrium")
  q <- oe$quantities
  sp <- 1L
  nage <- BS2017SS$nages[sp]
  M1 <- q$M1_at_age[sp, 1, 1:nage, 1]
  R1 <- q$R[sp, 1]                        # first-year recruitment = exp(rec_pars + rec_dev[1])
  N1 <- q$N_at_age[sp, 1, 1:nage, 1]
  cum <- c(0, head(cumsum(M1), -1))
  expect <- R1 * exp(-cum)                # 1-sex reference: sex_ratio = 1
  expect[nage] <- expect[nage] / (1 - exp(-M1[nage]))   # plus-group geometric closure
  testthat::expect_equal(as.numeric(N1), as.numeric(expect), tolerance = 1e-10)
})

testthat::test_that("OffsetEquilibrium carries no init_dev penalty (NonEquilibrium does)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  oe <- .oe_fit("OffsetEquilibrium")
  ne <- .oe_fit("NonEquilibrium")
  testthat::expect_equal(sum(oe$quantities$jnll_comp["Initial abundance deviates", ]), 0)
  testthat::expect_gt(sum(ne$quantities$jnll_comp["Initial abundance deviates", ]), 0)
  # init_dev is fully mapped out (fixed at 0) under OffsetEquilibrium.
  mp <- Rceattle:::build_map(oe$data_list, oe$estimated_params)$mapList
  testthat::expect_true(all(is.na(mp$init_dev)))
})

testthat::test_that("OffsetEquilibrium initial ages scale with the year-1 recruitment dev", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  data(BS2017SS, package = "Rceattle")
  oe <- .oe_fit("OffsetEquilibrium")
  ne <- .oe_fit("NonEquilibrium")
  # Bump the year-1 recruitment deviation by +0.5.
  p2 <- oe$estimated_params
  p2$rec_dev[1, 1] <- p2$rec_dev[1, 1] + 0.5
  oe_b <- .oe_fit("OffsetEquilibrium", inits = p2)
  ne_b <- .oe_fit("NonEquilibrium",    inits = p2)
  sp <- 1L; nage <- BS2017SS$nages[sp]; iage <- 2:nage
  rat_oe <- oe_b$quantities$N_at_age[sp, 1, iage, 1] / oe$quantities$N_at_age[sp, 1, iage, 1]
  rat_ne <- ne_b$quantities$N_at_age[sp, 1, iage, 1] / ne$quantities$N_at_age[sp, 1, iage, 1]
  # OffsetEquilibrium seeds off first-year recruitment -> all initial ages scale
  # by exp(0.5); NonEquilibrium seeds off R0 -> initial ages are unchanged.
  testthat::expect_equal(as.numeric(rat_oe), rep(exp(0.5), length(iage)), tolerance = 1e-10)
  testthat::expect_equal(as.numeric(rat_ne), rep(1, length(iage)), tolerance = 1e-10)
})
