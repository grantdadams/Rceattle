# Provenance: adversarial review of PR #117 (Tier 0 cleanup), defect 3.
#
# Reference-point recruitment is assigned by two arms in section 6.6. Option 1a
# required `proj_mean_rec == 1 & srr_pred_fun < 2`; Option 2 required
# `proj_mean_rec == 0`. A model with an estimated stock-recruit curve and
# proj_mean_rec left at the build_srr() default of TRUE matched NEITHER, so
# NByage0 and NByageF stayed 0 for every year after the first and SB0 was just
# the initial cohort decaying away.
#
# That is not a cosmetic gap: SB0(sp, nyrs - 1) is the depletion reference HCR 5
# and HCR 6 read, so a shrinking SB0 inflates perceived depletion and the
# resulting catch advice. Measured before the fix, on the fixture below:
#
#   NByage0[1,1,1,1:6]  1.212487 0 0 0 0 0
#   SB0[1,1:6]          3.344437 2.738193 2.241843 1.835466 1.502752 1.230350
#
# Option 2 now fires whenever a curve exists, not only when the projection runs
# on it. Reachable by an ordinary configuration -- an SRR with the default
# projection setting -- which is why it went unnoticed: build_srr()'s default
# srr_pred_fun is 0, so no bundled dataset or golden reference reaches it.

srr_default_proj_fit <- function() {
  d <- make_test_data(nyrs = 8, nages = 5, seed = 1)
  suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, estimateMode = 3, msmMode = 0,
    random_rec = FALSE,
    recFun = build_srr(srr_fun = "BevertonHolt", srr_pred_fun = "BevertonHolt",
                       srr_est_mode = 1, Bmsy_lim = 0),
    fit_control = fit_control(verbose = 0))))
}


test_that("an SRR sets reference-point recruitment even when the projection uses mean rec", {
  testthat::skip_on_cran()

  m <- srr_default_proj_fit()

  # The configuration that fell through both arms.
  expect_equal(as.integer(m$obj$env$data$proj_mean_rec), 1L)
  expect_gte(m$obj$env$data$srr_pred_fun, 2L)

  r <- m$obj$report(m$obj$env$last.par)

  # Unfished equilibrium recruitment is assigned for every year, not just year 1.
  expect_true(all(r$NByage0[1, 1, 1, ] > 0))
  expect_true(all(r$NByageF[1, 1, 1, ] > 0))

  # So SB0 is a flat equilibrium rather than a decaying cohort. It is the
  # terminal value that HCR 5 and 6 read as the depletion reference.
  sb0 <- r$SB0[1, ]
  expect_true(all(sb0 > 0))
  expect_equal(max(sb0) / min(sb0), 1, tolerance = 1e-8)
  expect_gt(sb0[length(sb0)], 0.9 * sb0[1])
})


test_that("the projection itself still runs on mean recruitment", {
  testthat::skip_on_cran()

  # The reference-point arm and the projection arm read proj_mean_rec
  # separately, so widening the first must not quietly switch the second onto
  # the curve. Projection recruitment is avg_R, flat across the forecast.
  m <- srr_default_proj_fit()
  r <- m$obj$report(m$obj$env$last.par)

  nyrs_hind <- length(m$data_list$styr:m$data_list$endyr)
  proj      <- r$R[1, seq(nyrs_hind + 1L, ncol(r$R))]

  expect_gt(length(proj), 1)
  expect_equal(max(proj) / min(proj), 1, tolerance = 1e-8)
  expect_equal(unname(proj[1]), unname(r$avg_R[1]), tolerance = 1e-8)
})
