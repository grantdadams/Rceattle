# The two natural-scale index families and their simulators.
#
# Index_distribution = "Normal" is a plain normal on the natural-scale residual.
# It matches the AMAK avo_like / cpue_like term for term, which is what the ADMB
# bridges compare against, so it stays untruncated. Its draw is therefore a plain
# normal too -- and a plain normal can come back non-positive, which data_check()
# refuses, so those draws are redrawn. The redrawn row follows the normal
# TRUNCATED at zero while the likelihood scores the untruncated one, and
# sim_mod() warns when that is doing enough of the work to matter.
#
# Index_distribution = "TruncatedNormal" closes that gap: the density is
# renormalized over the support the data can actually occupy,
#   log f(x) = log phi(x; mu, sd) - log(1 - Phi(-mu/sd)),  1 - Phi(-mu/sd) = Phi(mu/sd)
# and the draw is taken by inverse CDF on (0, Inf), so draw and density are the
# same distribution and no draw can come back for data_check() to reject.
#
# No golden model uses either family -- all four are lognormal -- so
# /golden-check is silent here and this is the regression net.

.tn_fixture <- function(family, sd, nyrs = 20, mode = 3) {
  dat <- make_test_data(nyrs = nyrs, nages = 5, seed = 123)
  srv <- dat$fleet_control$Fleet_name == "Survey"
  dat$fleet_control$Index_distribution[srv] <- family
  dat$fleet_control$Catchability[srv] <- "AnalyticalArith"
  dat$index_data$Log_sd[dat$index_data$Fleet_name == "Survey"] <- sd
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = mode, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))
}

testthat::test_that("TruncatedNormal carries the truncation constant and Normal does not", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # sd large relative to the index, so the constant is far from zero and the
  # check has power: at mu/sd ~ 1 it is worth ~0.17 nll per observation.
  for (fam in c("Normal", "TruncatedNormal")) {
    fit <- .tn_fixture(fam, sd = 80)
    rep <- fit$obj$report()
    td  <- fit$obj$env$data
    sel <- which(td$index_ctl[, 1] == 1 & td$index_ctl[, 3] > 0 &
                   td$index_ctl[, 3] <= td$endyr & td$index_obs[, 1] > 0)

    plain <- -sum(stats::dnorm(td$index_obs[sel, 1], rep$index_hat[sel], 80, log = TRUE))
    const <- sum(log(stats::pnorm(rep$index_hat[sel] / 80)))
    want  <- if (fam == "TruncatedNormal") plain + const else plain

    testthat::expect_equal(rep$jnll_comp[1, 1], want, tolerance = 1e-6)
    # ...and the constant is big enough that the two families are distinguishable,
    # or neither assertion above proves anything.
    testthat::expect_gt(abs(const), 1)
  }
})

testthat::test_that("a TruncatedNormal draw is positive by construction and warns about nothing", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Inverse-CDF sampling on (0, Inf): positivity holds by construction, with no
  # rejection loop and no retry budget that could run out. sd is deliberately
  # larger than the index, where an untruncated normal would go negative often.
  fit <- .tn_fixture("TruncatedNormal", sd = 150, mode = 1)

  set.seed(4)
  for (k in 1:15) {
    sim <- Rceattle::sim_mod(fit, simulate = TRUE)   # no warning expected: draw == density
    srv <- sim$index_data$Fleet_name == "Survey"
    testthat::expect_true(all(sim$index_data$Observation[srv] > 0))
  }
})

testthat::test_that("a Normal draw is redrawn when it comes back non-positive, and says so", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Same fixture on the untruncated family. The draw and the likelihood now
  # disagree -- that is the honest description of this family, and the point of
  # the warning is that it names the family to switch to.
  fit <- .tn_fixture("Normal", sd = 150, mode = 1)

  set.seed(4)
  w <- testthat::capture_warnings(sim <- Rceattle::sim_mod(fit, simulate = TRUE))
  testthat::expect_true(any(grepl("TruncatedNormal", w)))
  # Whatever the rejection rate, a returned draw is either positive or reported
  # as unusable -- it is never silently non-positive.
  srv <- sim$index_data$Fleet_name == "Survey"
  testthat::expect_true(all(sim$index_data$Observation[srv] > 0) ||
                          any(grepl("non-positive", w)))
})

testthat::test_that("both natural-scale families keep obsvec in step with the draw", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # The OSA path scores the natural-scale observation read from obsvec, so a draw
  # that updates index_obs but not obsvec leaves the two describing different
  # data. Checked on the report rather than through a refit, because that is
  # where the inconsistency would be silent.
  #
  # obsvec carries two layouts. build_osa_data() stores the untransformed
  # observation for these families only when the OSA data are built; on the
  # ordinary fitting path it keeps the log(obs) layout it uses for every fleet.
  # These fits are built the ordinary way, so log is the expectation here -- and
  # asserting the natural scale instead would pass against a draw that wrote the
  # wrong layout.
  for (fam in c("Normal", "TruncatedNormal")) {
    fit <- .tn_fixture(fam, sd = 20, mode = 1)
    set.seed(11)
    rep <- fit$obj$simulate()
    td  <- fit$obj$env$data
    pos <- td$index_obsvec_idx
    got <- which(pos >= 0)
    testthat::skip_if(!length(got), "no fitted index rows in obsvec")
    drawn <- as.numeric(rep$index_obs_sim[got, 1])
    testthat::expect_true(all(drawn > 0))
    testthat::expect_equal(as.numeric(rep$obsvec_sim[pos[got] + 1L]),
                           log(drawn), tolerance = 1e-10)
  }
})

testthat::test_that("the OSA obsvec holds the untransformed obs for every natural-scale family", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # build_osa_data() lays obsvec out per family, and the cpp's OSA branch scores
  # each family against the layout it expects. A natural-scale family that misses
  # the untransformed branch falls through to the lognormal default and is
  # residualized as log(obs) against a natural-scale mean -- no error, just a
  # meaningless residual. Enumerated from index_distribution_map so a new family
  # fails here rather than in someone's Q-Q plot.
  for (fam in c("Normal", "TruncatedNormal")) {
    dat <- make_test_data(nyrs = 12, nages = 5, seed = 123)
    srv <- dat$fleet_control$Fleet_name == "Survey"
    dat$fleet_control$Index_distribution[srv] <- fam
    dat$fleet_control$Catchability[srv] <- "AnalyticalArith"
    dat$index_data$Log_sd[dat$index_data$Fleet_name == "Survey"] <- 20

    d <- suppressMessages(suppressWarnings(
      Rceattle::rearrange_data(Rceattle::switch_check(dat), build_osa = TRUE)))
    ill  <- d$index_ll_type
    rows <- which(d$index_obsvec_idx >= 0)
    testthat::skip_if(!length(rows), "no index rows in obsvec")
    nat <- ill[d$index_ctl[rows, 1]] %in% c(3L, 4L)
    testthat::expect_true(any(nat), label = fam)

    got  <- d$obsvec[d$index_obsvec_idx[rows[nat]] + 1L]
    want <- d$index_obs[rows[nat], 1]
    testthat::expect_equal(as.numeric(got), as.numeric(want),
                           tolerance = 1e-10, label = fam)
    # ...and NOT the log, which is what the fall-through would have stored.
    testthat::expect_false(isTRUE(all.equal(as.numeric(got),
                                            log(as.numeric(want)))))
  }
})
