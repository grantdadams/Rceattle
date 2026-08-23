# Provenance: inst/dev/CLEANUP_BACKLOG.md Tier 0, the ceattle.cpp `caal_ll_type`
# entry. CAAL_distribution = "MultinomialAFSC" is in comp_loglike_map, which
# CAAL_distribution shares with Comp_distribution, so it passed
# validate_switches() and then died inside the template with
# "Invalid 'caal_ll_type'" -- the CAAL dispatch had only cases 0 and 1. The same
# shape as Estimate_catch_sd = "Analytical", fixed in 5.12.0 by implementing the
# missing case; this is that fix for CAAL.
#
# The form is the AFSC/AMAK multinomial pseudo-likelihood already implemented for
# age comps (comp_ll_type == -1), applied to the CAAL row's proportions with the
# CAAL sample size and weight:
#   -sum_a  w * N * (obs_a + off) * log((hat_a + off) / (obs_a + off))
# It drops the multinomial normalizing constant, so it is a pseudo-likelihood
# rather than a density -- the OSA branch residualizes it under the full
# multinomial and the simulation draw is multinomial, exactly as for comps.

caal_fit <- function(family) {
  d <- make_test_data(nyrs = 8, nages = 5, seed = 1, growth = "vonBertalanffy")
  d$fleet_control$CAAL_distribution <- family
  suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, estimateMode = 3, msmMode = 0,
    random_rec = FALSE, growthFun = build_growth(fun = "vonBertalanffy"),
    fit_control = fit_control(verbose = 0))))
}

jnll_caal_row <- function() {
  dir <- c("src/TMB", testthat::test_path("..", "..", "src", "TMB"))
  dir <- dir[dir.exists(dir)]
  testthat::skip_if(length(dir) == 0, "src/TMB not available")
  src <- readLines(file.path(dir[1], "ceattle.cpp"), warn = FALSE)
  hit <- grep("JNLL_CAAL[ ]*=[ ]*[0-9]+", src, value = TRUE)
  testthat::expect_length(hit, 1)
  as.integer(sub(".*JNLL_CAAL[ ]*=[ ]*([0-9]+).*", "\\1", hit)) + 1L
}


test_that("CAAL_distribution = 'MultinomialAFSC' fits instead of erroring", {
  testthat::skip_on_cran()

  m <- caal_fit("MultinomialAFSC")
  expect_true(is.finite(m$obj$fn(m$obj$par)))

  # A non-zero CAAL row: the branch has to be reached, not merely not crash.
  expect_gt(sum(m$obj$report(m$obj$env$last.par)$jnll_comp[jnll_caal_row(), ]), 0)
})


test_that("the CAAL AFSC value is the AFSC form, computed from its own inputs", {
  testthat::skip_on_cran()

  m   <- caal_fit("MultinomialAFSC")
  dat <- m$obj$env$data
  r   <- m$obj$report(m$obj$env$last.par)
  w   <- m$obj$env$parList()$caal_weights   # PARAMETER_VECTOR, not data
  off <- dat$comp_offset

  hand <- rep(0, ncol(r$jnll_comp))
  for (i in seq_len(nrow(dat$caal_ctl))) {
    flt <- dat$caal_ctl[i, 1]; sp <- dat$caal_ctl[i, 2]; yr <- dat$caal_ctl[i, 4]
    if (!(yr <= dat$endyr && yr > 0 && dat$caal_n[i, 1] > 0)) next
    o <- r$caal_obs[i, seq_len(dat$nages[sp])]
    h <- r$caal_hat[i, seq_len(dat$nages[sp])]
    hand[flt] <- hand[flt] -
      sum(w[flt] * dat$caal_n[i, 1] * (o + off) * log((h + off) / (o + off)))
  }

  expect_equal(unname(r$jnll_comp[jnll_caal_row(), ]), hand, tolerance = 1e-10)
})


test_that("the AFSC CAAL family is a different likelihood from the multinomial", {
  testthat::skip_on_cran()

  # Otherwise the new case could be silently falling through to case 0 and this
  # file would pass while testing nothing.
  afsc <- caal_fit("MultinomialAFSC")
  mult <- caal_fit("Multinomial")
  expect_false(isTRUE(all.equal(afsc$obj$fn(afsc$obj$par),
                                mult$obj$fn(mult$obj$par))))

  # And the multinomial itself is untouched by the new branch existing.
  expect_true(is.finite(mult$obj$fn(mult$obj$par)))
})
