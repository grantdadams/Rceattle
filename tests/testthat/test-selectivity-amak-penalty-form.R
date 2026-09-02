# The ADMB atka mackerel penalty forms on NonParametric (type 2) selectivity.
#
# Provenance. `Rceattle-models/BSAI atka mackerel` bridged to ADMB through a
# vendored fork of an old ceattle.cpp (`src/ceattle_v01_11_amak.cpp`). That fork
# writes two of the type-2 penalties differently from the package: the random
# walk is a bare Gaussian sum of squares, `0.5 * square(diff / sd)` with its
# `dnorm` line commented out, and the level term is `20 * square(avg_sel)`
# rather than a fixed 2. `Sel_penalty_form = "AMAK"` selects those forms, with
# the level weight read from `Sel_avgsel_pen`, so the bridge can run on the
# package instead of on a fork that no longer builds against the current data
# contract.
#
# The two differ by a constant only while `Time_varying_sel_sd` is fixed; once
# the sd is estimated they are different likelihoods, which is why this is a
# switch and not a reporting detail.
testthat::skip_on_cran()

testthat::test_that("the default path is unchanged, against a pinned objective", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOAcod")
  f <- function(d) suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, inits = NULL, msmMode = 0, estimateMode = "Hindcast",
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))$opt$objective

  # A literal, so this catches a rewrite of the default branch. Comparing the
  # column's absence against an explicit "Rceattle" on the same build would not:
  # both would move together. `test-golden-regression.R` carries the four
  # reference models; this is a fifth pin on the non-parametric path specifically.
  testthat::expect_equal(f(GOAcod), 11742.2844472541, tolerance = 1e-10)

  # Absent and explicitly "Rceattle" are the same setting.
  d <- GOAcod
  d$fleet_control$Sel_penalty_form <- "Rceattle"
  testthat::expect_equal(f(d), f(GOAcod))
})

testthat::test_that("the AMAK walk is the normalized walk minus its constant, exactly", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("Atka2022")
  mk <- function(form) {
    d <- Atka2022
    d$fleet_control$Time_varying_sel <- as.character(d$fleet_control$Time_varying_sel)
    d$fleet_control$Time_varying_sel[2]    <- "RandomWalk"
    d$fleet_control$Time_varying_sel_sd[2] <- 0.35
    d$fleet_control$Sel_penalty_form       <- form
    d$fleet_control$Sel_avgsel_pen         <- 20
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = d, inits = NULL, msmMode = 0, estimateMode = 3,
      fit_control = Rceattle::fit_control(verbose = 0))))
  }
  a <- mk("Rceattle"); b <- mk("AMAK")
  ja <- a$quantities$jnll_comp; jb <- b$quantities$jnll_comp
  dev <- grep("^Selectivity deviates", rownames(ja))

  # dnorm(x, mu, sd, log = TRUE) is the AMAK sum of squares minus
  # log(sd) + 0.5*log(2*pi) per increment, so the gap is that constant times the
  # number of increments SCORED. That count is the loop bound, so this also pins
  # which bins and years the walk covers: 45 increments over 10 of the 11 bins.
  # ADMB's norm2 runs over all 11, counting the plus-group increment twice; that
  # difference is recorded in inst/dev/TODO-nonparametric-penalty-unification.md.
  n_incr <- 45 * 10
  testthat::expect_equal(ja[dev, 2] - jb[dev, 2],
                         n_incr * (log(0.35) + 0.5 * log(2 * pi)))

  # Nothing outside the two selectivity rows moves: same parameters, same data.
  sel <- grep("Non-parametric selectivity|Selectivity deviates", rownames(ja))
  testthat::expect_equal(ja[-sel, ], jb[-sel, ])
})

testthat::test_that("the AMAK level weight is Sel_avgsel_pen, at its exact value", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("Atka2022")
  # Bottom_trawl is time-invariant (one curve, nsex 1), so its level term is
  # exactly weight * avg_sel^2 for a single year. Setting every estimated
  # coefficient to the same value c makes avg_sel = log(mean(exp(c))) = c, and
  # leaves the mean-centred curve flat -- so the decreasing and curvature
  # penalties are 0 and the fleet's whole row is the level term alone.
  cc <- 0.5
  lvl <- function(w) {
    d <- Atka2022
    d$fleet_control$Time_varying_sel <- "Off"
    d$fleet_control$Sel_penalty_form <- "AMAK"
    d$fleet_control$Sel_avgsel_pen   <- w
    m <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = d, inits = NULL, msmMode = 0, estimateMode = 3,
      fit_control = Rceattle::fit_control(verbose = 0))))
    p <- m$obj$par
    p[names(p) == "sel_coff"] <- cc
    r <- m$obj$report(p)
    r$jnll_comp[grep("Non-parametric selectivity", rownames(m$quantities$jnll_comp)), 1]
  }
  testthat::expect_equal(lvl(10), 10 * cc^2)
  testthat::expect_equal(lvl(20), 20 * cc^2)

  # Under the default form the same weight column is not read: the level is
  # pinned at 2 whatever Sel_avgsel_pen says.
  lvl_default <- function(w) {
    d <- Atka2022
    d$fleet_control$Time_varying_sel <- "Off"
    d$fleet_control$Sel_avgsel_pen   <- w
    m <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = d, inits = NULL, msmMode = 0, estimateMode = 3,
      fit_control = Rceattle::fit_control(verbose = 0))))
    p <- m$obj$par
    p[names(p) == "sel_coff"] <- cc
    m$obj$report(p)$jnll_comp[
      grep("Non-parametric selectivity", rownames(m$quantities$jnll_comp)), 1]
  }
  testthat::expect_equal(lvl_default(10), 2 * cc^2)
  testthat::expect_equal(lvl_default(20), 2 * cc^2)
})

testthat::test_that("AMAK with no level weight is refused, not quietly unidentified", {
  testthat::skip_if_not_installed("Rceattle")
  data("Atka2022")

  # Sel_avgsel_pen defaults to 0, and under "AMAK" it IS the level weight. Every
  # other likelihood term reads the per-year mean-centred curve, so at weight 0
  # adding a constant to all coefficients of a fleet changes nothing: the level
  # is a flat direction, and the fit would report a converged, non-identified
  # parameter rather than fail.
  d <- Atka2022
  d$fleet_control$Time_varying_sel <- "Off"   # else the IID guard fires too
  d$fleet_control$Sel_penalty_form <- "AMAK"
  d$fleet_control$Sel_avgsel_pen   <- 0
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = d, inits = NULL, msmMode = 0, estimateMode = 3,
      fit_control = Rceattle::fit_control(verbose = 0)))),
    "Sel_avgsel_pen")

  # The default value of the column reaches the same error, so it cannot be hit
  # by omission.
  d2 <- Atka2022
  d2$fleet_control$Time_varying_sel <- "Off"
  d2$fleet_control$Sel_avgsel_pen   <- NULL
  d2$fleet_control$Sel_penalty_form <- "AMAK"
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = d2, inits = NULL, msmMode = 0, estimateMode = 3,
      fit_control = Rceattle::fit_control(verbose = 0)))),
    "Sel_avgsel_pen")
})

testthat::test_that("an unrecognized penalty form is refused, naming the fleets", {
  testthat::skip_if_not_installed("Rceattle")
  data("GOAcod")
  d <- GOAcod
  d$fleet_control$Sel_penalty_form <- "Nonsense"
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = d, inits = NULL, msmMode = 0, estimateMode = 3,
      fit_control = Rceattle::fit_control(verbose = 0)))),
    "Sel_penalty_form")
})

testthat::test_that("AMAK is refused with IID, which has no ADMB form", {
  testthat::skip_if_not_installed("Rceattle")
  data("Atka2022")

  # The AMAK walk is written for Time_varying_sel = "RandomWalk" only. Left
  # unguarded, "IID" would take the AMAK level weight with the package's own
  # deviation density -- a third likelihood belonging to neither model.
  d <- Atka2022
  d$fleet_control$Time_varying_sel <- as.character(d$fleet_control$Time_varying_sel)
  d$fleet_control$Time_varying_sel[2] <- "IID"
  d$fleet_control$Sel_penalty_form    <- "AMAK"
  d$fleet_control$Sel_avgsel_pen      <- 20
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = d, inits = NULL, msmMode = 0, estimateMode = 3,
      fit_control = Rceattle::fit_control(verbose = 0)))),
    "Off' or 'RandomWalk")
})

testthat::test_that("a shared Selectivity_index resolves the penalty to the lead fleet", {
  testthat::skip_if_not_installed("Rceattle")
  data("Atka2022")

  # ceattle.cpp accumulates the non-parametric penalties once, on the lead of
  # each Selectivity_index group, so "AMAK" written on a mirror is never read.
  d <- Atka2022
  d$fleet_control$Selectivity_index <- c(1, 1)
  d$fleet_control$Time_varying_sel  <- "Off"
  d$fleet_control$Sel_penalty_form  <- c("Rceattle", "AMAK")
  d$fleet_control$Sel_avgsel_pen    <- c(20, 20)
  seen <- character(0)
  suppressMessages(withCallingHandlers(
    try(Rceattle:::data_check(suppressMessages(suppressWarnings(
      Rceattle:::switch_check(d)))), silent = TRUE),
    warning = function(x) { seen <<- c(seen, conditionMessage(x))
                            invokeRestart("muffleWarning") }))
  testthat::expect_true(any(grepl("Sel_penalty_form", seen, fixed = TRUE)))

  # And the level-weight guard reads the lead, not the mirror: "AMAK" on the
  # mirror alone is never accumulated, so a 0 weight there is not an error.
  d2 <- d
  d2$fleet_control$Sel_avgsel_pen <- c(20, 0)
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = d2, inits = NULL, msmMode = 0, estimateMode = 3,
      fit_control = Rceattle::fit_control(verbose = 0)))),
    NA)
})

testthat::test_that("a fleet sharing an index with a different selectivity type is its own lead", {
  testthat::skip_if_not_installed("Rceattle")
  data("Atka2022")

  # rearrange_data() keys flt_sel_lead on Selectivity_index AND selectivity type,
  # so two fleets sharing an index with different Selectivity values get one
  # penalty block EACH. A guard keyed on the index alone would call the second
  # one a mirror, skip it, and let the AMAK level term run at weight 0 -- the
  # unidentified level this guard exists to refuse, reached with no error.
  d <- Atka2022
  d$fleet_control$Selectivity_index <- c(1, 1)
  d$fleet_control$Selectivity       <- c("NonParametricPM", "NonParametric")
  d$fleet_control$Time_varying_sel  <- "Off"
  d$fleet_control$Sel_penalty_form  <- c("Rceattle", "AMAK")
  d$fleet_control$Sel_avgsel_pen    <- c(0, 0)
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = d, inits = NULL, msmMode = 0, estimateMode = 3,
      fit_control = Rceattle::fit_control(verbose = 0)))),
    "Sel_avgsel_pen")

  # And the guard still reads the lead within a genuine same-type mirror group.
  d2 <- Atka2022
  d2$fleet_control$Selectivity_index <- c(1, 1)
  d2$fleet_control$Time_varying_sel  <- "Off"
  d2$fleet_control$Sel_penalty_form  <- c("Rceattle", "AMAK")
  d2$fleet_control$Sel_avgsel_pen    <- c(20, 0)
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = d2, inits = NULL, msmMode = 0, estimateMode = 3,
      fit_control = Rceattle::fit_control(verbose = 0)))),
    NA)
})

testthat::test_that("the guards read 'amak' the same as 'AMAK'", {
  testthat::skip_if_not_installed("Rceattle")
  data("Atka2022")

  # validate_switches() and rearrange_data() both accept the lowercase spelling
  # and route it to the AMAK branch, so a guard comparing only the canonical
  # spelling would wave it through into exactly the states the other two refuse.
  d <- Atka2022
  d$fleet_control$Time_varying_sel <- "Off"
  d$fleet_control$Sel_penalty_form <- "amak"
  d$fleet_control$Sel_avgsel_pen   <- 0
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = d, inits = NULL, msmMode = 0, estimateMode = 3,
      fit_control = Rceattle::fit_control(verbose = 0)))),
    "Sel_avgsel_pen")

  d2 <- Atka2022
  d2$fleet_control$Time_varying_sel <- as.character(d2$fleet_control$Time_varying_sel)
  d2$fleet_control$Time_varying_sel[2] <- "IID"
  d2$fleet_control$Sel_penalty_form    <- "amak"
  d2$fleet_control$Sel_avgsel_pen      <- 20
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = d2, inits = NULL, msmMode = 0, estimateMode = 3,
      fit_control = Rceattle::fit_control(verbose = 0)))),
    "Off' or 'RandomWalk")
})
