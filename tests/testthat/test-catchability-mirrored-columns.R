# Fleets sharing a Catchability_index share one index_log_q, and the lead fleet
# -- the group's first non-Off fleet -- settles whether it is estimated. True for
# Fixed / Estimated / Estimated-with-prior.
#
# Not for Analytical / AnalyticalArith: those solve q from each fleet's own index
# OBSERVATIONS, bypassing the shared parameter, so a group containing one shares
# no catchability even when every fleet in it agrees. Environmental and AR1
# rebuild index_q too, but from the group's shared parameters, so they DO share.
# data_check() reports the analytical case.

testthat::skip_on_cran()

# Two survey fleets in one Catchability_index group. The second fleet's index is
# on a different scale, so a per-fleet solved q lands where the shared one cannot.
.q_group_dat <- function(lead_q = "Estimated", member_q = "Estimated",
                         lead_tvq = NULL, member_tvq = NULL) {
  d   <- make_test_data(nyrs = 15, nages = 6, seed = 42)
  fc  <- d$fleet_control
  srv <- which(fc$Fleet_type == "Survey")[1]
  testthat::skip_if(is.na(srv), "fixture lacks a survey fleet")

  m <- fc[srv, ]
  m$Fleet_name         <- "Survey2"
  m$Fleet_code         <- 3L
  m$Selectivity_index  <- fc$Selectivity_index[srv]
  m$Catchability_index <- fc$Catchability_index[srv]   # one q group
  m$Catchability       <- member_q
  m$Catchability_init  <- 1
  m$Estimate_index_sd  <- 0
  if (!is.null(member_tvq)) m$Time_varying_q <- member_tvq
  d$fleet_control <- rbind(fc, m)

  d$fleet_control$Catchability[srv]      <- lead_q
  d$fleet_control$Catchability_init[srv] <- 1
  if (!is.null(lead_tvq)) d$fleet_control$Time_varying_q[srv] <- lead_tvq

  add <- d$index_data[d$index_data$Fleet_code == fc$Fleet_code[srv], , drop = FALSE]
  add$Fleet_code   <- 3L
  add$Fleet_name   <- "Survey2"
  add$Observation  <- add$Observation * 3
  d$index_data <- rbind(d$index_data, add)
  d
}

# data_check() is internal, so call it bare -- the test env's parent is the
# package namespace. Rceattle::data_check() errors, and try() would swallow that
# into an empty warning list, passing every assertion for nothing.
.q_warnings <- function(d) {
  w <- character(0)
  withCallingHandlers(
    suppressMessages(try(
      data_check(suppressMessages(switch_check(d))), silent = TRUE)),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  grep("sharing Catchability_index", w, value = TRUE)
}


testthat::test_that("a consistent estimated q group raises nothing", {
  testthat::expect_length(.q_warnings(.q_group_dat()), 0)
})


testthat::test_that("a q setting that differs within the group is reported as overridden", {
  # Fixed against Estimated: the group shares the parameter and the lead decides.
  # Worth saying, but not divergence.
  w <- .q_warnings(.q_group_dat(lead_q = "Estimated", member_q = "Fixed"))
  testthat::expect_length(w, 1)
  testthat::expect_match(w[1], "different Catchability", fixed = TRUE)
  testthat::expect_match(w[1], "the others are ignored", fixed = TRUE)
})


testthat::test_that("a solved q form in a shared group is reported as not sharing", {
  for (form in c("Analytical", "AnalyticalArith")) {
    w <- .q_warnings(.q_group_dat(lead_q = "Estimated", member_q = form))
    testthat::expect_length(w, 1)
    testthat::expect_match(w[1], form, fixed = TRUE, info = form)
    testthat::expect_match(w[1], "does not share", fixed = TRUE, info = form)
  }
})


testthat::test_that("a solved q form is reported even when the whole group agrees", {
  # The case a mismatch check cannot see: two Analytical fleets each solve from
  # their own residuals, so the group shares the mapped parameter but not the
  # catchability the model uses.
  w <- .q_warnings(.q_group_dat(lead_q = "Analytical", member_q = "Analytical"))
  testthat::expect_length(w, 1)
  testthat::expect_match(w[1], "does not share", fixed = TRUE)
})


testthat::test_that("the reported divergence is real, not just a fleet_control reading", {
  # Compare the catchabilities the model actually used -- otherwise the warning
  # above is an assertion about a data frame, not about the model.
  d <- .q_group_dat(lead_q = "Estimated", member_q = "Analytical")
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                        phase = FALSE))))
  qq <- fit$quantities$index_q
  testthat::expect_false(isTRUE(all.equal(as.numeric(qq[1, 1]),
                                          as.numeric(qq[3, 1]))))

  # ...and an all-Estimated group of the same shape does share one, so the
  # comparison above is not just picking up two unrelated fleets.
  d2 <- .q_group_dat(lead_q = "Estimated", member_q = "Estimated")
  fit2 <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d2, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                        phase = FALSE))))
  q2 <- fit2$quantities$index_q
  testthat::expect_equal(as.numeric(q2[1, 1]), as.numeric(q2[3, 1]),
                         tolerance = 1e-10)
})


testthat::test_that("Time_varying_q that differs within the group is reported", {
  w <- .q_warnings(.q_group_dat(lead_tvq = "Off", member_tvq = "RandomWalk"))
  testthat::expect_length(w, 1)
  testthat::expect_match(w[1], "different Time_varying_q", fixed = TRUE)
})


testthat::test_that("Environmental and AR1 share, and are not reported as solved", {
  # These rebuild index_q per fleet, but from index_log_q / index_q_beta /
  # index_q_dev, all mapped to the lead fleet's, against an env_index row that is
  # not fleet specific. So the group does share, and must not be flagged as if it
  # did not. Only the Time_varying_q difference is worth saying.
  for (form in c("AR1", "Environmental")) {
    w <- .q_warnings(.q_group_dat(lead_q = form, member_q = form,
                                  lead_tvq = "1", member_tvq = "2"))
    testthat::expect_false(any(grepl("does not share", w, fixed = TRUE)),
                           info = form)
    testthat::expect_true(any(grepl("different Time_varying_q", w, fixed = TRUE)),
                          info = form)
  }
})


testthat::test_that("an Environmental group shares its q, on a coefficient that is not zero", {
  testthat::skip_if_not_installed("TMB")
  # The claim, measured rather than reasoned, and made discriminating on two
  # counts. At the parameter inits every env coefficient is 0, so index_q
  # collapses to exp(index_log_q) and the comparison would only re-prove that
  # index_log_q is map-tied -- equally true of a plain Estimated group. So
  # index_q_beta is perturbed off zero first, and the two fleets are given
  # DIFFERENT Time_varying_q, which is the documented claim: the lead's env
  # series is the one the group uses.
  d <- .q_group_dat(lead_q = "Environmental", member_q = "Environmental",
                    lead_tvq = "1", member_tvq = "2")
  d$fleet_control$Catchability_init <- c(2, 5, 5)[seq_len(nrow(d$fleet_control))]
  n <- d$endyr - d$styr + 1L
  d$env_data <- data.frame(Year  = d$styr:d$endyr,
                           BTemp = as.numeric(scale(sin(seq_len(n)))),
                           Ice   = as.numeric(scale(cos(seq_len(n)))))
  fit_at <- function(dl, mutate_pars = identity) {
    f <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
      dl, file = NULL, estimateMode = 3, msmMode = 0, random_rec = FALSE,
      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                          phase = FALSE))))
    p <- mutate_pars(f$estimated_params)
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      dl, inits = p, file = NULL, estimateMode = 3, msmMode = 0,
      random_rec = FALSE,
      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                          phase = FALSE))))
  }
  fit <- fit_at(d, function(p) { p$index_q_beta[] <- 0.4; p })
  qq  <- fit$quantities$index_q

  # The covariate is live, so this is not the flat exp(index_log_q) case.
  testthat::expect_gt(length(unique(round(as.numeric(qq[1, ]), 10))), 1L)
  # ...and the two fleets still agree, year by year, despite naming different
  # env series and carrying different Catchability_init.
  testthat::expect_equal(as.numeric(qq[1, ]), as.numeric(qq[3, ]),
                         tolerance = 1e-10)

  # Negative control: the same configuration with the member on its OWN
  # Catchability_index must NOT agree, or the check above proves nothing.
  d2 <- d
  d2$fleet_control$Catchability_index[3] <- 3L
  q2 <- fit_at(d2, function(p) { p$index_q_beta[] <- 0.4; p })$quantities$index_q
  testthat::expect_false(isTRUE(all.equal(as.numeric(q2[1, ]),
                                          as.numeric(q2[3, ]))))
})


testthat::test_that("a shared q group carrying AR1 catchability is refused, not shared", {
  testthat::skip_if_not_installed("TMB")
  # This used to pin the inert state: build_map() gates index_q_dev on
  # Time_varying_q %in% c("IID","AR1","RandomWalk"), but for
  # Catchability = "AR1" that column holds an env_data column index, so the
  # deviates were mapped out entirely and index_q was constant -- which made
  # "AR1 shares" indistinguishable from "AR1 has nothing to share".
  #
  # Settled in 5.12.0: the switch is refused outright, so the question does not
  # arise. Sharing under a q linkage is covered by test-linkage-catchability.R;
  # the refusal itself by test-catchability-qar1-deprecation.R. This asserts the
  # group case specifically, because a lead and a member both set to AR1 is how
  # the inert path was reached.
  d <- .q_group_dat(lead_q = "AR1", member_q = "AR1",
                    lead_tvq = "1", member_tvq = "1")
  n <- d$endyr - d$styr + 1L
  d$env_data <- data.frame(Year = d$styr:d$endyr,
                           BTemp = as.numeric(scale(sin(seq_len(n)))))

  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      d, file = NULL, estimateMode = 3, msmMode = 0, random_rec = FALSE,
      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                          phase = FALSE)))),
    "QAR1")
})


testthat::test_that("the bundled datasets are clean", {
  # GOA2018SS fleets 9 and 10 share selectivity and q -- the worked example -- so
  # a warning here would mean the check is wrong, not the data.
  for (nm in c("BS2017SS", "BS2017MS", "GOA2018SS")) {
    d <- get(nm, envir = asNamespace("Rceattle"))
    testthat::expect_length(.q_warnings(d), 0)
  }
})
