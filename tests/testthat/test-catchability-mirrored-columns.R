# Fleets sharing a Catchability_index share one index_log_q, and the lead fleet
# -- the group's first non-Off fleet -- settles whether it is estimated. True for
# Fixed / Estimated / Estimated-with-prior.
#
# Not for the solved forms: Analytical and AnalyticalArith solve q from each
# fleet's own residuals, Environmental and AR1 build it from each fleet's own
# covariate, each overwriting index_q per fleet. A group containing one shares no
# catchability, even when every fleet in it agrees. data_check() reports both.

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


testthat::test_that("Time_varying_q is not compared where it holds env indices", {
  # Under Environmental and AR1 the column carries env_data indices, not a mode,
  # so comparing them across the group is noise. Only the solved-form warning
  # should fire.
  w <- .q_warnings(.q_group_dat(lead_q = "AR1", member_q = "AR1",
                                lead_tvq = "1", member_tvq = "2"))
  testthat::expect_length(w, 1)
  testthat::expect_match(w[1], "does not share", fixed = TRUE)
  testthat::expect_false(any(grepl("Time_varying_q", w, fixed = TRUE)))
})


testthat::test_that("the bundled datasets are clean", {
  # GOA2018SS fleets 9 and 10 share selectivity and q -- the worked example -- so
  # a warning here would mean the check is wrong, not the data.
  for (nm in c("BS2017SS", "BS2017MS", "GOA2018SS")) {
    d <- get(nm, envir = asNamespace("Rceattle"))
    testthat::expect_length(.q_warnings(d), 0)
  }
})
