# Fleets sharing a Selectivity_index share ONE selectivity deviation sd
# (sel_dev_log_sd). TMB starts a shared parameter from the first group member's
# value, so a different Time_varying_sel_sd on any other fleet is discarded.
#
# The discard only happens when that sd is estimated, which is random_sel = TRUE
# with a random-effect Time_varying_sel. With random_sel = FALSE the map is NA
# and each fleet keeps its own value, so warning there would be a false alarm.
# Both directions are pinned here.

make_shared_sel_data <- function(sd1 = 0.3, sd2 = 0.7, tv_sel = "IID") {
  d <- make_test_data()
  fc <- d$fleet_control
  # Put both fleets in one selectivity group with the same curve, so the only
  # thing they disagree on is the deviation sd.
  fc$Selectivity_index   <- c(1, 1)
  fc$Selectivity         <- rep(fc$Selectivity[1], 2)
  fc$Time_varying_sel    <- rep(tv_sel, 2)
  fc$Time_varying_sel_sd <- c(sd1, sd2)
  d$fleet_control <- fc
  d
}

build_shared_sel_map <- function(d, random_sel) {
  suppressMessages(
    Rceattle::fit_mod(d, estimateMode = 3, msmMode = 0, random_sel = random_sel)
  )
}


test_that("a shared Selectivity_index warns when an estimated sd differs", {
  d <- make_shared_sel_data(sd1 = 0.3, sd2 = 0.7)

  expect_warning(
    build_shared_sel_map(d, random_sel = TRUE),
    "different Time_varying_sel_sd"
  )

  # The message names the group and the value actually used, so the assessor can
  # see which fleet won rather than only that someone lost. Collect every
  # warning: build_map raises others (the Time_varying_sel soft-deprecation)
  # ahead of this one.
  w <- character(0)
  withCallingHandlers(
    build_shared_sel_map(d, random_sel = TRUE),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond))
      invokeRestart("muffleWarning")
    }
  )
  hit <- grep("Time_varying_sel_sd", w, value = TRUE)
  expect_length(hit, 1)
  expect_match(hit, "Selectivity_index 1")
  expect_match(hit, "0.3", fixed = TRUE)
})


test_that("agreeing fleets are silent", {
  d <- make_shared_sel_data(sd1 = 0.5, sd2 = 0.5)
  expect_no_warning(
    withCallingHandlers(
      build_shared_sel_map(d, random_sel = TRUE),
      warning = function(w) {
        if (!grepl("Time_varying_sel_sd", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    )
  )
})


test_that("random_q estimates the catchability deviation sd, and random_sel the selectivity one", {
  # The two are symmetric: integrating the deviates out is what makes their sd a
  # marginal variance rather than one half of a degenerate joint mode. Without
  # the flag the sd stays at its fleet_control value.
  d <- make_test_data()
  d$fleet_control$Catchability[1]      <- "Estimated"
  d$fleet_control$Time_varying_q[1]    <- "IID"
  d$fleet_control$Time_varying_q_sd[1] <- 0.2

  est <- function(rq) {
    m <- suppressWarnings(suppressMessages(
      Rceattle::fit_mod(d, estimateMode = 3, msmMode = 0, random_q = rq)))
    m$map$mapFactor
  }

  expect_true(all(is.na(est(FALSE)$index_q_dev_log_sd)))
  expect_true(any(!is.na(est(TRUE)$index_q_dev_log_sd)))

  # index_q_log_sd is a prior sd the assessor sets; random_q must not free it.
  expect_true(all(is.na(est(TRUE)$index_q_log_sd)))
})


test_that("a fixed deviation sd is honoured per fleet, so it does not warn", {
  # random_sel = FALSE leaves sel_dev_log_sd mapped out. Nothing is shared and
  # nothing is discarded, so differing values are legitimate.
  d <- make_shared_sel_data(sd1 = 0.3, sd2 = 0.7)
  m <- build_shared_sel_map(d, random_sel = FALSE)
  expect_true(all(is.na(m$map$mapFactor$sel_dev_log_sd)))

  expect_no_warning(
    withCallingHandlers(
      build_shared_sel_map(d, random_sel = FALSE),
      warning = function(w) {
        if (!grepl("Time_varying_sel_sd", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    )
  )
})
