# Fleets sharing a Selectivity_index share one selectivity deviation block, so
# Sel_start_year must resolve to the group minimum -- the earliest year any
# sharing fleet has data. Otherwise adjust_map_shared_params() copies the lead
# fleet's mask over the mirror's and the fit depends on fleet_control row order.

testthat::skip_on_cran()

# Two survey fleets, LogisticPM + RandomWalk (a selectivity type whose penalty is
# anchored at Sel_start_year, so the map mask applies).
mirrored_dat <- function(y1, y2, mirror) {
  data("BS2017SS", envir = environment())
  d <- BS2017SS
  d$fleet_control$Sel_start_year <- NA_integer_
  d$fleet_control$Sel_curve_pen1 <- 0
  d$fleet_control$Sel_curve_pen2 <- 0
  d$fleet_control$Sel_curve_pen3 <- 0
  srv <- which(d$fleet_control$Fleet_type == 2)[1:2]
  d$fleet_control$Selectivity[srv]      <- "LogisticPM"
  d$fleet_control$Time_varying_sel[srv] <- "RandomWalk"
  d$fleet_control$Sel_curve_pen1[srv]   <- 50
  d$fleet_control$Sel_curve_pen2[srv]   <- 50
  d$fleet_control$Sel_curve_pen3[srv]   <- 8
  d$fleet_control$Sel_start_year[srv]   <- c(y1, y2)
  if (mirror) d$fleet_control$Selectivity_index[srv] <- 99
  list(data = d, rows = srv, codes = d$fleet_control$Fleet_code[srv])
}

free_devs <- function(y1, y2, mirror) {
  s <- mirrored_dat(y1, y2, mirror)
  fit <- suppressWarnings(suppressMessages(
    Rceattle::fit_mod(data_list = s$data, inits = NULL, estimateMode = 3, msmMode = 0,
                      fit_control = fit_control(phase = FALSE, verbose = 0))))
  mp <- fit$map$mapList$sel_inf_dev
  list(free  = vapply(s$codes, function(f) sum(!is.na(mp[1, f, 1, ])), integer(1)),
       start = fit$obj$env$data$flt_sel_start_yr[s$codes])
}

testthat::test_that("mirrored fleets resolve Sel_start_year to the group minimum", {
  testthat::skip_if_not_installed("TMB")

  a <- free_devs(1982, 1994, mirror = TRUE)
  b <- free_devs(1994, 1982, mirror = TRUE)

  # Both fleets in the group share one mask...
  testthat::expect_equal(a$free[[1]], a$free[[2]])
  testthat::expect_equal(b$free[[1]], b$free[[2]])
  # ...and the result does not depend on fleet_control row order.
  testthat::expect_equal(a$free, b$free)
  testthat::expect_equal(a$start, b$start)
  # The group starts at the EARLIER year (1982, styr 1979 -> 0-based index 3),
  # so deviations informed by the earlier-starting fleet are not masked off.
  testthat::expect_true(all(a$start == 1982L - 1979L))
})

testthat::test_that("unmirrored fleets keep their own Sel_start_year", {
  testthat::skip_if_not_installed("TMB")

  s <- free_devs(1982, 1994, mirror = FALSE)
  # Distinct start years -> distinct masks (the later start frees fewer deviates)
  testthat::expect_equal(s$start, c(1982L, 1994L) - 1979L)
  testthat::expect_gt(s$free[[1]], s$free[[2]])
})

testthat::test_that("data_check warns when a mirrored group has differing Sel_start_year", {
  testthat::skip_if_not_installed("TMB")

  s <- mirrored_dat(1982, 1994, mirror = TRUE)
  testthat::expect_warning(
    suppressMessages(Rceattle::fit_mod(data_list = s$data, inits = NULL, estimateMode = 3,
                                       msmMode = 0,
                                       fit_control = fit_control(phase = FALSE, verbose = 0))),
    regexp = "Selectivity_index"
  )
})
