# data_check() is callable on a list straight from read_data(), where the switch
# columns are still the integer codes the workbook stores -- and every bundled
# data set that carries a fleet_control stores them that way, so the integer is
# the shipped representation, not an edge case. `0 != "Off"` coerces to
# `"0" != "Off"`, which is TRUE, so an Off fleet read as ESTIMATED and was named
# among the fleets with estimated selectivity and no data to fit it.
# Found reviewing the 5.34.0-5.43.0 release (#158).

# GOA2018SS fleet 7 is Fleet_type 0 with an estimated Selectivity and no
# composition data. Given its own Selectivity_index it cannot borrow another
# fleet's data either, so it is the case the defect acted on.
raw_with_off_fleet <- function() {
  d <- Rceattle::GOA2018SS
  d$fleet_control$Selectivity_index[7] <- 99
  d
}
off_fleet_name <- function(d) d$fleet_control$Fleet_name[7]

msgs_of <- function(d) {
  out <- character(0)
  withCallingHandlers(
    tryCatch(Rceattle:::data_check(d), error = function(e) out <<- c(out, conditionMessage(e))),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage"))
  gsub("\n", " ", paste(out, collapse = " "))
}

testthat::test_that("data_check() does not call an integer-coded Off fleet estimated", {
  d <- raw_with_off_fleet()

  # Positive control: the fixture really does reach the check it is meant to.
  # Without this the negative assertion below could rot to green if the data
  # or the message ever changed and the condition stopped reproducing.
  testthat::expect_match(msgs_of(d), "estimated Selectivity but no comp_data")

  # And the Off fleet is not among the fleets it names.
  testthat::expect_equal(as.character(d$fleet_control$Fleet_type[7]), "0")
  testthat::expect_false(grepl(off_fleet_name(d), msgs_of(d), fixed = TRUE))
})

testthat::test_that("the old spelling is what made it estimated", {
  # Documents the defect rather than the fix: kept so that a future reader who
  # wonders why .canon_switch() is used here can see what `!=` did instead.
  ft <- Rceattle::GOA2018SS$fleet_control$Fleet_type
  off <- which(ft == 0)
  testthat::expect_gt(length(off), 0)
  testthat::expect_true(all((ft != "Off")[off]))                             # read as live
  testthat::expect_false(any((Rceattle:::.canon_switch(
    ft, Rceattle:::fleet_map) != "Off")[off]))                               # read as Off
})

testthat::test_that("the canonical path is untouched", {
  # fit_mod(), build_map() and build_params() canonicalize first, so the change
  # must be inert on a canonical fleet_control.
  for (nm in c("BS2017SS", "GOA2018SS", "GOApollock", "GeorgesBank3spp", "GOAatf")) {
    d <- suppressMessages(Rceattle:::switch_check(
      get(nm, envir = as.environment("package:Rceattle"))))
    ft <- d$fleet_control$Fleet_type
    testthat::expect_identical(Rceattle:::.canon_switch(ft, Rceattle:::fleet_map) != "Off",
                               ft != "Off")
  }
})
