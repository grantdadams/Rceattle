# The pre-4.4 non-parametric workbook format, and how switch_check() recognises it.
#
# Provenance. Those workbooks had no `Sel_curve_pen` COLUMNS and carried the two
# AMAK penalty weights in `Time_varying_sel` / `Time_varying_sel_sd` instead. The
# upgrade that moves them across used to be triggered by the VALUE in
# `Time_varying_sel`, which cannot tell a shape weight of 4 from the `RandomWalk`
# code 4 -- so `Time_varying_sel = 4` on a modern workbook silently turned the
# walk off and overwrote GOAcod's 20 / 12.5 penalty weights with 4 and the
# deviation sd, while the string "RandomWalk" worked. From 5.25.0 the trigger is
# the columns' absence, which is what actually distinguishes the two formats:
# `read_data()` leaves them out for a pre-4.4 workbook (verified against
# `atka_single_species_2022.xlsx`) and fills them for a modern one.
testthat::skip_on_cran()

testthat::test_that("an integer Time_varying_sel is a mode, not a legacy penalty weight", {
  testthat::skip_if_not_installed("Rceattle")

  data("GOAcod")
  d <- GOAcod
  d$fleet_control$Time_varying_sel <- as.character(d$fleet_control$Time_varying_sel)
  d$fleet_control$Time_varying_sel[2] <- "4"
  d$fleet_control$Time_varying_sel_sd[2] <- 0.2
  x <- suppressMessages(suppressWarnings(Rceattle:::switch_check(d)))
  testthat::expect_equal(x$fleet_control$Time_varying_sel[2], "RandomWalk")
  testthat::expect_equal(x$fleet_control$Sel_curve_pen1[2], 20)
  testthat::expect_equal(x$fleet_control$Sel_curve_pen2[2], 12.5)

  # A genuinely legacy workbook -- no Sel_curve_pen columns at all -- still
  # upgrades: the weights move across and both time-varying columns are cleared.
  d2 <- GOAcod
  d2$fleet_control$Sel_curve_pen1 <- NULL
  d2$fleet_control$Sel_curve_pen2 <- NULL
  d2$fleet_control$Sel_curve_pen3 <- NULL
  np <- which(d2$fleet_control$Selectivity == 2)
  d2$fleet_control$Time_varying_sel[np]    <- 20
  d2$fleet_control$Time_varying_sel_sd[np] <- 12.5
  x2 <- suppressMessages(suppressWarnings(Rceattle:::switch_check(d2)))
  testthat::expect_equal(x2$fleet_control$Sel_curve_pen1[np], rep(20, length(np)))
  testthat::expect_equal(x2$fleet_control$Sel_curve_pen2[np], rep(12.5, length(np)))
  testthat::expect_equal(x2$fleet_control$Time_varying_sel[np], rep("Off", length(np)))
  testthat::expect_equal(x2$fleet_control$Time_varying_sel_sd[np], rep(0, length(np)))

  # And a blank cell in a column that EXISTS is not a legacy workbook: it is a
  # modern one missing a required weight, and is refused rather than guessed at.
  d3 <- GOAcod
  d3$fleet_control$Sel_curve_pen1[2] <- NA
  d3$fleet_control$Time_varying_sel[2] <- 4
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle:::switch_check(d3))),
    "Sel_curve_pen1")
})
