# data_check() rejects a diet row whose predator or prey age exceeds that
# species' nages. Each species' maximum has to be paired with its OWN nages:
# group_by() drops a species absent from the column, so pairing by position
# puts species on the wrong limit -- in both directions, rejecting a valid
# table and passing an over-range age.

prep_ms <- function(x) {
  x$HCR <- 0                    # build_hcr() would set this
  x$msmMode <- 1                # fit_mod() copies its arguments onto the
  x$suitMode <- rep(0, x$nspp)  # data_list before calling data_check()
  suppressMessages(suppressWarnings(
    Rceattle::switch_check(Rceattle::clean_data(x))))
}

quiet_check <- function(x) suppressWarnings(suppressMessages(data_check(x)))


testthat::test_that("a species absent from the Pred column does not shift the age limits", {
  testthat::skip_on_cran()

  data("BS2017MS", package = "Rceattle", envir = environment())
  cleaned <- prep_ms(BS2017MS)
  # nages is 12 / 12 / 21, so dropping the middle species from Pred would pair
  # arrowtooth's maximum age (21) with cod's nages (12).
  testthat::expect_identical(as.numeric(cleaned$nages), c(12, 12, 21))

  dropped <- cleaned
  dropped$diet_data <- cleaned$diet_data[cleaned$diet_data$Pred != 2, ]
  dropped <- suppressMessages(suppressWarnings(Rceattle::clean_data(dropped)))

  testthat::expect_no_error(quiet_check(dropped))
})


testthat::test_that("a genuinely out-of-range age is still rejected", {
  testthat::skip_on_cran()

  data("BS2017MS", package = "Rceattle", envir = environment())
  cleaned <- prep_ms(BS2017MS)

  over_pred <- cleaned
  over_pred$diet_data$Pred_age[over_pred$diet_data$Pred == 1][1] <-
    cleaned$nages[1] + 1L
  testthat::expect_error(quiet_check(over_pred), "Pred ages in 'diet_data'")

  over_prey <- cleaned
  over_prey$diet_data$Prey_age[over_prey$diet_data$Prey == 1][1] <-
    cleaned$nages[1] + 1L
  testthat::expect_error(quiet_check(over_prey), "Prey ages in 'diet_data'")

  # The aggregated diet formats sit below minage and are not out of range.
  agg <- cleaned
  agg$diet_data$Pred_age[1] <- -999L
  agg$diet_data$Prey_age[1] <- -1L
  testthat::expect_no_error(quiet_check(agg))
})


testthat::test_that("an out-of-range species code neither recycles nor errors", {
  testthat::skip_on_cran()

  data("BS2017MS", package = "Rceattle", envir = environment())
  cleaned <- prep_ms(BS2017MS)

  # A code of 0 drops from `nages[Pred]` and leaves a short vector for the
  # comparison to recycle against; a negative one errors. Both must map to NA,
  # leaving the species-range checks to report the bad code.
  for (bad in c(0L, -1L)) {
    odd <- cleaned
    odd$diet_data$Pred[1] <- bad
    odd$diet_data$Prey[1] <- bad
    lim <- cleaned$nages[match(odd$diet_data$Pred, seq_along(cleaned$nages))]
    testthat::expect_length(lim, nrow(odd$diet_data))
    testthat::expect_true(is.na(lim[1]))
    testthat::expect_no_warning(
      any(odd$diet_data$Pred_age > lim, na.rm = TRUE))
  }
})
