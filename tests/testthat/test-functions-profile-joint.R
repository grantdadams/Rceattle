# Argument contract for profile()'s joint modes and the "q" alias.
#
# These assert the checks that fire BEFORE any fit, so they need no model and
# stay unguarded (fast). The joint modes exist because slots otherwise cross
# through expand.grid(): moving a ten-age M schedule together over 13 values is
# 13^10 fits under joint = "none", not 13. The behaviour of the fits themselves
# is asserted in test-functions-profile-param.R.

# Enough of a fitted object for profile() to validate against: it reads the
# class, `estimated_params`, and (for the q alias) `data_list$fleet_control`.
fake_fit <- function(qi = c(1L, 1L, 2L),
                     qform = c("Estimated", "Estimated", "Estimated")) {
  structure(list(
    sdrep = NULL,
    estimated_params = list(
      log_M1      = array(log(0.3), dim = c(1L, 1L, 3L)),
      index_log_q = stats::setNames(log(c(1, 0.8, 0.5)),
                                    c("Shelikof", "ADFG", "NMFS_BT"))),
    data_list = list(fleet_control = data.frame(
      Fleet_name = c("Shelikof", "ADFG", "NMFS_BT"),
      Catchability = qform,
      Catchability_index = qi, stringsAsFactors = FALSE))
  ), class = "Rceattle")
}

age_slots <- list(c(1, 1, 1), c(1, 1, 2), c(1, 1, 3))


testthat::test_that("a joint mode takes one grid, however many cells it moves", {
  # One vector covers every cell -- that is the whole point. Passing one per
  # slot is the crossing mode's contract and means the caller wants `none`.
  testthat::expect_error(
    profile(fake_fit(), param = "M1", slots = age_slots,
            values = list(1:3, 1:3, 1:3), joint = "multiply"),
    "must be a single vector")

  # And the crossing mode still demands one grid per slot.
  testthat::expect_error(
    profile(fake_fit(), param = "M1", slots = age_slots, values = list(1:3)),
    "same length as `slots`")
})


testthat::test_that("multiply and add refuse a transform they cannot invert", {
  # They act on the natural scale, so the stored value has to be recoverable.
  # Applying the multiplier to the stored value instead would silently scale a
  # log, which is a different model, not a different parameterisation.
  testthat::expect_error(
    profile(fake_fit(), param = "log_M1", slots = age_slots,
            values = list(seq(0.8, 1.2, by = 0.1)), joint = "multiply",
            transform = stats::qlogis),
    "natural scale")
  testthat::expect_error(
    profile(fake_fit(), param = "log_M1", slots = age_slots,
            values = list(seq(-0.1, 0.1, by = 0.05)), joint = "add",
            transform = stats::qlogis),
    "natural scale")
})


testthat::test_that("the q alias resolves a fleet by name", {
  testthat::expect_error(
    profile(fake_fit(), param = "q", slots = list("Sheliokf"),
            values = list(c(0.5, 1))),
    "No fleet named")
  # The message names the fleets, so the typo is fixable without opening the
  # workbook.
  testthat::expect_error(
    profile(fake_fit(), param = "q", slots = list("Sheliokf"),
            values = list(c(0.5, 1))),
    "Shelikof, ADFG, NMFS_BT")
})


testthat::test_that("an analytical catchability cannot be profiled", {
  # Catchability 3 / 7 solve q from the fleet's own index and overwrite
  # index_q, so index_log_q never reaches the likelihood. The grid would move
  # nothing and every point would return the same fit -- a flat curve that
  # reads as "the data do not inform q" when the profile did nothing.
  for (form in c("Analytical", "AnalyticalArith")) {
    testthat::expect_error(
      profile(fake_fit(qform = c(form, "Estimated", "Estimated")),
              param = "q", slots = list("Shelikof"), values = list(c(0.5, 1))),
      "ignores `index_log_q`")
  }

  # And an analytical fleet is not dragged into a group it does not share:
  # its q is its own, whatever the Catchability_index says.
  testthat::expect_silent(suppressWarnings(
    try(profile(fake_fit(qi = c(1L, 1L, 2L),
                         qform = c("Estimated", "Analytical", "Estimated")),
                param = "q", slots = list("Shelikof"),
                values = list(c(0.5, 1)), cores = 1),
        silent = TRUE)))
})


testthat::test_that("a shared Catchability_index group is profiled together", {
  # Fleets sharing a Catchability_index share ONE q parameter -- the map copies
  # the lead's slice across the group -- so fixing one member's cell alone would
  # leave the others estimating a common q, which is not that fleet's q
  # profiled. `try()` because the fake carries no model to fit once the checks
  # pass; the expansion happens before any fit.
  testthat::expect_message(
    try(profile(fake_fit(qi = c(1L, 1L, 2L)), param = "q",
                slots = list("Shelikof"), values = list(c(0.5, 1)), cores = 1),
        silent = TRUE),
    "share a Catchability_index")

  # Naming BOTH members of one group adds no fleet, so a count-based check
  # would see no growth and leave joint = "none" -- crossing two grids onto the
  # single q parameter that group shares.
  p <- try(profile(fake_fit(qi = c(1L, 1L, 2L)), param = "q",
                   slots = list("Shelikof", "ADFG"),
                   values = list(c(0.5, 1), c(0.5, 1)), cores = 1),
           silent = TRUE)
  testthat::expect_true(inherits(p, "try-error"))
  # It fails on the fake having no model to fit, NOT on values/slots mismatch:
  # reaching that error means joint was switched and one grid was accepted.
  testthat::expect_false(grepl("same length as", conditionMessage(attr(p, "condition"))))

  # A fleet in a group of its own is left alone, and so is one whose
  # Catchability_index is NA -- `qi == NA` selects nothing, which must be read
  # as "no group" rather than silently dropping the slot.
  testthat::expect_silent(suppressWarnings(
    try(profile(fake_fit(qi = c(1L, 2L, 3L)), param = "q",
                slots = list("Shelikof"), values = list(c(0.5, 1)), cores = 1),
        silent = TRUE)))
  testthat::expect_silent(suppressWarnings(
    try(profile(fake_fit(qi = c(NA, NA, NA)), param = "q",
                slots = list("Shelikof"), values = list(c(0.5, 1)), cores = 1),
        silent = TRUE)))
})


testthat::test_that("a joint grid that drives a log-scale cell to zero is an error", {
  # exp() of the result would be the fit silently continuing from a nonsense M.
  # The message names the cell and the grid value that did it.
  # Every grid value offends, so the check is what stops the run rather than
  # the fake's missing map further down.
  testthat::expect_error(
    profile(fake_fit(), param = "M1", slots = age_slots,
            values = list(-2), joint = "multiply", cores = 1),
    "not positive")
  # M1 here is 0.3, so an offset of -0.5 takes it below zero.
  testthat::expect_error(
    profile(fake_fit(), param = "M1", slots = age_slots,
            values = list(-0.5), joint = "add", cores = 1),
    "not positive")
  # The message says which cell and which grid value, so it is actionable.
  testthat::expect_error(
    profile(fake_fit(), param = "M1", slots = age_slots,
            values = list(-2), joint = "multiply", cores = 1),
    "c\\(1, 1, 1\\)")
})
