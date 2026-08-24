# Two ways of specifying an MSE that run_mse() reads off its arguments before
# any model is fitted: the years an assessment is completed, and the multiplier
# applied to the catch the control rule recommends.
#
# Provenance: both exist for a design in which the schedule is IRREGULAR -- one
# assessment missed inside an otherwise biennial cycle, with a precautionary
# reduction held only for the years advice is stale. A scalar period cannot
# express the first (it is a permanent policy, a different question) and a
# per-species multiplier cannot express the second (it applies in every
# projection year, which is a permanent harvest cut).
#
# Tested against the two internal helpers rather than a run_mse() call: these
# are bookkeeping rules, and fitting an operating and an estimation model to
# exercise them would take minutes.

test_that("a scalar assessment_period is the regular grid it always was", {
  # endyr 2020, projyr 2030: biennial assessments from 2022.
  expect_equal(.mse_assess_years(2, om_endyr = 2020, max_yr = 2030),
               seq(2022, 2030, by = 2))
  expect_equal(.mse_assess_years(1, om_endyr = 2020, max_yr = 2025),
               2021:2025)

  # The grid stops at max_yr rather than stepping past it.
  expect_equal(.mse_assess_years(4, om_endyr = 2020, max_yr = 2030),
               c(2024, 2028))
})


test_that("a vector assessment_period is the schedule itself", {
  biennial <- seq(2022, 2030, by = 2)
  skipped  <- setdiff(biennial, 2026)

  expect_equal(.mse_assess_years(skipped, om_endyr = 2020, max_yr = 2030),
               c(2022, 2024, 2028, 2030))

  # Sorted and de-duplicated, so the order the caller happened to build them in
  # does not change which year the loop treats as terminal.
  expect_equal(.mse_assess_years(c(2028, 2022, 2028, 2024),
                                 om_endyr = 2020, max_yr = 2030),
               c(2022, 2024, 2028))

  # The skipped schedule is NOT the triennial grid it superficially resembles.
  expect_false(identical(.mse_assess_years(skipped, 2020, 2030),
                         .mse_assess_years(3, 2020, 2030)))
})


test_that("assessment years outside the projection window are refused, not dropped", {
  # Before or at the operating model's terminal year: there is no projection
  # year for the assessment to advance the model to.
  expect_error(.mse_assess_years(c(2019, 2022), om_endyr = 2020, max_yr = 2030),
               "after the operating model's terminal year")
  expect_error(.mse_assess_years(c(2020, 2022), om_endyr = 2020, max_yr = 2030),
               "2020")

  # Past the horizon. Silently dropping it would run a shorter schedule than
  # the caller asked for, with nothing in the result to say so.
  expect_error(.mse_assess_years(c(2022, 2032), om_endyr = 2020, max_yr = 2030),
               "within the projection horizon")

  # A fractional year would set the operating model's endyr to a fraction, and
  # nothing downstream indexes on one.
  expect_error(.mse_assess_years(2.5, om_endyr = 2020, max_yr = 2030),
               "whole years")
  expect_error(.mse_assess_years(0, om_endyr = 2020, max_yr = 2030),
               "at least 1 year")
  expect_error(.mse_assess_years(c(2022, NA), om_endyr = 2020, max_yr = 2030),
               "must be a number of years")
})


test_that("one year on its own is refused rather than read as a period", {
  # The overload's only blind spot: length 1 is a period, so a schedule naming
  # a single assessment year cannot be told from one. Read as a period, 2025
  # would start the grid in 4045 and seq() would fail on the sign of `by`;
  # the point of the guard is that the message says which reading was taken.
  expect_error(.mse_assess_years(2025, om_endyr = 2020, max_yr = 2030),
               "number of years BETWEEN assessments")

  # The same guard catches a period simply longer than the projection, which
  # is the case that used to reach seq() and report on its `by` argument.
  expect_error(.mse_assess_years(20, om_endyr = 2020, max_yr = 2030),
               "no assessment would run")

  # A period that just reaches the last year is fine.
  expect_equal(.mse_assess_years(10, om_endyr = 2020, max_yr = 2030), 2030)
})


test_that("a vector catch_mult still applies to every year of every species", {
  # Three species, two projection years, one fleet each.
  catch   <- c(100, 200, 300, 100, 200, 300)
  year    <- rep(c(2021, 2022), each = 3)
  species <- rep(1:3, times = 2)

  out <- .mse_apply_catch_mult(catch, year, species, c(0.9, 1, 0.5))
  expect_equal(out, c(90, 200, 150, 90, 200, 150))
})


test_that("a data.frame catch_mult applies only in the years and species it lists", {
  catch   <- c(100, 200, 300, 100, 200, 300)
  year    <- rep(c(2021, 2022), each = 3)
  species <- rep(1:3, times = 2)

  # A 10% reduction in 2022 only, and on species 1 and 3 only.
  cm <- data.frame(Year = c(2022, 2022), Species = c(1, 3), mult = 0.9)

  out <- .mse_apply_catch_mult(catch, year, species, cm)

  # 2021 is untouched -- this is the whole point of the year-indexed form.
  expect_equal(out[1:3], c(100, 200, 300))
  # 2022 species 2 is untouched too: a pair the table omits is multiplied by 1.
  expect_equal(out[4:6], c(90, 200, 270))

  # And it is not the same as the equivalent per-species vector, which would
  # have cut 2021 as well.
  expect_false(isTRUE(all.equal(
    out, .mse_apply_catch_mult(catch, year, species, c(0.9, 1, 0.9)))))
})


test_that("a data.frame catch_mult with no matching row leaves catch unchanged", {
  catch   <- c(100, 200)
  year    <- c(2021, 2021)
  species <- c(1, 2)

  cm <- data.frame(Year = 2022, Species = 1, mult = 0.5)
  expect_equal(.mse_apply_catch_mult(catch, year, species, cm), catch)
})


test_that("catch_mult is validated at the call, in either form", {
  # Vector form: recycled from length 1, checked against nspp, unchanged
  # message so an existing script's error still reads the same.
  expect_equal(.mse_check_catch_mult(0.8, nspp = 3, 2021, 2030),
               rep(0.8, 3))
  expect_error(.mse_check_catch_mult(c(0.8, 0.9), nspp = 3, 2021, 2030),
               "not length 1 or length nspp")
  expect_error(.mse_check_catch_mult(-1, nspp = 3, 2021, 2030),
               "finite and non-negative")
  expect_null(.mse_check_catch_mult(NULL, nspp = 3, 2021, 2030))

  # A 3-column data.frame has length 3. Without a data.frame branch ahead of
  # the length test it would pass silently at nspp = 3 and then be indexed as
  # a vector -- the reason the check is on is.data.frame() first.
  cm <- data.frame(Year = 2022, Species = 1, mult = 0.9)
  expect_identical(.mse_check_catch_mult(cm, nspp = 3, 2021, 2030), cm)

  expect_error(.mse_check_catch_mult(data.frame(Year = 2022, mult = 0.9),
                                     nspp = 3, 2021, 2030),
               "missing Species")

  # A year outside the projection could never be applied, so a typo would
  # otherwise read as "no reduction" in a run that looks like it worked.
  expect_error(.mse_check_catch_mult(data.frame(Year = 2019, Species = 1, mult = 0.9),
                                     nspp = 3, 2021, 2030),
               "must be a projection year")

  expect_error(.mse_check_catch_mult(data.frame(Year = 2022, Species = 4, mult = 0.9),
                                     nspp = 3, 2021, 2030),
               "species number in 1:3")

  # match() takes the first of a duplicated key, so the second row would be
  # discarded without a word.
  expect_error(.mse_check_catch_mult(
    data.frame(Year = c(2022, 2022), Species = c(1, 1), mult = c(0.9, 0.5)),
    nspp = 3, 2021, 2030),
    "more than one multiplier")

  expect_error(.mse_check_catch_mult(data.frame(Year = 2022, Species = 1, mult = NA),
                                     nspp = 3, 2021, 2030),
               "finite and non-negative")

  expect_error(.mse_check_catch_mult(
    data.frame(Year = numeric(0), Species = numeric(0), mult = numeric(0)),
    nspp = 3, 2021, 2030),
    "no rows")
})
