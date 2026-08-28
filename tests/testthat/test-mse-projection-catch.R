# run_mse() fills a projection year's catch by selecting rows on is.na(Catch),
# so a projection row arriving with a number in it is skipped entirely: no catch
# is written, the assessment interval's total comes out 0, and the operating
# model takes the no-catch path -- rebuilt at estimateMode 3 with its map
# dropped, rather than advanced with the control rule's advice.
#
# Provenance: a Pacific hake MSE workbook (MSE_hake_yr24_final.xlsx) carried
# catch through 2023 with endyr still 2019. clean_data() only creates an NA
# projection row for a year with no row already, so 2020-2023 kept their
# recorded catch and only 2024-2025 were fillable; four of the six assessments
# never advanced the operating model. Nothing in the run said so.
#
# Tested against the internal helper rather than a run_mse() call: this is a
# setup rule, and fitting an operating and an estimation model to exercise it
# would take minutes.

test_that("catch past the terminal year is blanked, and the years are named", {
  catch <- data.frame(Fleet_code = 1L,
                      Year  = 2017:2025,
                      Catch = c(410, 412, 413, 380, 327, 323, 264, NA, NA))

  expect_warning(out <- .mse_blank_proj_catch(catch, endyr = 2019, "operating"),
                 "records catch in 2020, 2021, 2022, 2023")

  # The hindcast is untouched; every projection year is NA, which is the state
  # clean_data() creates a projection row in and the state the fill selects on.
  expect_equal(out$Catch[out$Year <= 2019], c(410, 412, 413))
  expect_true(all(is.na(out$Catch[out$Year > 2019])))
  expect_equal(nrow(out), nrow(catch))
})


test_that("the warning names the model and the terminal year it read", {
  # Either model can be the one carrying the stale years, and they can end in
  # different years, so both are named rather than reported as one condition.
  catch <- data.frame(Fleet_code = 1L, Year = 2019:2021, Catch = c(413, 380, 327))

  expect_warning(.mse_blank_proj_catch(catch, endyr = 2019, "operating"),
                 paste0("^The operating model records catch in 2020, 2021, ",
                        "past its terminal year \\(2019\\)"))
  expect_warning(.mse_blank_proj_catch(catch, endyr = 2020, "estimation"),
                 paste0("^The estimation model records catch in 2021, ",
                        "past its terminal year \\(2020\\)"))
})


test_that("a model whose projection rows are already NA is left alone", {
  # The ordinary case: clean_data() has created the projection rows and nothing
  # has filled them yet. No warning, and the frame comes back identical.
  catch <- data.frame(Fleet_code = 1L,
                      Year  = 2017:2021,
                      Catch = c(410, 412, 413, NA, NA))

  expect_silent(out <- .mse_blank_proj_catch(catch, endyr = 2019, "operating"))
  expect_identical(out, catch)

  # A model with no projection at all (endyr == projyr) is the same quiet case.
  hind <- catch[catch$Year <= 2019, ]
  expect_silent(out2 <- .mse_blank_proj_catch(hind, endyr = 2019, "operating"))
  expect_identical(out2, hind)
})


test_that("an empty or absent catch table is returned unchanged", {
  # Guarded so the setup cannot fail on a model that carries no catch table;
  # nrow(NULL) is NULL, which would error in the `if` rather than skip.
  empty <- data.frame(Fleet_code = integer(), Year = integer(),
                      Catch = numeric())
  expect_silent(out <- .mse_blank_proj_catch(empty, endyr = 2019, "operating"))
  expect_identical(out, empty)
  expect_null(.mse_blank_proj_catch(NULL, endyr = 2019, "operating"))
})


test_that("the MSE refuses models whose terminal years disagree", {
  # The assessment loop reads the row positions to fill off the ESTIMATION
  # model's catch table and indexes the OPERATING model's max_catch_hat with
  # them, so the two tables must describe the same rows in the same order. A
  # mismatch reads the exploitable-biomass limit off the wrong years or off NA,
  # which is catch advice that is wrong without looking wrong.
  # The check runs in run_mse()'s input block, before anything reads a fitted
  # quantity, so a minimal object carrying the two years is enough to reach it.
  mk <- function(endyr, projyr) {
    structure(list(data_list = list(endyr = endyr, projyr = projyr)),
              class = "Rceattle")
  }
  expect_error(run_mse(om = mk(2017, 2050), em = mk(2016, 2050)),
               "must share a terminal year and a projection horizon")
  expect_error(run_mse(om = mk(2017, 2050), em = mk(2017, 2040)),
               "projyr 2050 / 2040")
})


test_that("blanking makes every assessment year fillable", {
  # The property that matters downstream: after the blanking, every year in the
  # assessment schedule has a row the fill can select. This is the check that
  # would have caught the hake workbook.
  catch <- data.frame(Fleet_code = 1L,
                      Year  = 2019:2025,
                      Catch = c(413, 380, 327, 323, 264, NA, NA))
  endyr <- 2019

  fillable <- function(cd) {
    vapply(.mse_assess_years(1, om_endyr = endyr, max_yr = 2025),
           function(y) any(cd$Year == y & is.na(cd$Catch)), logical(1))
  }

  expect_equal(fillable(catch), c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE))

  out <- suppressWarnings(.mse_blank_proj_catch(catch, endyr, "operating"))
  expect_true(all(fillable(out)))
})
