# run_mse() carries `weight` (weight-at-age, kg) and `ration_data` (annual
# foraging days) into the projection at their terminal HINDCAST year. Two rules
# govern which rows it may add:
#
#   * A series supplied at Year 0 is time-invariant -- rearrange_data() fills
#     every hindcast year from the one row -- so it takes no projection rows.
#   * A projection year that already has a row keeps it. rearrange_data()
#     assigns row by row into the weight array, so a second row for the same
#     year wins on sort order and would replace an observation.
#
# Provenance: expanding a Year 0 series replaced one legal row with a series
# naming only the projection years. data_check() reads an index carrying more
# than one year as time-varying and requires it to span styr..endyr, so the
# first assessment (which advances endyr) failed with
#
#   Weight data for index = 4 & sex = 1 does not span all hindcast years
#
# raised from data_check() inside the operating-model refit and recorded as a
# bare "OM" failure with the message discarded. Hit on a Pacific hake MSE whose
# workbook supplies weight indices 4-6 at Year 0.
#
# Tested against the internal helper rather than a run_mse() call: these are
# setup rules, and fitting two models to exercise them would take minutes.

wt_frame <- function() {
  # Indices 1-2 time-varying over 2018:2020; index 3 time-invariant (Year 0),
  # the shape rearrange_data() fills across the whole hindcast.
  rbind(
    data.frame(Wt_index = 1L, Sex = 1L, Year = 2018:2020,
               Age1 = c(0.1, 0.2, 0.3)),
    data.frame(Wt_index = 2L, Sex = 1L, Year = 2018:2020,
               Age1 = c(1.1, 1.2, 1.3)),
    data.frame(Wt_index = 3L, Sex = 1L, Year = 0L, Age1 = 9.9)
  )
}


test_that("a time-invariant series is left at its single Year 0 row", {
  out <- .mse_expand_fixed_input(wt_frame(), "Wt_index", 2020, 2021:2023)

  inv <- out[out$Wt_index == 3L, ]
  expect_equal(nrow(inv), 1L)
  expect_equal(inv$Year, 0)
  expect_equal(inv$Age1, 9.9)

  # This is the property data_check() reads: an index with one year is
  # time-invariant and exempt from the hindcast-span rule. More than one year
  # and it must cover styr..endyr, which a Year 0 row plus projection years
  # never does.
  expect_equal(length(unique(inv$Year)), 1L)
})


test_that("a time-varying series is extended from its terminal year", {
  out <- .mse_expand_fixed_input(wt_frame(), "Wt_index", 2020, 2021:2023)

  v <- out[out$Wt_index == 1L, ]
  expect_equal(v$Year, c(2018:2020, 2021:2023))
  # The projection holds the terminal hindcast year, which is what run_mse()
  # documents: weight-at-age is not forecast.
  expect_equal(v$Age1, c(0.1, 0.2, 0.3, 0.3, 0.3, 0.3))
})


test_that("a year that already has a row is not given a second", {
  # A workbook whose weight series runs past its own endyr. The observed 2021
  # and 2022 rows stay as they are; only 2023 is added. A duplicate would not
  # be inert -- rearrange_data() assigns row by row, so the appended copy sorts
  # last and its value would replace the observation.
  wt <- data.frame(Wt_index = 1L, Sex = 1L, Year = 2018:2022,
                   Age1 = c(0.1, 0.2, 0.3, 0.4, 0.5))
  out <- .mse_expand_fixed_input(wt, "Wt_index", 2020, 2021:2023)

  expect_equal(out$Year, 2018:2023)
  expect_false(anyDuplicated(out$Year) > 0)
  # 2021 and 2022 keep their own values; 2023 carries the terminal HINDCAST
  # year (2020 = 0.3), not the last row in the table (2022 = 0.5).
  expect_equal(out$Age1, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.3))
})


test_that("the carried-forward value comes from endyr, not the last row", {
  # run_mse() documents holding weight-at-age at the operating model's terminal
  # hindcast year. Reading the table's last row instead would project a value
  # the hindcast never fitted.
  wt <- data.frame(Wt_index = 1L, Sex = 1L, Year = 2018:2021,
                   Age1 = c(0.1, 0.2, 0.3, 0.9))
  out <- .mse_expand_fixed_input(wt, "Wt_index", 2019, 2020:2022)

  expect_equal(out$Age1[out$Year == 2022], 0.2)   # 2019's value
  expect_equal(out$Age1[out$Year == 2021], 0.9)   # its own row, untouched
})


test_that("each sex is extended on its own", {
  # A two-sex series is two independent groups; expanding on the id alone would
  # carry one sex's terminal weight into the other.
  wt <- rbind(
    data.frame(Wt_index = 1L, Sex = 1L, Year = 2019:2020, Age1 = c(0.2, 0.3)),
    data.frame(Wt_index = 1L, Sex = 2L, Year = 2019:2020, Age1 = c(0.6, 0.7))
  )
  out <- .mse_expand_fixed_input(wt, "Wt_index", 2020, 2021:2022)

  expect_equal(out$Age1[out$Sex == 1L], c(0.2, 0.3, 0.3, 0.3))
  expect_equal(out$Age1[out$Sex == 2L], c(0.6, 0.7, 0.7, 0.7))
})


test_that("ration_data is expanded on Species, the column it is keyed by", {
  ration <- rbind(
    data.frame(Species = 1L, Sex = 1L, Year = 2019:2020, Age1 = c(200, 210)),
    data.frame(Species = 2L, Sex = 1L, Year = 0L, Age1 = 180)
  )
  out <- .mse_expand_fixed_input(ration, "Species", 2020, 2021:2022)

  expect_equal(out$Year[out$Species == 1L], c(2019, 2020, 2021, 2022))
  expect_equal(nrow(out[out$Species == 2L, ]), 1L)
})


test_that("a wholly time-invariant table comes back untouched", {
  wt <- data.frame(Wt_index = 1:3, Sex = 1L, Year = 0L, Age1 = c(1, 2, 3))
  expect_identical(.mse_expand_fixed_input(wt, "Wt_index", 2020, 2021:2023), wt)

  # And an absent or empty table is returned as given, so the setup cannot fail
  # on a model that carries no ration data.
  empty <- wt[0, ]
  expect_identical(.mse_expand_fixed_input(empty, "Wt_index", 2020, 2021:2023),
                   empty)
  expect_null(.mse_expand_fixed_input(NULL, "Wt_index", 2020, 2021:2023))
})


test_that("rows come back ordered by series then year", {
  # rearrange_data() reads these tables row by row, and the ordering is what
  # both call sites produced before.
  out <- .mse_expand_fixed_input(wt_frame(), "Wt_index", 2020, 2021:2023)
  expect_false(is.unsorted(out$Wt_index))
  for (i in unique(out$Wt_index)) {
    expect_false(is.unsorted(out$Year[out$Wt_index == i]))
  }
  expect_equal(rownames(out), as.character(seq_len(nrow(out))))
})
