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


test_that("the terminal year is found by year, not by table position", {
  # Nothing sorts `weight` on read, so a group's last ROW need not be its last
  # YEAR. Taking it by position carried 2018's weight-at-age through the whole
  # projection while 2019 sat above it in the table.
  wt <- data.frame(Wt_index = 1L, Sex = 1L, Year = c(2019L, 2017L, 2018L),
                   Age1 = c(0.9, 0.1, 0.2))
  out <- .mse_expand_fixed_input(wt, "Wt_index", 2019, 2020:2021)

  expect_equal(unique(out$Age1[out$Year %in% 2020:2021]), 0.9)

  # And the fallback, for a series that ends before the terminal year.
  old <- data.frame(Wt_index = 1L, Sex = 1L, Year = c(2016L, 2014L, 2015L),
                    Age1 = c(0.8, 0.1, 0.2))
  out2 <- .mse_expand_fixed_input(old, "Wt_index", 2019, 2020:2021)
  expect_equal(unique(out2$Age1[out2$Year %in% 2020:2021]), 0.8)
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


# `NByageFixed` is NOT one of these. Numbers-at-age for a species with
# `estDynamics > 0` are an input, so the projection years cannot be supplied by
# the MSE at all -- neither by holding the terminal year, which invents data and
# pins the population flat, nor by leaving them empty, which makes
# `rearrange_data()`'s array read 0 and the template multiply a zero population
# into `N_at_age`. `.mse_require_fixed_natage()` warns instead, naming the
# species and the uncovered years and saying what the consequence is, since the
# MSE still runs.
#
# Provenance: run_mse() used to expand it, with its own inline copy of the
# operation above, and that copy also dropped groups. It sliced one block of
# projection rows per (Species, Sex) group, then subset by
# `which(om_proj_yrs %!in% NByageFixed$Year)` -- positions within `om_proj_yrs`,
# not row positions in the grouped frame -- so only the FIRST group's block
# survived. Measured on `GOA2018SS` (4 groups, 32 projection years): 128 rows
# built, 32 kept, leaving species 2 (both sexes) and species 3 at the array's
# initialised zeros for every projection year. Both the dropping and the
# expansion itself are gone.

test_that("a fixed-numbers species whose input stops short is warned about", {
  nb <- rbind(
    data.frame(Species = 1L, Sex = 0L, Year = 2019:2020, Age1 = c(11, 12)),
    data.frame(Species = 2L, Sex = 1L, Year = 2019:2020, Age1 = c(21, 22))
  )
  d <- list(estDynamics = c(1, 1), NByageFixed = nb)
  # The MSE still runs, so the message has to state the CONSEQUENCE, not just
  # the gap: those years get numbers-at-age of 0.
  expect_warning(.mse_require_fixed_natage(d, 2021:2023),
                 "does not cover every projection year")
  expect_warning(.mse_require_fixed_natage(d, 2021:2023),
                 "numbers-at-age.*of 0")
  # And which species and years, so a multi-species MSE says where to look.
  expect_warning(.mse_require_fixed_natage(d, 2021:2023), "species 1, 2")
  expect_warning(.mse_require_fixed_natage(d, 2021:2023), "2021-2023")
  # Returns the short species invisibly, so a caller can act on it.
  expect_equal(suppressWarnings(.mse_require_fixed_natage(d, 2021:2023)),
               c(1L, 2L))
})


test_that("only the uncovered years are counted, not the whole horizon", {
  # A species covered for part of the projection is short by the remainder, and
  # the message must say so rather than implying the whole horizon is missing.
  # Covered for 2021-2022, short by 2023-2024.
  nb <- data.frame(Species = 1L, Sex = 0L, Year = 2019:2022, Age1 = 1:4)
  d <- list(estDynamics = 1, NByageFixed = nb)
  expect_warning(.mse_require_fixed_natage(d, 2021:2024),
                 "2 year\\(s\\) missing, 2023-2024")
})


test_that("an input that spans the horizon passes, and sex is not the key", {
  # Coverage is per SPECIES: rearrange_data() writes a Sex == 0 row to every
  # sex of a two-sex species, so requiring a row per sex would refuse a legal
  # table.
  nb <- rbind(
    data.frame(Species = 1L, Sex = 0L, Year = 2019:2023, Age1 = 1:5),
    data.frame(Species = 2L, Sex = 1L, Year = 2019:2023, Age1 = 6:10)
  )
  d <- list(estDynamics = c(1, 1), NByageFixed = nb)
  expect_silent(.mse_require_fixed_natage(d, 2021:2023))
  expect_null(.mse_require_fixed_natage(d, 2021:2023))
})


test_that("a species that estimates its own numbers is not asked for rows", {
  # estDynamics == 0 is the default and the usual case: the check must be inert
  # there, including when no NByageFixed table exists at all.
  d0 <- list(estDynamics = c(0, 0), NByageFixed = data.frame())
  expect_silent(.mse_require_fixed_natage(d0, 2021:2023))
  expect_silent(.mse_require_fixed_natage(list(estDynamics = NULL), 2021:2023))

  # And only the fixed-numbers species is required to cover the horizon.
  nb <- data.frame(Species = 2L, Sex = 1L, Year = 2019:2023, Age1 = 1:5)
  d <- list(estDynamics = c(0, 1), NByageFixed = nb)
  expect_silent(.mse_require_fixed_natage(d, 2021:2023))
})


test_that("a fixed-numbers species with no table at all is warned about", {
  # data_check() requires presence, but run_mse() can be handed a list directly.
  d <- list(estDynamics = c(1, 0), NByageFixed = data.frame())
  expect_warning(.mse_require_fixed_natage(d, 2021:2023), "species 1")
})
