# env_data rows are matched to model years by POSITION -- row r is model year
# styr + r - 1 -- so a row BEFORE styr shifts every later row and feeds the wrong
# covariate to consumption and, under multispecies, to predation mortality.
# rearrange_data() drops such rows from env_index, so a model without a linkage
# has always tolerated them; the linkage path refused the whole fit instead.
#
# Rows AFTER the model years are deliberately left alone. Alignment runs from the
# front, so they are inert -- but the fixed part of a linkage formula goes to
# model.matrix(), and cut() / poly() / scale() are computed on whatever rows
# env_data supplies, so dropping them moves the design matrix and the objective.

.env_tbl <- function(first, last) {
  data.frame(Year = first:last,
             temp = as.numeric(seq_along(first:last)))   # value == position
}

testthat::test_that("pre-styr rows are dropped, later rows keep their own year", {
  ed  <- .env_tbl(1975, 1990)
  out <- suppressWarnings(Rceattle:::.trim_env_data(ed, styr = 1980))

  testthat::expect_equal(out$Year[1], 1980)
  testthat::expect_equal(nrow(out), 11L)
  # 1980 keeps the value it had before the trim rather than inheriting 1975's.
  testthat::expect_equal(out$temp, ed$temp[ed$Year >= 1980])
})


testthat::test_that("rows after the model years are NOT dropped", {
  # Dropping them would rebuild a cut()/poly()/scale() design matrix on a
  # different row set and move the objective of a fit that was already correct.
  ed <- .env_tbl(1980, 2050)
  testthat::expect_silent(out <- Rceattle:::.trim_env_data(ed, styr = 1980))
  testthat::expect_identical(out, ed)
})


testthat::test_that("the drop is a warning, naming the count and styr", {
  testthat::expect_warning(
    Rceattle:::.trim_env_data(.env_tbl(1975, 1990), 1980),
    "dropped 5 row\\(s\\) before styr \\(1980\\)")
})


testthat::test_that("an NA year is kept, for the year check to reject", {
  # An NA year is unlabelled, not early. 5.27.0 rejected it by name; silently
  # dropping it here would change which covariates reach the model.
  ed  <- data.frame(Year = c(1980:1990, NA), temp = 1:12)
  out <- Rceattle:::.trim_env_data(ed, 1980)
  testthat::expect_equal(nrow(out), 12L)
  testthat::expect_error(Rceattle:::.check_env_data_years(out, 1980),
                         "sorted ascending with no duplicates or NA")
})


testthat::test_that("a table entirely before styr is an error, not an empty drop", {
  testthat::expect_error(
    Rceattle:::.trim_env_data(.env_tbl(1960, 1970), 1980),
    "env_data ends at 1970, before the model start year styr \\(1980\\)")
})


testthat::test_that("the trim leaves the extension and year check satisfiable", {
  out <- suppressWarnings(Rceattle:::.trim_env_data(.env_tbl(1975, 1985), 1980))
  out <- suppressMessages(Rceattle:::.extend_env_data(out, 1980))
  testthat::expect_equal(out$Year[1], 1980)
  testthat::expect_silent(Rceattle:::.check_env_data_years(out, 1980))
})


testthat::test_that("an empty or absent table is a no-op, not an error", {
  testthat::expect_null(Rceattle:::.trim_env_data(NULL, 1980))
  ed0 <- data.frame(Year = integer(0), temp = numeric(0))
  testthat::expect_equal(nrow(Rceattle:::.trim_env_data(ed0, 1980)), 0L)
  testthat::expect_equal(nrow(Rceattle:::.extend_env_data(ed0, 1980)), 0L)
  testthat::expect_silent(Rceattle:::.check_env_data_years(ed0, 1980))
})


testthat::test_that("a fit with a linkage and a pre-styr row runs, aligned", {
  testthat::skip_on_cran()
  # The case that refused the whole fit. Asserted on env_index -- the matrix the
  # model actually reads, with Year dropped and rows taken by position -- not on
  # env_data, which is the table fit_mod() just rewrote.
  d <- make_test_data(nyrs = 12, nprojyrs = 2, nages = 5)
  d$growth_model <- rep(0, d$nspp)
  srv <- which(as.character(d$fleet_control$Fleet_type) == "Survey")[1]
  testthat::skip_if(is.na(srv), "fixture has no survey fleet")
  d$fleet_control$Catchability[srv] <- "Estimated"
  yrs <- d$styr:d$projyr
  # temp == the model year, so a shift of one row is visible as an off-by-one.
  d$env_data <- data.frame(Year = c(d$styr - 3, d$styr - 1, yrs),
                           temp = as.numeric(c(d$styr - 3, d$styr - 1, yrs)))

  fit <- suppressWarnings(suppressMessages(fit_mod(
    data_list  = d, msmMode = 0, estimateMode = "DebugBuild",
    qFun       = build_catchability(linkages = list(
      q = linkage_spec(~ temp, by = ~ fleet, fleet = srv))),
    fit_control = fit_control(verbose = 0))))

  testthat::expect_s3_class(fit, "Rceattle")
  ei <- Rceattle:::rearrange_data(Rceattle::switch_check(fit$data_list))$env_index
  testthat::expect_equal(nrow(ei), length(yrs))
  # Row i of env_index is model year styr + i - 1, and temp was set to the year.
  testthat::expect_equal(as.numeric(ei[, "temp"]), as.numeric(yrs))
})
