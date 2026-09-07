# env_data rows are matched to model years by POSITION -- row r is model year
# styr + r - 1 -- so a year outside styr:projyr shifts every later row and feeds
# the wrong covariate to consumption and, under multispecies, to predation
# mortality. rearrange_data() drops such rows from env_index, so a model without
# a linkage has always tolerated them; the linkage path refused the whole fit
# instead. Both now drop them, so a workbook reads the same way either way.

.env_win <- function(styr = 1980, projyr = 1990, first = 1975, last = 1995) {
  data.frame(Year = first:last,
             temp = as.numeric(seq_along(first:last)))   # value == position, to spot a shift
}

testthat::test_that("rows outside the model years are dropped, not refused", {
  ed <- .env_win()
  out <- suppressMessages(Rceattle:::.trim_env_data(ed, styr = 1980, projyr = 1990))

  testthat::expect_equal(range(out$Year), c(1980, 1990))
  testthat::expect_equal(nrow(out), 11L)
  # The covariate still belongs to its own year: 1980 keeps the value it had
  # before the trim, rather than inheriting the first row's.
  testthat::expect_equal(out$temp, ed$temp[ed$Year %in% 1980:1990])
})


testthat::test_that("the drop is announced, with the count and the window", {
  testthat::expect_message(
    Rceattle:::.trim_env_data(.env_win(), 1980, 1990),
    "dropped 10 row\\(s\\) outside the model years 1980-1990")
})


testthat::test_that("a table already inside the window is untouched and silent", {
  ed <- data.frame(Year = 1980:1990, temp = 1:11)
  testthat::expect_silent(out <- Rceattle:::.trim_env_data(ed, 1980, 1990))
  testthat::expect_identical(out, ed)
})


testthat::test_that("no usable row is an error naming both ranges", {
  # Dropping everything would leave the linkage with no covariate at all, which
  # is a mis-specified model rather than a row to discard quietly.
  testthat::expect_error(
    Rceattle:::.trim_env_data(data.frame(Year = 1960:1970, temp = 1:11), 1980, 1990),
    "no row inside the model years 1980-1990; it covers 1960-1970")
})


testthat::test_that("the trim leaves the year check and the extension satisfiable", {
  # .extend_env_data() bails when the table starts before styr, and
  # .check_env_data_years() then errors -- which is how a pre-styr row refused
  # the fit. After the trim the pair completes.
  ed  <- .env_win(first = 1975, last = 1985)
  out <- suppressMessages(Rceattle:::.trim_env_data(ed, 1980, 1990))
  out <- suppressMessages(Rceattle:::.extend_env_data(out, 1980))
  testthat::expect_equal(out$Year[1], 1980)
  testthat::expect_silent(Rceattle:::.check_env_data_years(out, 1980))
})


testthat::test_that("an empty or absent table is a no-op, not an error", {
  testthat::expect_null(Rceattle:::.trim_env_data(NULL, 1980, 1990))
  ed0 <- data.frame(Year = integer(0), temp = numeric(0))
  testthat::expect_equal(nrow(Rceattle:::.trim_env_data(ed0, 1980, 1990)), 0L)
  # The pair downstream must survive it too.
  testthat::expect_equal(nrow(Rceattle:::.extend_env_data(ed0, 1980)), 0L)
  testthat::expect_silent(Rceattle:::.check_env_data_years(ed0, 1980))
})


testthat::test_that("a fit with a linkage and a pre-styr covariate row runs", {
  testthat::skip_on_cran()
  # The case that refused the whole fit: 15 workbooks in the sibling
  # repositories carry pre-styr environmental rows.
  d <- make_test_data(nyrs = 12, nprojyrs = 2, nages = 5)
  d$growth_model <- rep(0, d$nspp)
  srv <- which(as.character(d$fleet_control$Fleet_type) == "Survey")[1]
  testthat::skip_if(is.na(srv), "fixture has no survey fleet")
  d$fleet_control$Catchability[srv] <- "Estimated"   # a linkage needs an estimated q
  yrs <- d$styr:d$projyr
  d$env_data <- data.frame(Year = c(d$styr - 3, d$styr - 1, yrs),
                           temp = as.numeric(seq_len(length(yrs) + 2)))

  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list  = d, msmMode = 0, estimateMode = "DebugBuild",
    qFun       = build_catchability(linkages = list(
      q = linkage_spec(~ temp, by = ~ fleet, fleet = srv))),
    fit_control = fit_control(verbose = 0))))

  testthat::expect_s3_class(fit, "Rceattle")
  # The covariate kept its own years rather than being shifted by the two
  # dropped rows: temp for styr is the value the table gave styr.
  testthat::expect_equal(fit$data_list$env_data$Year[1], d$styr)
  testthat::expect_equal(fit$data_list$env_data$temp[1], 3)
})
