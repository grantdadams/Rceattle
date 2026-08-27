# env_data is read as numbers everywhere it is used: it becomes env_index, the
# DATA_MATRIX handed to TMB, and it is the data behind every linkage's
# model.matrix(). A table that carries text (a workbook column read back as
# character, a factor, a character matrix) used to survive all the way to
# MakeADFun, where it failed with "Only numeric matrices ... can be
# interfaced" -- after the mean-fill of missing years had already given up on
# the column ("argument is not numeric or logical: returning NA"). These tests
# pin the coercion that now happens on the way in.

testthat::test_that("clean_data() coerces env_data to numeric with an integer Year", {
  dat <- make_test_data(nyrs = 5, nages = 4, seed = 1)
  env_num <- dat$env_data$Index1
  dat$env_data$Year   <- as.character(dat$env_data$Year)
  dat$env_data$Index1 <- as.character(env_num)

  out <- Rceattle::clean_data(dat)

  testthat::expect_true(is.integer(out$env_data$Year))
  testthat::expect_true(is.numeric(out$env_data$Index1))
  testthat::expect_equal(out$env_data$Index1, env_num, tolerance = 1e-12)
})


testthat::test_that("a factor env_data column reaches TMB as a number", {
  dat <- make_test_data(nyrs = 5, nages = 4, seed = 1)
  dat$env_data$Index1 <- factor(format(dat$env_data$Index1))

  out <- Rceattle::clean_data(dat)
  testthat::expect_true(is.numeric(out$env_data$Index1))
})


testthat::test_that("env_data supplied as a matrix becomes a numeric data.frame", {
  env <- cbind(Year = 1980:1984, Qcov = c(0.1, 0.2, 0.3, 0.4, 0.5))
  storage.mode(env) <- "character"

  out <- Rceattle:::.rce_numeric_env_data(env)

  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_true(is.integer(out$Year))
  testthat::expect_equal(out$Qcov, c(0.1, 0.2, 0.3, 0.4, 0.5))
})


testthat::test_that("a one-column matrix column is flattened, a wider one refused", {
  # What data.frame(Year = yrs, Qcov = ecov) builds when `ecov` is a matrix --
  # as.matrix() would later expand it into several env_index columns, shifting
  # every index a Cindex / M1_indices / Time_varying_q refers to.
  d <- data.frame(Year = 1980:1982)
  d$Qcov <- matrix(c(1, 2, 3), ncol = 1)
  out <- Rceattle:::.rce_numeric_env_data(d)
  testthat::expect_null(dim(out$Qcov))
  testthat::expect_equal(out$Qcov, c(1, 2, 3))

  d2 <- data.frame(Year = 1980:1982)
  d2$Qcov <- matrix(1:6, ncol = 2)
  testthat::expect_error(Rceattle:::.rce_numeric_env_data(d2),
                         "2-column matrix")
})


testthat::test_that("a genuinely non-numeric env_data cell is named, not silently NA", {
  # Silently coercing to NA would hand the mean-fill a hole to fill with the
  # column mean -- a fabricated observation for that year.
  d <- data.frame(Year = 1980:1982, Qcov = c("0.1", "O.2", "0.3"))
  testthat::expect_error(Rceattle:::.rce_numeric_env_data(d),
                         "Non-numeric value\\(s\\) 'O.2' in env_data column 'Qcov'")

  # Blank cells and the literal NA / NaN stay missing, which is what the
  # mean-fill is for.
  d2 <- data.frame(Year = 1980:1982, Qcov = c("0.1", NA, "NaN"))
  out <- Rceattle:::.rce_numeric_env_data(d2)
  testthat::expect_equal(out$Qcov, c(0.1, NA, NaN))
})


testthat::test_that("a fractional Year is an error rather than a rounded year", {
  d <- data.frame(Year = c(1980, 1980.5), Qcov = c(1, 2))
  testthat::expect_error(Rceattle:::.rce_numeric_env_data(d),
                         "non-integer value")
})


testthat::test_that("rearrange_data() builds a numeric env_index from a character column", {
  dat <- make_test_data(nyrs = 5, nages = 4, seed = 1)
  env_num <- dat$env_data$Index1
  dat$env_data$Index1 <- as.character(env_num)

  out <- Rceattle::rearrange_data(
    Rceattle::switch_check(Rceattle::clean_data(dat)))

  testthat::expect_identical(storage.mode(out$env_index), "double")
  testthat::expect_equal(nrow(out$env_index),
                         length(dat$styr:dat$projyr))
  # Observed years keep their values; the projection years take the mean.
  testthat::expect_equal(unname(out$env_index[seq_along(env_num), 1]), env_num,
                         tolerance = 1e-12)
  testthat::expect_equal(unname(out$env_index[nrow(out$env_index), 1]),
                         mean(env_num), tolerance = 1e-12)
})


testthat::test_that("an env_data column with no values at all warns", {
  # An empty column has no mean to fill missing years with, so every year is
  # NaN and any index pointed at it makes the objective NaN.
  dat <- make_test_data(nyrs = 5, nages = 4, seed = 1)
  dat$env_data$Empty <- NA

  testthat::expect_warning(
    Rceattle::rearrange_data(Rceattle::switch_check(Rceattle::clean_data(dat))),
    "'Empty'")
})
