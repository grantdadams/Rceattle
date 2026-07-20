# The linkage formula grammar: what is a fixed effect, what is a random
# effect, and what must be rejected.

.p <- function(f) Rceattle:::.parse_linkage_formula(f)

testthat::test_that("terms outside a bar are fixed effects", {
  for (f in list(~ temp, ~ temp + PDO, ~ Year, ~ factor(Year))) {
    r <- .p(f)
    testthat::expect_length(r$re_terms, 0L)
    testthat::expect_length(r$re_structures, 0L)
    testthat::expect_equal(deparse(r$fixed), deparse(f))
  }
})


testthat::test_that("`~ Year` is a linear trend, `~ factor(Year)` is deviates", {
  # These must stay distinct: model.matrix() gives one column for the trend
  # and one per level (less the reference) for the deviates. Collapsing them
  # would silently change the model.
  env <- data.frame(Year = 2000:2004)
  x_trend <- stats::model.matrix(.p(~ Year)$fixed, env)
  x_devs  <- stats::model.matrix(.p(~ factor(Year))$fixed, env)
  testthat::expect_equal(ncol(x_trend), 2L)          # intercept + slope
  testthat::expect_equal(ncol(x_devs),  length(unique(env$Year)))
})


testthat::test_that("bar terms are random effects with a structure", {
  r <- .p(~ (1 | Year))
  testthat::expect_length(r$re_terms, 1L)
  testthat::expect_equal(r$re_structures, "us")      # plain bar

  r <- .p(~ ar1(Year + 0 | fleet))
  testthat::expect_equal(r$re_structures, "ar1")

  r <- .p(~ rw(Year + 0 | fleet))
  testthat::expect_equal(r$re_structures, "rw")
})


testthat::test_that("fixed and random parts separate correctly", {
  r <- .p(~ temp + PDO + ar1(Year + 0 | fleet) + (1 | species))
  testthat::expect_equal(deparse(r$fixed), "~temp + PDO")
  testthat::expect_equal(r$re_structures, c("ar1", "us"))
})


testthat::test_that("an unrecognised structure is an error, not a silent default", {
  # reformulas::splitForm() maps any unknown wrapper to the default
  # structure without warning, so `bogus(Year + 0 | fleet)` would come back
  # as "us" and the user would get an unstructured term while believing
  # they had asked for something else. Rceattle checks first.
  testthat::expect_error(.p(~ bogus(Year + 0 | fleet)),
                         "unknown covariance structure")
  testthat::expect_error(.p(~ temp + notAThing(Year + 0 | fleet)),
                         "notAThing")
  # The message must name the alternatives.
  testthat::expect_error(.p(~ bogus(Year + 0 | fleet)), "Available:")
})


testthat::test_that("every advertised structure parses", {
  for (s in Rceattle:::LINKAGE_STRUCTURES) {
    f <- stats::as.formula(sprintf("~ %s(Year + 0 | fleet)", s))
    testthat::expect_equal(.p(f)$re_structures, s, info = s)
  }
})


testthat::test_that("two-sided and non-formula input are rejected", {
  testthat::expect_error(.p(y ~ temp), "one-sided")
  testthat::expect_error(.p("temp"),   "one-sided")
})


testthat::test_that("bar terms are rejected by materialize_linkage, not mangled", {
  # model.matrix() evaluates `1 | Year` as a logical OR and produces a
  # column literally named "1 | YearTRUE". Splitting the formula first
  # turns that silent corruption into an error.
  env <- data.frame(Year = 1:6, temp = stats::rnorm(6))
  spec <- Rceattle::linkage_spec(~ (1 | Year), param = "M1")
  testthat::expect_error(
    Rceattle:::materialize_linkage(spec, "M", env, list(species = 1L)),
    "random-effect terms are not yet supported")
})


testthat::test_that("fixed formulas still materialize, with contrast coding", {
  env <- data.frame(Year = 1:6, temp = stats::rnorm(6))

  tbl <- Rceattle:::materialize_linkage(
    Rceattle::linkage_spec(~ temp, param = "M1"), "M", env, list(species = 1L))
  testthat::expect_equal(tbl$design_col, c("(Intercept)", "temp"))

  # `~ Year` is one column (a linear trend); `~ factor(Year)` is one per
  # level less the reference -- which is what makes the deviate vector
  # identifiable without a sum-to-zero constraint.
  tr <- Rceattle:::materialize_linkage(
    Rceattle::linkage_spec(~ Year, param = "M1"), "M", env, list(species = 1L))
  dv <- Rceattle:::materialize_linkage(
    Rceattle::linkage_spec(~ factor(Year), param = "M1"), "M", env,
    list(species = 1L))
  testthat::expect_equal(nrow(tr), 2L)
  testthat::expect_equal(nrow(dv), length(unique(env$Year)))
})
