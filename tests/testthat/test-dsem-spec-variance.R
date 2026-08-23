# check_dsem_spec()'s covariate_variance record.
#
# A variable whose lag-0 two-headed self-path is pinned at zero is DETERMINISTIC
# under dsem's projecting parameterizations: the model computes it from its
# incoming paths, so the env_data supplied for it is never read, it drops out of
# the precision sample_rec() draws from, and the recruitment bias correction
# cannot condition on it. Nothing errors -- the fit just quietly ignores the
# covariate -- so this check is the only thing that says so.
#
# The trigger is dsem's own rule (parameter and start both zero), NOT "has no
# two-headed path": make_dsem_ram() adds an ESTIMATED V[x] row for any variable
# that lacks one, so omitting the line is the ordinary way to ask for a free
# variance. Getting that backwards would flag most working sems.

.spec_of <- function(sem, vars = "BT") {
  d <- list(styr = 1990L, endyr = 2010L, projyr = 2015L, nspp = 1L,
            spnames = "a", sigma_rec = 1, random_rec = TRUE,
            proj_mean_rec = TRUE)
  ed <- data.frame(Year = 1990:2015)
  for (v in vars) ed[[v]] <- as.numeric(seq_along(ed$Year))
  d$env_data <- ed
  # Deliberately not clean_data(): build_dsem_objects() and check_dsem_spec()
  # read only the year span, nspp and env_data, so this stays a spec test rather
  # than needing a whole assessment's worth of inputs.
  built <- suppressWarnings(Rceattle::build_dsem_objects(
    Rceattle::build_DSEM(sem = sem, family = "fixed"), data_list = d))
  Rceattle::check_dsem_spec(d, built)
}

.flagged <- function(res) {
  v <- res$checks$covariate_variance
  if (is.null(v)) character(0) else as.character(v$data$variable)
}


testthat::test_that("a variance pinned at zero is reported", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("dsem")

  res <- .spec_of("recdevs1 <-> recdevs1, 0, sigmaR1, 0.6
BT <-> BT, 0, NA, 0
BT -> recdevs1, 0, bBT, 0.45")
  testthat::expect_equal(.flagged(res), "BT")
  testthat::expect_equal(res$status, "WARN")
})


testthat::test_that("an omitted variance line is not reported", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("dsem")

  # make_dsem_ram() supplies an estimated V[BT], so this sem is healthy. Getting
  # this wrong would warn on the shipped vignette, whose sem omits the line.
  res <- .spec_of("recdevs1 <-> recdevs1, 0, sigmaR1, 0.6
BT -> recdevs1, 0, bBT, 0.45")
  testthat::expect_equal(.flagged(res), character(0))

  # A variance FIXED at a non-zero value is also a real variance.
  res2 <- .spec_of("recdevs1 <-> recdevs1, 0, sigmaR1, 0.6
BT <-> BT, 0, NA, 1
BT -> recdevs1, 0, bBT, 0.45")
  testthat::expect_equal(.flagged(res2), character(0))
})


testthat::test_that("it reaches a variable that is not a direct predictor", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("dsem")

  # PDO drives BT drives recruitment. PDO is two steps upstream, so it is not in
  # the recruitment design the other Tier-1 checks look at -- but it removes the
  # same rows from the precision.
  res <- .spec_of("recdevs1 <-> recdevs1, 0, sigmaR1, 0.6
BT <-> BT, 0, sdBT, 1
PDO <-> PDO, 0, NA, 0
PDO -> BT, 0, bPB, 0.3
BT -> recdevs1, 0, bBT, 0.45", vars = c("BT", "PDO"))
  testthat::expect_equal(.flagged(res), "PDO")
})


testthat::test_that("the deterministic column really is ignored by the fit", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("dsem")
  testthat::skip_if_not_installed("TMB")

  # The claim the warning rests on: with the variance pinned at zero, scaling the
  # covariate leaves the objective bit-identical, so the series is not read.
  sem <- "recdevs1 <-> recdevs1, 0, sigmaR1, 0.6
recdevs2 <-> recdevs2, 0, sigmaR2, 0.6
recdevs3 <-> recdevs3, 0, sigmaR3, 0.6
BT <-> BT, 0, NA, 0
BT -> recdevs1, 0, bBT, 0.45"
  obj_at <- function(mult) {
    d <- Rceattle::BS2017SS
    yrs <- d$styr:d$projyr
    d$env_data <- data.frame(Year = yrs,
                             BT = mult * as.numeric(scale(sin(seq_along(yrs) / 3))))
    f <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
      data_list = d, inits = NULL, file = NULL, estimateMode = 3,
      random_rec = TRUE, msmMode = 0,
      dsem = Rceattle::build_DSEM(sem = sem, family = "fixed"),
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                          verbose = 0))))
    as.numeric(f$obj$fn())
  }
  testthat::expect_identical(obj_at(1), obj_at(100))
})
