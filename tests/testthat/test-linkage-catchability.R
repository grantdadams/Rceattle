# Environmental linkages on survey/index catchability.

.q_env_data <- function(d) {
  yrs <- d$styr:d$projyr
  set.seed(1)
  data.frame(Year = yrs, temp = as.numeric(scale(sin(seq_along(yrs) / 4))))
}

testthat::test_that("q is a wired linkage process", {
  testthat::expect_true("q" %in% Rceattle:::LINKAGE_PROCESSES_IMPLEMENTED)
})


testthat::test_that("build_catchability() validates its linkages", {
  testthat::expect_null(Rceattle::build_catchability()$linkages)

  ok <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ temp, by = ~ fleet)))
  testthat::expect_s3_class(ok$linkages$q, "Rceattle_linkage_spec")

  # Only `q` is a catchability parameter.
  testthat::expect_error(
    Rceattle::build_catchability(linkages = list(
      nope = Rceattle::linkage_spec(~ temp))),
    "unknown")
})


testthat::test_that("a q linkage produces one coefficient set per fleet", {
  dat <- Rceattle::switch_check(Rceattle::clean_data(Rceattle::BS2017SS))
  env <- .q_env_data(dat)
  spec <- Rceattle::linkage_spec(~ temp, param = "q", by = ~ fleet)
  n_flt <- nrow(dat$fleet_control)

  tbl <- Rceattle:::materialize_linkage(
    spec, "q", env, list(fleet = seq_len(n_flt)))

  # intercept + slope, per fleet
  testthat::expect_equal(nrow(tbl), 2L * n_flt)
  testthat::expect_true(all(tbl$process == "q"))
  testthat::expect_setequal(unique(tbl$fleet), seq_len(n_flt))
  # q is not stratified by species/sex/age
  testthat::expect_true(all(is.na(tbl$species)))
  testthat::expect_true(all(is.na(tbl$age_bin)))
})


testthat::test_that("the fleet filter restricts which fleets get coefficients", {
  dat <- Rceattle::switch_check(Rceattle::clean_data(Rceattle::BS2017SS))
  env <- .q_env_data(dat)
  n_flt <- nrow(dat$fleet_control)
  testthat::skip_if(n_flt < 3)

  tbl <- Rceattle:::materialize_linkage(
    Rceattle::linkage_spec(~ temp, param = "q", by = ~ fleet, fleet = c(1L, 3L)),
    "q", env, list(fleet = seq_len(n_flt)))

  testthat::expect_setequal(unique(tbl$fleet), c(1L, 3L))
})


testthat::test_that("a model without a q linkage is numerically unchanged", {
  testthat::skip_on_cran()
  fit <- Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, estimateMode = 1, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))
  # The accumulator runs unconditionally; with no q rows its offsets are
  # zero and the fit must match the golden reference exactly.
  testthat::expect_equal(fit$opt$objective, 10241.0304274978,
                         tolerance = 1e-8)
  testthat::expect_true(all(fit$quantities$q_linkage_offset == 0))
  testthat::expect_true(all(fit$quantities$q_linkage_offset_nat == 0))
})


testthat::test_that("a q linkage moves catchability and is estimated", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  d$env_data <- .q_env_data(d)

  fit <- Rceattle::fit_mod(
    data_list = d, estimateMode = 1, msmMode = 0,
    qFun = Rceattle::build_catchability(linkages = list(
      q = Rceattle::linkage_spec(~ temp, by = ~ fleet))),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))

  testthat::expect_true(any(fit$quantities$q_linkage_offset != 0))
  # The slope coefficients are free parameters, so the fit cannot be worse
  # than the no-linkage one.
  testthat::expect_lt(fit$opt$objective, 10241.0304274978)
})


testthat::test_that("Catchability = 'Environmental' warns but still fits", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  d$env_data <- .q_env_data(d)
  d$fleet_control$Catchability[1]   <- "Environmental"
  d$fleet_control$Time_varying_q[1] <- "1"

  warns <- character(0)
  fit <- withCallingHandlers(
    Rceattle::fit_mod(
      data_list = d, estimateMode = 1, msmMode = 0,
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                          verbose = 0)),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
    })

  # The deprecation is a redirect, not a removal: the legacy path still runs.
  testthat::expect_true(any(grepl("deprecated", warns)))
  testthat::expect_true(any(grepl("build_catchability", warns)))
  testthat::expect_true(is.finite(fit$opt$objective))
})


testthat::test_that("a q linkage intercept re-targets the base catchability", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  d$env_data <- .q_env_data(d)

  fit <- Rceattle::fit_mod(
    data_list = d, estimateMode = 1, msmMode = 0,
    qFun = Rceattle::build_catchability(linkages = list(
      q = Rceattle::linkage_spec(~ temp, by = ~ fleet))),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))

  tbl <- fit$data_list$linkage_table
  # An intercept-bearing formula fixes the intercept rows at 0 and leaves the
  # base index_log_q estimable to carry the level -- as growth/M/recruitment
  # linkages do with their base parameters.
  testthat::expect_true(
    all(is.na(fit$map$mapList$beta_linkage[tbl$design_col == "(Intercept)"])))
  testthat::expect_true(any(!is.na(fit$map$mapList$index_log_q)))
})
