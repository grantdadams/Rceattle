# Environmental linkages on selectivity parameters (logistic family).

.sel_env_data <- function(d) {
  yrs <- d$styr:d$projyr
  set.seed(1)
  data.frame(Year = yrs, temp = as.numeric(scale(sin(seq_along(yrs) / 4))))
}

.logistic_fleet <- function(d) {
  fc <- Rceattle::switch_check(Rceattle::clean_data(d))$fleet_control
  which(fc$Selectivity == "Logistic")[1]
}

testthat::test_that("sel is a wired linkage process", {
  testthat::expect_true("sel" %in% Rceattle:::LINKAGE_PROCESSES_IMPLEMENTED)
})


testthat::test_that("build_selectivity() validates its linkages", {
  testthat::expect_null(Rceattle::build_selectivity()$linkages)
  ok <- Rceattle::build_selectivity(linkages = list(
    inf_asc = Rceattle::linkage_spec(~ temp, by = ~ fleet)))
  testthat::expect_s3_class(ok$linkages$inf_asc, "Rceattle_linkage_spec")
  testthat::expect_error(
    Rceattle::build_selectivity(linkages = list(
      nope = Rceattle::linkage_spec(~ temp))),
    "unknown")
})


testthat::test_that("the logistic slots map to the right parameter arrays", {
  # slope slots on log_sel_slp, inflection slots on sel_inf.
  m <- Rceattle:::.SEL_PARAM_TO_SLOT
  testthat::expect_equal(m$slp_asc$arr,  "log_sel_slp")
  testthat::expect_equal(m$inf_desc$arr, "sel_inf")
  testthat::expect_equal(m$inf_desc$slot, 2L)
  # DoubleNormal aliases share the slots.
  testthat::expect_equal(m$sigma_asc$slot, m$slp_asc$slot)
  testthat::expect_equal(m$peak$slot,      m$inf_asc$slot)
})


testthat::test_that("a model without a sel linkage is numerically unchanged", {
  testthat::skip_on_cran()
  fit <- Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, estimateMode = 1, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))
  testthat::expect_equal(fit$opt$objective, 10241.0304274978,
                         tolerance = 1e-8)
})


testthat::test_that("a sel linkage on a logistic fleet is estimated and moves the fit", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  d$env_data <- .sel_env_data(d)
  flt <- .logistic_fleet(d)
  testthat::skip_if(is.na(flt))

  fit <- Rceattle::fit_mod(
    data_list = d, estimateMode = 1, msmMode = 0,
    selFun = Rceattle::build_selectivity(linkages = list(
      inf_asc = Rceattle::linkage_spec(~ temp, by = ~ fleet, fleet = flt))),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))

  tbl <- fit$data_list$linkage_table
  testthat::expect_equal(sum(tbl$process == "sel"), 2L)  # intercept + slope
  testthat::expect_lt(fit$opt$objective, 10241.0304274978)
})


testthat::test_that("unwired sel params and forms are rejected, not silently dropped", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  d$env_data <- .sel_env_data(d)

  # coff is reserved until the non-parametric consume site is wired.
  testthat::expect_error(
    Rceattle::fit_mod(
      data_list = d, estimateMode = 3, msmMode = 0,
      selFun = Rceattle::build_selectivity(linkages = list(
        coff = Rceattle::linkage_spec(~ temp, by = ~ fleet, fleet = 1))),
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                          verbose = 0)),
    "not yet wired")

  # A non-logistic fleet is rejected.
  fc <- Rceattle::switch_check(Rceattle::clean_data(d))$fleet_control
  nonlogi <- which(!fc$Selectivity %in%
                     c("Logistic", "DoubleLogistic", "DescendingLogistic") &
                     fc$Fleet_type != "Off")
  testthat::skip_if(length(nonlogi) == 0)
  testthat::expect_error(
    Rceattle::fit_mod(
      data_list = d, estimateMode = 3, msmMode = 0,
      selFun = Rceattle::build_selectivity(linkages = list(
        inf_asc = Rceattle::linkage_spec(~ temp, by = ~ fleet,
                                         fleet = nonlogi[1]))),
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                          verbose = 0)),
    "not yet wired")
})
