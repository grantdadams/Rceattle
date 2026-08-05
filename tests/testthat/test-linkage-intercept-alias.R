# =============================================================================
# Keying a prior / init / bound to the intercept.
#
# These lists are looked up by exact design-matrix column name, so a covariate
# keys by the name the user wrote (`temp`) but the intercept keys by the name
# model.matrix() gives it, `(Intercept)`, which has to be backtick-quoted.
# `intercept` is accepted as the same request, and a key matching no design
# column is an error rather than a silently unapplied prior -- the two are
# indistinguishable in a fitted model, and on a model where the prior is
# load-bearing that is the difference between converging and running to a bound.
# =============================================================================

.intercept_alias_data <- function() {
  d <- Rceattle::BS2017SS
  comp_flts <- sort(unique(d$comp_data$Fleet_code[d$comp_data$Fleet_code > 0]))
  d$fleet_control$Comp_distribution <- "MultinomialAFSC"
  d$fleet_control$Comp_distribution[comp_flts] <- "DirichletMultinomial"
  list(data = d, fleets = comp_flts)
}

# Prior contribution from a theta_comp prior keyed by `key`.
.intercept_alias_nll <- function(key) {
  s <- .intercept_alias_data()
  p <- list(Rceattle::prior_lognormal(0, 2))
  names(p) <- key
  cf <- Rceattle::build_composition(linkages = list(
    theta_comp = Rceattle::linkage_spec(~ 1, by = ~ fleet, fleet = s$fleets,
                                        priors = p)))
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = s$data, compFun = cf, estimateMode = 3, msmMode = 0,
    random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, phase = FALSE,
                                        verbose = 0))))
  sum(fit$quantities$jnll_comp["Linkage-table priors", ])
}


testthat::test_that("`intercept` keys the same prior as `(Intercept)`", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  canonical <- .intercept_alias_nll("(Intercept)")
  # Guard the fixture: a prior that was never applied would be 0, and then the
  # comparisons below would hold trivially.
  testthat::expect_gt(canonical, 0)
  testthat::expect_equal(.intercept_alias_nll("intercept"), canonical)
  testthat::expect_equal(.intercept_alias_nll("Intercept"), canonical)
})


testthat::test_that("`intercept` keys init and bounds, not just priors", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # init and bounds go through the same translation, and an intercept init is
  # pushed into the base parameter (2-build_params.R), so a dropped one moves
  # where the fit starts. Composition linkages are prior-only, so this uses a
  # catchability linkage. Assert against the built design table rather than the
  # spec: keying the spec is not what matters, reaching the row is.
  d <- Rceattle::BS2017SS
  yrs <- d$styr:d$projyr
  d$env_data <- data.frame(Year = yrs,
                           temp = as.numeric(scale(sin(seq_along(yrs) / 4))))
  built <- function(arg, value) {
    args <- list(~ temp, by = ~ fleet, fleet = 7L)   # fleet 7 = the estimated-q survey
    args[[arg]] <- stats::setNames(list(value), "intercept")
    fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
      data_list = d,
      qFun = Rceattle::build_catchability(linkages = list(
        q = do.call(Rceattle::linkage_spec, args))),
      estimateMode = 3, msmMode = 0, random_rec = FALSE,
      fit_control = Rceattle::fit_control(getsd = FALSE, phase = FALSE,
                                          verbose = 0))))
    as.data.frame(fit$data_list$linkage_table)
  }

  tb <- built("init", 0.25)
  is_int <- tb$design_col == "(Intercept)"
  testthat::expect_true(any(is_int))
  testthat::expect_true(all(tb$init[is_int] == 0.25))
  testthat::expect_true(all(tb$init_supplied[is_int]))

  tb <- built("bounds", c(-2, 2))
  is_int <- tb$design_col == "(Intercept)"
  testthat::expect_true(all(tb$lower[is_int] == -2))
  testthat::expect_true(all(tb$upper[is_int] ==  2))
})


testthat::test_that("a covariate named `intercept` keeps its own prior", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # `intercept` here names a real covariate, so the key means that column --
  # translating it would silently move the prior onto the intercept instead.
  # Both spellings of a real covariate name must survive: a key that already
  # names a design column always means that column.
  testthat::expect_identical(
    names(Rceattle:::.canonical_intercept_names(
      list(intercept = 1), c("(Intercept)", "intercept"))), "intercept")
  testthat::expect_identical(
    names(Rceattle:::.canonical_intercept_names(
      list(Intercept = 1), c("(Intercept)", "Intercept"))), "Intercept")
  # And with no intercept in the design at all, nothing is translated -- the
  # key check then reports it rather than silently dropping it.
  testthat::expect_identical(
    names(Rceattle:::.canonical_intercept_names(list(intercept = 1), "temp")),
    "intercept")

  d <- Rceattle::BS2017SS
  yrs <- d$styr:d$projyr
  d$env_data <- data.frame(Year = yrs,
                           intercept = as.numeric(scale(sin(seq_along(yrs) / 4))))
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d,
    qFun = Rceattle::build_catchability(linkages = list(
      q = Rceattle::linkage_spec(~ intercept, by = ~ fleet, fleet = 7L,
                                 priors = list(intercept = Rceattle::prior_normal(0, 3))))),
    estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, phase = FALSE,
                                        verbose = 0))))
  tb <- as.data.frame(fit$data_list$linkage_table)
  fam <- stats::setNames(tb$prior_family, tb$design_col)
  testthat::expect_equal(unname(fam["intercept"]), "normal")   # the covariate
  testthat::expect_equal(unname(fam["(Intercept)"]), "none")   # not the intercept
})


testthat::test_that("supplying both spellings is an error, not a silent pick", {
  testthat::expect_error(
    Rceattle:::.canonical_intercept_names(
      list(`(Intercept)` = 1, intercept = 2), c("(Intercept)", "temp")),
    "once")
})


testthat::test_that("a key matching no design column is rejected", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  s <- .intercept_alias_data()
  bad <- function(arg, value) {
    args <- list(~ 1, by = ~ fleet, fleet = s$fleets)
    args[[arg]] <- stats::setNames(list(value), "typo")
    cf <- Rceattle::build_composition(linkages = list(
      theta_comp = do.call(Rceattle::linkage_spec, args)))
    suppressWarnings(suppressMessages(Rceattle::fit_mod(
      data_list = s$data, compFun = cf, estimateMode = 3, msmMode = 0,
      random_rec = FALSE,
      fit_control = Rceattle::fit_control(getsd = FALSE, phase = FALSE,
                                          verbose = 0))))
  }
  testthat::expect_error(bad("priors", Rceattle::prior_lognormal(0, 2)),
                         "unknown `priors` name")
  testthat::expect_error(bad("init", 0.5),   "unknown `init` name")
  testthat::expect_error(bad("bounds", c(-1, 1)), "unknown `bounds` name")
})
