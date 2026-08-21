# The shared arguments on the predation plotters, asserted against the BUILT
# plot rather than the spec.
#
# These five declared line_col, lwd, lty, minyr, incl_mean and add_ci and read
# none of them, so a smoke test that only checks "it runs" cannot tell the fix
# from the bug. Every assertion here reads ggplot_build() output or the layer
# data, which is what the argument is supposed to change.
#
# Two cases need a fixture the bundled data cannot supply:
#   * every bundled model is sex-combined (nsex = 1), so the line-type-keys-on-
#     sex path is invisible without a two-sex fit;
#   * models sharing a figure normally share an endyr, so the per-model hindcast
#     mean is invisible without one that ends earlier.
testthat::skip_on_cran()

ms_fit <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    testthat::skip_if_not_installed("TMB")
    data("BS2017MS", package = "Rceattle")
    cached <<- suppressMessages(suppressWarnings(
      Rceattle::fit_mod(data_list = BS2017MS, estimateMode = 3, msmMode = 1,
                        fit_control = Rceattle::fit_control(verbose = 0))))
    cached
  }
})

# Give a sex-combined fit a second sex, so the line-type-keys-on-sex path is
# reachable. The bundled multispecies data is sex-combined, so the sex slot has
# to be grown, not just filled: sex 2 is sex 1 scaled, which keeps the series
# distinguishable without pretending to be a real sexed fit.
as_two_sex <- function(fit, scale = 0.5) {
  d <- fit
  d$data_list$nsex <- rep(2L, d$data_list$nspp)
  for (q in c("M_at_age", "consumption_at_age", "biomass_at_age", "N_at_age")) {
    a <- d$quantities[[q]]
    if (is.null(a) || length(dim(a)) != 4L) next
    if (dim(a)[2] < 2L) {
      dm <- dim(a); dm[2] <- 2L
      grown <- array(NA_real_, dim = dm)
      grown[, 1, , ] <- a[, 1, , ]
      a <- grown
    }
    a[, 2, , ] <- a[, 1, , ] * scale
    d$quantities[[q]] <- a
  }
  d
}

built <- function(p, layer = 1L) ggplot2::ggplot_build(p)$data[[layer]]

# The five, with the column each maps to colour and the value column.
predation_plotters <- list(
  list(fn = "plot_b_eaten",        colour = "models",    value = "value"),
  list(fn = "plot_b_eaten_prop",   colour = "predators", value = "value"),
  list(fn = "plot_m_at_age",       colour = "models",    value = "M"),
  list(fn = "plot_m2_at_age_prop", colour = "predators", value = "Proportion"),
  list(fn = "plot_ration",         colour = "models",    value = "value")
)


testthat::test_that("lwd, lty and line_col reach the drawn line", {
  fit <- ms_fit()
  for (spec in predation_plotters) {
    f <- getExportedValue("Rceattle", spec$fn)
    testthat::expect_equal(unique(built(f(fit))$linewidth), 1,
                           info = paste(spec$fn, "default lwd = 3"))
    testthat::expect_equal(unique(built(f(fit, lwd = 6))$linewidth), 2,
                           info = paste(spec$fn, "lwd = 6"))
    # Recycling one colour over several predators warns; that is checked below.
    testthat::expect_true(
      "red" %in% suppressWarnings(built(f(fit, line_col = "red"))$colour),
      info = paste(spec$fn, "line_col"))
    # A base-graphics palette index must not error.
    testthat::expect_true(
      "black" %in% suppressWarnings(built(f(fit, line_col = 1))$colour),
      info = paste(spec$fn, "line_col = 1"))
  }
})

testthat::test_that("too few colours for the mapped variable warns", {
  # plot_b_eaten_prop colours predators, not models, so `line_col = 1` -- which
  # the ESR scripts pass -- would silently draw every predator black.
  fit <- ms_fit()
  testthat::expect_warning(Rceattle::plot_b_eaten_prop(fit, line_col = 1),
                           "colour separates")
  testthat::expect_warning(Rceattle::plot_m2_at_age_prop(fit, line_col = 1),
                           "colour separates")
  # Enough colours, no warning.
  cols <- .okabe_ito[seq_len(fit$data_list$nspp)]
  testthat::expect_silent(Rceattle::plot_b_eaten_prop(fit, line_col = cols))
})

testthat::test_that("line type keys on sex, and both sexes are drawn", {
  # The regression this guards: passing a scalar `lty` as a fixed geom parameter
  # overrode the linetype = Sex mapping and drew both sexes as one line.
  fit <- as_two_sex(ms_fit())
  for (fn in c("plot_ration", "plot_m_at_age")) {
    f <- getExportedValue("Rceattle", fn)
    p <- f(fit)
    testthat::expect_setequal(unique(as.character(p$data$Sex)),
                              c("Female", "Male"))
    testthat::expect_length(unique(built(p)$linetype), 2L)
    # A vector overrides the defaults, one value per sex.
    testthat::expect_setequal(unique(built(f(fit, lty = c(1, 3)))$linetype),
                              c(1, 3))
  }
})

testthat::test_that("line type keys on model where the figure separates models", {
  # plot_b_eaten_prop() and plot_m2_at_age_prop() use colour for predators, so
  # models are told apart by line type. A scalar lty must not flatten that.
  fit <- ms_fit()
  two <- list(fit, fit)
  for (fn in c("plot_b_eaten_prop", "plot_m2_at_age_prop")) {
    f <- getExportedValue("Rceattle", fn)
    p <- f(two, model_names = c("A", "B"))
    testthat::expect_length(unique(built(p)$linetype), 2L)
    testthat::expect_setequal(
      unique(built(f(two, model_names = c("A", "B"), lty = c(1, 3)))$linetype),
      c(1, 3))
  }
  # plot_b_eaten maps no line type of its own, so a scalar applies to the line.
  testthat::expect_equal(
    unique(built(Rceattle::plot_b_eaten(fit, lty = 2))$linetype), 2)
})

testthat::test_that("a varying lty warns when its key has one level", {
  # Sex-combined models are the common case, so silently dropping the extra
  # values would leave `lty` looking functional while doing nothing.
  fit <- ms_fit()
  testthat::expect_equal(unique(fit$data_list$nsex), 1)
  testthat::expect_warning(Rceattle::plot_ration(fit, lty = c(1, 2)),
                           "which has 1 level")
  testthat::expect_warning(Rceattle::plot_m_at_age(fit, lty = c(1, 2)),
                           "which has 1 level")
})

testthat::test_that("minyr and maxyr narrow the data so the panel rescales", {
  # Clipping the axis alone left the y scale trained on the hidden years.
  fit <- ms_fit()
  yr <- fit$data_list$styr + 10L
  for (spec in predation_plotters) {
    f <- getExportedValue("Rceattle", spec$fn)
    p <- f(fit, minyr = yr)
    testthat::expect_gte(min(p$data$Year), yr, label = spec$fn)
    testthat::expect_lte(max(f(fit, maxyr = yr)$data$Year), yr, label = spec$fn)
  }
  # The y range must follow the window, not the whole series.
  full <- ggplot2::ggplot_build(Rceattle::plot_b_eaten(fit))
  cut  <- ggplot2::ggplot_build(Rceattle::plot_b_eaten(fit, minyr = yr))
  testthat::expect_false(identical(full$layout$panel_params[[1]]$y.range,
                                   cut$layout$panel_params[[1]]$y.range))
})

testthat::test_that("incl_mean draws each model's own hindcast mean", {
  # With models of different endyr, cutting every one at the first model's
  # endyr averages a peel over years it never fitted -- and the answer then
  # depends on the order of the list.
  a <- ms_fit()
  b <- a
  b$data_list$endyr <- b$data_list$endyr - 7L
  sp1 <- a$data_list$spnames[1]

  means_for <- function(models, names) {
    p <- Rceattle::plot_b_eaten(models, model_names = names,
                                incl_proj = TRUE, incl_mean = TRUE)
    hl <- p$layers[[which(vapply(p$layers,
      function(l) inherits(l$geom, "GeomHline"), logical(1)))]]
    d <- hl$data[as.character(hl$data$Species) == sp1, ]
    stats::setNames(d$.mean, as.character(d$Model))
  }
  ab <- means_for(list(a, b), c("A", "B"))
  ba <- means_for(list(b, a), c("B", "A"))
  testthat::expect_equal(ab[["A"]], ba[["A"]])
  testthat::expect_equal(ab[["B"]], ba[["B"]])
  testthat::expect_false(isTRUE(all.equal(ab[["A"]], ab[["B"]])))

  # And the value is the mean over that model's own hindcast.
  d <- Rceattle::plot_b_eaten(list(a, b), model_names = c("A", "B"),
                              incl_proj = TRUE)$data
  hind_b <- d[d$Model == "B" & as.character(d$Species) == sp1 &
                d$Year <= b$data_list$endyr, ]
  testthat::expect_equal(ab[["B"]], mean(hind_b$value))
})

testthat::test_that("add_ci warns rather than silently drawing nothing", {
  # None of these quantities are ADREPORTed, so no interval is available.
  fit <- ms_fit()
  for (spec in predation_plotters) {
    f <- getExportedValue("Rceattle", spec$fn)
    testthat::expect_warning(f(fit, add_ci = TRUE), "not available",
                             info = spec$fn)
  }
})

testthat::test_that("species selects and orders panels the same way everywhere", {
  fit <- ms_fit()
  panel_of <- list(plot_b_eaten = "Species", plot_b_eaten_prop = "Prey",
                   plot_m_at_age = "Species", plot_m2_at_age_prop = "Prey",
                   plot_ration = "Species")
  nm <- fit$data_list$spnames
  for (fn in names(panel_of)) {
    f <- getExportedValue("Rceattle", fn)
    col <- panel_of[[fn]]
    testthat::expect_equal(levels(droplevels(f(fit, species = 2)$data[[col]])),
                           nm[2], info = fn)
    # Panels follow the order asked for, not model order.
    testthat::expect_equal(levels(droplevels(f(fit, species = c(2, 1))$data[[col]])),
                           nm[c(2, 1)], info = fn)
    # Selection by name matches selection by index.
    testthat::expect_equal(f(fit, species = nm[2])$data[[col]],
                           f(fit, species = 2)$data[[col]], info = fn)
  }
})

testthat::test_that("plot_m2_at_age_prop draws shares that sum to one", {
  # M2_prop holds each predator's contribution to M2, not a share of it --
  # summing it over predators reproduces M2_at_age. Plotting that sum gave a
  # "proportion" reaching 1500.
  fit <- ms_fit()
  d <- Rceattle::plot_m2_at_age_prop(fit)$data
  testthat::expect_true(all(d$Proportion >= 0 & d$Proportion <= 1,
                            na.rm = TRUE))

  # The shares of one prey panel in one year account for all of its predation.
  tot <- tapply(d$Proportion, list(d$Prey, d$Year), sum, na.rm = TRUE)
  tot <- tot[!is.na(tot) & tot > 0]
  testthat::expect_equal(as.numeric(tot), rep(1, length(tot)))

  # The underlying contributions still reproduce M2_at_age, i.e. the
  # normalisation is what changed, not the quantity.
  q <- fit$quantities
  testthat::expect_equal(sum(q$M2_prop[, 1, , 1, 1]), q$M2_at_age[1, 1, 1, 1])
})

testthat::test_that("plot_ration multiplies the ration by numbers, not biomass", {
  # consumption_at_age is one fish's annual ration in kg; numbers-at-age are in
  # thousands, so the product is mt. Multiplying by biomass-at-age weights the
  # age-sum by weight-at-age and is not a quantity in any unit.
  fit <- ms_fit()
  q <- fit$quantities
  sp <- 1L
  ages <- seq_len(fit$data_list$nages[sp])
  expected <- sum(q$consumption_at_age[sp, 1, ages, 1] *
                    q$N_at_age[sp, 1, ages, 1]) / 1e6

  d <- Rceattle::plot_ration(fit)$data
  got <- d$value[as.character(d$Species) == fit$data_list$spnames[sp] &
                   d$Year == min(d$Year)]
  testthat::expect_equal(got, expected)

  # And it is not the old biomass-weighted number.
  wrong <- sum(q$consumption_at_age[sp, 1, ages, 1] *
                 q$biomass_at_age[sp, 1, ages, 1]) / 1e6
  testthat::expect_false(isTRUE(all.equal(got, wrong)))
})

testthat::test_that("the MSE path validates and honours its arguments", {
  fit <- ms_fit()
  mse <- list(list(OM = fit), list(OM = fit), list(OM = fit))

  # Ribbon transparency follows alpha rather than a hardcoded value.
  testthat::expect_equal(
    unique(built(Rceattle::plot_b_eaten(mse, mse = TRUE, alpha = 0.2))$alpha),
    0.2)
  # incl_mean is no longer skipped by the early return.
  geoms <- vapply(Rceattle::plot_b_eaten(mse, mse = TRUE,
                                         incl_mean = TRUE)$layers,
                  function(l) class(l$geom)[1], character(1))
  testthat::expect_true("GeomHline" %in% geoms)
  # lwd is validated here too, not just on the non-MSE path.
  testthat::expect_error(Rceattle::plot_b_eaten(mse, mse = TRUE, lwd = -1),
                         "non-negative")
  testthat::expect_equal(
    unique(built(Rceattle::plot_b_eaten(mse, mse = TRUE, lwd = 6), 2L)$linewidth),
    2)
  # line_col has nothing to colour once the simulations are pooled.
  testthat::expect_warning(
    Rceattle::plot_b_eaten(mse, mse = TRUE, line_col = "red"),
    "not used with")
})
