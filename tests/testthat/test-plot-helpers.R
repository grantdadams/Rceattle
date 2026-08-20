# Unit tests for the shared plotting-argument helpers (R/7-plot_helpers.R).
#
# These are what make `line_col`, `lwd`, `lty`, `alpha` and `species` behave the
# same way in every plot_*() function. They need no fitted model, so they stay
# unguarded (fast) while the plotters that use them are skip_on_cran.
#
# The load-bearing assertion is the `lwd / 3` bridge: every plotter converted to
# .rce_line_params() previously hardcoded `linewidth = 1`, so the default
# `lwd = 3` must keep producing exactly that or every existing figure changes.

# A minimal stand-in for a fitted model: .resolve_species() only reads
# `data_list$nspp` and `data_list$spnames`.
fake_models <- function(nspp = 3L, spnames = c("Pollock", "Cod", "ATF")) {
  list(list(data_list = list(nspp = nspp, spnames = spnames)))
}

# Colours of the first layer of a built plot.
built_colours <- function(p) {
  unique(ggplot2::ggplot_build(p)$data[[1]]$colour)
}

demo_df <- data.frame(
  x     = 1:4,
  y     = 1:4,
  Model = factor(rep(c("a", "b"), each = 2)),
  Year  = 1:4
)


# --- .as_colour ---------------------------------------------------------------
testthat::test_that(".as_colour resolves base-graphics colour specifications", {
  # `line_col = 1` and `line_col = c(1, 2, 4)` are how the real assessment
  # scripts pass colours; ggplot2's manual scales reject bare integers.
  testthat::expect_equal(.as_colour(1), "black")
  testthat::expect_length(.as_colour(c(1, 2, 4)), 3L)
  testthat::expect_equal(.as_colour("black"), "black")
  testthat::expect_equal(.as_colour("#0072B2"), "#0072B2")
  testthat::expect_null(.as_colour(NULL))
  # Base treats 0 as transparent rather than as palette entry 0.
  testthat::expect_true(is.na(.as_colour(0)))
  # Recycling matches base: index 9 wraps around an 8-colour palette.
  testthat::expect_equal(.as_colour(length(grDevices::palette()) + 1L),
                         .as_colour(1))
})


# --- .resolve_species ---------------------------------------------------------
testthat::test_that(".resolve_species accepts every spelling of a selection", {
  m <- fake_models()

  testthat::expect_equal(.resolve_species(m)$index, 1:3)
  testthat::expect_equal(.resolve_species(m, "all")$index, 1:3)
  testthat::expect_equal(.resolve_species(m, c(TRUE, FALSE, TRUE))$index, c(1L, 3L))
  testthat::expect_equal(.resolve_species(m, c("Cod", "ATF"))$index, c(2L, 3L))

  # The caller's order is preserved, so `species = c(3, 1)` faceted in that
  # order rather than silently sorting back to model order.
  res <- .resolve_species(m, c(3, 1))
  testthat::expect_equal(res$index, c(3L, 1L))
  testthat::expect_equal(res$labels, c("ATF", "Pollock"))
})

testthat::test_that(".resolve_species reads unmatched names as labels", {
  # Back-compatibility: plot_selectivity() and plot_maturity() used `species`
  # for display labels, and the assessment scripts still call them that way.
  m <- fake_models()
  ebs <- c("Walleye pollock", "Pacific cod", "Arrowtooth flounder")

  testthat::expect_message(res <- .resolve_species(m, ebs), "display labels")
  testthat::expect_equal(res$index, 1:3)      # nothing is dropped
  testthat::expect_equal(res$labels, ebs)     # the strings become the labels

  # An explicit `spnames` wins: the labels are not overwritten.
  testthat::expect_message(
    res2 <- .resolve_species(m, ebs, spnames = c("A", "B", "C")))
  testthat::expect_equal(res2$labels, c("A", "B", "C"))
})

testthat::test_that(".resolve_species warns on a partial match and errors on bad input", {
  m <- fake_models()

  # Some names matched: a selection with a typo, not a set of labels.
  testthat::expect_warning(res <- .resolve_species(m, c("Cod", "Halibut")),
                           "did not match")
  testthat::expect_equal(res$index, 2L)

  testthat::expect_error(.resolve_species(m, 5), "out of range")
  testthat::expect_error(.resolve_species(m, 0), "out of range")
  testthat::expect_error(.resolve_species(m, integer(0)), "selected no species")
  testthat::expect_error(.resolve_species(m, list(1)), "must be species indices")
})

testthat::test_that(".resolve_species falls back to positional labels", {
  m <- list(list(data_list = list(nspp = 2L, spnames = NULL)))
  testthat::expect_equal(.resolve_species(m)$spnames,
                         c("Species 1", "Species 2"))
})


# --- .rce_line_params ---------------------------------------------------------
testthat::test_that("the default lwd renders as linewidth 1", {
  # The invariant that keeps every converted plotter's default figure identical.
  lp <- .rce_line_params(lwd = 3)
  testthat::expect_equal(lp$args$linewidth, 1)
  testthat::expect_null(lp$mapping)
  testthat::expect_length(lp$scales, 0L)

  testthat::expect_equal(.rce_line_params(lwd = 6)$args$linewidth, 2)
  testthat::expect_equal(.rce_line_params(lwd = 2)$args$linewidth, 2 / 3)
})

testthat::test_that("a varying lwd/lty is mapped when there is a key for it", {
  lp <- .rce_line_params(lwd = c(3, 6), lwd_by = "Model",
                         lwd_levels = c("a", "b"))
  testthat::expect_false(is.null(lp$mapping))
  testthat::expect_length(lp$scales, 1L)
  testthat::expect_null(lp$args$linewidth)   # mapped, not fixed

  lp2 <- .rce_line_params(lty = c(1, 2), lty_by = "Sex",
                          lty_levels = c("Female", "Male"))
  testthat::expect_false(is.null(lp2$mapping))
  testthat::expect_null(lp2$args$linetype)
})

testthat::test_that("a varying lwd/lty without a key warns and uses the first", {
  testthat::expect_warning(lp <- .rce_line_params(lwd = c(3, 6)),
                           "no per-series line width")
  testthat::expect_equal(lp$args$linewidth, 1)

  testthat::expect_warning(lp2 <- .rce_line_params(lty = c(1, 2)),
                           "no per-series line type")
  testthat::expect_equal(lp2$args$linetype, 1)
})

testthat::test_that("alpha is passed through only when supplied", {
  testthat::expect_null(.rce_line_params()$args$alpha)
  testthat::expect_equal(.rce_line_params(alpha = 0.25)$args$alpha, 0.25)
})

testthat::test_that(".rce_add_line reaches the built plot", {
  p <- ggplot2::ggplot(demo_df, ggplot2::aes(x = .data$x, y = .data$y,
                                             colour = .data$Model))
  b <- ggplot2::ggplot_build(.rce_add_line(p, .rce_line_params(lwd = 6, lty = 2)))
  testthat::expect_equal(unique(b$data[[1]]$linewidth), 2)
  testthat::expect_equal(unique(b$data[[1]]$linetype), 2)
})


# --- .rceattle_scale ----------------------------------------------------------
testthat::test_that("line_col = NULL keeps the package default palettes", {
  # Guards against the extension changing any existing figure.
  disc <- ggplot2::ggplot(demo_df, ggplot2::aes(.data$x, .data$y,
                                                colour = .data$Model)) +
    ggplot2::geom_line()
  testthat::expect_equal(built_colours(.rceattle_scale(disc, aesthetics = "colour")),
                         .okabe_ito[1:2])

  cont <- ggplot2::ggplot(demo_df, ggplot2::aes(.data$x, .data$y,
                                                colour = .data$Year,
                                                group = .data$Model)) +
    ggplot2::geom_line()
  testthat::expect_true(all(grepl("^#",
    built_colours(.rceattle_scale(cont, discrete = FALSE, aesthetics = "colour")))))
})

testthat::test_that("line_col overrides a discrete colour scale", {
  p <- ggplot2::ggplot(demo_df, ggplot2::aes(.data$x, .data$y,
                                             colour = .data$Model)) +
    ggplot2::geom_line()
  testthat::expect_equal(
    built_colours(.rceattle_scale(p, aesthetics = "colour",
                                  line_col = c("red", "blue"),
                                  levels = c("a", "b"))),
    c("red", "blue"))

  # Integer colours must not error.
  testthat::expect_equal(
    built_colours(.rceattle_scale(p, aesthetics = "colour", line_col = 1,
                                  levels = c("a", "b"))),
    "black")

  # Without levels the colours are still applied, in level order.
  testthat::expect_equal(
    built_colours(.rceattle_scale(p, aesthetics = "colour",
                                  line_col = c("red", "blue"))),
    c("red", "blue"))
})

testthat::test_that("line_col supplies ramp anchors on a continuous colour scale", {
  p <- ggplot2::ggplot(demo_df, ggplot2::aes(.data$x, .data$y,
                                             colour = .data$Year,
                                             group = .data$Model)) +
    ggplot2::geom_line()
  # One colour: the whole year fan drawn in it -- "line_col = 'black'" meaning
  # what a base-graphics user expects it to mean.
  testthat::expect_equal(
    built_colours(.rceattle_scale(p, discrete = FALSE, aesthetics = "colour",
                                  line_col = "black")),
    "#000000")
  # Two or more: interpolated between them, so the fan still reads as a ramp.
  many <- built_colours(.rceattle_scale(p, discrete = FALSE,
                                        aesthetics = "colour",
                                        line_col = c("red", "blue")))
  testthat::expect_gt(length(many), 1L)
  testthat::expect_true("#FF0000" %in% many)
})


# --- year window --------------------------------------------------------------
testthat::test_that(".rce_year_limits clips without filtering", {
  # Year on the x axis: p$data must stay complete, because the accuracy tests
  # compare it against the model's own quantities.
  testthat::expect_null(.rce_year_limits(NULL, NULL))
  lim <- .rce_year_limits(1990, 2000)
  testthat::expect_s3_class(lim, "CoordCartesian")
})

testthat::test_that(".rce_year_filter drops rows outside the window", {
  df <- data.frame(Year = 1990:1999, v = 1:10)
  testthat::expect_equal(nrow(.rce_year_filter(df, 1995, 1997)), 3L)
  testthat::expect_equal(nrow(.rce_year_filter(df, minyr = 1998)), 2L)
  testthat::expect_equal(nrow(.rce_year_filter(df)), 10L)
  testthat::expect_error(.rce_year_filter(df, 2050), "No years left")
})


# --- projection divider / mean line -------------------------------------------
testthat::test_that(".rce_proj_divider only draws when projecting", {
  m <- list(list(data_list = list(endyr = 2017)))
  testthat::expect_null(.rce_proj_divider(m, incl_proj = FALSE))
  testthat::expect_s3_class(.rce_proj_divider(m, incl_proj = TRUE), "ggproto")
  # A model without an endyr must not produce a vline at NA.
  testthat::expect_null(
    .rce_proj_divider(list(list(data_list = list(endyr = NULL))), TRUE))
})

testthat::test_that(".rce_mean_line averages hindcast years only", {
  df <- data.frame(Species = "a", Model = "m", Year = 1:10,
                   value = c(rep(1, 5), rep(100, 5)))
  testthat::expect_null(.rce_mean_line(df, incl_mean = FALSE))

  lyr <- .rce_mean_line(df, incl_mean = TRUE, by = c("Species", "Model"),
                        hind_endyr = 5)
  # Mean over years 1:5 only; the projection's 100s must not drag it up.
  testthat::expect_equal(lyr$data$.mean, 1)

  lyr_all <- .rce_mean_line(df, incl_mean = TRUE, by = c("Species", "Model"))
  testthat::expect_equal(lyr_all$data$.mean, mean(df$value))
})


# --- confidence-interval plumbing ---------------------------------------------
testthat::test_that(".rce_series_sd returns NULL rather than a short vector", {
  testthat::expect_null(.rce_series_sd(list(sdrep = NULL), "ssb", 5))
  fit <- list(sdrep = list(value = stats::setNames(1:3, rep("ssb", 3)),
                           sd = c(0.1, 0.2, 0.3)))
  testthat::expect_null(.rce_series_sd(fit, "ssb", 5))          # too few
  testthat::expect_equal(.rce_series_sd(fit, "ssb", 2), c(0.1, 0.2))
  testthat::expect_null(.rce_series_sd(fit, "biomass", 1))      # wrong name
})

testthat::test_that(".rce_no_ci warns only when the interval was asked for", {
  testthat::expect_silent(res <- .rce_no_ci(FALSE, "ration", "no standard errors"))
  testthat::expect_false(res)
  testthat::expect_warning(res2 <- .rce_no_ci(TRUE, "ration", "no standard errors"),
                           "not available for ration")
  testthat::expect_false(res2)
})
