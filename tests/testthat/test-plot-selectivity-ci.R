# plot_selectivity(add_ci = TRUE) draws the interval that
# fit_control(selectivity_se = TRUE) makes available (issue #107).
#
# A figure that assembles is not a pass -- these force ggplot_build(), per
# inst/dev/TRAPS.md.

layer_geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], character(1))

testthat::test_that("add_ci draws a band that brackets the curve", {
  testthat::skip_on_cran()

  m <- suppressMessages(suppressWarnings(fit_mod(
    data_list = Atka2022, msmMode = 0, estimateMode = "Hindcast",
    fit_control = fit_control(selectivity_se = TRUE, verbose = 0))))

  plain <- plot_selectivity(m)
  testthat::expect_false("GeomRibbon" %in% layer_geoms(plain))

  p <- plot_selectivity(m, add_ci = TRUE)
  testthat::expect_true("GeomRibbon" %in% layer_geoms(p))

  g <- ggplot2::ggplot_build(p)
  rib <- g$data[[which(layer_geoms(p) == "GeomRibbon")[1]]]
  lin <- g$data[[which(layer_geoms(p) == "GeomLine")[1]]]
  testthat::expect_gt(nrow(rib), 0)
  testthat::expect_true(all(is.finite(rib$ymin)))
  testthat::expect_true(all(is.finite(rib$ymax)))
  testthat::expect_true(all(rib$ymin <= rib$ymax))
  testthat::expect_true(all(rib$ymin > 0))   # exp() of a log interval

  # The band has to bracket the curve it belongs to, not sit beside it.
  testthat::expect_equal(nrow(rib), nrow(lin))
  testthat::expect_true(all(rib$ymin <= lin$y + 1e-8))
  testthat::expect_true(all(rib$ymax >= lin$y - 1e-8))
})


testthat::test_that("add_ci is declined, with a reason, when the fit has no SEs", {
  testthat::skip_on_cran()

  m <- suppressMessages(suppressWarnings(fit_mod(
    data_list = Atka2022, msmMode = 0, estimateMode = "Hindcast",
    fit_control = fit_control(verbose = 0))))

  testthat::expect_warning(p <- plot_selectivity(m, add_ci = TRUE),
                           "selectivity_se")
  # Declined, not failed: the curve is still drawn.
  testthat::expect_false("GeomRibbon" %in% layer_geoms(p))
  testthat::expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})


testthat::test_that("a Fixed fleet draws no band while the rest of the figure does", {
  testthat::skip_on_cran()

  # BS2017SS's EIT_Pollock is Selectivity = "Fixed". It estimates nothing, so it
  # has no standard error to draw, but it must not suppress the other six.
  m <- suppressMessages(suppressWarnings(fit_mod(
    data_list = BS2017SS, msmMode = 0, estimateMode = "Hindcast",
    fit_control = fit_control(selectivity_se = TRUE, verbose = 0))))

  se <- Rceattle:::.rce_sel_se(m)
  fc <- m$data_list$fleet_control
  fixed <- fc$Fleet_code[fc$Selectivity == "Fixed"]
  testthat::expect_gt(length(fixed), 0)
  testthat::expect_length(intersect(se$Fleet, fixed), 0)
  testthat::expect_gt(length(setdiff(unique(se$Fleet), fixed)), 0)

  p <- plot_selectivity(m, add_ci = TRUE)
  testthat::expect_true("GeomRibbon" %in% layer_geoms(p))
  testthat::expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})


testthat::test_that("a mirrored fleet borrows its lead's interval", {
  # The template reports the lead once, so a mirror has no rows of its own. It
  # shares the lead's parameter block, so a bare line beside a banded lead would
  # read as the more certain of the two. No bundled dataset mirrors selectivity,
  # so the fallback is exercised directly.
  sel <- c(0.1, 0.4, 0.8, 1, 1)
  se <- data.frame(Fleet = 1L, Sex = 1L, Age = 1:5, Year = 2000L,
                   log_sel = log(sel), sd = 0.2)

  own <- Rceattle:::.rce_sel_ci(se, flt = 1L, sex = 1L, bins = 1:5, year = 2000L)
  testthat::expect_equal(own$lower, exp(se$log_sel - 1.96 * se$sd))

  # Fleet 2 mirrors fleet 1 and has no rows: it gets fleet 1's band.
  mirror <- Rceattle:::.rce_sel_ci(se, flt = 2L, sex = 1L, bins = 1:5,
                                   year = 2000L, lead = c(1L, 2L), curve = sel)
  testthat::expect_equal(mirror$lower, own$lower)
  testthat::expect_equal(mirror$upper, own$upper)

  # A fleet that mirrors nothing and has no rows still gets nothing.
  testthat::expect_null(
    Rceattle:::.rce_sel_ci(se, flt = 3L, sex = 1L, bins = 1:5, year = 2000L))

  # An age the fit did not report comes back NA rather than shifting the band
  # onto the wrong age.
  gap <- Rceattle:::.rce_sel_ci(se, flt = 1L, sex = 1L, bins = 1:7, year = 2000L)
  testthat::expect_length(gap$lower, 7L)
  testthat::expect_true(all(is.na(gap$lower[6:7])))
})


testthat::test_that("a mirror whose curve differs from its lead's borrows nothing", {
  # Sharing a Selectivity_index shares the parameter block, not necessarily the
  # curve: data_check() only WARNS when a group differs in a shaping column, and
  # a differing Sel_norm_bin alone rescales the whole curve. Borrowing then puts
  # the lead's band around a curve it does not belong to.
  sel <- c(0.1, 0.4, 0.8, 1, 1)
  se <- data.frame(Fleet = 1L, Sex = 1L, Age = 1:5, Year = 2000L,
                   log_sel = log(sel), sd = 0.2)

  testthat::expect_null(
    Rceattle:::.rce_sel_ci(se, flt = 2L, sex = 1L, bins = 1:5, year = 2000L,
                           lead = c(1L, 2L), curve = sel * 2))
  # The same mirror, normalized the same way, still gets the band.
  testthat::expect_false(is.null(
    Rceattle:::.rce_sel_ci(se, flt = 2L, sex = 1L, bins = 1:5, year = 2000L,
                           lead = c(1L, 2L), curve = sel)))
})


testthat::test_that("the lead is the group's first estimated fleet, not its index value", {
  # rearrange_data() builds flt_sel_lead with .group_lead(): the key is
  # Selectivity_index AND Selectivity together, and an Off fleet never leads. So
  # the group key need not be any member's Fleet_code, and reading it as one
  # drops the band on every mirror in such a group.
  fc <- data.frame(
    Fleet_code = 1:4,
    Fleet_type = c("Off", "Fishery", "Survey", "Survey"),
    Selectivity = c("Logistic", "Logistic", "Logistic", "NonParametric"),
    Selectivity_index = c(7L, 7L, 7L, 7L),
    stringsAsFactors = FALSE)

  # Fleet 3 mirrors fleet 2. Fleet 1 shares the key but is Off and never leads;
  # no fleet has Fleet_code 7 at all.
  testthat::expect_equal(Rceattle:::.rce_sel_lead(fc, 3L), c(2L, 3L))
  testthat::expect_equal(Rceattle:::.rce_sel_lead(fc, 2L), c(2L, 3L))
  # Fleet 4 shares the index but not the form, so it is its own group.
  testthat::expect_equal(Rceattle:::.rce_sel_lead(fc, 4L), 4L)
  # No Selectivity_index column at all: a fleet leads itself.
  testthat::expect_equal(
    Rceattle:::.rce_sel_lead(fc[, setdiff(names(fc), "Selectivity_index")], 3L), 3L)
})
