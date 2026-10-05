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


testthat::test_that("a mirror on a FITTED object still resolves its lead", {
  testthat::skip_on_cran()
  # The one that matters: `fit$data_list` is the PRE-rearrange_data() list, so it
  # carries no flt_sel_lead and no flt_sel_type. A lead resolved from those
  # fields is NULL on every real fit, and the mirror silently loses its band
  # beside a banded lead -- which reads as the more certain of the two. Asserted
  # on a fit rather than a hand-built list, because a hand-built list is exactly
  # what hid this.
  d <- Atka2022
  testthat::skip_if(nrow(d$fleet_control) < 2, "need two fleets to mirror")
  d$fleet_control[2, ] <- d$fleet_control[1, ]          # fleet 2 mirrors fleet 1
  d$fleet_control$Fleet_code[2] <- 2
  d$fleet_control$Fleet_name[2] <- "mirror_of_1"
  d$fleet_control$Selectivity_index[2] <- d$fleet_control$Selectivity_index[1]

  m <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, msmMode = 0, estimateMode = "Hindcast",
    fit_control = fit_control(selectivity_se = TRUE, verbose = 0))))

  testthat::expect_null(m$data_list$flt_sel_lead)       # the field that is absent
  fc <- m$data_list$fleet_control
  grp <- which(fc$Selectivity_index == fc$Selectivity_index[2])
  testthat::skip_if(length(grp) < 2, "fleets did not end up sharing selectivity")

  # The template reports the lead only, and the mirror resolves to it.
  se <- Rceattle:::.rce_sel_se(m)
  lead <- Rceattle:::.rce_sel_lead(m$data_list, 2L)
  testthat::expect_true(any(setdiff(lead, 2L) %in% se$Fleet))

  # And the band actually comes back, rather than NULL.
  yr  <- m$data_list$styr
  sel <- m$quantities$sel_at_age
  bins <- seq_len(m$data_list$nages[fc$Species[2]]) - 1L +
    m$data_list$minage[fc$Species[2]]
  ci <- Rceattle:::.rce_sel_ci(se, flt = 2L, sex = 1L, bins = bins, year = yr,
                               lead = lead,
                               curve = as.numeric(sel[2, 1, seq_along(bins), 1]))
  testthat::expect_false(is.null(ci))
  testthat::expect_true(any(is.finite(ci$lower)))
})


testthat::test_that("the lead is keyed the way .group_lead() keys it", {
  # The key is Selectivity_index AND the Selectivity code together, and an Off
  # fleet never leads -- so the group key need not be any member's Fleet_code,
  # and reading it as one drops the band on every mirror in such a group.
  fc <- data.frame(
    Fleet_code = 1:4,
    Fleet_type = c("Off", "Fishery", "Survey", "Survey"),
    Selectivity = c(1L, 1L, 1L, 2L),        # Logistic, Logistic, Logistic, NonParametric
    Selectivity_index = c(7L, 7L, 7L, 7L),
    stringsAsFactors = FALSE)
  dl <- list(fleet_control = fc)

  # Fleet 2 leads; fleet 1 shares the key but is Off, and no fleet has code 7.
  testthat::expect_equal(Rceattle:::.rce_sel_lead(dl, 3L), c(2L, 3L))
  # A lead borrows from nobody: it has rows of its own.
  testthat::expect_equal(Rceattle:::.rce_sel_lead(dl, 2L), 2L)
  # Fleet 4 shares the index but not the form, so it is its own group.
  testthat::expect_equal(Rceattle:::.rce_sel_lead(dl, 4L), 4L)
  # No Selectivity_index column at all: a fleet leads itself.
  dl2 <- dl; dl2$fleet_control$Selectivity_index <- NULL
  testthat::expect_equal(Rceattle:::.rce_sel_lead(dl2, 3L), 3L)
})


# The sdreport a fit ends on is not always the one that estimated the curve.
# Under an estimating HCR it is the projection's, where build_hcr_map() has fixed
# every selectivity parameter, so the delta method returns 0 at every age. Drawn,
# that is a hairline band reading as certainty -- worse than the bare line.
testthat::test_that("an all-zero standard error is declined, not drawn as a hairline", {
  testthat::skip_on_cran()

  d <- Atka2022
  fsh <- which(as.character(d$fleet_control$Fleet_type) %in% c("Fishery", "1"))
  d$fleet_control$Proj_F_proportion[fsh] <- 1 / length(fsh)

  # fit_mod() says so at the point the setting is made, before the fit is spent.
  testthat::expect_warning(
    m <- suppressMessages(fit_mod(
      data_list = d, msmMode = 0, estimateMode = "Estimate",
      HCR = build_hcr(HCR = 5, Ftarget = 0.4, Flimit = 0.35,
                      Plimit = 0.2, Alpha = 0.05),
      fit_control = fit_control(selectivity_se = TRUE, verbose = 0))),
    "standard error of 0 for every age")

  # The errors are reported, so the plotter cannot decline on their absence.
  se <- Rceattle:::.rce_sel_se(m)
  testthat::expect_false(is.null(se))
  testthat::expect_true(all(se$sd == 0))

  testthat::expect_warning(p <- plot_selectivity(m, add_ci = TRUE),
                           "standard error of 0 for every age")
  # Declined, not failed: the curve is still drawn, just without a band.
  testthat::expect_false("GeomRibbon" %in% layer_geoms(p))
  testthat::expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})
