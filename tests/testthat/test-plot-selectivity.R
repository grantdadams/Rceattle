# plot_selectivity(): the arguments it takes, the models it draws, and the
# dimension it draws them on.
#
# This function carried three defects at once: line_col/lwd/model_names/species
# were declared and unread, a list of models silently lost all but the first,
# and a length-based fleet was drawn as the growth-matrix conversion of its
# curve on an axis labelled "Age". Each has a test here.
testthat::skip_on_cran()

ss_fit <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    testthat::skip_if_not_installed("TMB")
    data("BS2017SS", package = "Rceattle")
    cached <<- suppressMessages(suppressWarnings(
      Rceattle::fit_mod(data_list = BS2017SS, estimateMode = 3, msmMode = 0,
                        fit_control = Rceattle::fit_control(verbose = 0))))
    cached
  }
})

# A model whose fleets are length-based. No bundled dataset sets
# Selectivity_dimension, which is why the wrong-axis bug survived.
len_fit <- local({
  cached <- NULL
  function(dimension = "Length") {
    key <- paste(dimension, collapse = "|")
    if (!is.null(cached[[key]])) return(cached[[key]])
    testthat::skip_if_not_installed("TMB")
    set.seed(17)
    # An even number of years: the fixture's default Fmort halves it.
    sim <- make_msm_test_data(years = 1:24, use_size_sel = TRUE,
                              log_phi = matrix(-Inf, 2, 2, byrow = TRUE))
    d <- sim$data_list
    d$fleet_control$Selectivity_dimension <-
      rep(dimension, length.out = nrow(d$fleet_control))
    fit <- suppressMessages(suppressWarnings(
      Rceattle::fit_mod(data_list = d, estimateMode = 3, msmMode = 0,
                        fit_control = Rceattle::fit_control(verbose = 0))))
    cached[[key]] <<- list(fit = fit, data_list = d)
    cached[[key]]
  }
})

built <- function(p) ggplot2::ggplot_build(p)$data[[1]]


testthat::test_that("a list of models draws every model, not just the first", {
  # The silent-drop bug: only Rceattle[[1]] was ever read.
  fit <- ss_fit()
  p <- Rceattle::plot_selectivity(list(fit, fit), model_names = c("A", "B"))
  testthat::expect_equal(levels(p$data$Model), c("A", "B"))
  testthat::expect_equal(nrow(p$data),
                         2L * nrow(Rceattle::plot_selectivity(fit)$data))
})

testthat::test_that("colour follows the model count, and colour_by overrides it", {
  fit <- ss_fit()
  # One model: the year fan, as before.
  testthat::expect_equal(
    rlang::as_label(Rceattle::plot_selectivity(fit)$mapping$colour), "Year")
  # Several: colour separates them and the fan moves to transparency.
  p2 <- Rceattle::plot_selectivity(list(fit, fit), model_names = c("A", "B"))
  testthat::expect_equal(rlang::as_label(p2$mapping$colour), "Model")
  testthat::expect_equal(rlang::as_label(p2$mapping$alpha), "Year")
  testthat::expect_equal(range(built(p2)$alpha), c(0.25, 1))

  testthat::expect_equal(
    rlang::as_label(Rceattle::plot_selectivity(fit, colour_by = "model")$mapping$colour),
    "Model")
  testthat::expect_equal(
    rlang::as_label(Rceattle::plot_selectivity(list(fit, fit),
                                               model_names = c("A", "B"),
                                               colour_by = "year")$mapping$colour),
    "Year")
})

testthat::test_that("line_col and lwd reach the drawn line", {
  # The user report: neither passed through.
  fit <- ss_fit()
  # One model: colour is the year, so line_col gives the ramp anchors. A single
  # colour draws the whole fan in it.
  testthat::expect_equal(
    unique(built(Rceattle::plot_selectivity(fit, line_col = "black"))$colour),
    "#000000")
  # Several models: line_col is the model palette.
  testthat::expect_setequal(
    unique(built(Rceattle::plot_selectivity(list(fit, fit),
                                            model_names = c("A", "B"),
                                            line_col = c("red", "blue")))$colour),
    c("red", "blue"))
  testthat::expect_equal(unique(built(Rceattle::plot_selectivity(fit, lwd = 6))$linewidth), 2)
  testthat::expect_equal(unique(built(Rceattle::plot_selectivity(fit, lwd = 3))$linewidth), 1)
})

testthat::test_that("species selects fleets and spnames relabels", {
  fit <- ss_fit()
  fc <- fit$data_list$fleet_control
  n2 <- length(unique(fc$Fleet_name[fc$Species == 2]))
  testthat::expect_equal(nlevels(Rceattle::plot_selectivity(fit, species = 2)$data$Fleet), n2)
  # By name matches by index.
  nm <- fit$data_list$spnames[2]
  testthat::expect_equal(
    Rceattle::plot_selectivity(fit, species = nm)$data$Fleet,
    Rceattle::plot_selectivity(fit, species = 2)$data$Fleet)
  # Selecting no fleets is an error, not an empty panel.
  testthat::expect_error(Rceattle::plot_selectivity(fit, species = 99), "out of range")
})

testthat::test_that("minyr and maxyr thin the year fan", {
  fit <- ss_fit()
  all_yrs <- unique(Rceattle::plot_selectivity(fit)$data$Year)
  cut <- unique(Rceattle::plot_selectivity(fit, minyr = max(all_yrs) - 5)$data$Year)
  testthat::expect_lt(length(cut), length(all_yrs))
  testthat::expect_gte(min(cut), max(all_yrs) - 5)
})

testthat::test_that("a length-based fleet is drawn on length bins, not ages", {
  # sel_at_age for a length-based fleet is the growth-matrix conversion of the
  # estimated curve, so plotting it under an "Age" label showed neither the
  # fitted quantity nor its dimension.
  lf <- len_fit("Length")
  fit <- lf$fit
  dl  <- lf$data_list
  testthat::expect_false(identical(dl$nlengths[1], dl$nages[1]))  # tells them apart

  p <- Rceattle::plot_selectivity(fit)
  testthat::expect_equal(p$labels$x, "Length bin")
  testthat::expect_equal(length(unique(p$data$Bin)), dl$nlengths[1])
  testthat::expect_setequal(unique(p$data$Dimension), "Length")

  # The values are sel_at_length, not sel_at_age.
  f1 <- p$data[p$data$Fleet == levels(p$data$Fleet)[1] &
                 p$data$Year == min(p$data$Year), ]
  testthat::expect_equal(
    f1$Selectivity[order(f1$Bin)],
    as.numeric(fit$quantities$sel_at_length[1, 1, seq_len(dl$nlengths[1]), 1]))
})

testthat::test_that("an age-based fleet is still drawn on ages", {
  lf <- len_fit("Age")
  fit <- lf$fit
  dl  <- lf$data_list
  p <- Rceattle::plot_selectivity(fit)
  testthat::expect_equal(p$labels$x, "Age")
  testthat::expect_equal(length(unique(p$data$Bin)), dl$nages[1])
  f1 <- p$data[p$data$Fleet == levels(p$data$Fleet)[1] &
                 p$data$Year == min(p$data$Year), ]
  testthat::expect_equal(
    f1$Selectivity[order(f1$Bin)],
    as.numeric(fit$quantities$sel_at_age[1, 1, seq_len(dl$nages[1]), 1]))
})

testthat::test_that("a model mixing dimensions returns one figure per dimension", {
  # Ages and length bins share no axis, so they must not land in one panel grid.
  lf <- len_fit(c("Length", "Age", "Length", "Age"))
  out <- Rceattle::plot_selectivity(lf$fit)
  testthat::expect_type(out, "list")
  testthat::expect_setequal(names(out), c("length", "age"))
  testthat::expect_equal(out$length$labels$x, "Length bin")
  testthat::expect_equal(out$age$labels$x, "Age")
  testthat::expect_equal(length(unique(out$length$data$Bin)), lf$data_list$nlengths[1])
  testthat::expect_equal(length(unique(out$age$data$Bin)), lf$data_list$nages[1])
})

testthat::test_that("the plotted curve is the model's own selectivity", {
  # Selectivity is not bounded at 1 unless the fleet normalizes (Sel_norm_bin),
  # so the check is against the reported array rather than against a range.
  fit <- ss_fit()
  d <- Rceattle::plot_selectivity(fit)$data
  testthat::expect_true(all(d$Selectivity >= 0, na.rm = TRUE))

  fc <- fit$data_list$fleet_control
  flt <- fc$Fleet_code[1]
  sp  <- fc$Species[1]
  nag <- fit$data_list$nages[sp]
  one <- d[d$Fleet == levels(d$Fleet)[1] & d$Year == min(d$Year), ]
  testthat::expect_equal(
    one$Selectivity[order(one$Bin)],
    as.numeric(fit$quantities$sel_at_age[flt, 1, seq_len(nag), 1]))
})
