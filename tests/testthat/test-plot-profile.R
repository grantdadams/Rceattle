# Unit tests for profile_components() (R/9-profile.R) and plot_profile()
# (R/7-plot_profile.R).
#
# The profiles here are constructed rather than fitted. What these functions do
# is read `quantities$jnll_comp` and label its cells, and the labelling is the
# part that can be silently wrong: the matrix's columns count fleets on some
# rows and species on others, so a mislabelled cell would attribute a survey's
# likelihood to a species. Building the matrix by hand is the only way to know
# what each cell should come back as. A live profile's class and shape are
# asserted in test-functions-profile-param.R.

# Two species, three fleets: jnll_comp is dimensioned max(n_fleets, nspp), and
# the two axes have different lengths so a row read off the wrong one shows up.
FLEETS <- c("Fishery", "Shelikof acoustic", "NMFS bottom trawl")
SPP    <- c("Pollock", "Cod")

# One 21-row jnll_comp, laid out by component name so the rows stay tied to
# R/6-rename_output.R rather than to integer positions.
mk_jnll <- function(cells) {
  m <- matrix(0, nrow = length(.JNLL_ROW_AXIS), ncol = 3,
              dimnames = list(names(.JNLL_ROW_AXIS), NULL))
  for (nm in names(cells)) m[nm, ] <- cells[[nm]]
  m
}

mk_fit <- function(m, fleets = FLEETS, spp = SPP) {
  list(quantities = list(jnll_comp = m, unweighted_jnll_comp = m / 2),
       data_list = list(fleet_control = data.frame(Fleet_name = fleets,
                                                   stringsAsFactors = FALSE),
                        spnames = spp))
}

# A profile whose components disagree on purpose: the Shelikof index prefers
# x = 1, the bottom trawl prefers x = 5, and the fishery composition is flat.
# That is the conflict the figure exists to show.
mk_profile <- function(x = 1:5, drop = integer(0), nll = NULL,
                       fleets = FLEETS, spp = SPP) {
  fits <- lapply(x, function(v) {
    mk_fit(mk_jnll(list(
      "Index data"                 = c(0, (v - 1)^2, (v - 5)^2),
      "Composition data"           = c(1, 0, 0),
      "Recruitment deviates"       = c(0.01 * v, 0, 0)
    )), fleets, spp)
  })
  tot <- vapply(fits, function(f) sum(f$quantities$jnll_comp), numeric(1))
  for (i in drop) { fits[[i]] <- NULL; tot[i] <- NA_real_ }
  # Dropping by NULL assignment shortens the list, so rebuild it aligned.
  full <- vector("list", length(x))
  keep <- setdiff(seq_along(x), drop)
  full[keep] <- lapply(keep, function(i) {
    mk_fit(mk_jnll(list(
      "Index data"                 = c(0, (x[i] - 1)^2, (x[i] - 5)^2),
      "Composition data"           = c(1, 0, 0),
      "Recruitment deviates"       = c(0.01 * x[i], 0, 0)
    )), fleets, spp)
  })
  structure(list(Rceattle_list = full,
                 grid  = data.frame(slot_1 = x),
                 nll   = if (is.null(nll)) tot else nll,
                 param = "log_M1", slots = list(c(1, 1, 1)), alias = "M1"),
            class = "Rceattle_profile")
}


# --- profile_components: cell labelling ---------------------------------------
testthat::test_that("profile_components labels each cell from its own axis", {
  d <- profile_components(mk_profile())

  # Index data is a fleet row, so column 2 is the second FLEET.
  idx <- unique(d[d$component == "Index data", c("unit", "axis", "series")])
  testthat::expect_equal(sort(idx$unit),
                         sort(c("Shelikof acoustic", "NMFS bottom trawl")))
  testthat::expect_true(all(idx$axis == "fleet"))
  testthat::expect_true("Shelikof acoustic: Index data" %in% idx$series)

  # Recruitment deviations are a species row, so column 1 is the first SPECIES.
  rec <- unique(d[d$component == "Recruitment deviates",
                  c("unit", "axis", "series")])
  testthat::expect_equal(rec$unit, "Pollock")
  testthat::expect_equal(rec$axis, "species")
  testthat::expect_equal(as.character(rec$series), "Pollock: Recruitment deviates")

  # Fleet 1 of Index data is zero at every grid point: a component the model
  # does not fit, not a flat one, so it is dropped rather than drawn at zero.
  testthat::expect_false("Fishery: Index data" %in% levels(d$series))

  # Every retained cell is reported at every grid point.
  testthat::expect_equal(as.integer(table(d$series)[table(d$series) > 0]),
                         rep(5L, length(unique(d$series))))
})


testthat::test_that("profile_components drops the unit when there is only one", {
  # "Recruitment deviates" reads better than "Pollock: Recruitment deviates"
  # when Pollock is the only species there is. One fleet and one species means
  # jnll_comp has a single column, since it is dimensioned max(n_fleets, nspp).
  fits <- lapply(1:5, function(v) {
    m <- matrix(0, nrow = length(.JNLL_ROW_AXIS), ncol = 1,
                dimnames = list(names(.JNLL_ROW_AXIS), NULL))
    m["Index data", 1] <- (v - 3)^2
    m["Recruitment deviates", 1] <- 0.01 * v
    mk_fit(m, fleets = "Fishery", spp = "Pollock")
  })
  prof <- structure(
    list(Rceattle_list = fits, grid = data.frame(slot_1 = 1:5),
         nll = vapply(fits, function(f) sum(f$quantities$jnll_comp), numeric(1)),
         param = "log_M1", slots = list(c(1, 1, 1)), alias = "M1"),
    class = "Rceattle_profile")

  d <- profile_components(prof)
  testthat::expect_true("Recruitment deviates" %in% levels(d$series))
  testthat::expect_true("Index data" %in% levels(d$series))
  testthat::expect_false(any(grepl(":", levels(d$series), fixed = TRUE)))
})


testthat::test_that("profile_components names a row it does not recognise", {
  prof <- mk_profile()
  prof$Rceattle_list <- lapply(prof$Rceattle_list, function(f) {
    rn <- rownames(f$quantities$jnll_comp)
    rownames(f$quantities$jnll_comp) <- replace(
      rn, rn == "Recruitment deviates", "Some new likelihood")
    f
  })
  # A row added to the JnllRow enum and not to .JNLL_ROW_AXIS still plots, but
  # the caller is told its fleet or species is not named.
  testthat::expect_warning(d <- profile_components(prof),
                           "Some new likelihood")
  testthat::expect_true("Some new likelihood" %in% levels(d$series))
  testthat::expect_true(is.na(unique(d$unit[d$component == "Some new likelihood"])))
})


# --- profile_components: re-zeroing -------------------------------------------
testthat::test_that("relative = 'own' puts every series at its own minimum", {
  d <- profile_components(mk_profile())
  mins <- tapply(d$value, d$series, min, na.rm = TRUE)
  testthat::expect_true(all(abs(mins[is.finite(mins)]) < 1e-12))

  # The x at which each series bottoms out is the whole point of the figure:
  # the Shelikof index wants 1 and the bottom trawl wants 5.
  best_x <- function(s) d$slot_1[d$series == s][which.min(d$value[d$series == s])]
  testthat::expect_equal(best_x("Shelikof acoustic: Index data"), 1)
  testthat::expect_equal(best_x("NMFS bottom trawl: Index data"), 5)
})


testthat::test_that("relative = 'minimum' zeroes every series at the total's minimum", {
  prof <- mk_profile()
  d <- profile_components(prof, relative = "minimum")
  best <- which.min(prof$nll)
  at_best <- vapply(split(d, d$series, drop = TRUE),
                    function(g) g$value[g$fit == best], numeric(1))
  testthat::expect_true(all(abs(at_best) < 1e-12))
})


testthat::test_that("relative = 'none' returns the raw likelihoods", {
  prof <- mk_profile()
  d <- profile_components(prof, relative = "none")
  testthat::expect_equal(d$value[d$series == "Total"], prof$nll)
  # With no random effects the components do sum to the objective.
  comp <- d[d$series != "Total", ]
  testthat::expect_equal(as.numeric(tapply(comp$value, comp$fit, sum)),
                         prof$nll)
})


testthat::test_that("Total is the objective, and a Laplace gap is reported", {
  # Under random_rec the objective is the marginal likelihood while the
  # components are the inner joint one, so they do not sum. Reported rather
  # than reconciled: the shapes are comparable, the sums are not.
  prof <- mk_profile()
  prof$nll <- prof$nll + 5
  testthat::expect_warning(profile_components(prof),
                           "Laplace-approximated marginal")
})


# --- profile_components: filtering, ordering, failed fits ---------------------
testthat::test_that("minfraction drops flat components but never the total", {
  prof <- mk_profile()
  all_comp <- profile_components(prof)
  # Recruitment deviations move 0.04 over the grid against a total that moves
  # 16, i.e. 0.25% -- below a 1% threshold and above a 0.1% one.
  kept_1  <- levels(droplevels(profile_components(prof, minfraction = 0.01)$series))
  kept_01 <- levels(droplevels(profile_components(prof, minfraction = 0.001)$series))
  testthat::expect_false("Pollock: Recruitment deviates" %in% kept_1)
  testthat::expect_true("Pollock: Recruitment deviates" %in% kept_01)
  testthat::expect_true("Total" %in% kept_1)
  testthat::expect_lt(length(kept_1), length(levels(all_comp$series)))
})


testthat::test_that("series are ordered by how much they move, total first", {
  # Legend order is read as importance, so it is the order the components move
  # in -- with the total pinned first, as the reference the rest are read
  # against. Last is the fishery composition, which is constant over the grid.
  lv <- levels(profile_components(mk_profile())$series)
  testthat::expect_equal(lv[1], "Total")
  # The two index series move by the same amount, so their order between
  # themselves is a tie and not asserted.
  testthat::expect_setequal(lv[2:3], c("NMFS bottom trawl: Index data",
                                       "Shelikof acoustic: Index data"))
  testthat::expect_equal(lv[length(lv)], "Fishery: Composition data")
})


testthat::test_that("a failed grid point becomes NA rather than a closed gap", {
  # Dropping the row instead would draw a straight segment across the failure,
  # which reads as a likelihood that is smooth there.
  d <- profile_components(mk_profile(drop = 3L))
  testthat::expect_equal(sum(is.na(d$value[d$series != "Total"])),
                         length(unique(d$series[d$series != "Total"])))
  testthat::expect_true(all(is.na(d$value[d$fit == 3L])))
  # Rows are still present, so `fit` stays aligned with the grid.
  testthat::expect_true(all(d$slot_1[d$fit == 3L] == 3))
})


testthat::test_that("profile_components rejects what it cannot read", {
  testthat::expect_error(profile_components(list()), "Rceattle_profile")
  testthat::expect_error(profile_components(mk_profile(), minfraction = -1),
                         "non-negative")
  dead <- mk_profile(drop = 1:5)
  testthat::expect_error(profile_components(dead), "grid point in this profile")
})


testthat::test_that("weighted = FALSE reads the unweighted matrix", {
  prof <- mk_profile()
  w  <- profile_components(prof, relative = "none", include_total = FALSE)
  uw <- profile_components(prof, relative = "none", include_total = FALSE,
                           weighted = FALSE)
  testthat::expect_equal(uw$value, w$value / 2)
})


# --- plot_profile -------------------------------------------------------------
testthat::test_that("plot_profile draws the components coloured and the total black", {
  p <- plot_profile(mk_profile())

  # The total is its own black layer, not a level of the colour scale, so the
  # component palette does not shift when minfraction changes what is shown.
  testthat::expect_false("Total" %in% as.character(p$data$series))
  b <- ggplot2::ggplot_build(p)
  cols <- unlist(lapply(b$data, function(d) unique(d$colour)))
  testthat::expect_true("black" %in% cols)

  # The axis is named in the units the grid is in: the alias, not the log-scale
  # slot the model estimates.
  testthat::expect_equal(p$labels$x, "M1[1, 1, 1]")
  testthat::expect_equal(p$labels$y, "Change in negative log-likelihood")
})


testthat::test_that("plot_profile marks each series' minimum", {
  p <- plot_profile(mk_profile(), minfraction = 0)
  b <- ggplot2::ggplot_build(p)
  pts <- do.call(rbind, lapply(b$data, function(d) {
    if (!is.null(d$shape) && nrow(d)) d[, c("x", "y")] else NULL
  }))
  # One point per component plus one for the total, each at that series' own
  # minimum, which relative = "own" has put at zero.
  testthat::expect_true(all(abs(pts$y) < 1e-8))
  testthat::expect_equal(nrow(pts),
                         length(unique(profile_components(mk_profile())$series)))
})


testthat::test_that("plot_profile facets a list of profiles and checks they match", {
  p <- plot_profile(list(mk_profile(), mk_profile()),
                    model_names = c("Base", "Alternative"))
  testthat::expect_equal(levels(p$data$Model), c("Base", "Alternative"))

  other <- mk_profile()
  other$param <- "R_log_sd"
  testthat::expect_error(plot_profile(list(mk_profile(), other)),
                         "different parameters")
})


testthat::test_that("two fleets sharing a name stay two series", {
  # Only Fleet_code is required to be unique, so a workbook can name two fleets
  # alike. Merged, their two curves would be drawn as one line through
  # interleaved values -- a figure that is wrong without looking wrong.
  prof <- mk_profile(fleets = c("Survey", "Survey", "Survey"))
  d <- profile_components(prof)
  idx <- unique(as.character(d$series[d$component == "Index data"]))
  testthat::expect_length(idx, 2L)
  testthat::expect_true(all(grepl("fleet [23]", idx)))
  # Each still carries its own five grid points, not ten interleaved ones.
  testthat::expect_equal(as.integer(table(droplevels(d$series[
    d$component == "Index data"]))), c(5L, 5L))
})


testthat::test_that("plot_profile refuses to overlay two scales", {
  # profile(param = "M1") and profile(param = "log_M1") both report
  # param = "log_M1", so the parameter check alone would pass while the grids
  # hold M and log M.
  natural <- mk_profile()
  logged  <- mk_profile()
  logged$alias <- NA_character_
  logged$grid$slot_1 <- log(logged$grid$slot_1)
  testthat::expect_error(plot_profile(list(natural, logged)),
                         "not on the same scale")
})


testthat::test_that("plot_profile says when the profiles are of different cells", {
  # Species 1 against species 2 is a fair comparison, but only the first one's
  # cell reaches the axis label.
  other <- mk_profile()
  other$slots <- list(c(2, 1, 1))
  testthat::expect_warning(plot_profile(list(mk_profile(), other)),
                           "different cells")
})


testthat::test_that("plot_profile refuses a cross-profile", {
  # Two profiled cells is a surface, not a line, and silently plotting one of
  # them against the pooled NLL would draw a curve that is not a profile.
  cross <- mk_profile()
  cross$grid <- expand.grid(slot_1 = 1:5, slot_2 = 1:5)[1:5, ]
  testthat::expect_error(plot_profile(cross), "cross-profile")
})


testthat::test_that("plot_profile honours line_col, lwd and the cutoff line", {
  p <- plot_profile(mk_profile(), minfraction = 0,
                    line_col = c("red", "blue", "green", "purple"),
                    add_cutoff = TRUE)
  b <- ggplot2::ggplot_build(p)
  cols <- unlist(lapply(b$data, function(d) unique(d$colour)))
  testthat::expect_true(all(c("red", "blue") %in% cols))
  # The cutoff is drawn where it was asked for.
  ys <- unlist(lapply(b$data, function(d) d$yintercept))
  testthat::expect_true(1.92 %in% ys)

  # lwd is on the base-graphics scale, as everywhere else in the package.
  p2 <- plot_profile(mk_profile(), lwd = 3)
  lw <- unlist(lapply(ggplot2::ggplot_build(p2)$data,
                      function(d) unique(d$linewidth)))
  testthat::expect_true(1 %in% lw)
})


testthat::test_that("plot_profile says when minfraction has emptied the figure", {
  # Drawing the total alone is a blank answer to the question the figure asks,
  # so the argument that emptied it is named rather than left to be guessed.
  testthat::expect_warning(plot_profile(mk_profile(), minfraction = 2),
                           "minfraction")
})


testthat::test_that("relative = 'minimum' needs a minimum to re-zero at", {
  prof <- mk_profile()
  prof$nll <- rep(NA_real_, 5)
  testthat::expect_warning(d <- profile_components(prof, relative = "minimum"),
                           "no minimum to re-zero at")
  # Raw values, not a crash: the components are still readable against
  # each other.
  testthat::expect_false(all(is.na(d$value[d$series != "Total"])))
})


testthat::test_that("plot() on a profile is plot_profile()", {
  testthat::expect_equal(plot(mk_profile())$labels$x,
                         plot_profile(mk_profile())$labels$x)
})
