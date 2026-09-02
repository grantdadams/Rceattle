# parameter_index() locates every estimated parameter in the model's own
# coordinates, so a convergence diagnostic can name the fleet, bin and year it
# means rather than repeating the parameter block.
#
# Provenance: the estimability check reported "(sel_coff_dev, log_sel_slp,
# sel_inf)" for 49 flagged parameters, which does not say where to look. The
# recovery goes through TMB's own parList(), so it tracks build_map() exactly;
# these tests pin that against obj$par rather than against a re-derivation.

testthat::test_that("the dictionary's declared rank matches the built array", {
  testthat::skip_on_cran()
  # The dictionary names what each dimension MEANS; parameter_index() reads the
  # extent off the array. That only works if the two agree on how many
  # dimensions a block has. M1_dev_log_sd was declared [nspp, nsex, 2] against a
  # [nspp, nsex] array, which is the drift this catches.
  fit <- Rceattle::fit_mod(data_list = Rceattle::BS2017SS, estimateMode = 3,
                           msmMode = 0,
                           fit_control = Rceattle::fit_control(getsd = FALSE,
                                                               verbose = 0))
  pl   <- fit$obj$env$parList(x = seq_along(fit$obj$par))
  dict <- Rceattle::parameter_dictionary()

  strip <- function(s) gsub("]", "", gsub("[", "", s, fixed = TRUE), fixed = TRUE)
  for (block in names(pl)) {
    dd <- dict$dims[match(block, dict$internal)]
    testthat::expect_false(is.na(dd),
      info = paste(block, "is not in parameter_dictionary()"))
    toks <- trimws(strsplit(strip(dd), ",")[[1]])
    toks <- toks[nzchar(toks)]
    actual <- dim(pl[[block]])
    if (is.null(actual)) actual <- length(pl[[block]])
    testthat::expect_equal(length(toks), length(actual),
      info = paste0(block, ": dims '", dd, "' declares ", length(toks),
                    " dimension(s), array has ", length(actual)))
  }
})


testthat::test_that("every estimated parameter is located exactly once", {
  testthat::skip_on_cran()
  fit <- Rceattle::fit_mod(data_list = Rceattle::BS2017SS, estimateMode = 1,
                           msmMode = 0,
                           fit_control = Rceattle::fit_control(getsd = FALSE,
                                                               verbose = 0,
                                                               phase = FALSE))
  idx <- Rceattle::parameter_index(fit)

  # One row per element of obj$par, and the blocks agree with TMB's own names.
  testthat::expect_equal(nrow(idx), length(fit$obj$par))
  testthat::expect_equal(sort(idx$par_index), seq_along(fit$obj$par))
  testthat::expect_equal(idx$block[order(idx$par_index)],
                         unname(names(fit$obj$par)))
})


testthat::test_that("coordinates are the model's own labels, not indices", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  fit <- Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
                           fit_control = Rceattle::fit_control(getsd = FALSE,
                                                               verbose = 0,
                                                               phase = FALSE))
  idx <- Rceattle::parameter_index(fit)

  # Fleets by name, years as calendar years -- an index would be no better than
  # the block name this replaces.
  f <- idx[idx$block == "log_F", ]
  testthat::expect_true(all(f$fleet %in% d$fleet_control$Fleet_name))
  yr <- as.numeric(stats::na.omit(f$year))
  testthat::expect_true(all(yr >= d$styr & yr <= d$projyr))

  # rec_dev is indexed by species, and BS2017SS names three.
  r <- idx[idx$block == "rec_dev", ]
  testthat::expect_true(all(r$species %in% d$spnames))
})


testthat::test_that("a mirrored parameter is one row naming every fleet it drives", {
  testthat::skip_on_cran()
  # Fleets sharing a Selectivity_index share ONE parameter block, so reporting
  # them as several parameters would overstate how many things are unidentified.
  d <- Rceattle::BS2017SS
  fit <- Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
                           fit_control = Rceattle::fit_control(getsd = FALSE,
                                                               verbose = 0,
                                                               phase = FALSE))
  idx <- Rceattle::parameter_index(fit)

  testthat::expect_true(all(idx$n_cells >= 1))
  # However many cells a parameter drives, it is still one parameter.
  testthat::expect_equal(sum(duplicated(idx$par_index)), 0L)
  mirrored <- idx[idx$n_cells > 1, ]
  if (nrow(mirrored) > 0) {
    testthat::expect_true(all(grepl(" \\+ ", mirrored$fleet) |
                              grepl(" \\+ ", mirrored$species)))
  }
})


testthat::test_that("an axis the model does not distinguish is not printed", {
  testthat::skip_on_cran()
  # GOApollock is one species and one sex, so naming either on every line says
  # nothing. The structured column still carries the species.
  d <- Rceattle::GOApollock
  fit <- Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
                           fit_control = Rceattle::fit_control(getsd = FALSE,
                                                               verbose = 0,
                                                               phase = FALSE))
  idx <- Rceattle::parameter_index(fit)

  testthat::expect_equal(d$nsex, 1)
  testthat::expect_true(all(is.na(idx$sex)))          # single-sex: no sex exists
  testthat::expect_false(any(grepl("combined", idx$label)))
  testthat::expect_false(any(grepl(d$spnames[1], idx$label)))  # constant, so unprinted
  testthat::expect_true(any(idx$species == d$spnames[1], na.rm = TRUE))
})


testthat::test_that("the summary collapses ordinal axes and counts the group", {
  testthat::skip_on_cran()
  fit <- Rceattle::fit_mod(data_list = Rceattle::GOApollock, estimateMode = 1,
                           msmMode = 0,
                           fit_control = Rceattle::fit_control(getsd = FALSE,
                                                               verbose = 0,
                                                               phase = FALSE))
  idx <- Rceattle::parameter_index(fit)
  lf  <- idx$par_index[idx$block == "log_F"]
  out <- Rceattle:::.rce_par_summary(lf, idx)

  testthat::expect_length(out, 1L)
  # Years are ordinal, so they read as a range rather than 49 values.
  testthat::expect_match(out, "log_F")
  testthat::expect_match(out, "\\d{4}-\\d{4}")
  testthat::expect_match(out, sprintf("\\(%d\\)$", length(lf)))
})
