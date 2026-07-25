# Composition young/old tail-accumulation (AFSC ac_yng / ac_old). The
# Comp_accum_young / Comp_accum_old fleet_control columns fold the age/length
# composition tails into a boundary bin before the composition likelihood.
#
# These tests pin, for BOTH composition likelihoods that read the folded comps:
#   * MultinomialAFSC (comp_ll_type == -1, the package default) -- proportions
#     read directly; verified against a hand-computed folded NLL (the oracle
#     below), for combined-sex AND joint-sex (Sex = 3) blocks;
#   * Multinomial (comp_ll_type == 0) -- folding changes the likelihood;
# plus the invariant that keeps every existing model unchanged: the
# no-accumulation configuration is bit-identical to omitting the columns.

# Composition-data likelihood component for a build-only (estimateMode = 3) fit.
.comp_jnll <- function(dl) {
  sum(.build_accum(dl)$quantities$jnll_comp["Composition data", ], na.rm = TRUE)
}

.build_accum <- function(dl) {
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    dl, file = NULL, estimateMode = 3, random_rec = FALSE, msmMode = 0,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))
}

# Oracle: the MultinomialAFSC (case -1) composition NLL for one fleet, folding
# the young/old tails into [young, old] PER SEX BLOCK, computed in R straight
# from the model's own comp_obs / comp_hat. This mirrors the template exactly
# (jnll -= comp_weights * n * (obs + off) * log((hat + off)/(obs + off))), so an
# equal value proves the C++ fold folds the right bins for case -1 -- the path
# that reads proportions directly rather than the numbers vector.
.afsc_fold_nll <- function(f, fleet, young, old = 0L) {
  dat <- f$obj$env$data
  obs <- f$quantities$comp_obs; hat <- f$quantities$comp_hat
  ctl <- dat$comp_ctl; cn <- dat$comp_n[, 2]
  cw  <- f$estimated_params$comp_weights[fleet]
  off <- dat$comp_offset; endyr <- f$data_list$endyr; flt_type <- dat$flt_type
  nll <- 0
  for (i in seq_len(nrow(ctl))) {
    if (ctl[i, 1] != fleet) next
    yr <- ctl[i, 5]
    if (!(yr <= endyr && yr > 0 && flt_type[fleet] > 0 && cn[i] > 0)) next
    sp <- ctl[i, 2]; sex <- ctl[i, 3]; ctype <- ctl[i, 4]
    nblk <- if (sex == 3) 2L else 1L                      # joint sex = two blocks
    nb   <- if (ctype == 0) f$data_list$nages[sp] else f$data_list$nlengths[sp]
    yng <- young; if (is.na(yng) || yng < 1) yng <- 1L
    o_  <- old;   if (is.na(o_) || o_ < 1 || o_ > nb) o_ <- nb
    if (yng > o_) yng <- o_
    tgt <- pmin(pmax(seq_len(nb), yng), o_)               # <yng -> yng, >old -> old
    for (b in 0:(nblk - 1)) {
      o <- obs[i, b * nb + seq_len(nb)]; h <- hat[i, b * nb + seq_len(nb)]
      of <- tapply(o, tgt, sum); hf <- tapply(h, tgt, sum)
      nll <- nll - sum(cw * cn[i] * (of + off) * log((hf + off) / (of + off)))
    }
  }
  nll
}

testthat::test_that("Multinomial (case 0) young/old accumulation folds the likelihood", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data()                      # 5 ages, Multinomial comps
  base <- .comp_jnll(d)
  testthat::expect_true(is.finite(base))

  d_y <- d; d_y$fleet_control$Comp_accum_young <- 3L
  testthat::expect_false(isTRUE(all.equal(.comp_jnll(d_y), base)))

  d_o <- d; d_o$fleet_control$Comp_accum_old <- 3L
  testthat::expect_false(isTRUE(all.equal(.comp_jnll(d_o), base)))
})

testthat::test_that("MultinomialAFSC (case -1) folds the exact hand-computed bins", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # The default likelihood (MultinomialAFSC) reads proportions directly, so it
  # must fold comp_obs/comp_hat the same as cases 0/1. Compare to the oracle.
  d <- make_test_data()
  d$fleet_control$Comp_distribution <- "MultinomialAFSC"
  d$fleet_control$Comp_accum_young  <- 3L
  f <- .build_accum(d)
  testthat::expect_true(all(f$obj$env$data$comp_ll_type == -1L))

  for (flt in unique(f$obj$env$data$comp_ctl[, 1])) {
    reported <- f$quantities$jnll_comp["Composition data", flt]
    testthat::expect_equal(reported, .afsc_fold_nll(f, flt, young = 3L),
                           tolerance = 1e-8, info = paste("fleet", flt))
  }
})

testthat::test_that("MultinomialAFSC (case -1) folds per sex block for joint-sex (Sex 3) comps", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not(exists("GOA2018SS"))

  # GOA2018SS fleet 9 is a joint-sex (Sex = 3) age composition for species 2
  # (21 ages -> a 42-entry, two-block joint row). Force it to the AFSC pseudo-
  # likelihood and accumulate the young tail; the fold must stay within each sex
  # block. One build supplies comp_obs/comp_hat for both the folded and the
  # unfolded oracle (comp_hat is param-driven, identical either way).
  d <- GOA2018SS
  is9 <- d$fleet_control$Fleet_code == 9
  d$fleet_control$Comp_distribution[is9] <- "MultinomialAFSC"
  d$fleet_control$Comp_accum_young[is9]  <- 5L
  f <- .build_accum(d)
  testthat::expect_equal(unname(f$obj$env$data$comp_ll_type[9]), -1)

  reported <- f$quantities$jnll_comp["Composition data", 9]
  folded   <- .afsc_fold_nll(f, 9, young = 5L)   # per-block fold, young 1:4 -> 5
  unfolded <- .afsc_fold_nll(f, 9, young = 1L)   # no fold

  testthat::expect_equal(reported, folded, tolerance = 1e-6)     # block-correct fold
  testthat::expect_false(isTRUE(all.equal(folded, unfolded)))    # fold is non-trivial
})

testthat::test_that("no accumulation is bit-identical to omitting the columns (both likelihoods)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  for (dist in c("Multinomial", "MultinomialAFSC")) {
    d <- make_test_data()
    d$fleet_control$Comp_distribution <- dist
    base <- .comp_jnll(d)                       # columns absent -> defaulted no-op

    d_n <- d                                    # explicit no-op: young = 1, old = nbins
    d_n$fleet_control$Comp_accum_young <- 1L
    d_n$fleet_control$Comp_accum_old   <- max(d$nages)
    testthat::expect_identical(.comp_jnll(d_n), base, info = dist)

    d_na <- d                                   # all-NA (what switch_check materializes)
    d_na$fleet_control$Comp_accum_young <- NA_integer_
    d_na$fleet_control$Comp_accum_old   <- NA_integer_
    testthat::expect_identical(.comp_jnll(d_na), base, info = dist)
  }
})

testthat::test_that("data_check rejects out-of-range or inverted accumulation bins", {
  # Fast unit check (no fit): the template lower/upper-clamps defensively, but a
  # nonsensical value must be caught up front rather than silently folded into a
  # boundary bin. 5-age age comps -> valid bins are 1:5. data_check() runs on the
  # data list after fit_mod() has stamped the run-level switches, so set the few
  # scalars its earlier switch checks read.
  d <- make_test_data()
  d$msmMode <- 0L
  d$HCR <- "NoFishing"

  d1 <- d; d1$fleet_control$Comp_accum_young <- 99L
  testthat::expect_error(data_check(d1), "Comp_accum_young")

  d2 <- d; d2$fleet_control$Comp_accum_old <- 99L
  testthat::expect_error(data_check(d2), "Comp_accum_old")

  d3 <- d                                       # young > old (with a real old bin)
  d3$fleet_control$Comp_accum_young <- 4L
  d3$fleet_control$Comp_accum_old   <- 2L
  testthat::expect_error(data_check(d3), "must be <= Comp_accum_old")

  # Valid values and the no-accum sentinels pass (data_check may message; the
  # point is that it does not raise).
  quiet_check <- function(x) suppressMessages(suppressWarnings(data_check(x)))
  d4 <- d; d4$fleet_control$Comp_accum_young <- 3L; d4$fleet_control$Comp_accum_old <- 5L
  testthat::expect_no_error(quiet_check(d4))
  d5 <- d; d5$fleet_control$Comp_accum_old <- 0L    # 0 = no old accumulation
  testthat::expect_no_error(quiet_check(d5))
})

testthat::test_that("OSA residuals refuse an active composition accumulation", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # OSA residuals are built on the full, un-accumulated bins, so combining them
  # with an active fold would residualize a different model than was fit --
  # build_osa_data() must refuse rather than silently disagree.
  d <- make_test_data()
  d$fleet_control$Comp_accum_young <- 3L
  f <- .build_accum(d)
  testthat::expect_error(build_osa_data(f$obj$env$data, build_osa = TRUE),
                         "tail accumulation")

  # Without accumulation the OSA build succeeds.
  d0 <- make_test_data()
  f0 <- .build_accum(d0)
  testthat::expect_no_error(build_osa_data(f0$obj$env$data, build_osa = TRUE))
})
