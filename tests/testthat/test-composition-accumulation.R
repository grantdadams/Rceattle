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
#
# The final block pins that OSA residuals now SUPPORT an active accumulation:
# build_osa_data() folds the obsvec identically to the fit, so the OSA
# decomposition reproduces the folded composition likelihood.

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

# Oracle for the Multinomial (case 0) and Dirichlet-multinomial (case 1)
# composition NLL for one fleet, folding the tails per sex block from the model's
# own comp_obs / comp_hat. dmultinom_osa() with all bins kept is algebraically the
# full multinomial (chain rule), so the case-0 kernel is the plain lgamma
# multinomial; case 1 is ddirmultinom() with alpha = Neff * (hat + off) *
# exp(comp_weights). Independent of the C++ loop.
.comp_fold_nll <- function(f, fleet, young, kind) {
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
    nblk <- if (sex == 3) 2L else 1L
    nb   <- if (ctype == 0) f$data_list$nages[sp] else f$data_list$nlengths[sp]
    yng <- young; if (is.na(yng) || yng < 1) yng <- 1L
    tgt <- pmin(pmax(seq_len(nb), yng), nb)
    of <- numeric(0); hf <- numeric(0)
    for (b in 0:(nblk - 1)) {
      of <- c(of, tapply(obs[i, b * nb + seq_len(nb)], tgt, sum))
      hf <- c(hf, tapply(hat[i, b * nb + seq_len(nb)], tgt, sum))
    }
    if (kind == "mult") {
      x <- (of + off) * cn[i]; p <- (hf + off) / sum(hf + off)
      nll <- nll - cw * (lgamma(sum(x) + 1) - sum(lgamma(x + 1)) + sum(x * log(p)))
    } else {                                   # Dirichlet-multinomial
      # alpha's leading N is the offset-inflated observed total sum(ob), the same N
      # the density uses, matching the C++ `comp_obs_tmp.sum()` (not raw cn[i]).
      ob <- (of + off) * cn[i]; al <- sum(ob) * (hf + off) * exp(cw)
      nll <- nll - (lgamma(sum(ob) + 1) + lgamma(sum(al)) - lgamma(sum(ob) + sum(al)) +
                    sum(-lgamma(ob + 1) + lgamma(ob + al) - lgamma(al)))
    }
  }
  nll
}

testthat::test_that("Multinomial (case 0) and Dirichlet-multinomial (case 1) fold the exact hand-computed likelihood", {
  testthat::skip_on_cran()
  # DM fold drifts ~1e-5 under covr's -O0 instrumentation vs the -O2 build.
  testthat::skip_on_covr()
  testthat::skip_if_not_installed("TMB")

  # Pin the folded value (not just "it changed") against an independent oracle,
  # for both densities that read the folded numbers vector.
  d0 <- make_test_data(); d0$fleet_control$Comp_accum_young <- 3L
  f0 <- .build_accum(d0)
  for (flt in unique(f0$obj$env$data$comp_ctl[, 1]))
    testthat::expect_equal(f0$quantities$jnll_comp["Composition data", flt],
                           .comp_fold_nll(f0, flt, 3L, "mult"),
                           tolerance = 1e-8, info = paste("mult fleet", flt))

  d1 <- make_test_data()
  d1$fleet_control$Comp_distribution <- "DirichletMultinomial"
  d1$fleet_control$Comp_accum_young  <- 3L
  f1 <- .build_accum(d1)
  for (flt in unique(f1$obj$env$data$comp_ctl[, 1]))
    testthat::expect_equal(f1$quantities$jnll_comp["Composition data", flt],
                           .comp_fold_nll(f1, flt, 3L, "dm"),
                           tolerance = 1e-8, info = paste("DM fleet", flt))
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
  testthat::expect_error(data_check(d3), "must be < Comp_accum_old")

  d3b <- d                                      # young == old folds to a single bin
  d3b$fleet_control$Comp_accum_young <- 3L
  d3b$fleet_control$Comp_accum_old   <- 3L
  testthat::expect_error(data_check(d3b), "must be < Comp_accum_old")

  # Valid values and the no-accum sentinels pass (data_check may message; the
  # point is that it does not raise).
  quiet_check <- function(x) suppressMessages(suppressWarnings(data_check(x)))
  d4 <- d; d4$fleet_control$Comp_accum_young <- 3L; d4$fleet_control$Comp_accum_old <- 5L
  testthat::expect_no_error(quiet_check(d4))
  d5 <- d; d5$fleet_control$Comp_accum_old <- 0L    # 0 = no old accumulation
  testthat::expect_no_error(quiet_check(d5))
})

# Independent oracle: the folded obsvec bin counts for one comp row, built the
# same way the C++ fitting path does -- fold the RAW proportions per sex block
# into [young, old], then (folded + offset) * Neff. Returns a list per comp
# data_row of the expected obsvec counts, so a test can compare against the real
# obsvec segments build_osa_data() lays down. tapply here is independent of the
# package's .fold_comp_bins() helper.
.osa_fold_counts <- function(dat, fleet, young, old = 0L) {
  ctl <- dat$comp_ctl; cn <- dat$comp_n[, 2]; off <- dat$comp_offset
  endyr <- dat$endyr; flt_type <- dat$flt_type
  out <- list()
  for (i in seq_len(nrow(ctl))) {
    if (ctl[i, 1] != fleet) next
    yr <- ctl[i, 5]
    if (!(yr <= endyr && yr > 0 && flt_type[fleet] > 0 && cn[i] > 0)) next
    sp <- ctl[i, 2]; sex <- ctl[i, 3]; ctype <- ctl[i, 4]
    nblk <- if (sex == 3) 2L else 1L
    nb   <- if (ctype == 0) dat$nages[sp] else dat$nlengths[sp]
    yng <- young; if (is.na(yng) || yng < 1) yng <- 1L
    o_  <- old;   if (is.na(o_) || o_ < 1 || o_ > nb) o_ <- nb
    if (yng > o_) yng <- o_
    tgt <- pmin(pmax(seq_len(nb), yng), o_)             # <yng -> yng, >old -> old
    counts <- numeric(0)
    for (b in 0:(nblk - 1)) {
      folded <- tapply(dat$comp_obs[i, b * nb + seq_len(nb)], tgt, sum)  # raw props
      counts <- c(counts, (folded + off) * cn[i])       # single offset per folded bin
    }
    out[[as.character(i)]] <- unname(counts)
  }
  out
}

testthat::test_that("OSA build folds the comp obsvec for an active accumulation (case 0)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # build_osa_data() now supports an active fold: the obsvec is laid down over
  # the folded (nblk * nkeep) bins, not the full bins, and matches the C++
  # fitting fold count-for-count. 5-age Multinomial comps, young = 3 -> keep bins
  # 3:5 (nkeep = 3).
  d <- make_test_data()
  d$fleet_control$Comp_accum_young <- 3L
  f <- .build_accum(d)
  dat <- f$obj$env$data
  osa <- suppressWarnings(build_osa_data(dat, build_osa = TRUE))   # no longer errors
  oc  <- osa$obs_ctl

  for (flt in unique(dat$comp_ctl[, 1])) {
    expected <- .osa_fold_counts(dat, flt, young = 3L)
    for (i in names(expected)) {
      seg <- oc[oc$source == "comp" & oc$data_row == as.integer(i), ]
      seg <- seg[order(seg$obs_pos), ]
      got <- osa$obsvec[seg$obs_pos + 1L]
      testthat::expect_equal(length(got), length(expected[[i]]),
                             info = paste("fleet", flt, "row", i))   # folded count
      testthat::expect_equal(got, expected[[i]], tolerance = 1e-10,
                             info = paste("fleet", flt, "row", i))   # folded values
      # Exactly one sum-to-N last bin per composition group.
      testthat::expect_equal(sum(seg$is_last_bin), 1L)
    }
  }
})

testthat::test_that("OSA mode reproduces the folded comp fitting likelihood (decomposition invariant)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # The crux: with keep == 1 the conditional-binomial decomposition (dmultinom_osa)
  # sums back to the joint multinomial, so building the model in OSA mode
  # (osa_mode = 1, reading the FOLDED obsvec) must reproduce the composition
  # likelihood slot of ordinary fitting (osa_mode = 0, folding comp_obs) WITH an
  # active accumulation. This proves the R obsvec fold and the un-gated C++ fold
  # agree bin-for-bin. Multinomial (case 0), which is a proper density; the AFSC
  # pseudo-likelihood (-1) has no density and is residualized under the multinomial.
  d <- make_test_data()
  d$fleet_control$Comp_distribution <- "Multinomial"
  d$fleet_control$Comp_accum_young  <- 3L
  fit <- .build_accum(d)
  testthat::expect_true(all(fit$obj$env$data$comp_ll_type == 0L))

  comp_row <- 3L                       # jnll_comp C++ slot 2 (composition) -> R row 3
  jc0  <- fit$obj$report(fit$obj$par)$jnll_comp                 # osa_mode = 0 (folded fit)
  obj1 <- Rceattle:::.osa_build_obj(fit)                        # rebuild osa_mode = 1
  jc1  <- obj1$report(obj1$par)$jnll_comp

  testthat::expect_true(sum(jc0[comp_row, ]) != 0)             # composition is active
  testthat::expect_equal(sum(jc1[comp_row, ]), sum(jc0[comp_row, ]),
                         tolerance = 1e-8)                      # OSA reproduces the fold
})

testthat::test_that("OSA build folds an active OLD tail (young + old window)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Pin the old-tail fold on the OSA path (the young-only tests above leave an
  # active Comp_accum_old unverified for value). 5-age Multinomial, keep bins 2:4
  # (young = 2, old = 4 -> nkeep = 3): the obsvec matches a hand fold and the
  # decomposition invariant holds.
  d <- make_test_data()
  d$fleet_control$Comp_distribution <- "Multinomial"
  d$fleet_control$Comp_accum_young  <- 2L
  d$fleet_control$Comp_accum_old    <- 4L
  fit <- .build_accum(d)
  dat <- fit$obj$env$data
  osa <- suppressWarnings(build_osa_data(dat, build_osa = TRUE))
  oc  <- osa$obs_ctl

  for (flt in unique(dat$comp_ctl[, 1])) {
    expected <- .osa_fold_counts(dat, flt, young = 2L, old = 4L)
    for (i in names(expected)) {
      seg <- oc[oc$source == "comp" & oc$data_row == as.integer(i), ]
      seg <- seg[order(seg$obs_pos), ]
      testthat::expect_equal(nrow(seg), 3L, info = paste("fleet", flt, "row", i))
      testthat::expect_equal(osa$obsvec[seg$obs_pos + 1L], expected[[i]],
                             tolerance = 1e-10, info = paste("fleet", flt, "row", i))
    }
  }

  comp_row <- 3L
  jc0  <- fit$obj$report(fit$obj$par)$jnll_comp
  obj1 <- Rceattle:::.osa_build_obj(fit)
  jc1  <- obj1$report(obj1$par)$jnll_comp
  testthat::expect_equal(sum(jc1[comp_row, ]), sum(jc0[comp_row, ]), tolerance = 1e-8)
})

testthat::test_that("OSA mode reproduces the folded Dirichlet-multinomial (case 1)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # The DM OSA path (ddirmultinom_osa on the folded alphas) must reproduce the
  # folded DM fitting likelihood. Neither branch multiplies by comp_weights (the
  # DM par is exp(comp_weights) inside alphas), so the invariant holds directly.
  d <- make_test_data()
  d$fleet_control$Comp_distribution <- "DirichletMultinomial"
  d$fleet_control$Comp_accum_young  <- 3L
  fit <- .build_accum(d)
  testthat::expect_true(all(fit$obj$env$data$comp_ll_type == 1L))

  comp_row <- 3L
  jc0  <- fit$obj$report(fit$obj$par)$jnll_comp
  obj1 <- Rceattle:::.osa_build_obj(fit)
  jc1  <- obj1$report(obj1$par)$jnll_comp
  testthat::expect_true(sum(jc0[comp_row, ]) != 0)
  testthat::expect_equal(sum(jc1[comp_row, ]), sum(jc0[comp_row, ]), tolerance = 1e-8)
})

testthat::test_that("OSA build folds per sex block for joint-sex (Sex 3) comps", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not(exists("GOA2018SS"))

  # GOA2018SS fleet 9 (ATF, species 2) carries a joint-sex (Sex = 3) age comp: a
  # 42-entry, two-block row (21 ages x 2 sexes). With young = 5 each sex block
  # keeps bins 5:21 (nkeep = 17), so the folded obsvec is 2 * 17 = 34 bins, folded
  # within each block. Force Multinomial (case 0) so the decomposition invariant
  # also holds; the obsvec fold is likelihood-type-independent, so this covers the
  # MultinomialAFSC block-major layout too.
  d <- GOA2018SS
  is9 <- d$fleet_control$Fleet_code == 9
  d$fleet_control$Comp_distribution[is9] <- "Multinomial"
  d$fleet_control$Comp_accum_young[is9]  <- 5L
  d$fleet_control$Comp_weights[is9]      <- 1   # OSA density is unweighted; unit
                                                # weight makes it equal the fit slot
  fit <- .build_accum(d)
  dat <- fit$obj$env$data
  osa <- suppressWarnings(build_osa_data(dat, build_osa = TRUE))
  oc  <- osa$obs_ctl

  expected <- .osa_fold_counts(dat, 9L, young = 5L)
  testthat::skip_if(length(expected) == 0, "no included fleet-9 comps")
  for (i in names(expected)) {
    seg <- oc[oc$source == "comp" & oc$data_row == as.integer(i), ]
    seg <- seg[order(seg$obs_pos), ]
    testthat::expect_equal(nrow(seg), 34L, info = paste("row", i))   # 2 blocks x 17
    testthat::expect_equal(osa$obsvec[seg$obs_pos + 1L], expected[[i]],
                           tolerance = 1e-8, info = paste("row", i)) # per-block fold
  }

  # Decomposition invariant on the joint-sex fold. Compare only fleet 9's column
  # (jnll_comp columns are fleets): the other GOA fleets keep the default
  # MultinomialAFSC (case -1), whose OSA value (residualized under the multinomial)
  # legitimately differs from its AFSC pseudo-likelihood, so only the Multinomial
  # fleet 9 satisfies osa_mode 1 == osa_mode 0.
  comp_row <- 3L
  jc0  <- fit$obj$report(fit$obj$par)$jnll_comp
  obj1 <- Rceattle:::.osa_build_obj(fit)
  jc1  <- obj1$report(obj1$par)$jnll_comp
  testthat::expect_true(jc0[comp_row, 9] != 0)                       # fleet 9 active
  testthat::expect_equal(jc1[comp_row, 9], jc0[comp_row, 9], tolerance = 1e-8)
})

testthat::test_that("osa_residuals() runs end-to-end with an active accumulation", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # End-to-end: a converging fit (estimateMode = 1) with an active fold yields
  # finite, correctly-counted comp residuals -- one fewer than the folded bin
  # count per group (the sum-to-N last bin is dropped, residual NA).
  d <- make_test_data()
  d$fleet_control$Comp_distribution <- "Multinomial"
  d$fleet_control$Comp_accum_young  <- 3L
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, file = NULL, estimateMode = 1, random_rec = FALSE, msmMode = 0,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))

  osa <- osa_residuals(fit, source = "comp")
  testthat::expect_s3_class(osa, "rceattle_osa")
  testthat::expect_true(all(is.finite(osa$residual)))

  ctl <- Rceattle:::build_osa_data(fit$obj$env$data, build_osa = TRUE)$obs_ctl
  n_comp_bins   <- sum(ctl$source == "comp")
  n_comp_groups <- length(unique(ctl$group_id[ctl$source == "comp"]))
  testthat::expect_equal(nrow(osa), n_comp_bins - n_comp_groups)   # one dropped / group
})
