# sim_mod(process = ) -- redrawing process error alongside the observations.
# The fast tests cover the switch translation; the fits cover the draw itself.
# Moment checks (does an ar1 draw have the right marginal sd?) live in
# tools/verify/verify-sim-linkage-re.R, which can afford hundreds of replicates.

# ---- switch translation (fast, no fit) --------------------------------------

testthat::test_that("sim_mod() process switches translate to state codes", {
  none <- rep(0L, 5); all1 <- rep(1L, 5)
  testthat::expect_identical(Rceattle:::.sim_state_codes(FALSE),  none)
  testthat::expect_identical(Rceattle:::.sim_state_codes("none"), none)
  testthat::expect_identical(Rceattle:::.sim_state_codes(NULL),   none)
  testthat::expect_identical(Rceattle:::.sim_state_codes(TRUE),   all1)
  testthat::expect_identical(Rceattle:::.sim_state_codes("all"),  all1)

  # dynamics = recruitment, M, growth; observation = catchability, selectivity
  testthat::expect_identical(Rceattle:::.sim_state_codes("dynamics"),
                             c(1L, 1L, 1L, 0L, 0L))
  testthat::expect_identical(Rceattle:::.sim_state_codes("observation"),
                             c(0L, 0L, 0L, 1L, 1L))
  testthat::expect_identical(Rceattle:::.sim_state_codes(c("recruitment", "growth")),
                             c(1L, 0L, 1L, 0L, 0L))
  testthat::expect_error(Rceattle:::.sim_state_codes("mortality"),
                         "unknown process")
})

testthat::test_that("simulate_state slots are the linkage process codes", {
  # ceattle.cpp gates a random linkage with simulate_state(linkage_process(i)),
  # with no translation between the two, so the orders have to agree. They were
  # crossed at 3/4 (q against selectivity) while this was being written.
  codes <- Rceattle:::LINKAGE_PROCESS_CODES
  slot_of <- function(nm) which(Rceattle:::.sim_state_codes(nm) == 1L) - 1L

  testthat::expect_identical(slot_of("recruitment"),  unname(codes[["recruitment"]]))
  testthat::expect_identical(slot_of("M"),            unname(codes[["M"]]))
  testthat::expect_identical(slot_of("growth"),       unname(codes[["growth"]]))
  testthat::expect_identical(slot_of("catchability"), unname(codes[["q"]]))
  testthat::expect_identical(slot_of("selectivity"),  unname(codes[["sel"]]))

  # Composition is a linkage process too, but carries priors rather than random
  # effects, so it has no simulate_state slot. The template must not index one.
  testthat::expect_length(Rceattle:::.sim_state_codes(TRUE), 5L)
  testthat::expect_true(codes[["comp"]] >= 5L)
})

# ---- the draw (needs a built model) -----------------------------------------

# Draw the linkage deviations n times with sigma pinned, returning a
# draw x slot matrix. sigma is pinned because these fixtures hold little
# information about a variance, so a fitted one collapses toward zero.
.re_draws <- function(fit, state, sigma, n, period = c(1, 0)) {
  obj <- Rceattle:::.sim_obj(fit)
  obj$fn(obj$par)
  par <- obj$env$last.par.best
  par[names(par) == "log_sigma_linkage"] <- log(sigma)
  obj$env$last.par.best <- par
  # Sized from the first draw, not from the parameter count: a rw group's first
  # deviation is mapped off, so it has one more slot than free parameter.
  out <- NULL
  for (i in seq_len(n)) {
    v <- as.numeric(Rceattle:::.sim_draw(obj, state = state, period = period)$beta_linkage_re_sim)
    if (is.null(out)) out <- matrix(NA_real_, n, length(v))
    out[i, ] <- v
  }
  out
}

.m_re_fit <- function(form) {
  suppressWarnings(Rceattle::fit_mod(
    data_list = make_test_data(), estimateMode = 3, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0),
    M1Fun = Rceattle::build_M1(
      M1_model = 1,
      linkages = list(M1 = Rceattle::linkage_spec(form, by = ~ species)))))
}

testthat::test_that("a random linkage is drawn only for the process asked for", {
  testthat::skip_on_cran()
  fit <- .m_re_fit(~ (1 | Year))
  ref <- {
    obj <- Rceattle:::.sim_obj(fit); obj$fn(obj$par)
    p <- obj$env$last.par.best
    p[names(p) == "log_sigma_linkage"] <- log(0.4)
    as.numeric(p[names(p) == "beta_linkage_re"])
  }

  # Every process except M leaves the M linkage where it was.
  for (p in list(FALSE, "recruitment", "growth", "catchability", "selectivity")) {
    d <- .re_draws(fit, Rceattle:::.sim_state_codes(p), 0.4, n = 5L)
    testthat::expect_equal(max(abs(sweep(d, 2, ref))), 0,
                           info = paste("process =", paste(p, collapse = ",")))
  }

  # M does draw, and differently each time.
  dm <- .re_draws(fit, Rceattle:::.sim_state_codes("M"), 0.4, n = 5L)
  testthat::expect_gt(max(abs(sweep(dm, 2, ref))), 0)
  testthat::expect_gt(nrow(unique(dm)), 1L)

  # The density covers the fitted window, so outside it nothing is drawn.
  dp <- .re_draws(fit, Rceattle:::.sim_state_codes("M"), 0.4, n = 5L, period = c(0, 0))
  testthat::expect_equal(max(abs(sweep(dp, 2, ref))), 0)
})

testthat::test_that("a random-walk linkage leaves its first deviation alone", {
  testthat::skip_on_cran()
  # The rw density scores first differences only, so the level of slot 0 is not
  # identified and the map pins it. Drawing it would move a parameter the
  # likelihood cannot see.
  d <- .re_draws(.m_re_fit(~ rw(1 | Year)), Rceattle:::.sim_state_codes("M"),
                 sigma = 0.3, n = 10L)
  testthat::expect_equal(stats::sd(d[, 1]), 0)
  testthat::expect_gt(stats::sd(as.numeric(d[, -1])), 0)
})

testthat::test_that("sim_mod() attaches the deviations it drew, and only those", {
  testthat::skip_on_cran()
  fit <- Rceattle::fit_mod(
    data_list = make_test_data(), estimateMode = 1, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))

  # Nothing redrawn -> nothing to recover, so no truth is attached.
  testthat::expect_null(attr(Rceattle::sim_mod(fit, simulate = TRUE), "process_sim"))

  tr <- attr(Rceattle::sim_mod(fit, simulate = TRUE, process = "recruitment"),
             "process_sim")
  testthat::expect_true(all(c("rec_dev", "init_dev") %in% names(tr)))
  # M was not drawn, so its fitted values are not offered as a truth.
  testthat::expect_null(tr$log_M1_dev)
  # The truth is the draw, not the fit it started from.
  testthat::expect_gt(max(abs(tr$rec_dev - fit$estimated_params$rec_dev)), 0)
})

testthat::test_that("a process the model cannot draw is reported, not fabricated", {
  testthat::skip_on_cran()
  # make_test_data() has M1_re = 0 and no growth linkage, so neither process has
  # a distribution to draw from. Returning the fitted log_M1_dev as a "truth"
  # would make self_test() report perfect recovery of a process that was never
  # simulated -- the same trap as the growth_re switch this release removes.
  fit <- Rceattle::fit_mod(
    data_list = make_test_data(), estimateMode = 1, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))

  testthat::expect_warning(
    m <- Rceattle::sim_mod(fit, simulate = TRUE, process = "M"),
    "no distribution to draw from")
  testthat::expect_null(attr(m, "process_sim")$log_M1_dev)
  testthat::expect_null(attr(m, "process_sim"))

  testthat::expect_warning(
    g <- Rceattle::sim_mod(fit, simulate = TRUE, process = "growth"),
    "no distribution to draw from")
  testthat::expect_null(attr(g, "process_sim"))

  # Recruitment always has a density, so it is drawn and reported as usual.
  r <- Rceattle::sim_mod(fit, simulate = TRUE, process = "recruitment")
  testthat::expect_false(is.null(attr(r, "process_sim")$rec_dev))
  # ...and carries no linkage vector, since this model has no random linkage.
  testthat::expect_null(attr(r, "process_sim")$beta_linkage_re)
})

testthat::test_that("self_test(process = ) returns each replicate's truth", {
  testthat::skip_on_cran()
  fit <- Rceattle::fit_mod(
    data_list = make_test_data(), estimateMode = 1, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))

  sims <- suppressWarnings(suppressMessages(
    Rceattle::self_test(fit, nsim = 2L, process = "recruitment", cores = 1L,
                        seed = 7L, debug = TRUE)))
  truth <- attr(sims, "process_sim")
  testthat::expect_length(truth, length(sims))
  testthat::expect_named(truth, names(sims))
  # Each replicate is its own draw, so no two share a truth.
  testthat::expect_false(identical(truth[[1]]$rec_dev, truth[[2]]$rec_dev))

  # Without process error there is no truth to carry.
  plain <- suppressWarnings(suppressMessages(
    Rceattle::self_test(fit, nsim = 1L, cores = 1L, seed = 7L, debug = TRUE)))
  testthat::expect_null(attr(plain, "process_sim"))
})

testthat::test_that("growth process error comes from the linkage, not the legacy array", {
  testthat::skip_on_cran()
  # A random linkage on a growth parameter is the only time-varying growth there
  # is. The legacy log_growth_par_devs array carried no density and was never
  # un-mapped, so it and the growth_re switch that advertised it were removed.
  # The fixture's survey is age-dimensioned, so data_check() warns that its CAAL
  # rows cannot be predicted. Unrelated to what is being tested here.
  fit <- suppressWarnings(Rceattle::fit_mod(
    data_list = make_test_data(growth = "vonBertalanffy"),
    estimateMode = 3, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0),
    growthFun = Rceattle::build_growth(
      fun = "vonBertalanffy",
      linkages = list(K = Rceattle::linkage_spec(~ (1 | Year), by = ~ species)))))

  # Assert absence directly. `fit$map` is list(mapFactor, mapList), so
  # `fit$map$log_growth_par_devs` is NULL whatever happens and
  # `all(is.na(NULL))` is TRUE -- a check that could never fail.
  testthat::expect_false("log_growth_par_devs" %in% names(fit$estimated_params))
  testthat::expect_false("log_growth_par_devs" %in% names(fit$map$mapList))
  testthat::expect_false("log_growth_par_devs" %in%
                           names(Rceattle::build_params(make_test_data())))

  d <- .re_draws(fit, Rceattle:::.sim_state_codes("growth"), sigma = 0.35, n = 10L)
  testthat::expect_gt(stats::sd(as.numeric(d)), 0)
})
