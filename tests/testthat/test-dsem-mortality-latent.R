# A DSEM natural-mortality latent must land on the right M_linkage_offset cell
# -- right species, right sex, right year -- and must be FLAT over ages. The
# flatness is the whole distinction from `M1_re`: one latent says "mortality was
# high for this stock in this year", while an age-structured deviation is what
# M1_re estimates. Nothing but the reported tensor can check that: the registry
# tests next door pin the latent NAMES, and the index map built from them is
# pinned here, but neither says where a latent arrives in the tensor or how far
# over ages it spreads.
#
# The latent is a log-scale multiplier on residual natural mortality M1, so
# section 6.1 reads it inside the exponent:
#
#   M1_at_age = exp(log_M1 + log_M1_dev + M1_mult.sum() + M_linkage_offset)
#               + M_linkage_offset_nat
#
# which the predation test below asserts rather than assumes.

testthat::skip_on_cran()

.m_fit <- function(sem, nyrs = 12L, nages = 8L) {
  e <- new.env(parent = asNamespace("Rceattle"))
  sys.source(testthat::test_path("helpers.R"), envir = e)
  d <- e$make_test_data(nyrs = nyrs, nprojyrs = 3L, nages = nages, minage = 1L,
                        growth = "vonBertalanffy", seed = 1)
  d$random_rec <- TRUE
  set.seed(4)
  d$env_data <- data.frame(Year = seq(d$styr, d$endyr), BT = stats::rnorm(nyrs))
  suppressWarnings(Rceattle::fit_mod(
    data_list = d, estimateMode = 3, msmMode = 0,
    growthFun = Rceattle::build_growth(fun = "vonBertalanffy"),
    dsem = Rceattle::build_DSEM(sem = sem),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0)))
}

# dsem_x_tj is column-major, so cell (t, j) is parameter (j-1)*nrow + t. A spike
# is how the age and year spread get tested: at the starting values the latent
# field is all zeros, so the offset is the same constant in every year and
# flatness would hold trivially.
.spike <- function(fit, column, t, delta = 0.5) {
  obj <- fit$obj
  p <- obj$env$last.par.best
  if (is.null(p)) p <- obj$env$par
  xi <- which(names(p) == "dsem_x_tj")
  nr <- nrow(fit$dsem$tmb_inputs$data$y_tj)
  j  <- which(colnames(fit$dsem$tmb_inputs$data$y_tj) == column)
  testthat::expect_length(j, 1L)
  k <- xi[(j - 1L) * nr + t]
  p[k] <- p[k] + delta
  obj$report(p)
}

REC1 <- "recdevs1 <-> recdevs1, 0, sigmaR1, 1\n"
SEM_M <- paste0(REC1, "mdevs1 <-> mdevs1, 0, sigM1, 0.2\n")


testthat::test_that("an M spike moves every age of its own year and no other year", {
  testthat::skip_if_not_installed("dsem")

  # The acceptance criterion that separates M from growth: M1-at-age is read
  # year by year with no cohort recursion, so a spike at year t moves ALL ages
  # in year t by exactly the spike and leaves every other year untouched.
  # Growth's L1 and Linf instead propagate forward through lagkappa/laglinf.
  fit <- .m_fit(SEM_M)
  nyd <- fit$dsem$tmb_inputs$data$nyrs_dsem
  nag <- fit$data_list$nages[1]
  TS  <- 5L
  testthat::expect_gt(nyd, TS)

  r0 <- fit$obj$report(fit$obj$env$par)
  r1 <- .spike(fit, "mdevs1", TS, delta = 0.5)
  d  <- (r1$M_linkage_offset - r0$M_linkage_offset)[1, 1, seq_len(nag),
                                                    seq_len(nyd)]

  # Flat over ages, and exactly the spike -- the same scalar added nages times.
  testthat::expect_equal(as.numeric(d[, TS]), rep(0.5, nag), tolerance = 1e-12)
  # Every other year is untouched: no cohort carry-over, no leak backwards.
  testthat::expect_equal(sum(abs(d[, -TS])), 0, tolerance = 1e-14)

  # The absolute offset is flat over ages too, not only the difference.
  off <- r1$M_linkage_offset[1, 1, seq_len(nag), seq_len(nyd)]
  for (yr in seq_len(nyd)) {
    testthat::expect_equal(diff(range(off[, yr])), 0, tolerance = 1e-14)
  }
  # And the spike made the offset year-varying, so that is not trivial.
  testthat::expect_gt(diff(range(off[1, ])), 1e-6)
})


testthat::test_that("the M latent lands on the cell the index map names", {
  testthat::skip_if_not_installed("dsem")

  # The identity section 5.5b-M claims, asserted absolutely rather than as a
  # difference: a difference cancels the bias correction, which review showed
  # passes identically with the -margvar/2 term present, absent or sign-flipped.
  fit <- .m_fit(SEM_M)
  r   <- fit$obj$report(fit$obj$env$par)
  dt  <- fit$dsem$tmb_inputs$data
  nyd <- dt$nyrs_dsem
  ba  <- fit$data_list$bias_adjust_proc

  testthat::expect_length(dt$mdev_col, 1L)
  for (m in seq_along(dt$mdev_col)) {
    col <- dt$mdev_col[m] + 1L
    sp  <- dt$mdev_sp[m]  + 1L
    sx  <- dt$mdev_sex[m] + 1L
    want <- r$dsem_z_tj[seq_len(nyd), col] -
      ba * r$dsem_margvar_tj[seq_len(nyd), col] / 2
    for (age in seq_len(fit$data_list$nages[sp])) {
      got <- r$M_linkage_offset[sp, sx, age, seq_len(nyd)]
      testthat::expect_equal(as.numeric(got), as.numeric(want),
                             tolerance = 1e-12)
    }
  }

  # The bias correction is doing something here, so the identity is not trivial.
  mv <- r$dsem_margvar_tj[seq_len(nyd), dt$mdev_col + 1L, drop = FALSE]
  testthat::expect_gt(ba, 0)
  testthat::expect_gt(max(abs(mv)), 0)
})


testthat::test_that("the latent survives the niter loop and reaches M2 under predation", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("dsem")

  # Both halves need msmMode > 0. `niter` is forced to 1 when msmMode == 0
  # (ceattle.cpp:189), so no single-species fixture can see whether the
  # multispecies loop re-zeros or re-accumulates M_linkage_offset -- the tensor
  # is written once, before the loop opens. And M2 is predator demand over
  # available prey, which depends on survivorship, which depends on M1, so the
  # latent DOES move M2 even though no code writes it there.
  e <- new.env(parent = asNamespace("Rceattle"))
  sys.source(testthat::test_path("helpers.R"), envir = e)
  sys.source(testthat::test_path("helpers-make-msm-data.R"), envir = e)
  d <- e$make_msm_test_data()$data_list
  d$random_rec <- TRUE
  set.seed(4)
  d$env_data <- data.frame(Year = seq(d$styr, d$endyr),
                           BT = stats::rnorm(d$endyr - d$styr + 1L))
  fit <- suppressWarnings(Rceattle::fit_mod(
    data_list = d, estimateMode = 3, msmMode = 1, niter = 3,
    dsem = Rceattle::build_DSEM(sem = paste0(
      "recdevs1 <-> recdevs1, 0, sR1, 1\nrecdevs2 <-> recdevs2, 0, sR2, 1\n",
      "mdevs1 <-> mdevs1, 0, sM1, 0.2\n")),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0)))
  dt <- fit$dsem$tmb_inputs$data
  testthat::expect_length(dt$mdev_col, 1L)
  nyd <- dt$nyrs_dsem
  nag <- fit$data_list$nages[1]

  r0 <- fit$obj$report(fit$obj$env$par)
  r1 <- .spike(fit, "mdevs1", 10L, delta = 0.5)

  # The offset is still exactly the spike after three predation iterations.
  d_off <- (r1$M_linkage_offset - r0$M_linkage_offset)[1, 1, seq_len(nag), ]
  testthat::expect_equal(as.numeric(d_off[, 10L]), rep(0.5, nag),
                         tolerance = 1e-12)
  testthat::expect_equal(sum(abs(d_off[, -10L])), 0, tolerance = 1e-14)
  # And M1-at-age carries it, so no pass re-zeroed or doubled the tensor.
  testthat::expect_equal(
    as.numeric(r1$M1_at_age[1, 1, seq_len(nag), 10L] /
                 r0$M1_at_age[1, 1, seq_len(nag), 10L]),
    rep(exp(0.5), nag), tolerance = 1e-10)

  # M2 moves, which is what makes "the latent does not touch M2" a statement
  # about the code and not about the model.
  testthat::expect_gt(max(abs(r1$M2_at_age - r0$M2_at_age)), 1e-8)
})


testthat::test_that("a recruitment-only DSEM leaves the M offset at zero", {
  testthat::skip_if_not_installed("dsem")
  # The no-op guarantee for every existing DSEM fit: with no M latent the index
  # map is empty, so section 5.5b-M's loop never executes.
  fit <- .m_fit(REC1)
  testthat::expect_length(fit$dsem$tmb_inputs$data$mdev_col, 0L)
  testthat::expect_equal(
    sum(abs(fit$obj$report(fit$obj$env$par)$M_linkage_offset)), 0)
})


testthat::test_that("an M latent beside a growth latent reaches only its own tensor", {
  testthat::skip_if_not_installed("dsem")

  # The two writers share section 5.5b and both use `+=`, so a crossed index
  # would put mortality into growth or the reverse without changing either
  # tensor's magnitude enough to notice.
  fit <- .m_fit(paste0(SEM_M, "l1devs1 <-> l1devs1, 0, sigL1, 0.3\n"))
  r   <- fit$obj$report(fit$obj$env$par)
  dt  <- fit$dsem$tmb_inputs$data
  nyd <- dt$nyrs_dsem

  testthat::expect_length(dt$mdev_col, 1L)
  testthat::expect_length(dt$growth_dev_col, 1L)
  testthat::expect_false(dt$mdev_col == dt$growth_dev_col)

  # Growth's L1 slot moved; no other growth slot did.
  g_par <- dt$growth_dev_par + 1L
  testthat::expect_gt(max(abs(r$growth_linkage_offset[1, 1, seq_len(nyd), g_par])), 1e-6)
  other <- setdiff(seq_len(dim(r$growth_linkage_offset)[4]), g_par)
  testthat::expect_equal(
    sum(abs(r$growth_linkage_offset[1, 1, seq_len(nyd), other])), 0)

  # M moved, and by its OWN column, not growth's.
  off <- r$M_linkage_offset[1, 1, 1, seq_len(nyd)]
  ba  <- fit$data_list$bias_adjust_proc
  want <- r$dsem_z_tj[seq_len(nyd), dt$mdev_col + 1L] -
    ba * r$dsem_margvar_tj[seq_len(nyd), dt$mdev_col + 1L] / 2
  testthat::expect_equal(as.numeric(off), as.numeric(want), tolerance = 1e-12)
})


# ---- process error: the M field is drawn, selectively ----------------------
# These pin what adversarial review found missing when M joined the registry.
# `.dsem_draw_process()` partitioned columns into recruitment, growth and "latent
# covariate", and a covariate column is drawn UNCONDITIONALLY -- so an M column
# fell through and `process = "recruitment"` silently redrew M as well. Every
# replicate was then generated under M process error the caller did not ask for.

testthat::test_that("sim_mod() draws the M field, and only what was asked for", {
  testthat::skip_if_not_installed("dsem")

  fit <- .m_fit(SEM_M)
  if (is.null(fit$obj$env$last.par.best))
    fit$obj$env$last.par.best <- fit$obj$env$par
  dt  <- fit$dsem$tmb_inputs$data
  mcol <- dt$mdev_col + 1L
  rcol <- dt$rec_dev_col + 1L

  sd_of <- function(d, col) stats::sd(attr(d, "process_sim")$dsem_x_tj[, col])

  set.seed(11)
  d_m <- suppressWarnings(Rceattle::sim_mod(fit, simulate = TRUE, process = "M"))
  set.seed(11)
  d_r <- suppressWarnings(Rceattle::sim_mod(fit, simulate = TRUE,
                                            process = "recruitment"))

  # Asking for M moves M and holds recruitment; asking for recruitment does the
  # reverse. The second half is the regression: it used to redraw M too.
  testthat::expect_gt(sd_of(d_m, mcol), 1e-8)
  testthat::expect_equal(sd_of(d_m, rcol), 0, tolerance = 1e-12)
  testthat::expect_gt(sd_of(d_r, rcol), 1e-8)
  testthat::expect_equal(sd_of(d_r, mcol), 0, tolerance = 1e-12)
})


testthat::test_that("an M-only sem can draw its process", {
  testthat::skip_if_not_installed("dsem")

  # An M-only sem reached the "this fit records no latent process columns"
  # refusal, which counted recruitment and growth only -- the twin of the
  # drives-nothing bug in the build path.
  fit <- .m_fit("mdevs1 <-> mdevs1, 0, sigM1, 0.2\n")
  if (is.null(fit$obj$env$last.par.best))
    fit$obj$env$last.par.best <- fit$obj$env$par
  mcol <- fit$dsem$tmb_inputs$data$mdev_col + 1L

  set.seed(12)
  d <- suppressWarnings(Rceattle::sim_mod(fit, simulate = TRUE, process = TRUE))
  testthat::expect_gt(
    stats::sd(attr(d, "process_sim")$dsem_x_tj[, mcol]), 1e-8)
})


testthat::test_that("a DSEM M latent is not reported as having no distribution", {
  testthat::skip_if_not_installed("dsem")

  # .check_dsem_m_support() REQUIRES M1_re = 0 to allow the latent, and the
  # warning's M arm tested exactly `M1_re = 0 and no random linkage on M` -- so
  # it fired on every M-driving DSEM, saying nothing was redrawn while the GMRF
  # draw was happening. Same false alarm the growth arm had.
  fit <- .m_fit(SEM_M)
  if (is.null(fit$obj$env$last.par.best))
    fit$obj$env$last.par.best <- fit$obj$env$par
  testthat::expect_equal(fit$data_list$M1_re, 0)

  w <- testthat::capture_warnings(
    Rceattle::sim_mod(fit, simulate = TRUE, process = "M"))
  testthat::expect_false(any(grepl("no distribution to draw from", w)))

  # The arm still fires for a model that really has no M distribution.
  f_rec <- .m_fit(REC1)
  if (is.null(f_rec$obj$env$last.par.best))
    f_rec$obj$env$last.par.best <- f_rec$obj$env$par
  w2 <- testthat::capture_warnings(
    Rceattle::sim_mod(f_rec, simulate = TRUE, process = "M"))
  testthat::expect_true(any(grepl("no distribution to draw from", w2)))
})


testthat::test_that("an M latent is refused beside a RANDOM M1 linkage, not a fixed one", {
  testthat::skip_if_not_installed("dsem")

  # A random M1 linkage writes the SAME tensor as the latent, each scored by its
  # own density. data_check() soft-deprecates M1_re 1-5 in favour of that
  # grammar, so a user following the advice lands here. A FIXED linkage is the
  # intended combination and composes through `+=`.
  e <- new.env(parent = asNamespace("Rceattle"))
  sys.source(testthat::test_path("helpers.R"), envir = e)
  mk <- function(frm) {
    d <- e$make_test_data(nyrs = 12L, nprojyrs = 3L, nages = 8L, minage = 1L,
                          growth = "vonBertalanffy", seed = 1)
    d$random_rec <- TRUE
    set.seed(4)
    d$env_data <- data.frame(Year = seq(d$styr, d$endyr), BT = stats::rnorm(12))
    d$M1_linkages <- Rceattle::build_M1(
      M1_model = 0, M1_re = 0,
      linkages = list(M1 = Rceattle::linkage_spec(formula = frm,
                                                  param = "M1")))$linkages
    d
  }
  sem <- Rceattle::build_DSEM(sem = SEM_M)

  testthat::expect_error(
    build_dsem_objects(dsem_settings = sem, data_list = mk(~ (1 | Year))),
    "RANDOM M1 linkage")
  # The fixed one is allowed and the index map is built.
  g <- build_dsem_objects(dsem_settings = sem, data_list = mk(~ BT))
  testthat::expect_length(g$tmb_inputs$data$mdev_col, 1L)
})


testthat::test_that("the index maps are refused at disagreeing lengths", {
  testthat::skip_if_not_installed("dsem")

  # Provenance: measured on GOA2018SS with a sem naming `mdevs1` + `mdevs2f`,
  # whose correct sum|M_linkage_offset| is 26.0400. A one-entry `mdev_sp`
  # against a two-entry `mdev_col` RAN and gave 16.8000 -- both latents on
  # species 1, one species' offset halved and the other's doubled -- and an
  # empty `mdev_sp`/`mdev_sex` segfaulted. The maps are parallel vectors
  # supplied as independent fields, so only a length check catches this; the
  # value guards cannot see it. `growth_dev_*` had the same hole.
  fit <- .m_fit(SEM_M)
  d <- fit$data_list
  refit <- function(dsem_obj) suppressWarnings(Rceattle::fit_mod(
    data_list = d, estimateMode = 3, msmMode = 0, dsem = dsem_obj,
    growthFun = Rceattle::build_growth(fun = "vonBertalanffy"),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0)))

  g <- fit$dsem
  testthat::expect_length(g$tmb_inputs$data$mdev_col, 1L)

  bad_m <- g; bad_m$tmb_inputs$data$mdev_sp <- integer(0)
  testthat::expect_error(refit(bad_m), "mdev_col/sp/sex have 1/0/1 entries")

  bad_m2 <- g; bad_m2$tmb_inputs$data$mdev_col <- c(0L, 0L)
  testthat::expect_error(refit(bad_m2), "parallel and must agree")

  # Same guard covers the growth maps.
  gg <- .m_fit(paste0(SEM_M, "l1devs1 <-> l1devs1, 0, sigL1, 0.3\n"))$dsem
  bad_g <- gg; bad_g$tmb_inputs$data$growth_dev_sex <- integer(0)
  testthat::expect_error(refit(bad_g),
                         "growth_dev_col/sp/sex/par have 1/1/0/1 entries")
})
