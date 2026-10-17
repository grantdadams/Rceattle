# `dsem_on` used to mean two things at once: "a SEM is attached" and "that SEM
# owns the recruitment deviations". They are the same question only while
# recruitment is the one process a SEM can drive. The moment one can drive
# another, attaching it would have:
#
#   * mapped out rec_dev and R_log_sd            (R/6-fit_mod.R)
#   * dropped rec_dev from `random`              (R/6-fit_mod.R)
#   * skipped the JNLL_REC_DEV density           (ceattle.cpp, section 13)
#   * gated off the rec_dev simulate draw        (ceattle.cpp)
#
# i.e. deleted the recruitment process. It would ALSO still have overwritten
# rec_dev from the latent states, so the density would have scored a value that
# is not its own parameter while R_log_sd sat free with zero gradient.
#
# None of that is reachable today: a sem naming no recdevs latent is refused in
# build_dsem_objects() and again in the template, so this is groundwork rather
# than a fix to a live defect. The point is that the refusal cannot be relaxed
# -- which is what driving growth or M needs -- until the two are separated.
#
# `dsem_drives_rec` separates them. These tests pin both halves: that a
# recruitment SEM still takes over recruitment exactly as before, and that a SEM
# which does NOT drive recruitment leaves the recruitment process intact while
# keeping its own GMRF density.

.pa_fixture <- function() {
  d <- Rceattle::BS2017SS
  d$env_data <- data.frame(Year = d$styr:d$endyr, BT = 0)
  d
}

.pa_fc <- function() {
  Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
}

# Which parameter BLOCKS TMB integrated out. `obj$env$random` is an index
# vector into `obj$env$par`, not a vector of names, and `fit$random_vars` is
# recorded only at estimateMode 0 or 1 -- these fits build without optimizing.
.pa_random_blocks <- function(fit) {
  ran <- fit$obj$env$random
  if (!length(ran)) return(character(0))
  unique(names(fit$obj$env$par)[ran])
}


testthat::test_that("a recruitment DSEM takes over recruitment, and says so", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = .pa_fixture(), inits = NULL, file = NULL, estimateMode = 3,
    random_rec = TRUE, msmMode = 0, dsem = Rceattle::build_DSEM(),
    fit_control = .pa_fc())))

  # The flag the template read.
  testthat::expect_equal(fit$dsem$tmb_inputs$data$drives_rec, 1L)
  testthat::expect_true(Rceattle:::.dsem_drives(fit, "recruitment"))

  # The GMRF scores the deviations, so the standalone recruitment density must
  # not: scoring both would count them twice.
  jc <- fit$quantities$jnll_comp
  testthat::expect_equal(sum(jc["Recruitment deviates", ]), 0)
  testthat::expect_true(sum(jc["DSEM", ]) != 0)

  # rec_dev is derived from the latent states, so it is pinned and the field is
  # the random block instead.
  testthat::expect_true(all(is.na(fit$map$mapList$rec_dev)))
  ran <- .pa_random_blocks(fit)
  testthat::expect_false("rec_dev" %in% ran)
  testthat::expect_true("dsem_x_tj" %in% ran)
})


testthat::test_that(
  "a DSEM that does not drive recruitment leaves recruitment intact", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d <- .pa_fixture()

  # Fit normally, then refit from THAT fit's own built objects with one field
  # flipped, so the two runs differ in `drives_rec` and nothing else. Building
  # the objects separately is not equivalent: build_dsem_objects() reads
  # `data_list$random_rec`, which fit_mod() sets from its own argument but a raw
  # bundled data_list carries as 0, and pins each species' sigmaR entry in the
  # map when it is false. On the default IID sem those entries are the only free
  # paths, so the SEM comes back with nothing to estimate and this would be a
  # weaker test than it looks.
  ref <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 3,
    random_rec = TRUE, msmMode = 0, dsem = Rceattle::build_DSEM(),
    fit_control = .pa_fc())))

  # Both the template (dsem_drives_rec) and R (.dsem_drives) read this one
  # field, so there is no second place for them to disagree.
  #
  # Reachable only by setting it: a SEM that names recdevs latents while
  # declaring it does not drive recruitment is incoherent as a model, and a
  # real non-recruitment SEM would carry no recdevs latents at all. What is
  # under test is the gate, not this configuration -- and because the gate now
  # also stops the recruitment takeover, the two densities here score different
  # quantities: the GMRF scores the latent field, JNLL_REC_DEV scores the free
  # rec_dev parameter.
  b_dsem <- ref$dsem
  b_dsem$tmb_inputs$data$drives_rec <- 0L

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 3,
    random_rec = TRUE, msmMode = 0, dsem = b_dsem, fit_control = .pa_fc())))

  # The control: the free path coefficients survived the round trip, so this is
  # a SEM with something to estimate rather than one with every path pinned.
  testthat::expect_equal(fit$map$mapList$dsem_beta_z,
                         ref$map$mapList$dsem_beta_z)

  testthat::expect_false(Rceattle:::.dsem_drives(fit, "recruitment"))

  jc <- fit$quantities$jnll_comp

  # The headline regression: recruitment is still a scored process.
  testthat::expect_true(sum(jc["Recruitment deviates", ]) != 0)

  # And the SEM's own density is untouched -- compared against the unflipped
  # fit, not merely asserted non-zero. This is the other half of the same
  # defect: gating the GMRF on dsem_drives_rec instead of dsem_on would drop
  # `jnll_dsem` from the objective entirely, and that scalar carries the
  # covariate measurement density as well as the GMRF, so the field would be
  # declared random with nothing scoring it at all.
  testthat::expect_true(sum(jc["DSEM", ]) != 0)
  testthat::expect_equal(jc["DSEM", ],
                         ref$quantities$jnll_comp["DSEM", ])

  # rec_dev is an ordinary estimated deviation again: free in the map, and
  # integrated out under random_rec = TRUE.
  testthat::expect_false(all(is.na(fit$map$mapList$rec_dev)))
  testthat::expect_true("rec_dev" %in% .pa_random_blocks(fit))

  # The assertion that actually discriminates. The first version of this change
  # left the recruitment takeover in section 5.5b gated on dsem_on, so rec_dev
  # was scored by JNLL_REC_DEV and THEN overwritten from the latent states, and
  # R_sd was overridden from the SEM's variance path while R_log_sd stayed
  # estimated. Every assertion above passed in that state: the density row was
  # non-zero, the map was free, and the block was declared random. What gives it
  # away is the gradient -- a parameter that is overwritten before it is used
  # carries exactly zero -- so a free parameter with no gradient is the defect,
  # and these two coordinates are where it showed.
  # The JOINT gradient, over fixed and random parameters together: rec_dev is a
  # random effect here, so it is in obj$env$par and not in obj$par, and
  # obj$gr() would not see it at all.
  g   <- fit$obj$env$f(fit$obj$env$last.par, order = 1)
  nms <- names(fit$obj$env$par)
  testthat::expect_true(all(is.finite(g)))
  for (blk in c("rec_dev", "R_log_sd")) {
    idx <- which(nms == blk)
    testthat::expect_true(length(idx) > 0)
    testthat::expect_true(any(g[idx] != 0))
  }
  # dsem_x_tj is deliberately NOT asserted here. With recruitment detached the
  # field's only term is the GMRF, whose gradient at the all-zero starting
  # values is zero by construction -- correct, and the reason this arm says
  # nothing about a field that really drives a process.

  # summary() still reports the SEM path estimates: they are the SEM's own
  # coefficients whatever process it drives, so they are gated on the SEM
  # existing, not on what it drives.
  s <- summary(fit)
  testthat::expect_false(is.null(s$dsem_coefficients))
  testthat::expect_true(nrow(s$dsem_coefficients) > 0)
})


# `drives_rec` postdates every saved $dsem object in the assessment repos, and
# those refit through .dsem_tmb_inputs() on retrospective(), jitter(), project()
# and plain refits. A SEM built before the flag existed could only have driven
# recruitment, so its absence must read as 1 -- in BOTH the template inputs and
# .dsem_drives(), or a saved fit would refit as a different model than it
# reports.
testthat::test_that("a dsem object with no drives_rec flag reads as recruitment", {
  testthat::skip_if_not_installed("dsem")

  d <- .pa_fixture()
  built <- Rceattle::build_dsem_objects(
    dsem_settings = Rceattle::build_DSEM(), data_list = d)

  # Exactly what a pre-flag saved object looks like.
  built$tmb_inputs$data$drives_rec <- NULL

  testthat::expect_equal(
    Rceattle:::.dsem_tmb_inputs(built)$data$dsem_drives_rec, 1L)
  testthat::expect_true(
    Rceattle:::.dsem_drives(list(dsem = built), "recruitment"))
})


testthat::test_that(
  ".dsem_drives() is FALSE without a DSEM and refuses an unknown process", {
  testthat::expect_false(Rceattle:::.dsem_drives(NULL, "recruitment"))
  testthat::expect_false(Rceattle:::.dsem_drives(list(), "recruitment"))
  # A process name that is not in the registry is a caller mistake, not FALSE:
  # returning FALSE would read as "the SEM does not drive growth" on a tree
  # where growth is not a process a SEM can drive at all.
  testthat::expect_error(Rceattle:::.dsem_drives(list(), "growth"), "must be one of")
})
