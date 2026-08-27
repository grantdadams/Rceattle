# retrospective() on a DSEM.
#
# A peel does not shorten the model: retrospective() sets endyr_peel and turns
# off DATA after it, so the latent states still span every year. They stay free
# and in `random`, and the Laplace approximation integrates the peeled-year
# states out against the GMRF prior. That is the peeled marginal likelihood.
#
# Mirroring what rec_dev does -- zero the tail, map it out -- would be wrong.
# Pinning is inert for an INDEPENDENT deviate, but DSEM states are coupled
# through the RAM, so pinned zeros stay in the quadratic form and shrink the
# terminal retained state (1 + rho^2 inflated precision). Terminal recruitment
# drives terminal SSB, so Mohn's rho would measure the peel's own artifact.
#
# The sem below is LAGGED and CARRIES A COVARIATE on purpose. A default
# build_DSEM() sem is IID and has no covariate column, and can detect neither
# failure: pinning is inert without coupling, and a covariate bug needs one.

dsem_retro_data <- function() {
  d <- Rceattle::BS2017SS
  yrs <- d$styr:d$endyr
  set.seed(1)
  d$env_data <- data.frame(
    Year = yrs, temp = as.numeric(scale(cumsum(stats::rnorm(length(yrs))))))
  d
}

# No `temp <-> temp`: the builder adds V[temp] itself, and giving an exogenous
# variable two variance terms makes the precision singular and the objective NaN.
DSEM_RETRO_SEM <- "
  recdevs1 -> recdevs1, 1, rho_R,     0.3
  temp     -> recdevs1, 0, temp_to_R, 0.2
  recdevs1 <-> recdevs1, 0, sigmaR1,  0.6
  recdevs2 <-> recdevs2, 0, sigmaR2,  0.6
  recdevs3 <-> recdevs3, 0, sigmaR3,  0.6
"

testthat::test_that("a DSEM peel reproduces a model that genuinely ends at endyr_peel", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d  <- dsem_retro_data()
  fc <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  mk <- function() Rceattle::build_DSEM(sem = DSEM_RETRO_SEM, family = "fixed")

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = TRUE, msmMode = 0, dsem = mk(), fit_control = fc)))

  retro <- suppressWarnings(suppressMessages(
    Rceattle::retrospective(fit, peels = 1, cores = 1, getsd = FALSE)))
  ml <- retro$Rceattle_list
  # rev(c(list(Rceattle), peels)): the LAST element is the input object, so an
  # "unpeeled peel matches the parent" assertion would compare a fit to itself.
  testthat::expect_equal(length(ml), 2L)
  pk <- ml[[1]]
  endyr_peel <- pk$data_list$endyr_peel
  testthat::expect_equal(endyr_peel, d$endyr - 1L)

  # A model that genuinely ends at endyr_peel: the peeled year does not exist.
  # Integrating a state out must equal never having had it, for a Gaussian
  # process -- so the retained years must match. A peel that PINNED its
  # peeled-year state would inject information this model does not have.
  dt <- d
  dt$index_data <- dt$index_data[dt$index_data$Year <= endyr_peel, ]
  dt$comp_data  <- dt$comp_data[dt$comp_data$Year <= endyr_peel, ]
  dt$catch_data$Catch[dt$catch_data$Year > endyr_peel] <- 0
  dt$endyr <- endyr_peel
  dt$projyr <- endyr_peel
  dt$env_data <- dt$env_data[dt$env_data$Year <= endyr_peel, ]

  ends_at_peel <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = dt, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = TRUE, msmMode = 0, dsem = mk(), fit_control = fc)))

  keep <- seq_len(endyr_peel - d$styr + 1L)
  sp <- 1L   # the species carrying the lag and the covariate
  R_peel <- pk$quantities$R[sp, keep]
  R_ends <- ends_at_peel$quantities$R[sp, keep]
  R_par  <- fit$quantities$R[sp, keep]

  # Marginalized: matches the shorter model to optimizer tolerance.
  testthat::expect_equal(as.numeric(R_peel), as.numeric(R_ends), tolerance = 1e-3)
  testthat::expect_equal(pk$quantities$ssb[sp, keep],
                         ends_at_peel$quantities$ssb[sp, keep], tolerance = 1e-3)

  # ...and the comparison is not vacuous: the peel must DIFFER from the parent,
  # which is the retrospective signal itself.
  testthat::expect_gt(max(abs(R_peel - R_par)) / max(abs(R_par)), 0.01)
})

testthat::test_that("a DSEM peel keeps its covariate and reports a finite Mohn's rho", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d  <- dsem_retro_data()
  fc <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = TRUE, msmMode = 0,
    dsem = Rceattle::build_DSEM(sem = DSEM_RETRO_SEM, family = "fixed"),
    fit_control = fc)))

  retro <- suppressWarnings(suppressMessages(
    Rceattle::retrospective(fit, peels = 1, cores = 1, getsd = FALSE)))
  pk <- retro$Rceattle_list[[1]]

  # Under family = "fixed" the covariate column of x_tj IS the environmental
  # data, held fixed by the map. Zeroing x_tj to "peel" it would delete the
  # covariate outright rather than withhold it.
  xt <- as.matrix(pk$estimated_params$dsem_x_tj)
  testthat::expect_equal(ncol(xt), 4L)              # 3 recdev columns + temp
  testthat::expect_gt(max(abs(xt[, 4])), 1e-6)

  rho <- retro$mohns
  ssb_rho <- rho[rho$Object == "ssb", -(1:3), drop = FALSE]
  testthat::expect_true(any(is.finite(as.matrix(ssb_rho))))
})
